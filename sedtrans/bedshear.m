% Calculate bed shear stress

% properties = [g, rho, nu, kappa, s, d50]

function tau = bedshear(h, uav, vav, Hw, Tw, theta, properties, varargin)

JJ = size(h,1);
if size(uav,1) ~= JJ, disp('Input vectors must be same length'), return, end
if size(vav,1) ~= JJ, disp('Input vectors must be same length'), return, end
if size(Hw,1) ~= JJ, disp('Input vectors must be same length'), return, end
if size(Tw,1) ~= JJ, disp('Input vectors must be same length'), return, end
if size(theta,1) ~= JJ, disp('Input vectors must be same length'), return, end

if ~isempty(varargin)
    if length(varargin) == 1
        output = varargin{1};
    else
        disp('Incorrect varargin specification')
    end
else
    output = 'rms';
end

g = properties(1);
rho = properties(2);
nu = properties(3);
kappa = properties(4);
s = properties(5);
d50 = properties(6);
ks = 2.5 * d50;
U = sqrt(uav.^2 + vav.^2);
Uvec = [ones(JJ,1), zeros(JJ,1)];
tauc = zeros(JJ,1);
ind = find(U>0&h>ks);
if ~isempty(ind)
    Uvec(ind,:) = [uav(ind)./U(ind), vav(ind)./U(ind)];
    ufcr = kappa * U ./ log(11 * h / ks);
    ufcs = ufcr;
    error = 9999;
    n = 1;
    while n < 10
        ufcsnew = real(kappa * U(ind) ./ log(3.8 * h(ind) .* ufcs(ind) / nu));
        error = sqrt((ufcsnew - ufcs(ind)).^2./ufcsnew.^2) .* (ufcsnew > 0.001);
        ufcs(ind) = ufcsnew;
        n = n + 1;
        if max(error) < 0.001, break, end
    end
    if n == 10, disp('Smooth regime calcs failed to converge'), end
    ufc = max(ufcr,ufcs);
    tauc = rho*ufc.^2;
end

tauw_max = zeros(JJ,1);
Tw(Tw<0.1) = 0.1;
k0h = 4*pi^2./(g*Tw.^2) .* h;
kh =kh_calc(k0h);
Aw = 0.5*Hw./sinh(kh);
ind = find(Aw>0.001);
if ~isempty(ind)
    Uw = 2*pi*Aw(ind)./Tw(ind);
    Rew = Uw.*Aw(ind)/nu;
    fws = 0.035*Rew.^(-0.16);
    fwr = exp(5.21*(ks./Aw(ind)).^0.194-5.98);
    fw = max(fws,fwr);
    tauw_max(ind) = 0.5*rho*fw.*Uw.^2;
end

taucw_mean = zeros(JJ,1);
ind = find(U>0);
taucw_mean(ind) = tauc(ind).*(1+1.2*(tauw_max(ind)./(tauc(ind)+tauw_max(ind))).^3.2);
taucw_max = sqrt((taucw_mean+tauw_max.*cos(theta/180*pi)).^2 + ...
                (tauw_max.*sin(theta/180*pi)).^2);
taucw_rms = sqrt(taucw_mean.^2+0.5*tauw_max.^2);

switch lower(output)
    case 'mean'
        tau = [taucw_mean, taucw_mean] .* Uvec;
    case 'max'
        tau = [taucw_max, taucw_max] .* Uvec;
    case 'rms'
        tau = [taucw_rms, taucw_rms] .* Uvec;
    otherwise
        disp('Incorrect output specification'), return
end