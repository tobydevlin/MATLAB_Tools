% [out in] = bijker(in)
%
% in = bijker([]) (empty input data structure)
%
% Calculates sediment transport using method of Bijker (1967, 1971).
% Input and output are via data structures 'in' and 'out'.
% The input data can be entered via command line prompts by passing '[]'
% to the function, or alternatively an empty input data structure can be
% generated.
% The water depth field 'in.h' must be present, and can be a scalar,
% a vector or a matrix.  The other inputs can be omitted and will assume
% default values.
% Sediment transport rates returned in the 'out' data structure are in
% kg/s.

function varargout = bijker(in)

% Initialise output structure
out = struct('Qsx',[],'Qsy',[],'Qs_tot',[],'Qs_susp',[],'Qs_bed',[]);

% If input is empty ask for manual input
if isempty(in)
    in = struct('h',[],'Vx',[],'Vy',[],'Hrms',[],'T',[],...
        'd50',[],'d90',[],'ws',[],'ks',[],'rhow',[],'rhos',[],'b',[]);
    in.h = input(['Input depth vector, h\n',...
        '(or hit return for empty input structure), :']);
    if isempty(in.h)
        varargout{1} = out;
        varargout{2} = in;
        return
    end
    in.Vx = input('Input x-direction current, Vx: ');   
    in.Vy = input('Input x-direction current, Vy: ');
    in.Hrms = input('Input rms wave height, Hrms: ');
    in.T = input('Input wave period, T: ');
    in.d50 = input('Input sediment size, d50: ');
    in.d90 = input('Input sediment size, d90: ');
    in.ws = input('Input sediment fall velocity, ws: ');
    in.ks = input('Input roughness, ks: '); 
    in.rhow = input('Input fluid density, rhow: '); 
    in.rhos = input('Input sediment density, rhos: ');
    in.b = input('Input empirical coefficient, b (approx. 5): ');
end

% Check input and provide default values where necessary
if (~isfield(in,'h')||isempty(in.h)), disp('Must have depth input'), return, end
if (~isfield(in,'Vx')||isempty(in.Vx)), in.Vx = zeros(size(in.h)); end
if all(size(in.Vx)~=size(in.h)), disp('size(Vx) must equal size(h)'), return, end
if (~isfield(in,'Vy')||isempty(in.Vy)), in.Vy = zeros(size(in.h)); end
if all(size(in.Vy)~=size(in.h)), disp('size(Vy) must equal size(h)'), return, end
if (~isfield(in,'Hrms')||isempty(in.Hrms)), in.Hrms = zeros(size(in.h)); end
if all(size(in.Hrms)~=size(in.h)), disp('size(Hrms) must equal size(h)'), return, end
if (~isfield(in,'T')||isempty(in.T)), in.T = zeros(size(in.h)); end
if all(size(in.T)~=size(in.h)), disp('size(T) must equal size(h)'), return, end
if (~isfield(in,'d50')||isempty(in.d50)), in.d50 = 0.00025; end
if (~isfield(in,'d90')||isempty(in.d90)), in.d90 = 0.00050; end
if (~isfield(in,'ws')||isempty(in.ws)), in.ws = 0.030; end
if (~isfield(in,'ks')||isempty(in.ks)), in.ks = 0.05; end
if (~isfield(in,'rhow')||isempty(in.rhow)), in.rhow = 1025; end
if (~isfield(in,'rhos')||isempty(in.rhos)), in.rhos = 2650; end
if (~isfield(in,'b')||isempty(in.b)), in.b = 5; end

% Calculate velocity magnitude and unit vector
[M,N] = size(in.h);
Vx = in.Vx(:); Vy = in.Vy(:);
V = sqrt(Vx.^2+Vy.^2);
Vvec = [ones(size(V)), zeros(size(V))];
ind = find(V>0);
if ~isempty(ind), Vvec(ind,:) = [Vx(ind)./V(ind), Vy(ind)./V(ind)]; end

% Set lower limits
in.h(in.h<in.ks) = in.ks+0.01;
in.Hrms(in.Hrms<0.01) = 0.01;

% Call bijker_calc function
[Qbc, Qsc] = bijker_calc(in.h(:),V(:),in.Hrms(:),in.T(:),...
    in.ks,in.ws,in.d50,in.d90,in.rhow,in.rhos,in.b);

% Transfer results to output data structure
out.Qs_bed = reshape(Qbc,M,N);
out.Qs_susp = reshape(Qsc,M,N);
out.Qs_tot = reshape(Qbc+Qsc,M,N);
out.Qsx = reshape((Qbc+Qsc).*Vvec(:,1),M,N);
out.Qsy = reshape((Qbc+Qsc).*Vvec(:,2),M,N);

varargout{1} = out;
varargout{2} = in;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Qbc, Qsc, Z_break, rel_rough] = bijker_calc(h,Vc,Hw,Tw,ks,ws,d50,d90,rhow,rhos,b)

NH = length(h);
NZ = 20;

kappa = 0.40;
g = 9.81;

Uw = zeros(NH,1);
Tw(Tw<0.1) = 0.1;
k0h = 4*pi^2./(g*Tw.^2) .* h;
kh =kh_calc(k0h);
Aw = 0.5*Hw./sinh(kh);
ind = Aw>0.001;
if ~isempty(ind)
    Uw(ind) = 2*pi*Aw(ind)./Tw(ind);
end

A = ks./h;

C = 18*log10(12*h/ks);
Cprm = 18*log10(12*h/d90);
mu = (C./Cprm).^1.5;

fc = 8*g./C.^2;
fw = zeros(NH,1);
fw(Aw>0.001) = exp(-6+5.2*(Aw/ks).^(-0.19));

taubc = 0.125 * rhow * fc .* Vc.^2;
taubw = zeros(NH,1);
taubw(fw>0) = 0.25 * rhow * fw .* Uw.^2;
taubcw = taubc + taubw;

taubcw(taubcw<0.001) = 0.001;

ustc = (taubc/rhow).^0.5;
ustcw = (taubcw/rhow).^0.5;

Z = ws./(kappa*ustcw);

Qbc = rhos * b * d50 * ustc .* exp( (-0.27*(rhos-rhow)*g*d50) ./ (mu.*taubcw));

I1 = zeros(NH,1);
I2 = zeros(NH,1);
i1 = zeros(1,NZ);
i2 = zeros(1,NZ);
for n = 1 : NH
    zprm = logspace(log10(A(n)),0,NZ);
    i1(:) = 0.216 * A(n).^((Z(n)-1)/(1-A(n)).^Z(n)) .* (((1-zprm)./zprm)).^Z(n);
    I1(n) = trapz(zprm,i1);
    i2(:) = i1 .* log(zprm);
    I2(n) = trapz(zprm,i2);
end

Qsc = 1.83 * Qbc .* (I2 + I1 .* log(33./A));

% Z_break = Z(50);
% rel_rough = ks./h(50);
% ratio = Qsc(60)./Qbc(50);
