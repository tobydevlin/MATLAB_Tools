% Calculate cyclone wind field using Holland model
%
% [p, vm_x, vm_y, vg] = holland(x,y,x0,y0,...
%                                        p0,pn,rmax,B,...
%                                        vfm,deltafm,thetamax,...
%                                        rhoa,lat)
%
% The inputs:
% x and y are the points (equal size scalars, vectors or matrices)for which
% calculations are to be performed.
% x0 and y0 are the coordinates of the cyclone eye.
% p0 is the cyclone central pressure.
% pn is the ambient pressure.
% rmax is the radius to maximum winds.
% B is the "peakedness" coefficient.
% vfm is a vector [vfm_x, vfm_y] of the cyclone forward movement velocity.
% thetamax is the line of maximum wind angle.
% rhoa is the density of air.
% lat is the latitude (-ve for Southern Hemisphere).
%
% The outputs:
% p is a pressure field.
% vm_x and vm_y are the x- and y-components of the 10m wind field.
% vg is the gradient wind speed field.
%
% Ian Teakle WBM Pty Ltd 17/05/06

function [p, vm_x, vm_y, vg] = holland(x,y,x0,y0,...
                                        p0,pn,rmax,B,...
                                        vfm,deltafm,thetamax,...
                                        rhoa,lat)
                                    
% coriolis parameter, rad/s
f = 2 * 2*pi/(24*3600) * sin(lat/180*pi);

% wind asymmetry angle, rad
thetamax = thetamax/180*pi;

% forward motion speed and direction
vfm_mag = sqrt(vfm(1)^2 + vfm(2)^2);
vfm_dir = atan2(vfm(2),vfm(1));

% radius (m)
r = sqrt((x-x0).^2+(y-y0).^2);
% Set min radius to 10 m to avoid numerical issues
r = r .* (r>10) + 10 * (r<=10);
% angle from centre (rad)
thetar = atan2(y-y0,x-x0);

% pressure field (hPa)
p = p0 + (pn-p0) .* exp(-(rmax./r).^B); 

% gradient wind speed field (m/s)
vg = sqrt((pn-p0) * 100 / rhoa * B .* (rmax./r).^B ...
    .* exp(-(rmax./r).^B) + r.^2*f^2/4) - sqrt(r.^2*f^2/4);

vgmax = sqrt((pn-p0) * 100 / (exp(1) * rhoa) * B + rmax.^2*f^2/4) ...
    - sqrt(rmax.^2*f^2/4);

% wind boundary layer coefficient, vm(10) = km * vg
km = 0.81 * (vg<6) ...
    + (0.81 - 2.96e-3 * (vg - 6)) .* (vg>=6 & vg<19.5) ...
    + (0.77 - 4.31e-3 * (vg - 19.5)) .* (vg>=19.5 & vg<45) ...
    + 0.66 * (vg>=45);

% km = 0.67;

deltavg = 1.0 * (r <= rmax) ...
    + (vg ./ vgmax) .* (r > rmax);

% 10 m wind speed (m/s)
vm_mag = km .* vg + ...
    deltafm .* deltavg .* vfm_mag .* cos(thetar - (vfm_dir + thetamax));

% Inflow angle (deg)
beta = (10 * (r / rmax)) .* (r<rmax) ...
    + (10 + 75 * (r / rmax - 1)) .* (r>=rmax & r<(1.2*rmax)) ...
    + (25) .* (r>=(1.2*rmax));
beta = beta / 180 * pi; % (rad)

% 10 m wind direction (rad)
if lat < 0 % Southern Hemisphere
    vm_dir = 3*pi/2 + thetar - beta;
else % Northern Hemisphere
    vm_dir = pi/2 + thetar + beta;
end

% 10 m wind vector (m/s)
vm_x =  vm_mag .* cos(vm_dir);
vm_y = vm_mag .* sin(vm_dir);


    



