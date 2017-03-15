% /////// circles.m ///////
% [faces,verts] = circles(x,y,x_orig,y_orig,varargin)
%
% Outputs the faces and verticees used for creating 2D circles from patch objects
% The circles are placed at x_orig and y_orig
%
% inputs
%   x      = x vector components, can be empty if optional input fix is defined
%   y      = y vector components, can be empty if optional input fix is defined
%   x_orig = x coordinates of circles' origins
%   y_orig = y coordinates of circles' origins
%
% optional inputs
%   'scale' / scale      multiplier for circle size
%   'fix'   / magnitude  set the magnitude and hence fix the size of all circles
%
% outputs
%   faces = refer to "Faces" property of patch objects
%   verts = refer to "Vertices" property of patch objects
%   nseg  = number of segments used to define the circle, the
%           higher the number the smoother the circle. Default: 16
%
% JN November 2014

function [faces,verts] = circles(x,y,ox,oy,varargin)

% defaults
scale = 1000;
fix = [];
nseg = 16;

% optional inputs
noi = length(varargin);
if mod(noi,2) ~= 0
    error('expecting optional inputs as property / value pairs')
end
for aa = 1:2:noi
    switch lower(varargin{aa})
        case 'scale'
            scale = varargin{aa+1};
        case 'fix'
            fix = varargin{aa+1};
        case 'nseg'
            nseg = varargin{aa+1};
        otherwise
            error('unexpected variable argument type')
    end
end

% faces
nc = length(ox);
tmp = 1:nseg:nc*nseg;
faces = bsxfun(@plus,tmp',0:nseg-1);

% verticees
% -- convert to axes units
if isempty(fix)
    mag = hypot(x,y) * scale;
else
    mag = fix * scale;
end

% define the circles around the origin
tmp = 360/nseg:360/nseg:360;
verts = zeros(nc*nseg,2);
phi = repmat(tmp',1,nc);
if size(mag,1) == 1
    x = bsxfun(@times,cosd(phi),mag);
    y = bsxfun(@times,sind(phi),mag);
else
    x = bsxfun(@times,cosd(phi),mag');
    y = bsxfun(@times,sind(phi),mag');
end

% translate the circles
if size(ox,1) == 1
    x = bsxfun(@plus,x,ox);
    y = bsxfun(@plus,y,oy);
else
    x = bsxfun(@plus,x,ox');
    y = bsxfun(@plus,y,oy');
end

verts(:,1) = x(:);
verts(:,2) = y(:);