% [out in] = vanrijn(in)
%
% in = vanrijn([]) (to generate empty input data structure)
%
% Calculates sediment transport using method of van Rijn (TRANSPOR).
% Input and output are via data structures 'in' and 'out'.
% The input data can be entered via command line prompts by passing '[]'
% to the function, or alternatively an empty input data structure can be
% generated.
% The water depth field 'in.h' must be present, and can be a scalar,
% a vector or a matrix.  The other inputs can be omitted and will assume
% default values.
% Sediment transport rates returned in the 'out' data structure are in
% kg/s.

function varargout = vanrijn(in)

% Initialise output structure
out = struct('Qsx',[],'Qsy',[],'Qs_tot',[],'Qs_susp1',[],'Qs_susp2',[],'Qs_bed',[]);

% If input is empty ask for manual input
if isempty(in)
    in = struct('h',[],'Vx',[],'Vy',[],'Hs',[],'Tp',[],'Phi',[],...
        'd50',[],'d90',[],'ws',[],'rc',[],'rw',[],'ds',[],...
        'zref',[],'nu',[],'rhow',[],'rhos',[],'beta',[]);
    in.h = input(['Input depth vector, h\n',...
        '(or hit return for empty input structure), :']);
    if isempty(in.h), varargout{1} = in; return, end
    in.Vx = input('Input x-direction current, Vx: ');   
    in.Vy = input('Input x-direction current, Vy: ');
    in.Hs = input('Input significant wave height, Hs: ');
    in.Tp = input('Input wave peak period, Tp: ');
    in.Phi = input('Input wave direction relative to current: ');
    in.d50 = input('Input sediment size, d50: ');
    in.d90 = input('Input sediment size, d90: ');
    in.ws = input('Input sediment fall velocity, ws: ');
    in.rc = input('Input current related roughness, rc: ');
    in.rw = input('Input wave related roughness, rw: ');
    in.ds = input('Input mixing layer thickness, ds: ');  
    in.zref = input('Input reference level, zref: ');  
    in.nu = input('Input kinematic viscosity, nu: ');  
    in.rhow = input('Input fluid density, rhow: '); 
    in.rhos = input('Input sediment desnity, rhos: '); 
    in.beta = input('Input ratio sediment/fluid mixing, beta: '); 
end

% Check input and provide default values where necessary
if (~isfield(in,'h')||isempty(in.h)), disp('Must have depth input'), return, end
if (~isfield(in,'Vx')||isempty(in.Vx)), in.Vx = zeros(size(in.h)); end
if all(size(in.Vx)~=size(in.h)), disp('size(Vx) must equal size(h)'), return, end
if (~isfield(in,'Vy')||isempty(in.Vy)), in.Vy = zeros(size(in.h)); end
if all(size(in.Vy)~=size(in.h)), disp('size(Vy) must equal size(h)'), return, end
if (~isfield(in,'Hs')||isempty(in.Hs)), in.Hs = zeros(size(in.h)); end
if all(size(in.Hs)~=size(in.h)), disp('size(Hs) must equal size(h)'), return, end
if (~isfield(in,'Tp')||isempty(in.Tp)), in.Tp = zeros(size(in.h)); end
if all(size(in.Tp)~=size(in.h)), disp('size(Tp) must equal size(h)'), return, end
if (~isfield(in,'Phi')||isempty(in.Phi)), in.Phi = zeros(size(in.h)); end
if all(size(in.Phi)~=size(in.h)), disp('size(Phi) must equal size(h)'), return, end
if (~isfield(in,'d50')||isempty(in.d50)), in.d50 = 0.00025; end
if (~isfield(in,'d90')||isempty(in.d90)), in.d90 = 0.00050; end
if (~isfield(in,'ws')||isempty(in.ws)), in.ws = 0.025; end
if (~isfield(in,'rc')||isempty(in.rc)), in.rc = 0.05; end
if (~isfield(in,'rw')||isempty(in.rw)), in.rw = in.rc; end
if (~isfield(in,'ds')||isempty(in.ds)), in.ds = in.rc; end
if (~isfield(in,'zref')||isempty(in.zref)), in.zref = in.rc/6; end
if (~isfield(in,'nu')||isempty(in.nu)), in.nu = 1e-6; end
if (~isfield(in,'rhow')||isempty(in.rhow)), in.rhow = 1025; end
if (~isfield(in,'rhos')||isempty(in.rhos)), in.rhos = 2650; end
if (~isfield(in,'beta')||isempty(in.beta)), in.beta = 1; end

names = fieldnames(in);
% Convert fields to double precision if necessary
for i = 1 : length(names)
    if ~isa(in.(names{i}),'double')
        in.(names{i}) = double(in.(names{i}));
    end
end  

% Calculate velocity magnitude and unit vector
[M,N] = size(in.h);
Vx = in.Vx(:); Vy = in.Vy(:);
V = sqrt(Vx.^2+Vy.^2);
Vvec = [ones(size(V)), zeros(size(V))];
ind = V>0.0001;
V(~ind) = 0.0001;
Vvec(ind,1) = Vx(ind)./V(ind);
Vvec(ind,2) = Vy(ind)./V(ind);

% Call transpor mex function
par = [in.d50,in.d90,in.ws,in.rc,in.rw,in.ds,in.zref,in.nu,...
    in.rhow,in.rhos,in.beta];
[SSR, SSF, SBF] = transpor(in.h(:),V(:),in.Hs(:),in.Tp(:),in.Phi(:),par);

% Transfer results to output data structure
out.Qs_susp1 = reshape(SSR,M,N);
out.Qs_susp2 = reshape(SSF,M,N);
out.Qs_bed = reshape(SBF,M,N);
out.Qs_tot = reshape(SSR+SBF,M,N);
out.Qsx = reshape((SSR+SBF).*Vvec(:,1),M,N);
out.Qsy = reshape((SSR+SBF).*Vvec(:,2),M,N);

varargout{1} = out;
varargout{2} = in;