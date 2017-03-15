% /////// fvcurtain ///////
% C = fvcurtain(C,resfil,geofil,pline,it,variables)
% construct / update a curtain object (profile view of x-section through 3D or 2D TUFLOW-FV model)
% 2 options for viewing curtain
%   opt1, [default]          : elevation against chainage (x-z plane)
%   opt2, C.chainage = false    : elevation against coordinates (x-y-z) plane
%
% inputs
%   C  = (structure), can be empty on first call and populated thereafter
%   resfil = .nc output file containing 3D results
%   geofil = _geo.nc file containing information on the model mesh
%   pline = [npts 2] xy matrix defining the points making up the cross section
%   it = integer timestep
%   variables = scalar variable or pair of variables with vector components ie {'TSS'} or {'V_x'; 'V_y'}
%
% outputs
%   C = stucture with information regarding your curtain object.
%       C.p is the handle to the patch object
%
% Jesper Nielsen October 2012


function C = fvcurtain(C,resfil,geofil,pline,it,variable)

% feeding in data
if isfield(C,'data') && ~isempty(C.data)
    if ~isempty(variable)
        error('input "variable" must be empty when inputting data')
    end
else
    C.data = [];
end

% Defaults
if ~isfield(C,'f')
    if isfield(C,'h')
        C.f = get(C.h,'parent');
    else
        C.f = figure;
    end
end
if ~isfield(C,'h'), C.h = axes('parent',C.f); end
if ~isfield(C,'shading'),C.shading = 'flat'; end
if ~isfield(C,'edgecolor'),C.edgecolor = 'none'; end
if ~isfield(C,'chainage'),C.chainage = true; end

% does curtain object already exist
if ~isfield(C,'p')
    new  = true;
    set(C.h,'NextPLot','add')
    set(C.f,'DoubleBuffer','on')
    %     set(C.f,'Renderer','painters')
    if C.chainage
        view(C.h,0,0)
    else
        view(3);
    end
else
    new = false;
end

% recognised TUFLOW-FV variables
if ~isempty(C.data)
    nv = size(C.data,3);
    variable = {};
    if nv > 2
        error('length of 3rd dimension for inputted data must <= 2')
    end
    for aa = 1:nv
        variable = cat(1,variable,['data_' num2str(aa)]);
    end
else
    variable = fv_variables(variable);
    nv = length(variable);
end

if new
    if isfield(C,'fv_get_curtain')
        C = rmfield(C,{'fv_get_curtain';'face_centre'}); % face_centre is only variable created using indexing
    end
end

% call fv_get_curtain
C = fv_get_curtain(C,resfil,geofil,pline,it,'variables',variable,'smooth',false,'chainage',C.chainage,'data',C.data);

% build patch object
if new
    C.p = patch('Faces',C.face,'Vertices',C.vert,'Parent',C.h,'EdgeColor',C.edgecolor,'FaceColor',C.shading);
    
    % -- don't want your curtain object to appear in legend
    hAnnotation = get(C.p,'Annotation');
    hLegendEntry = get(hAnnotation','LegendInformation');
    set(hLegendEntry,'IconDisplayStyle','off')
end

% update layerface_Z
set(C.p,'Vertices',C.vert);

% update scalar
if nv == 3
    mag = sqrt(C.(variable{1}).^2 + C.(variable{2}).^2 + C.(variable{3}).^2);
elseif nv == 2
    mag = sqrt(C.(variable{1}).^2 + C.(variable{2}).^2);
elseif nv == 1
    mag = C.(variable{1});
else
    error('a maximum of 3 variables can be specified')
end

set(C.p,'FaceVertexCData',mag)

% calculate fluxes
if nv == 2 || nv == 3
    C.flux = dot([C.V_x C.V_y],C.unorm,2) .* C.area;
end