% /////// fvtransect ///////
% C = fvtransect(C,resfil,geofil,pline,it)
% construct / update a transect object (plan view of velocity vectors across transect)
%
% inputs
%   C  = (structure), can be empty on first call and populated thereafter
%   resfil = .nc output file containing 2D or 3D TUFLOW-FV results
%   geofil = _geo.nc file containing information on the model mesh
%   pline = [npts 2] xy matrix defining the points making up the cross section
%   it = integer timestep
%
% outputs
%   C = stucture with information regarding your curtain object.
%       C.p is the handle to the patch object
%
% Jesper Nielsen November 2012

function C = fvtransect(C,resfil,geofil,pline,it)

% feeding in data
if ~isfield(C,'data_vector')
    C.data_vector = [];
end

if ~isfield(C,'data_scalar')
    C.data_scalar = [];
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

if ~isfield(C,'variables_vector'), C.variables_vector = {'V_x';'V_y'}; end
if ~isfield(C,'variables_scalar'), C.variables_scalar = {}; end % can be 'magnitude' (compatible when feeding data) or TUFLOW-FV variables
if ~isempty(C.variables_scalar)
    fancy = true;
else
    fancy = false;
    if ~isfield(C,'facecolor'),C.facecolor = 'k'; end
    if ~isfield(C,'edgecolor'),C.edgecolor = 'k'; end
end
if ~isfield(C,'spacing'),C.spacing = 10; end        % spacing (m) between arrows
if ~isfield(C,'shift'),C.shift = 0; end             % distance from start of transect to 1st arrow
if ~isfield(C,'spherical'),C.spherical = false; end  % true when modelling in spherical coordinates (pline must then also be in spherical coordinates)

if ~isfield(C,'v_multi'), C.v_multi = 100; end           % 1 m/s or Kg/m/s etc etc = 100 pixels on screen
if ~isfield(C,'v_shape_head'), C.v_shape_head = [0.15 0.30 1]; end  % [width head, legth head, logical (true if frac of arrow length, false if abs values)]
if ~isfield(C,'v_shape_tail'), C.v_shape_tail = [0.2 0]; end % [width of stem, logical (true if frac of arrow length, false if abs value)];
if ~isfield(C,'linewidth'),C.linewidth = 0.5; end

if ~isfield(C,'ref')
    C.ref = 'sigma';
    C.range = [0 1];
end
C.zcell = false; % specified for fv_get_dave etc. NOT YET IMPLEMENTED IN fvtransect
C.stamp = [];

% does transect object already exist
if ~isfield(C,'p')
    new  = true;
    set(C.f,'CurrentAxes',C.h)
    %     axis fill
    %     axis equal
    set(C.h,'NextPLot','add')
    set(C.f,'DoubleBuffer','on')
    set(C.f,'RendererMode','auto')
else
    new = false;
end

if new
    % segment pline by C.spacing (metres)
    x = pline(:,1);
    y = pline(:,2);
    if C.spherical
        [x,y,gz] = ll2utm(x,y);
        if length(unique(gz(:,1))) ~= 1
            error('results extend beyond a single gridzone')
        end
    end
    
    dx = diff(x);
    dy = diff(y);
    dx = vertcat(0,dx);
    dy = vertcat(0,dy);
    d = sqrt(dx.^2 + dy.^2);
    chain = cumsum(d);
    
    % -- arrows spaced at C.spacing along transect
    chain_plot = C.shift:C.spacing:chain(end);
    
    % -- coords of arrows tails
    xp = interp1(chain,x,chain_plot);
    yp = interp1(chain,y,chain_plot);
    
    if C.spherical
        [xp,yp] = utm2ll(xp,yp,gz(1,:));
    end
    
    % index into model results
    if C.spherical
        ic2 = fv_get_ids([xp' yp'],resfil,'cell',true);
    else
        ic2 = fv_get_ids([xp' yp'],resfil,'cell',true);
    end
    % -- remove points which fv_get_ids coluld not assign a cell
    if any(isnan(ic2))
        xp(isnan(ic2)) = [];
        yp(isnan(ic2)) = [];
        ic2(isnan(ic2)) = [];
    end
    
    C.xp = xp;
    C.yp = yp;
    C.ic2 = ic2;
    
    % vector variables
    if ~isempty(C.data_vector)
        nv_vector = size(C.data,3);
        variables_vector = {};
        if nv_vector > 2
            error('length of 3rd dimension for inputted data must <= 2')
        end
        for aa = 1:nv_vector
            C.variables_vector = cat(1,variables_vector,['data_' num2str(aa)]);
        end
    else
        C.variables_vector = fv_variables(C.variables_vector);
        nv_vector = length(C.variables_vector);
    end
end

% scaling
[fx fy] = getscale(C.h);
C.v_scale = sqrt(fx.^2 + fy.^2) * C.v_multi;


% retrieve and process vector variables into 2D results if needed
switch lower(C.ref)
    case {'top','bot'}
        C = fv_get_layer(C,resfil,it,'variable',C.variables_vector,'ref',C.ref,'range',C.range,'data',C.data_vector,'stat',C.stamp,'zcell',C.zcell);
    case {'sigma','elevation','depth','height'}
        C = fv_get_dave(C,resfil,it,'variable',C.variables_vector,'ref',C.ref,'range',C.range,'data',C.data_vector,'stat',C.stamp,'zcell',C.zcell);
end

% when you want coloured arrows ===========================================
if fancy
    % scalar variables
    if ~isempty(C.data_scalar)
        nv_scalar = size(C.data,3);
        variables_scalar = {};
        if nv_scalar > 2
            error('length of 3rd dimension for inputted data must <= 2')
        end
        for aa = 1:nv_scalar
            C.variables_scalar = cat(1,variables_scalar,['data_' num2str(aa)]);
        end
    else
        C.variables_scalar = fv_variables(C.variables_scalar);
        nv_scalar = length(C.variables_scalar);
    end
    
    % retrieve and process scalar variables in 2D results if needed
    switch lower(C.variables_scalar)
        case {'magnitude','mag'}
            C.mag = bsxfun(@hypot,C.(C.variables_vector{1}),C.(C.variables_vector{2}));
        otherwise
            switch lower(S.ref)
                case {'top','bot'}
                    C = fv_get_layer(C,resfil,it,'variable',C.variables_scalar,'ref',C.ref,'range',C.range,'data',C.data_scalar,'stat',C.stamp,'zcell',C.zcell);
                case {'sigma','elevation','depth','height'}
                    C = fv_get_dave(C,resfil,it,'variable',C.variables_scalar,'ref',C.ref,'range',C.range,'data',C.data_scalar,'stat',C.stamp,'zcell',C.zcell);
            end
            if nv_scalar == 2
                C.mag = bsxfun(@hypot,C.(C.variables_scalar{1}),C.(C.variables_scalar{2}));
            else
                C.mag = C.variables_scalar;
            end
    end
end
% =========================================================================

% verticees making patch objects for vectors;
[C.v_pvx C.v_pvy] = arrow2(C.(C.variables_vector{1})(C.ic2),C.(C.variables_vector{2})(C.ic2),C.xp',C.yp','scale',C.v_scale,'shape_head',C.v_shape_head,'shape_tail',C.v_shape_tail);

% create / update patch object
if new
    if fancy
        C.p = patch(C.v_pvx,C.v_pvy,C.mag(C.ic2)','EdgeColor',C.edgecolor,'parent',C.h,'linewidth',C.linewidth);
    else
        C.p = patch('XData',C.v_pvx,'YData',C.v_pvy,'FaceColor',C.facecolor,'EdgeColor',C.edgecolor,'parent',C.h,'linewidth',C.linewidth);
    end
else
    if fancy
        set(C.p,'XData',C.v_pvx,'YData',C.v_pvy,'CData',C.mag(C.ic2));
    else
        set(C.p,'XData',C.v_pvx,'YData',C.v_pvy);
    end
end