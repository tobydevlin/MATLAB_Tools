% /////// fv_get_curtain ///////
% function C = fv_get_curtain(C,resfil,geofil,pline,it,varargin)
% Extracts 3D TUFLOW-FV results at the cells which are intersected by pline
% Designed as the engine for fvcurtain
%
% inputs
%   C = (structure), empty on first call but on further calls contains information which greatly boosts performance
%   resfil = TUFLOW-FV .nc results file - from 3D simulation
%   geofil = TUFLOW-FV _geo.nc file containing information on model mesh
%   pline = [npts,2], matrix or points (x,y) which define the x-section
%   it = integer timestep
%
% optional inputs
%   'variable' / {'var1';'var2';'var3'....}, variable/s to extract    default: all variables with 1st dimension #3D cells and 2nd dimension time
%   'smooth' / logical,                      true for interpolated shading and smooth bed (NOT YET SUPPORTED)
%   'data' / [nc3,1,nd],                     matrix of 3D data, must correspond to resfil & it, ie. data = SED_1 + SED_2 (built from 3D variables)
%   'chainage' / logical,                    true to plot z against chainage, points defining verticees in patch object (chain,0,z), false to plot z against coordinates (x,y,z)
% outputs
%   C.fv_get_curtain = information stored on first call
%   C.coords = coords of intersections between pline and 2D cell edges
%   C.chain = chainage along pline at the above coords
%   C.vert = verticees (patch object property) for curatin produced by fvcurtain
%   C.face = indexing of faces to vertices (patch object property) for curatin produced by fvcurtain
%
% Jesper Nielsen, October 2012

function C = fv_get_curtain(C,resfil,geofil,pline,it,varargin)

% do hard work once
if ~isfield(C,'fv_get_curtain')
    
    % avoid repeating the slow netcdf.open routine
    if ischar(resfil)
        nci = netcdf.open(resfil,'NC_NOWRITE');
    elseif isnumeric(resfil)
        nci = resfil;
    else
        error('unexpected resfil (should be nci or filename)')
    end
    
    % defaults
    smooth = false;  % NOT YET SUPPORTED
    feed = false;
    chainage = true;
    spherical = false;
    [variables,~] = netcdf_variables_unlimited(nci);
    variables = setxor(variables,{'ResTime';'layerface_Z';'stat'});
    
    % variable arguments
    nva = length(varargin);
    if mod(nva,2) > 0
        error('Expecting variable arguments as descriptor/value pairs')
    end
    
    for i = 1 : 2 : nva
        varargtyp{i} = varargin{i};
        varargval{i} = varargin{i+1};
        switch lower(varargtyp{i})
            case {'variable','variables','var','vars'}
                variables = varargval{i};
            case 'smooth'
                smooth = varargval{i};
            case 'data'
                feed = true;
                data = varargval{i};
            case 'chainage'
                chainage = varargval{i};
            otherwise
                error('unexpected variable argument type')
        end
    end
    
    % basic checks
    if ~iscell(variables)
        error('expecting cell array for optional input variables')
    end
    
    if ~islogical(smooth)
        error('expecting logical input for optional input smooth')
    end
    
    if ~islogical(chainage)
        error('expecting logical input for optional input chainage')
    end
    
    if isempty(data) % fvcurtain always passes data thorugh to fv_get_curtain
        feed = false;
    end
    
    if ~isscalar(it)
        error('expecting scalar input for input it')
    end
    
    % check for standard TUFLOW-FV variables
    if feed
        variables = {};
        nv = size(data,3);
        for aa = 1:nv
            variables = cat(1,variables,['data_' num2str(aa)]);
        end
    else
        variables = fv_variables(variables);
        nv = length(variables);
    end
    
    % timesteps
    TMP = netcdf_get_var(nci,'names',{'ResTime'});
    t = TMP.ResTime;
    nt = length(t);
    
    % call fv_get_curtain_ids
    [ic2, coords, chain] = fv_get_curtain_ids(pline,geofil);
    nic2 = length(ic2);
    
    % was modelling performed in spherical coordinates
    tmp = ncreadatt(geofil,'/','spherical');
    switch tmp
        case 'true'
            spherical = true;
        case 'false'
            spherical = false;
    end
    % switch coordinates to meters
    if spherical
        if chainage
            %coords = coords * 180 / pi; % revert to degrees from radians NO LONGER NEEDED SINCE IAN UPDATED THE _GEO.NC file
            [x,y,gz] = ll2utm(coords(:,1),coords(:,2));
            if sum(diff(gz,1,1)) ~= [0 0];
                error('cannot convert to cartesian coordinates as points extend beyond a single grid zone')
            end
            dx = cat(1,0,diff(x));
            dy = cat(1,0,diff(y));
            chain = cumsum(sqrt(dx.^2 + dy.^2));
        end
    end
    
    % info info info
    names = {'layerface_Z';'NL';'idx2';'idx3'};
    TMP = netcdf_get_var(nci,'names',names,'timestep',1);
    layerface_Z = TMP.layerface_Z;
    idx2 = TMP.idx2;
    idx3 = TMP.idx3;
    nl = TMP.NL;
    
    nlf = nl + 1;
    nc2 = length(idx3);
    nc3 = length(idx2);
    
    % check for 3D TUFLOW-FV variables
    if feed
        if size(data,1) ~= nc3
            error('data is inconsistent with TUFLOW-FV results file or is 2D')
        end
    else
        i = true(nv,1);
        TMP = netcdf_get_var(nci,'names',variables,'timestep',1);
        for aa = 1:nv
            v_name = variables{aa};
%             if length(TMP.(v_name)) ~= nc3
%                 i(aa) = false;
%                 display(['variable ' v_name 'is not a 3D variable and will not be processed by fv_get_curtain'])
%             end
        end
        variables = variables(i);
        nv = length(variables);
    end
    
    % unit normal vectors
    norm = [-diff(coords(:,2)) diff(coords(:,1))];
    unorm_tmp(:,1) = norm(:,1) ./ hypot(norm(:,1),norm(:,2));
    unorm_tmp(:,2) = norm(:,2) ./ hypot(norm(:,1),norm(:,2));
    unorm = [];
    for aa = 1:nic2
        i = ic2(aa);
        unorm = cat(1,unorm,repmat(unorm_tmp(aa,:),nl(i),1));
    end
  
    % unit tangent vectors
    tang = [diff(coords(:,1)) diff(coords(:,2))];
    utang_tmp(:,1) = tang(:,1) ./ hypot(tang(:,1),tang(:,2));
    utang_tmp(:,2) = tang(:,2) ./ hypot(tang(:,1),tang(:,2));
    utang = [];
    for aa = 1:nic2
        i = ic2(aa);
        utang = cat(1,utang,repmat(utang_tmp(aa,:),nl(i),1));
    end
    
    % logical indexing into variables (faces values)
    %     ir = false(nc3,1);
    ir = [];
    for aa = 1:nic2
        i = ic2(aa);
        k = idx3(i);
        kk = k + nl(i) - 1;
        ir = cat(1,ir,(k:kk)');
        %         ir(k:kk) = true;
    end
     
    % indexing into layerfaces_Z (vert values)
    itop = double(idx3) + (0:nc2-1)';
    il = [];
    for aa = 1:nic2
        i = ic2(aa);
        tmp = itop(i):(itop(i) + nl(i));
        il = cat(1,il,tmp',tmp');
    end
    
    % patches
    % -- verticees
    nvert = length(il);
    vert = zeros(nvert,3);
    kk = 0;
    for aa = 1:nic2
        i = ic2(aa);
        % -- -- left side of rectangular patches
        j = kk + 1;
        jj = j + nlf(i) - 1;
        % -- -- right side of rectangular patches
        k = jj + 1;
        kk = k + nlf(i) - 1;
        if chainage
            vert(j:jj,1) = chain(aa);
            vert(k:kk,1) = chain(aa+1);
        else
            vert(j:jj,1:2) = repmat(coords(aa,:),nlf(i),1);
            vert(k:kk,1:2) = repmat(coords(aa+1,:),nlf(i),1);
        end
    end
    vert(:,3) = layerface_Z(il);
    
    % -- face to vert indexing
    nface = length(ir);
    face = zeros(nface,4);
    k = 1;
    for aa = 1:nic2
        i = ic2(aa);
        for bb = 1:nl(i)
            if bb == 1
                if aa == 1
                    face(k,1) = 1;
                else
                    face(k,1) = face(k-1,3) + 1;
                end
            else
                face(k,1) = face(k-1,4);
            end
            face(k,2) = face(k,1) + nlf(i);
            face(k,3) = face(k,2) + 1;
            face(k,4) = face(k,1) + 1;
            k = k+1;
        end
    end  
     
    % if you like it smooth
    % -- update verts
    
    % -- update face to vert indexing
    
    % preallocate
    mod3 = zeros(nc3,1,nv);
    lay = zeros(length(layerface_Z),1);
    
    % ready yourself for netcdf.getVar
    if ~feed
        varid = zeros(nv,1);
        dim1len = zeros(nv,1);
        for aa = 1:nv
            varid(aa) = netcdf.inqVarID(nci,variables{aa});
            [~, ~, dimids, ~] = netcdf.inqVar(nci,varid(aa));
            [~, dim1len(aa)] = netcdf.inqDim(nci,dimids(1));
        end
    end
    % -- z_Layerface
    varid_lf = netcdf.inqVarID(nci,'layerface_Z');
    [~, ~, dimids, ~] = netcdf.inqVar(nci,varid_lf);
    [~, dim1len_lf] = netcdf.inqDim(nci,dimids(1));
    
    % clear variables you don't need
    C.coords = coords;
    C.chain = chain;
    C.unorm = unorm;
    C.utang = utang;
    C.vert = vert;
    C.face = face;
    
    clear('TMP','aa','bb','chain','vert','face','i','ic2','idx2','idx3','itop','j','jj','k','kk','layerface_Z','names','nc2','nc3','nface','nic2','nl','nlf','nva','nvert','tmp','v_name','varargtype','varargval')
    
    % store variables you do need (don't redo the hard work)
    f_name = 'fv_get_curtain';
    C.(f_name) = v2struct;
    C.(f_name).avoidOverWrite = ''; % don't overwrite 'it';
    
else
    f_name = 'fv_get_curtain';
    checkme(C.(f_name),f_name,resfil,varargin) % function which checks whether you are using fv_get_curtain correctly
    v2struct(C.(f_name));
    % get your data
    if feed
        [~,loc] = ismember('data',varargin(1:2:end));
        data = varargin{2*loc};
    end
end

% extract at specified timestep
if it > nt
    it = nt;
elseif it < 1;
    it = 1;
end


% -- face values
if feed
    mod3 = data;
else
    for aa = 1:nv
        mod3(:,1,aa) = netcdf.getVar(nci,varid(aa),[0 it-1],[dim1len(aa) 1]);
    end
end

% -- z layerfaces
lay(:) = netcdf.getVar(nci,varid_lf,[0 it-1],[dim1len_lf 1]);
C.vert(:,3) = lay(il);

% -- face areas
dx = C.vert(C.face(:,2),1) - C.vert(C.face(:,1),1);
dy = C.vert(C.face(:,2),2) - C.vert(C.face(:,1),2);
dz = C.vert(C.face(:,2),3) - C.vert(C.face(:,3),3);
C.area = hypot(dx,dy) .* dz;

% -- face centres
C.face_centre(:,1) = (C.vert(C.face(:,2),1) + C.vert(C.face(:,1),1)) / 2;
C.face_centre(:,2) = (C.vert(C.face(:,2),2) + C.vert(C.face(:,1),2)) / 2;
C.face_centre(:,3) = (C.vert(C.face(:,2),3) + C.vert(C.face(:,3),3)) / 2;

% store away variables
for aa = 1:nv
    v_name = variables{aa};
    C.(v_name) = mod3(ir,:,aa);
end





