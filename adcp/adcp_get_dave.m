% /////// adcp_get_dave ///////
% Performs what fv_get_sheet_points does on TUFLOW-FV model
% results on data from bottom mounted ADCPs. The user needs to be aware
% that the reference when using the height and depth options for depth
% averaging do not respectively refer to the sea bed and the water surface
% but the respective limits of the bin faces for which data is
% returned by the ADCP.
%
% Processes 3D variables (V_x, V_y, W, TSS etc.) by depth averaging as specified.
% The specified vaiables are processed and assigned to the output stucture
% in a field named after the means of depth averaging.
% 2D variables recieve no processing and are passed into fields named 'twod'.
%
% inputs
%   C = input structure (usually loaded from preformatted .mat file)
%   variables = variable/s within C which are to be processed (if 3D) and passed to the output structure.
%
% optional inputs as descriptor / value pairs
%   'ref'      / {'ref1';'ref2';'ref3'....} one of 'sigma','elevation','height','depth' or any combination, default: sigma
%   'range'    / {range1;range2;range3....} corresponding to above ie, [s1 s2], [e1 e2], etc.               default: [0 1] (matching the above default)
% outputs
%   OUT = Stucture with fields for the depth averaged variables.
%
% Jesper Nielsen, April 2012, August 2014

function OUT = adcp_get_dave(C,variables,varargin)

% defaults
ref_all = {'sigma'};
range_all = {[0 1]};

nva = length(varargin);
if mod(nva,2) > 0
    error('Expecting variable arguments as descriptor/value pairs')
end

for aa = 1 : 2 : nva
    switch lower(varargin{aa})
        case 'ref'
            ref_all = lower(varargin{aa+1});
        case 'range'
            range_all = varargin{aa+1};
        otherwise
            error('unexpected variable argument type')
    end
end

% top and bot are not yet supported as ref inputs
if any(strcmpi(ref_all,'top')) || any(strcmpi(ref_all,'bot'))
    error('ref inputs of top or bot are not yet supported in adcp_get_dave')
end

% check means of depth averaging
fv_check_dave(ref_all,range_all)
no = length(ref_all);

% my 3D variables
nv = length(variables);
is_3d = false(nv,1);
if isfield(C,'zface')
    nl = size(C.zface,1) - 1;
    for aa = 1:nv
        v_name = variables{aa};
        if size(C.(v_name),1) == 1 % variable is 2D
            continue
        elseif size(C.(v_name),1) ~= nl
            error(['1st dimension of variable ' v_name ' inconsistent with variable zcell in data file'])
        else
            is_3d(aa) = true;
        end
    end
    
    variables_3D = variables(is_3d);
    nv3 = length(variables_3D);
    
    % check 2nd dimension, time dimension
    nt = length(C.TIME_hydro);  % 3D variables can only correspond to TIME_hydro
    for aa = 1:nv3
        v_name = variables_3D{aa};
        if size(C.(v_name),2) ~= nt
            error(['2nd dimension of variable ' v_name ' inconsistent with variable TIME_hydro in data file'])
        end
    end
else
    nv3 = 0;
end

% process 3D variables
if nv3 > 0;
    % -- layer faces
    lfz = C.zface;
    if size(lfz,2) ~= nt
        lfz = repmat(lfz(:,1),1,nt);
    end
    
    % -- depths defining limits to average between
    for aa = 1:no
        
        top = lfz(1:end-1,:);
        bot = lfz(2:end,:);
        
        ref = lower(ref_all{aa});
        range = range_all{aa};
        
        % index into the top layer face
        v_name = variables_3D{1};
        [row, col] = find(~isnan(C.(v_name)));
        [col, i, ~] = unique(col);
        row = row(i);
        
        % -- ensembles (time steps) where all bins were NaN
        ens = (1:nt)';
        i = ~ismember(ens,col);
        col_fil = ens(i);
        row_fil = ones(length(col_fil),1); % the indexing does not matter as these ensembles will be filled with NaNs
        col_all = [col; col_fil];
        row_all = [row; row_fil];
        [col_all,i] = sort(col_all);
        row_all = row_all(i); % currently indexing into top most bin for each ensemble (time step). This is equivalent to indexing into the top faces.
        i = sub2ind([nl nt],row_all,col_all);
        
        switch ref % d1 is below d2
            case 'sigma'
                depth = top(i)' - bot(end,1);
                d1 = bot(end,1) + range(1) * depth;
                d2 = bot(end,1) + range(2) * depth;
            case 'elevation'
                d1 = bot(end,:);
                d2 = top(i)';
                d1 = max(d1,range(1));
                d2 = min(d2,range(2));
                if any(d2 < d1)
                    error('data does not always exist between the depth averaging limits') % Could be changed so results are set to NaN when this happens (low tide etc.)
                end
            case 'height'
                d1 = bot(end,:) + range(1);
                d2 = bot(end,:) + range(2);
                d2 = min(d2,top(i)');
            case 'depth'
                d1 = top(i)' - range(2);
                d2 = top(i)' - range(1);
                d1 = max(d1,bot(end,:));
        end
        
        % -- engine room
        bot = bsxfun(@max,bot,d1);
        top = bsxfun(@min,top,d2);
        frc = bsxfun(@rdivide,(top-bot),(d2-d1));
        frc = max(frc,0);
        
        % -- process 3D results into 2D results
        res = zeros(nl,nt,nv3);
        for bb = 1:nv3
            v_name = variables_3D{bb};
            res(:,:,bb) = C.(v_name);
        end
        out = bsxfun(@times,res,frc);
        i = isnan(res);
        j = all(i,1);
        out(i) = 0; % empty bins set to 0. Except for small holes in data these should correspond to values beyond the depth averaging limits. Where the corresponding frc is 0
        out(:,j) = NaN; % empty ensembles are assigned NaN
        dav = sum(out,1);
        
        % -- store the processed results with a tag detailing how it was processed
        o_name = dave_ref_names(ref,range);
        for bb = 1:nv3
            v_name = variables_3D{bb};
            OUT.(v_name).(o_name) = dav(:,:,bb);
        end
    end
end

% add in your 2D variables
for aa = 1:nv
    v_name = variables{aa};
    if ~is_3d(aa)
        o_name = 'twod';
        OUT.(v_name).(o_name) = C.(v_name);
    end
end

% take coordinatess and time vectors with you
OUT.coordinates = C.coordinates;
time_variables = {'TIME_hydro';'TIME_sedi';'TIME_wave';'TIME_wind'};
for aa = 1:4
    t_name = time_variables{aa};
    if isfield(C,t_name)
        OUT.(t_name) = C.(t_name);
    end
end
