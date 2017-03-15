% /////// fv_dave_profile ///////
%
% function OUT = fv_dave_profile(resfil,varargin)
%
% fv_dave_profile is designed to extract TUFLOW-FV results from PROFILES
% output.
% 3D results are processed into 2D results.
% 3D variables can be processed in different ways with a single call to
% this function.
% 2D variables undergo no processing.
%
% options for depth averaging:
%   sigma:     [s1 s2]     - average from s1*depth above the bed up to s2*depth above the bed
%   elevation: [e1 e2]     - average from e1 metres up to e2 metres (refereced to model datum)
%   height:    [h1 h2]     - average from h1 metres above the bed up to h2 metres above the bed
%   depth:     [d1 d2]     - average from d1 metres below the surface down to d2 metres below the surface
%   top:       [t1 t2]     - average from t1 layer down to t2 layer. 1 = top layer
%   bot:       [b1 b2]     - average from b1 layer up to b2 layer. 1 = bottom layer
%
% inputs
%   resfil = .nc file containing all OUTputs from FV run (2D or 3D)
%
% optional inputs as descriptor / value pairs
%   'ref'      / {'ref1';'ref2';'ref3'....}, one or a combination of 'sigma','elevation','height','depth'         default: 'sigma'
%   'range'    / {range1;range2;range3....}, corresponding to above ie [s1 s2], [e1 e2], h, d or any combination  default: [0 1] (corresponding to above default)
%   'variable' / {'var1';'var2';'var3'....}, variable/s to extract                                                default: all variables with 2nd dimension of time excluding "stat" & "layerface_Z"
%   'locations'/ {'loc1';'loc2';'loc3'....}, location/s in file to extract                                        default: all locations within file
%
% outputs (as example)
%   OUT.SITE_01.depth_5.V_x = variable V_x depth averaged from the
%       surface to a depth of 5 metres at the location 'SITE_01' corresponding to the
%       time vector OUT.time
%
% NOTE:
% fv_get_sheet_points is designed to depth average the 3D results from an idividual cell for all timesteps in one hit.
% fv_get_sheet is designed to simultaneously depth average 3D results for all cells in the model mesh timestep by timestep.
% fv_dave_profile is designed to depth average the 3D results from cell information stored in a TUFLOWFV 'Profile' output, or extracted profile output from fv_create_profiles
%
% Toby Devlin, Jesper Nielsen, Copyright (C) BMTWBM 2014% BMTWBM
% 


function OUT = fv_dave_profile(modfil,varargin)

% defaults
ref_all = {'sigma'};                                   % default to depth average entire water column
range_all = {[0 1]};
setv = false;
locations = {};

% variables arguments
if mod(length(varargin),2) > 0
    error('Expecting variable arguments as descriptor/value pairs')
end

for aa = 1 : 2 : length(varargin)
    varargtyp = varargin{aa};
    varargval = varargin{aa+1};
    switch lower(varargtyp)
        case 'ref'
            ref_all = varargval;
        case 'range'
            range_all = varargval;
        case {'variable','variables','var','vars'}
            variables = varargval;
            setv = true;
        case {'location';'locations';'loc';'locs'}
            locations = varargval;
        otherwise
            error('unexpected variable argument type')
    end
end


% open netcdf
nci = netcdf.open(modfil);
grid = netcdf.inqGrps(nci);

% Check the locations
locsin = cell(size(grid));
for aa = 1 : length(grid)
    locsin{aa} = netcdf.inqGrpName(grid(aa));
end
if isempty(locations)
    locations = locsin;
end

if ~setv
    variables = {};
    nv = netcdf.inqVarIDs(grid(1));
    % loop through variables
    for aa = nv
        [v_name, ~, dimids, ~] = netcdf.inqVar(grid(1),aa);
        name = netcdf.inqDim(grid(1),max(dimids));
        if strcmp(name,'Time')
            if ~strcmp(v_name,'layerface_Z')
                variables = cat(1,variables,v_name);
            end
        end
    end
end

% check for standard TUFLOW-FV variables
variables = fv_variables(variables);
variables = [variables ; 'layerface_Z'];
% check means of depth averaging
fv_check_dave(ref_all,range_all)

% Check Relevant variables
nv = length(variables);
is3D = false(nv,1);
for aa = 1:nv
    v_name = variables{aa};
    TMP = ncinfo(modfil,['/',locsin{1}, '/', v_name]);
    dim = TMP.Dimensions(1).Name;
    len = TMP.Dimensions(1).Length;
    switch dim
        case 'NumLayers'
            if len ~= 1
                is3D(aa) = true;
            end
    end
end
variables_3D = variables(is3D);
nv3 = length(variables_3D);

% Do all locations
for aa = 1 : length(locations)
    c_name = locations{aa};
	i = ismember(locsin,c_name);
   if isempty(i)
        disp(['WARNING: Location ' c_name ' was not found in the file']);
    else
        TMP = netcdf_get_var(grid(i),'names',variables);
        if nv3>0
            % -- faces of 3D cells
            lfz = TMP.layerface_Z;
            [nlf,nt] = size(lfz);
            nl = nlf-1;
            
            for bb = 1:length(range_all)
                range = range_all{bb};
                ref = ref_all{bb};
                nores = false;
                if ismember(ref,{'top';'bot'})
                    if range(2) > nl
                        range(2) = nl;
                    end
                    if range(1) > nl;
                        nores = true;
                        range(1) = 1; % dummy assignment
                    end
                end
                
                % -- depths defining limits to average between d1 is below d2
                top = lfz(1:end-1,:);
                bot = lfz(2:end,:);
                
                switch ref % d1 is below d2
                    case 'sigma'
                        depth = top(1,:) - bot(end,:);
                        d1 = bot(end,:) + range(1) * depth;
                        d2 = bot(end,:) + range(2) * depth;
                    case 'elevation'
                        d1 = max(bot(end,:),range(1));
                        d2 = min(top(1,:),range(2));
                    case 'height'
                        d1 = bot(end,:) + range(1);
                        d2 = min(bot(end,:) + range(2),top(1,:));
                    case 'depth'
                        d1 = max(top(1,:)-range(2),bot(end,:));
                        d2 = top(1,:) - range(1);
                    case 'top'
                        d1 = lfz(range(2)+1,:);
                        d2 = lfz(range(1),:);
                    case 'bot'
                        d1 = lfz(nl-range(1) + 2,:);
                        d2 = lfz(nl-range(2) + 1,:);
                end
                
                if ismember(ref,{'elevation';'height';'depth'})
                    if any(d1 > top(1,:)) || any(d2 < bot(end,:))
                        nores = true;
                    end
                end
                
                if ~nores
                    % -- Engine
                    bot = bsxfun(@max,bot,d1);
                    top = bsxfun(@min,top,d2);
                    frc = bsxfun(@rdivide,(top-bot),(d2-d1));
                    frc = max(frc,0);
                    
                    % -- process 3D results into 2D results
                    res = zeros(nl,nt,nv3);
                    for cc = 1:nv3
                        v_name = variables_3D{cc};
                        res(:,:,cc) = TMP.(v_name);
                    end
                    res = bsxfun(@times,res,frc);
                    res_2D = sum(res,1);
                else
                    display(['no 3D results exist in ' c_name ' for ref of ' ref ', range of ' num2str(range)])
                    res_2D = NaN(1,nt,nv3);
                end
                
                % -- store the processed results with a tag detailing how it was processed
                o_name = dave_ref_names(ref,range);
                for cc = 1:nv3
                    v_name = variables_3D{cc};
                    OUT.(c_name).(v_name).(o_name) = res_2D(:,:,cc);
                end
            end
        end
    end
    
    % add in your 2D variables
    for bb = 1:nv
        v_name = variables{bb};
        if ~is3D(bb)
            o_name = 'twod';
            OUT.(c_name).(v_name).(o_name) = TMP.(v_name);
        end
    end
    
end

% keep in time
OUT.time = ncread(modfil,'ResTime');

% clean up
netcdf.close(nci);

end