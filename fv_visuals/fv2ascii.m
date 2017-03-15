% /////// fv2ascii ///////
% fv2ascii(resObj,grd,varargin)
%
% Creates an ascii grid file (.asc) which can be read by GIS software
% packages such as MapInfo. The model results are not interpolated onto the
% grid points. Grid points are assigned the cell-centred output
% corresponding to the cell which ecompasses the grid point.
%
% When two fields are stored within resObj.ResultsCell for example V_x and V_y
% then the magnitude is calculated.
%
% inputs
%   resObj = fvres_sheet object created by fvres_sheet.m
%   grd    = grid spacing, metres | degrees depending on whether the simulation was performed in spherical
%
% optional inputs
%   'bounds' / [x_BotLeft y_BotLeft; x_TopRight y_TopRight] output grid limits ie. create a subset of the model results default: no subsetting
%   'outfil' / name of the outfil default: a name based on the simulation name and the grid spacing
%   'grd2cell' / grid to cell indexing: cell array created on a previous call to this function.
%       If you are creating a set of .asc files within a loop, on your first
%       call to fv2ascii you can create the indexing and then save yourself the
%       trouble for the following calls to fv2ascii. eg.
%           for aa = 1:n_sed
%            ......
%            if aa == 1
%               [grd2cell{1},grd2cell{2}] = fv2ascii(fvresObj,grd,'outfil',outfil);
%            else
%               fv2ascii(fvresObj,grd,'outfil',outfil,'grd2cell',grd2cell);
%            end
%           end
%
% optional outputs
%   'grd2cell' see above optional inputs.
%
% Jesper Nielsen, December 2013
% Overhauled in November 2014 to comply with the object based programming

function varargout = fv2ascii(resObj,grd,varargin)

% defaults
bounds = [];
outfil = []; % default outfil specified later
grd2cell = {};

% variable arguments
nva = length(varargin);
if mod(nva,2) > 0
    error('Expecting variable arguments as descriptor / value pairs')
end

for aa = 1 : 2 : nva
    varargtyp{aa} = varargin{aa};
    varargval{aa} = varargin{aa+1};
    switch lower(varargtyp{aa})
        case {'bound','bounds'}
            bounds = varargval{aa};
        case 'outfil'
            outfil = varargval{aa};
        case 'grd2cell'
            grd2cell = varargval{aa};
        otherwise
            error('unexpected variable argument type')
    end
end

% check bounds
if ~isempty(bounds) && size(bounds,1) ~= 2 && size(bounds,2) ~= 2
    error('input bounds must be of size [2,2]')
end

% results file
resfil = char(resObj.ResFil);

% are the TUFLOW-FV results in spherical coordinates
tmp = ncreadatt(resfil,'/','spherical');
switch tmp
    case 'true'
        spherical = true;
    case 'false'
        spherical = false;
end

% generate the grid
% -- extents
if isempty(bounds)
    TMP = netcdf_get_var(resfil,'names',{'node_X';'node_Y'});
    x1 = min(TMP.node_X);
    y1 = min(TMP.node_Y);
    x2 = max(TMP.node_X);
    y2 = max(TMP.node_Y);
else
    x1 = bounds(1,1);
    y1 = bounds(1,2);
    x2 = bounds(2,1);
    y2 = bounds(2,2);
end
xvec = x1:grd:x2+grd;
yvec = y1:grd:y2+grd;
yvec = fliplr(yvec);
nx = length(xvec);
ny = length(yvec);
[xgrd,ygrd] = meshgrid(xvec,yvec);

% index grid points into TUFLOW-FV 2D results
if isempty(grd2cell)
    tmp = fv_get_ids_2([xgrd(:) ygrd(:)],resfil,'cell',true);
    i_in = ~isnan(tmp);
    i_fv = tmp(i_in);
else
    i_in = grd2cell{1};
    i_fv = grd2cell{2};
end

% did you want to store the indexing
if nargout == 2
    varargout{1} = i_in;
    varargout{2} = i_fv;
end

% variable/s
vars = fieldnames(resObj.ResultsCell);
vars = setxor(vars,'stat');
nv = length(vars);
if nv == 2
    i = strfind(vars,'_x');
    if isempty(i)
        error('when two variables exist within the resObj they must be the x and y components')
    elseif isempty(i{1})
        j = 2;
    else
        j = 1;
    end
    v_name = strrep(vars{j},'_x','');
elseif nv == 1
    v_name = vars;
elseif nv > 2
    error('the resObj can contain a maximum of 2 variables')
end

% preallocate
res = -9999 * ones(nx*ny,1);

% assign to grid
if nv == 2
    res(i_in) = hypot(resObj.ResultsCell.(vars{1})(i_fv),resObj.ResultsCell.(vars{2})(i_fv));
else
    res(i_in) = resObj.ResultsCell.(vars{1})(i_fv);
end

% -- reshape back onto a grid
res_grd = reshape(res,ny,[]);

% -- set NaNs to -9999
res_grd(isnan(res_grd)) = -9999;

% output file
if isempty(outfil)
    if spherical
        ext = ['_' v_name '_' num2str(grd,'%.4f') 'deg'];
    else
        ext = ['_' v_name '_' num2str(grd,'%.0f') 'm'];
    end
    ext = strrep(ext,'.','pnt');
    outfil = strrep(resfil,'.nc',[ext '.asc']);
end

% -- write the headers
fid = fopen(outfil,'w');
fprintf(fid,'%s %d\n','ncols',nx);
fprintf(fid,'%s %d\n','nrows',ny);
fprintf(fid,'%s %f\n','xllcorner',x1);
fprintf(fid,'%s %f\n','yllcorner',y1);
fprintf(fid,'%s %f\n','cellsize',grd);
fprintf(fid,'%s %d\n','NODATA_value',-9999);

% -- write the body
fmat = ['%f' repmat(' %f',1,nx-1) '\n'];
inc = 0;
tic
for aa = 1:ny
    fprintf(fid,fmat,res_grd(aa,:));
    inc = mytimer(aa,[1 ny],inc);
end
fclose(fid);



