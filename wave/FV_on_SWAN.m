% /////// FV_on_SWAN ///////
% Function used to write current and water level boundary condition files for SWAN.
% Current and water levels are extracted from TUFLOW-FV netcdf output files.
% No interpolation occurs, grid points are assigned the value from the 2D cell which they occur within.
% Depth averaging when required is performed by fv_get_dave.
% Writes the text insert required to reference the boundary condition files in your .swn SWAN input file.
% Ouput files are written to the folder holding TUFLOW-FV output .nc file.
%
% As SWAN has no option to flag / identify data to ignore in its input grids other than the bed
% input, currents are assigned a value of zero when no TUFLOW-FV model
% exists and water levels are assigned a value of -50 to ensure SWAN turns
% the elements dry when the TUFLOW-FV cell has gone dry
%
% inputs
%   C = structure with fields xp,yp,alp,mx,my,dx & dy =  as in SWAN manual
%   ncfil                                             =  TUFLOW-FV netcdf output file
%
% optional inputs as descriptor / value pairs
%   tlim  = [t1 t2] (matlab convention) start and end of bc file     default: starts at beginning of TUFLOW-FV output, end with the output
%   dt = timestep in bc file (hours)                                 default: TUFLOW-FV output timestep
%   ref = option for depth averaging, see fv_get_dave.m              default: 'sigma'
%   range = corresponds to above                                     default: [0 1]
% calls on fv_get_ids, fv_get_dave
%
% Jesper Nielsen July 2012, minor touch ups to use fvres_sheet

function FV_on_SWAN(C,ncfil,varargin)

% values for dry cells
h_dry = -10; % SWAN4091 won't allow outputs WATLEV outputs below -15
v_dry = 0;

% defaults
tlim = [-Inf Inf];
dt = 'auto';
ref = 'sigma';
range = [0 1];

% variables arguments
for i = 1 : 2 : length(varargin)
    varargtyp{i} = varargin{i};
    varargval{i} = varargin{i+1};
    switch lower(varargtyp{i})
        case 'tlim'
            tlim = varargval{i};
        case 'dt'
            dt = varargval{i};
        case 'ref'
            ref = varargval{i};
        case 'range'
            range = varargval{i};
        otherwise
            error('unexpected variable argument type')
    end
end

% check means of depth averaging
fv_check_dave(ref,range)            % throws an error if inputs are non compliant

% generate SWAN input grid
[x_grid, y_grid] = gridme(C);

% indexing from TUFLOW-FV 2D results onto SWAN bc input grid
tmp = fv_get_ids_2([x_grid y_grid],ncfil,'cell',true);
i_in = ~isnan(tmp);
i_fv = tmp(i_in);

% time vector for SWAN  boundary condition files
% this is not necessarily exactly what you specify but shifted to match
% TUFLOW-FV timesteps to avoid timely interpolation (let SWAN do the work)
TMP = netcdf_get_var(ncfil,'names',{'ResTime'});
t = TMP.ResTime;
t = convtime(t);
dt_tmp = round(mean(diff(t) * 24 * 60)) / 60; % hours
switch dt
    case 'auto'
        dt = dt_tmp;
        skip = 1;
    otherwise
        skip = round(dt/dt_tmp);
        if skip < 1
            error('cannot specify a timestep for SWAN bc file smaller than TUFLOW-FV output timestep')
        end
end

its = find(t >= tlim(1),1,'first');
ite = find(t <= tlim(2),1,'last');
if isempty(its) || isempty(ite)
    error('TUFLOW-FV output does not exist within specified timelimits')
end

% it = its:skip:ite;
% nt = length(it);

% names of your output files
[pat,fil,ext] = fileparts(ncfil);
fil_h = [fil '_waterlevels_grid.txt'];
fil_v = [fil '_currents_grid.txt'];
fil_h_control = [fil '_waterlevels_grid_contol.txt'];
fil_v_control = [fil '_currents_grid_contol.txt'];
hfil = fullfile(pat,fil_h);
vfil = fullfile(pat,fil_v);
hfil_control = fullfile(pat,fil_h_control);
vfil_control = fullfile(pat,fil_v_control);

h_fid = fopen(hfil,'w');
v_fid = fopen(vfil,'w');
h_fid_control = fopen(hfil_control,'w');
v_fid_control = fopen(vfil_control,'w');

% intitialise resObj
resObj = fvres_sheet(ncfil,'Variables',{'H';'V_x';'V_y'},'Ref',ref,'Range',range);

% loop through time concurrently extacting TUFLOW-FV results and writing SWAN boundary condition file
fmat = repmat('%8.1f',1,C.mx);
fmat = strcat(fmat,'%8.1f\n');
fmat = repmat(fmat,1,C.my+1);

np = length(x_grid);

display('loading model results and writing SWAN bc files')
inc = 0;
tic;
for aa = its:skip:ite
    
    h_tmp = h_dry * ones(np,1);
    vx_tmp = v_dry * ones(np,1);
    vy_tmp = v_dry * ones(np,1);
    
    % load model results
    set(resObj,'TimeCurrent',t(aa))
%     D = fv_get_dave(D,ncfil,it(aa),'variables',{'H';'V_x';'V_y'},'ref',ref,'range',range);
    
    % index onto grid
    h_tmp(i_in) = resObj.ResultsCell.H(i_fv);
    vx_tmp(i_in) = resObj.ResultsCell.V_x(i_fv);
    vy_tmp(i_in) = resObj.ResultsCell.V_y(i_fv);
    
    % replace NaN's from fv_get_dave flagging dry values
    dry = isnan(h_tmp);
    h_tmp(dry) = h_dry;
    vx_tmp(dry) = v_dry;
    vy_tmp(dry) = v_dry;
    
    % write SWAN boundary condition file
    fprintf(h_fid,fmat,h_tmp);
    fprintf(v_fid,fmat,vx_tmp);
    fprintf(v_fid,fmat,vy_tmp);
    
    inc = mytimer(aa,[its skip ite],inc);
    
end

% write control text
display('writing inserts for SWAN control files')

tmp = floor(t(its) * 24 * 60);
ts_str = datestr(tmp/60/24,'yyyymmdd.HHMMSS');
tmp = floor(t(ite) * 24 * 60);
te_str = datestr(tmp/60/24,'yyyymmdd.HHMMSS');

back = '..\bc\hydros\';
fil_ref_h = ['''',back fil_h,''''];
fil_ref_v = ['''',back fil_v,''''];

fprintf(h_fid_control,'INPGRID WLEV REGULAR %.7f %.7f %d %d %d %.4f %.4f & \n',C.xp,C.yp,C.alp,C.mx,C.my,C.dx,C.dy);
fprintf(h_fid_control,'NONSTATIONARY %s %.3f HR %s\n',ts_str,dt,te_str);
fprintf(h_fid_control,'READINP WLEV +1 %s %d %d FREE',fil_ref_h,1,0);

fprintf(v_fid_control,'INPGRID CURRENT REGULAR %.7f %.7f %d %d %d %.4f %.4f & \n',C.xp,C.yp,C.alp,C.mx,C.my,C.dx,C.dy);
fprintf(v_fid_control,'NONSTATIONARY %s %.3f HR %s\n',ts_str,dt,te_str);
fprintf(v_fid_control,'READINP CURRENT +1 %s %d %d FREE',fil_ref_v,1,0);

fclose('all');


display('done & done :-)')

% /////// nested functions ///////
%
% gridme(C)
% create grid points in order which best suits mother function
% inputs: C = structure containing fields xp,yp,alp,mx,my,dx & dy

function [x, y] = gridme(C)
xp = C.xp; yp = C.yp; alp = C.alp; mx = C.mx; my = C.my; dx = C.dx; dy = C.dy;

np = (mx+1)*(my+1);
x = zeros(np,1);
y = zeros(np,1);

y_vec = [yp+my*dy:-dy:yp];
x_vec = [xp:dx:xp+mx*dx];


k = 1;
for aa = 1:my+1
    kk = k + mx;
    x(k:kk) = x_vec';
    y(k:kk) = y_vec(aa) * ones(mx+1,1);
    k = kk+1;
end

if alp ~= 0
    cosa = cosd(alp);
    sina = sind(alp);
    rot = [cosa, -sina; sina, cosa]';
    newxy = [x-xp, y-yp];
    newxy = newxy*rot;
    x = xp + newxy(:,1);
    y = yp + newxy(:,2);
end





