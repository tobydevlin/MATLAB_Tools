% --------- fv2grid ------------
% Takes TUFLOW FV data and inspects it onto a grid.
% The grid specified is computed from the models extent - unless 'bounds'
% is specified - and the grid spacing 'grd'
% Functionality is much the same as fv2ascii though outputting the gridded
% results instead of writing a file.
% Uses ray-casting like graphics algorithm to inspect grid.
%
% Inputs:
%         resfil - TUFLOW FV results file with geo data
%         data   - data to inspect off, must match above files dimensions
%         xv     - vector of x values for grid
%         yv     - vector of y values for grid
%         
% Optional:
%         idmap  - mapping output from previous call, to speed things up
%
%
% Output:
%         res_grd- data inspected onto the grid
%         idmap  - mapping of original scatter to the grid
%
% TDEVLIN Dec 2016


function [res_grd,idmap] = fv2grid(resfil,data,xv,yv,varargin)

% defaults
idmap = [];

% variable arguments
if nargin>4
    idmap = varargin{1};
end

% generate the grid
% -- extents
TMP = netcdf_get_var(resfil,'names',{'node_X';'node_Y';'cell_node'});
vertx = TMP.node_X; verty = TMP.node_Y;
cell_node = TMP.cell_node';

xvec = sort(xv);
yvec = sort(yv,'descend');
nx = length(xvec);
ny = length(yvec);

% Get idmaps
if isempty(idmap)
    % Loop through and find total faces
    nc3 = sum(cell_node(:,4)==0);  nc4 = length(cell_node)-nc3;
    nfmax = nc4*4+nc3*3;
    faces = zeros(nfmax,2);  cells = zeros(nfmax,1);  typ = false(nfmax,1);
    ntmp=0;
    for aa = 1 : length(cell_node)
        if cell_node(aa,4)==0,ne=3;else ne=4;end
        for bb = 1:ne
            ntmp = ntmp+1;
            [faces(ntmp, :)] = cell_node(aa,[bb,mod(bb,ne)+1]);
            cells(ntmp) = aa;
        end
    end
   
    faces = faces(1:ntmp,:);
    cells = cells(1:ntmp);
    typ = verty(faces(:,1))<=verty(faces(:,2)); % is right faces
    faces = sort(faces,2);
    
    
    [~,ind ] = sortrows(faces);
    faces = faces(ind,:);  cells = cells(ind);   typ = typ(ind);
    stat = [true;any(faces(1:end-1,1:2)~=faces(2:end,1:2),2)];
    typ2 = [(stat(1:ntmp-1) & stat(2:ntmp)); stat(end)]; % Ext
    
    cells(typ) = NaN;    
    logs = find(typ2 | ~typ);
    
    faces = faces(logs,:);
    cells = cells(logs);
    
    x1 = vertx(faces(:,1));
    x2 = vertx(faces(:,2));
    y1 = verty(faces(:,1));
    y2 = verty(faces(:,2));
    M = (x2-x1)./(y2-y1);
    
    idmap = NaN(ny,nx);
    for aa = 1:ny
        y = yvec(aa);
        yd1 = y-y1;
        yd2 = y-y2;
        isl = find(xor(yd1>=0,yd2>=0));
        x = x1(isl)+M(isl).*yd1(isl);
        ctmp = cells(isl);
        [x,ii] = sort(x);
        ctmp = ctmp(ii);
        for bb = 1 : length(x)
            idmap(aa,xvec>=x(bb)) = ctmp(bb);
        end
    end
else
    % Check if idmap matches size 
    if max(idmap)>size(cell_node,1);
        error('idmap doesn''t match number of cells in model')
    end
    if ~isequal(size(idmap),[ny,nx])
        error('idmap of wrong size')
    end
end

% preallocate
res_grd = NaN(ny,nx);
log = ~isnan(idmap);
res_grd(log) = data(idmap(log));

end