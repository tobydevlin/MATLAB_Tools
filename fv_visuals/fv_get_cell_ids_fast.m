% --- fv_get_cell_ids_fast --- %
% Nice and quick version of the old favourite
% this uses a face intersect algorithm to get the cell ids
% should be very quick at a large number of points
% 
% Input:
%        modfil - file 
%        xvec   - x values of points
%        yvec   - y values of points
% Output:
%        id     - list of ids corresponding to the points (NaN for missing)
%
%
% TDevlin Nov 2016
%


function id = fv_get_cell_ids_fast(modfil, xvec, yvec)
    tmp = netcdf_get_var(modfil,'names',{'cell_node','node_X','node_Y'});
    vertx = tmp.node_X;   verty = tmp.node_Y;
    cell_node = tmp.cell_node';
    nc3 = sum(cell_node(:,4)==0);  nc4 = length(cell_node)-nc3;
    nfmax = nc4*4+nc3*3;
    faces = zeros(nfmax,2);  cells = zeros(nfmax,1);  ntmp=0;
    np = length(cell_node);
    for aa = 1 : np
        if cell_node(aa,4)==0,ne=3;else ne=4;end
        for bb = 1:ne
            ntmp = ntmp+1;
            [faces(ntmp, :)] = cell_node(aa,[bb,mod(bb,ne)+1]);
            cells(ntmp) = aa;
        end
    end
    
    faces = faces(1:ntmp,:); cells = cells(1:ntmp);
    typ = verty(faces(:,1))<=verty(faces(:,2)); % is right faces
    faces = sort(faces,2);
    
    [~,ind ] = sortrows(faces);
    faces = faces(ind,:);  cells = cells(ind);   typ = typ(ind);
    stat = [true;any(faces(1:end-1,1:2)~=faces(2:end,1:2),2)];
    typ2 = [(stat(1:ntmp-1) & stat(2:ntmp)); stat(end)]; % Ext
    
    cells(typ) = NaN;      logs = find(typ2 | ~typ);
    faces = faces(logs,:); cells = cells(logs);
    
    x1 = vertx(faces(:,1));  x2 = vertx(faces(:,2));
    y1 = verty(faces(:,1));  y2 = verty(faces(:,2));
    M = (x2-x1)./(y2-y1);
    
    id = NaN(np,1);
    for aa = 1:np
        y = yvec(aa);
        yd1 = y-y1;
        yd2 = y-y2;
        isl = find(xor(yd1>=0,yd2>=0));
        x = x1(isl)+M(isl).*yd1(isl);
        [x,ii] = sort(x);
        ithis = find(xvec(aa)<x,1,'last')-1;
        id(aa) = cells(isl(ii(ithis)));
    end
    
end