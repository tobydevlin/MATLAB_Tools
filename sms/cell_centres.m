% /////// cell_centres ///////
%
% generates the .csv file of cell centres corresponding to the inputted 2DM
% ready for inspection in MI or use as a scatter set in SMS
%
% inputs
%   2dm_fil
%   outfil = .csv file of cell centres corresponding to 2dm_fil
%
% JN

function cell_centres(infil,outfil)

% read 2DM file
MESH = RD2DM(infil);

% get cell centres
if (isfield(MESH,'E3T'))
    ne3 = length(MESH.E3T);
else
    ne3 = 0;
end
if (isfield(MESH,'E4Q'))
    ne4 = length(MESH.E4Q);
else
    ne4 = 0;
end

ne = ne3 + ne4;
ctrd = zeros(ne,2);
id = zeros(ne,1);
z = zeros(ne,1);
k = 1;

if (isfield(MESH,'E3T'))
    for i = 1:length(MESH.E3T)
        pts = MESH.E3T(i,2:4);
        x = MESH.ND(pts,2);
        y = MESH.ND(pts,3);
        z(k) = mean(MESH.ND(pts,4));
        ctrd(k,:) = polycentre(x,y);
        id(k) = MESH.E3T(i,1);
        k = k+1;
    end
end

if (isfield(MESH,'E4Q'))
    for i = 1:length(MESH.E4Q)
        pts = MESH.E4Q(i,2:5);
        x = MESH.ND(pts,2);
        y = MESH.ND(pts,3);
        z(k) = mean(MESH.ND(pts,4));
        ctrd(k,:) = polycentre(x,y);
        id(k) = MESH.E4Q(i,1);
        k = k+1;
    end
end

% order cell ids
[id i] = sort(id,'ascend');
ctrd = ctrd(i,:);
z = z(i);

% write .csv file
fid = fopen(outfil,'w');
fprintf(fid,'%s\n','ID,X,Y,Z');

for aa = 1:ne
    fprintf(fid,'%i,%.7f,%.7f,%.7f\n',id(aa),ctrd(aa,1),ctrd(aa,2),z(aa));
end

fclose(fid);

display('done & done :-)')