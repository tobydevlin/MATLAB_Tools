% /////// FV_islands ///////
%
% B = fv_islands(ptsfil)
%
% B = structure with a field defining the mesh boundaries.
% The 1st field defines the outer mesh boundary and the following fields define
% any islands within the mesh. This function is usually used when an image
% of the mesh needs to be drawn (as patch objects) onto an axes or when you
% want to know if your points are within your mesh.
%
% ptsfil = .csv file of points (x,y) created by exporting points from MI
% When working with spherical coordinates it is often necessary to ensure
% that the (x,y) are "decimal" numbers and not "float" when inside MI
%
% The function works because in MI the last point defining a polygon is the
% same as the 1st. ie 5 points define a square.
%
% JN overhauled in July 2013

function B = fv_islands(ptsfil)

% read in .csv
pts = csvread(ptsfil);
tmp = pts;

% identify the boundaries
k = 1;
p1 = tmp(1,:);
while 1+1 == 2       % a dummy condition
    i = ismember(tmp,p1,'rows');
    i = find(i,1,'last');
    names{k} = ['b_',num2str(k)];
    B.(names{k}) = tmp(1:i,:);
    if i == length(tmp)
        break
    else
        tmp = tmp(i+1:end,:);
        p1 = tmp(1,:);
        k = k+1;
    end
end

% ensure the 1st boundary is the mesh boundary
nb = k;
mini = zeros(nb,1);
for aa = 1:nb
    mini(aa) = min(B.(names{aa})(:,1));
end

i = find(mini == min(mini));
if i > 1
    mbound = ['b_' num2str(i)]; % mesh boundary
    tmp = B.b_1;
    B.b_1 = B.(mbound);
    B.(mbound) = tmp;
end