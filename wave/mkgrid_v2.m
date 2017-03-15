% /////// mkgrid ///////
% same as Ian's code except the accuracy has been increased to 8 decimals




% Make a uniform grid for input into MAPINFO
%
% mkgrid(outname,xp,yp,alp,mx,my,dx,dy)

function outdata = mkgrid(outname,xp,yp,alp,mx,my,dx,dy)

N = (mx+1)*(my+1);
outdata = zeros(N,5);
% X = xp;
% Y = yp;
% i = 1;
% j = 1;
% outdata(1,:) = [X,Y,i,j,0];
% for n = 2:N
%     if i < (mx + 1)
%         i = i + 1;
%         X = X + dx;
%         outdata(n,:) = [X,Y,i,j,0];
%     else
%         i = 1;
%         j = j + 1;
%         X = xp;
%         Y = Y + dy;
%         outdata(n,:) = [X,Y,i,j,0];
%     end
% end

n = 1;
for j = 1:my+1
    Y = yp + (j-1)*dy;
    for i = 1:mx+1
        X = xp + (i-1)*dx;
        outdata(n,:) = [X,Y,i,j,0];
        n = n + 1;
    end
end        

if alp ~= 0
    cosa = cos(alp/180*pi);
    sina = sin(alp/180*pi);
    rot = [cosa, -sina; sina, cosa]';
    newxy = [outdata(:,1)-xp, outdata(:,2)-yp];
    newxy = newxy*rot;
    outdata(:,1) = xp + newxy(:,1);
    outdata(:,2) = yp + newxy(:,2);
end

fid = fopen(outname,'w');
fprintf(fid, '%s,%s,%s,%s,%s\n', 'X','Y', 'i', 'j', 'depth');
fprintf(fid, '%13.8f,%13.8f,%5.0f,%5.0f,%5.0f\n', outdata');
fclose(fid);



