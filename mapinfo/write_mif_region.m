
function write_mif_region(fid,x,y)

N = size(x,1);
if size(y,1)~=N, error('x & y must be same size'); end
fprintf(fid,'Region  1\n');
fprintf(fid,'%3d\n',N);
for i = 1 : N
    fprintf(fid,'%15f %15f\n',x(i),y(i));
end
ctrd = polycentre(x,y);
fprintf(fid,'    Pen (1,2,16711935)\n');
fprintf(fid,'    Brush (1,0,16777215)\n');
fprintf(fid,'    Center %0.2f %0.2f\n',ctrd(1),ctrd(2));