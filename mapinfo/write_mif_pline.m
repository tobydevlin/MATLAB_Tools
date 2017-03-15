
function write_mif_pline(fid,x,y)

N = size(x,1);
if size(y,1)~=N, error('x & y must be same size'); end
fprintf(fid,'Pline %d\n',N);
for i = 1 : N
    fprintf(fid,'%15f %15f\n',x(i),y(i));
end
fprintf(fid,'    Pen (1,2,0)\n');