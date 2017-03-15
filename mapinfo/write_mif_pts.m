
function write_mif_pts(fid,x,y)

fmt = '(35,0,12)';
N = size(x,1);
if size(y,1)~=N, error('x & y must be same size'); end
for i = 1 : N
    fprintf(fid,'Point %15f %15f\n',x(i),y(i));
    fprintf(fid,'    Symbol %s\n',fmt);
end
    