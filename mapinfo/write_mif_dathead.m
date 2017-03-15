
function write_mif_dathead(fid,dathead)

if ~iscell(dathead), error('dathead should be a cell array'), end
N = size(dathead,1);
if size(dathead,2)~=2
    error('dathead should have two columns; data_name and data_type')
end
fprintf(fid,'Columns %d\n',N);
for i = 1 : N
    fprintf(fid,' %s %s\n',dathead{i,1},dathead{i,2});
end
fprintf(fid,'Data\n\n');