

function add_endcard(datfil)

fid = fopen(datfil,'r+');
fseek(fid, 0 ,'eof');
fwrite(fid, 210, 'int32');
fclose(fid);