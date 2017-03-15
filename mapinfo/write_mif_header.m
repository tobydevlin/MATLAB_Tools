
function write_mif_header(fid)

hid = fopen('HEADER.MIF');
if hid ~= -1
    display('Found HEADER.MIF')
    for i = 1 : 4
        line = fgetl(hid);
        fprintf(fid,'%s\n',line);
    end
    fclose(hid);
else
    display('No HEADER.MIF found.  Using non-earth coordinates')
    fprintf(fid,'Version 300\n');
    fprintf(fid,'Charset "Windowsatin1"\n');
    fprintf(fid,'Delimiter ","\n');
    fprintf(fid,'CoordSys NonEarth Units "m" Bounds (-10000000, -10000000) (10000000, 1000000)\n');
end
