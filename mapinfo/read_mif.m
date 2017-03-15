function mif = read_mif(fnam)

mif = struct();

fid = fopen(fnam);
if fid == -1, error('ERROR:read_mif:fopen', 'Unable to open %s', fnam), end
disp(['Reading ', fnam])

i = 0;
j = 0;
k = 0;
l = 0;
while feof(fid)~=1
    lin = fgetl(fid);
    if length(lin)>2
        switch lower(deblank(lin(1:3)))
            case 'pli' 
                i=i+1;
                tmp = textscan(fid,'%15f %15f');
                mif.poly(i).coords = cell2mat(tmp);
                case 'reg'
                j=j+1;
                lin = fgetl(fid);
                tmp = textscan(fid,'%15f %15f');
                mif.region(j).coords = cell2mat(tmp);
            case 'poi'
                k = k+1;
                tmp = textscan(fid,'%15f %15f');
                mif.point(k).coords = cell2mat(tmp);
            case 'tex'
                l = l + 1;
                lin = fgetl(fid);
                mif.text(l).string = strrep(deblank(lin),'"','');
                tmp = textscan(fid,'%15f %15f');
                mif.text(l).coords = cell2mat(tmp);
      end
    end  
end
fclose(fid);