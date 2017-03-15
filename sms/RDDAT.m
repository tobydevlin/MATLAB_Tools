% Read SMS binary .dat file into the MATLAB workspace

function OUT = RDDAT(datfil)

tmax = 1500;

OUT = struct();
scl = false;
vec = false;

fid = fopen(datfil);
if fid == -1, disp(['Unable to open ', datfil]), return, end

vers = fread(fid, 1, 'int32');
if vers~=3000
    disp('File is not a correct SMS binary .dat format'), return
end

while feof(fid)~=1
    icard  = fread(fid, 1, 'int32');
    if isempty(icard), break, end
    
    if icard == 100, objtype = fread(fid, 1, 'int32');
        
    elseif icard == 110
        sflt = fread(fid, 1, 'int32');
        if sflt == 4, sflt = 'float32';,
        elseif sflt == 8, sflt = 'float64';;
        else disp('Incorrect floating point precision');
        end
        
    elseif icard == 120
        sflg = fread(fid, 1, 'int32');
        if sflg == 1, sflg = 'int8';,
        elseif sflg == 2, sflg = 'int16';;
        elseif slfg == 4, sflg = 'int32';;
        else disp('Incorrect integer precision');
        end
        
    elseif icard == 130, scl = true;
        
    elseif icard == 140, vec = true;
        
    elseif icard == 150, vectype = fread(fid, 1, 'int32');
        
    elseif icard == 160, objid = fread(fid, 1, 'int32');
        
    elseif icard == 170, numdata = fread(fid, 1, 'int32');
        
    elseif icard == 180, numcells = fread(fid, 1, 'int32');
        
    elseif icard == 190
        name = fread(fid, [1 40], 'char');
        name = char(name(1:max(find(name))));
        
        t = 0;
        OUT.(name).time = zeros(tmax,1);
        if scl, OUT.(name).data = zeros(numdata,1,tmax);
        elseif vec, OUT.(name).data = zeros(numdata,2,tmax);
        end
        
    elseif icard == 200
        istat = fread(fid, 1, sflg);
        t = t + 1
        if t > tmax, disp('Timestep array bounds exceeded'), end
        OUT.(name).time(t) = fread(fid, 1, sflt);
        if istat == 1, stat = fread(fid, numcells, sflg);, end
        if scl
            OUT.(name).data(:,:,t) = fread(fid, numdata, sflt);
        elseif vec
            OUT.(name).data(:,:,t) = fread(fid, [numdata 2], sflt);
        end
               
    elseif icard == 210
        OUT.(name).time = OUT.(name).time(1:t);
        OUT.(name).data = OUT.(name).data(:,:,1:t);
        scl = false;, vec = false;
        
    elseif icard == 250
        fseek(fid,4,'cof');
        
    else icard, disp(['Error encountered reading ', datfil]), fclose(fid), return
        
    end
end

fclose(fid);