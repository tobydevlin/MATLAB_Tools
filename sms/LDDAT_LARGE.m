% Load a map of an SMS binary .dat file into the MATLAB workspace.
%
% OUT = LDDAT_LARGE(datfil)
%
% Where datfil is the filename to be loaded.
% OUT is a structure array with the following fields;
% OUT.names contains the dataset names
% OUT.sflt and OUT.sflg contains info about the numeric precision
% OUT.numdata and OUT.numcells contains the number of nodes and number of
% cells
% OUT.(name).time is a vector of times for that dataset
%           .map is a vector of locations (in bytes from bof) of the data
%               for each timestep of that dataset
%           .vec is logical variable specifying whether the dataset is
%               vector or scalar
%           .data is the data for that dataset
% OUT.t contains the currently loaded timestep (final timestep immediately
% following execution of LDDAT).
%
% After running LDDAT_LARGE the final timestep of each dataset is active within
% memory.  The RDTDAT_LARGE function can subsequently be used to access data for
% any timestep of any dataset.
%
% Ian Teakle WBM Pty Ltd

function OUT = LDDAT_LARGE(datfil)

% Initialise some variables
tmax = 50000;
OUT = struct();
OUT.datfil = datfil;
OUT.names = {};
scl = false;
vec = false;

% Open datfil
fid = fopen(datfil);
if fid == -1, disp(['Unable to open ', datfil]), return, end
fseek(fid, 0 ,'eof');
fsize = ftell(fid);
fseek(fid, 0, 'bof');
prog = 0;
disp(['Loading ', datfil])
disp([num2str(prog), '%'])

% Read filetype
vers = fread(fid, 1, 'int32');
if vers~=3000
    disp('File is not a correct SMS binary .dat format'), return
end

% Loop through data identifying card id and reading data as appropriate
while feof(fid)~=1
    icard  = fread(fid, 1, 'int32');
    if isempty(icard), break, end
    
    if icard == 100, objtype = fread(fid, 1, 'int32');
        
    elseif icard == 110
        sflt = fread(fid, 1, 'int32');
        if sflt == 4, sflt = 'float32';
        elseif sflt == 8, sflt = 'float64';
        else disp('Incorrect floating point precision');
        end
        OUT.sflt = sflt;
        
    elseif icard == 120
        sflg = fread(fid, 1, 'int32');
        if sflg == 1, sflg = 'int8';
        elseif sflg == 2, sflg = 'int16';
        elseif sflg == 4, sflg = 'int32';
        else disp('Incorrect integer precision');
        end
        OUT.sflg = sflg;
        
    elseif icard == 130, scl = true;
        
    elseif icard == 140, vec = true;
        
    elseif icard == 150, vectype = fread(fid, 1, 'int32');
        
    elseif icard == 160, objid = fread(fid, 1, 'int32');
        
    elseif icard == 170
        numdata = fread(fid, 1, 'int32');
        OUT.numdata = numdata;
        
    elseif icard == 180
        numcells = fread(fid, 1, 'int32');
        OUT.numcells = numcells;
    
    % Read new dataset name
    elseif icard == 190
        name = fread(fid, [1 40], 'char');
        name = char(name(1:max(find(name))));
        name = strrep(name,' ','');
        name = strrep(name,'(','');
        name = strrep(name,')','');
        name = strrep(name,'+','_');
        name = strrep(name,'%','');
        name = strrep(name,'-','_');
        
        % Initialise output fields for this dataset
        t = 0;
        if isempty(OUT.names), OUT.names{1} = name;
        else, OUT.names = {OUT.names; name};
        end
        OUT.(name).time = zeros(tmax,1);
        OUT.(name).map = zeros(tmax,1);
        if scl
            OUT.(name).vec = false;
            OUT.(name).data = zeros(numdata,1);
        elseif vec
            OUT.(name).vec = true;
            OUT.(name).data = zeros(numdata,2);
        end
        
    % Read data
    elseif icard == 200
        t = t + 1;
        pos = ftell(fid);
        if (pos/fsize*100-prog) > 5
            prog = prog + 5;
            disp([num2str(prog), '%'])
        end
        % Return location (in bytes) of start of data for this timestep of
        % this dataset
        OUT.(name).map(t) = ftell(fid);
        istat = fread(fid, 1, sflg);
        if t > tmax, disp('Timestep array bounds exceeded'), end
        OUT.(name).time(t) = fread(fid, 1, sflt);
        if istat == 1, 
            OUT.(name).stat = fread(fid, numcells, sflg);
        else
            OUT.(name).stat = true(numcells,1);
        end
        if scl
            OUT.(name).data = fread(fid, numdata, sflt);
        elseif vec
            OUT.(name).data = fread(fid, [2 numdata], sflt)';
        end
           
    % End of dataset
    elseif icard == 210
        OUT.(name).map = OUT.(name).map(1:t);
        OUT.(name).time = OUT.(name).time(1:t);
        scl = false;, vec = false;
        
    elseif icard == 250
        fseek(fid,4,'cof');
        
    else icard, disp(['Error encountered reading ', datfil]), return
        
    end
end

% Close datfil
fclose(fid);