% Load a map of an SMS binary .dat file into the MATLAB workspace.
%
% OUT = LDDAT(datfil)
%
% Where datfil is the filename to be loaded.
% OUT is a structure array with the following fields;
% OUT.datfil is the .dat filename
% OUT.names contains the dataset names
% OUT.sflt and OUT.sflg contains info about the numeric precision
% OUT.numdata and OUT.numcells contains the number of nodes and number of
% cells
% OUT.(name).time is a vector of times for that dataset
%           .map is a memory map of the data
%           .vec is logical variable specifying whether the dataset is
%               vector or scalar
%           .data is the data for that dataset
% OUT.t contains the currently loaded timestep (final timestep immediately
% following execution of LDDAT).
%
% After running LDDAT the final timestep of each dataset is active within
% memory.  The RDTDAT function can subsequently be used to access data for
% any timestep of any dataset.
%
% Ian Teakle WBM Pty Ltd

function OUT = LDDAT(datfil)

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
        if sflt == 4, sfltstr = 'single';
        elseif sflt == 8, sfltstr = 'double';
        else disp('Incorrect floating point precision');
        end
        OUT.sflt = sflt;
        
    elseif icard == 120
        sflg = fread(fid, 1, 'int32');
        if sflg == 1, sflgstr = 'int8';
        elseif sflg == 2, sflgstr = 'int16';
        elseif sflg == 4, sflgstr = 'int32';
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
        name = strrep(name,'.','_');
        
        % Initialise output fields for this dataset
        t = 0;
        if isempty(OUT.names), OUT.names{1} = name;
        else, OUT.names = {OUT.names; name};
        end
        OUT.(name).time = zeros(tmax,1);
        if scl
            OUT.(name).vec = false;
            n = 1;
        elseif vec
            OUT.(name).vec = true;
            n = 2;
        end
        % Return location (in bytes) of start of this dataset
        pos = ftell(fid);
        OUT.(name).pos = pos;
        
    % Read data
    elseif icard == 200
        t = t + 1;
        pos = ftell(fid);
        if (pos/fsize*100-prog) > 5
            prog = prog + 5;
            disp([num2str(prog), '%'])
        end
        OUT.(name).istat = fread(fid, 1, sflgstr);
        if t > tmax, disp('Timestep array bounds exceeded'), end
        OUT.(name).time(t) = fread(fid, 1, sfltstr);
        if OUT.(name).istat == 1, fseek(fid, numcells * sflg, 'cof');, end
        fseek(fid, numdata * n * sflt, 'cof');
           
    % End of dataset
    elseif icard == 210
        OUT.(name).t = t;
        OUT.(name).time = OUT.(name).time(1:t);                           
        scl = false;, vec = false;
        
    elseif icard == 250
        fseek(fid,4,'cof');
        
    else icard, disp(['Error encountered reading ', datfil]), return
        
    end
end
%%%%%% Work around files for simulations that did not finish 
if ~isfield(OUT,'t')
        OUT.(name).t = t;
        OUT.(name).time = OUT.(name).time(1:t);                           
        scl = false;, vec = false;
end
% Close datfil
fseek(fid, 0, 'bof');
fclose(fid);

% Create memory map/s
OUT = MEMMAPDAT(OUT);
