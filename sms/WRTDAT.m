% Write an SMS binary .dat file
%
% WRTDAT(datfil,dataname,data,time,numcells)
% WRTDAT(datfil,dataname,data,time,numcells,'append')
% WRTDAT(datfil,dataname,data,time,numcells,'overwrite')
% WRTDAT(datfil,dataname,data,time,numcells,...,stat)
%
% Where;
% datfil is the filename to be created.
% dataname is the dataset name e.g. "velocity"
% data is the data in a Nx1xT or Nx2xT size matrix depending on whether the
% dataset is scalar or vector.
% time is an Tx1 vector of Matlab datenumbers.
% numcells is the number of elements in the model.
% an optional sixth argument can be provided to specify whether the write
% should be appended to an existing file of the same name or whether this
% file should be overwritten (default).
% an optional seventh argument can be input to specify the cell wet/dry
% status.  stat should be a vector of size numcellsx1xT populated
% by 0s (dry) and 1s (wet).
%
% Ian Teakle WBM Pty Ltd

function WRTDAT(datfil,dataname,data,time,numcells,varargin)

append = false;

N = size(data,1);
T = length(time);
if size(data,3)~=T, disp('Inconsistent number of timesteps'), return, end
stat = ones(numcells,1,T);

% Deal with the input arguments
narg = nargin;
if narg == 5
elseif (narg == 6) || (narg == 7)
    if strcmpi(varargin{1},'append')
        append = true;
    elseif strcmpi(varargin{1},'overwrite')
    else
        disp('Sixth argument must be either ''append'' or ''overwrite''')
        return
    end
    if narg == 7
        stat = varargin{2};
        if numel(stat) ~= numcells*1*T;
            disp('Seventh argument must be a numeric array of size Nx1xT')
            return
        end
    end
else
    error('Incorrect number of arguments');
end

if append == false
    disp(['Writing ',datfil])
    fid = fopen(datfil,'w');
else
    fid = fopen(datfil,'r+');
end

if fid == -1
    append = false;
    [fid, message] = fopen(datfil,'w');
    disp(message)
end

if append == false
    % VERSION
    fwrite(fid, 3000, 'int32');
    % OBJTYPE
    fwrite(fid, [100 3], 'int32');
    % SFLT
    fwrite(fid, [110 4], 'int32');
    % SFLG
    fwrite(fid, [120 1], 'int32');
    % TIME UNITS
    fwrite(fid,[250 0], 'int32');
    % BEGSCL / BEGVEC
    if size(data, 2) == 1
        fwrite(fid, 130, 'int32');
    elseif size(data, 2) == 2
        fwrite(fid, 140, 'int32');
    else disp('Too many columns in input data');
    end
    % NUMDATA
    fwrite(fid, [170 N], 'int32');
    % NUMCELLS
    fwrite(fid, [180 numcells], 'int32');
    % NAME
    fwrite(fid, 190, 'int32');
    dataname = dataname(1:min(40,length(dataname)));
    temp = [dataname, zeros(1, 40 - length(dataname))];
    fwrite(fid, temp, 'char');
else
    fseek(fid, -4, 'eof');
    if fread(fid, 1, 'int32') == 210;
        fseek(fid, -4, 'eof');
    else
        disp('Error appending to file')
    end
end


% TS
for t = 1:T

    fwrite(fid, 200, 'int32');
    % istat
    fwrite(fid, 1, 'int8');
    % time
    fwrite(fid, time(t), 'float32');
    % stat
    fwrite(fid, stat(:,1,t), 'int8');
    % val
    fwrite(fid, data(:,:,t)', 'float32');
end

% ENDDS
fwrite(fid, 210, 'int32');
    
fclose(fid);

        



