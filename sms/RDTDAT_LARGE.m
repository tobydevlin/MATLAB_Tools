% Read a single dataset and timestep from a binary SMS .dat file.
% The file must first have been loaded using the LDDAT.m function.
%
% OUT = RDTDAT_LARGE(OUT,name,tout)
%
% Where datfil is the name of the datfil to be read from.
% OUT is the structure array returned by the LDDAT function.
% name is the name of the dataset to be read.
% tout is the timestep to be read (units of time, not an integer timestep)
% the 'nearest neighbour' timestep to tout is returned.
% Refer to function LDDAT_LARGE for info on the OUT fields.
%
% Ian Teakle WBM Pty Ltd

function OUT = RDTDAT_LARGE(OUT,name,tout)

datfil = OUT.datfil;
fid = fopen(datfil);
if fid == -1, disp(['Unable to open ', datfil]), return, end

T = length(OUT.(name).time);
if tout <= OUT.(name).time(1), t = 1;
elseif tout >= OUT.(name).time(T), t = T;
else, t = interp1(OUT.(name).time,1:T,tout,'nearest');
end
OUT.(name).t = OUT.(name).time(t);
fseek(fid, OUT.(name).map(t), 'bof');
istat = fread(fid, 1, OUT.sflg);
temp = fread(fid, 1, OUT.sflt);
if istat == 1, 
    OUT.(name).stat = fread(fid, OUT.numcells, OUT.sflg);
else
    OUT.(name).stat = true(OUT.numcells,1);
end
if ~OUT.(name).vec
    OUT.(name).data = fread(fid, OUT.numdata, OUT.sflt);
else
    OUT.(name).data = fread(fid, [2 OUT.numdata], OUT.sflt)';
end

fclose(fid);
