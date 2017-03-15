% Read a single dataset and timestep from a binary SMS .dat file.
% The file must first have been loaded using the LDDAT.m function.
%
% OUT = RDTDAT(OUT,name,tout)
%
% OUT is the structure array returned by the LDDAT function.
% name is the name of the dataset to be read.
% tout is the timestep to be read (units of time, not an integer timestep)
% the 'nearest neighbour' timestep to tout is returned.
% Refer to function LDDAT for info on the OUT fields.
%
% Ian Teakle WBM Pty Ltd

function OUT = RDTDAT(OUT,name,tout)

T = length(OUT.(name).time);
if tout <= OUT.(name).time(1), t = 1;
elseif tout >= OUT.(name).time(T), t = T;
else, t = interp1(OUT.(name).time,1:T,tout,'nearest');
end
OUT.(name).t = t;
if OUT.(name).istat == 1
    OUT.(name).stat = OUT.(name).map.data(t).stat~=0;
else
    out.(name).stat = true(OUT.numcells,1);
end
OUT.(name).data = OUT.(name).map.data(t).dat';
