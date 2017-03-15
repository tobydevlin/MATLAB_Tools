% Create a memory map for an SMS .dat file

function OUT = MEMMAPDAT(OUT)

J = length(OUT.names);
for j = 1 : J

    datfil = OUT.datfil;
    numdata = OUT.numdata;
    numcells = OUT.numcells;
    name = OUT.names{j};
    sflt = OUT.sflt;
    sflg = OUT.sflg;
    t = OUT.(name).t;

    if sflt == 4, sfltstr = 'single';
    elseif sflt == 8, sfltstr = 'double';
    else disp('Incorrect floating point precision');
    end

    if sflg == 1, sflgstr = 'int8';
    elseif sflg == 2, sflgstr = 'int16';
    elseif sflg == 4, sflgstr = 'int32';
    else disp('Incorrect integer precision');
    end

    if OUT.(name).vec == false, n = 1; else n = 2; end

    if OUT.(name).istat == 1
        Format = {'int32' 1 'icard'; ...
            sflgstr 1 'istat'; ...
            sfltstr 1 'time'; ...
            sflgstr numcells 'stat'; ...
            'single' [n numdata] 'dat'};
    else
        Format = {'int32' 1 'icard'; ...
            sflgstr 1 'istat'; ...
            sfltstr 1 'time'; ...
            sfltstr [n numdata] 'dat'};
    end
    
    OUT.(name).map = memmapfile(datfil, ...
        'Offset', OUT.(name).pos, ...
        'Format', Format, ...
        'Repeat', length(OUT.(name).time));
    
    if OUT.(name).istat == 1
        OUT.(name).stat = OUT.(name).map.data(t).stat~=0;
    else
        out.(name).stat = true(numcells,1);
    end
    
    OUT.(name).data = OUT.(name).map.data(t).dat';

end