% Calculate mean, min and max from an SMS binary .dat file or an RMA10S
% .sol file.  The file must first be loaded using LDDAT or LDSOL.

function RES = MMM(OUT,timespan)

RES = struct();
if isfield(OUT,'solfil')
    sol = true;
    infil = OUT.solfil;
elseif isfield(OUT,'datfil')
    sol = false;
    infil = OUT.datfil;
else
    disp('Unrecognised file type')
    return
end

% .sol file
if sol
    disp(['Calculating mean, max and min of ',OUT.solfil,...
    ' between ',num2str(timespan(1)),' and ',num2str(timespan(2)),'.'])
    T = length(OUT.time);
    RES.data_head = OUT.data_head;
    meanval = zeros(size(OUT.data));
    cumsum = meanval;
    maxval = ones(size(OUT.data))*-9999999;
    minval = ones(size(OUT.data))*9999999;
    count = 0;
    for t = 1 : T
        disp(['t = ',num2str(t)])
        if OUT.time(t) >= timespan(1) ...
            && OUT.time(t) <= timespan(2)
            %OUT = RDTSOL(OUT,OUT.time(t));
            OUT = RDTSOL2(OUT,OUT.time(t));
            count = count + 1;
            cumsum = cumsum + double(OUT.data);
            maxval = maxval .* (OUT.data <= maxval) + ...
                double(OUT.data) .* (OUT.data > maxval);
            minval = minval .* (OUT.data >= minval) + ...
                double(OUT.data) .* (OUT.data < minval);
        end     
    end
    if count > 0, meanval = cumsum / count;, end
    RES.mean = meanval;
    RES.min = minval;
    RES.max = maxval;

% .dat file
else
    disp(['Calculating mean, max and min of ',OUT.datfil,...
    ' between ',timespan(1),' and ',timespan(2),'.'])
    for n = 1 : length(OUT.names)
        name = OUT.names{n};
        T = length(OUT.(name).time);
        meanval = zeros(size(OUT.(name).data));
        cumsum = meanval;
        maxval = ones(size(OUT.(name).data))*-9999999;
    minval = ones(size(OUT.(name).data))*9999999;
        count = 0;
        for t = 1 : T
            disp(['t = ',num2str(t)])
            if OUT.(name).time(t) >= timespan(1) ...
                    && OUT.(name).time(t) <= timespan(2)
                %OUT = RDTDAT(OUT,name,OUT.(name).time(t));
                OUT = RDTDAT(OUT,name,OUT.(name).time(t));
                count = count + 1;
                cumsum = cumsum + double(OUT.(name).data);
                maxval = maxval .* (OUT.(name).data <= maxval) + ...
                    double(OUT.(name).data) .* (OUT.(name).data > maxval);
                minval = minval .* (OUT.(name).data >= minval) + ...
                    double(OUT.(name).data) .* (OUT.(name).data < minval);
            end
        end
        if count > 0, meanval = cumsum / count;, end
        RES.(name).mean = meanval;
        RES.(name).min = minval;
        RES.(name).max = maxval;
    end
end