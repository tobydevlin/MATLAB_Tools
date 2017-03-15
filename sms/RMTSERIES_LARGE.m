% Extract timeseires from an SMS binary .dat file or an RMA10S
% .sol file.  The file must first be loaded using LDDAT_LARGE or LDSOL_LARGE.
%
% TSERIES = RMTSERIES(OUT,nodes,timespan)

function TSERIES = RMTSERIES_LARGE(OUT,nodes,timespan)

J = length(nodes);

if isfield(OUT,'solfil')
    disp(['Extracting timeseries from ',OUT.solfil])
    T = length(OUT.time);
    tstat = OUT.time>=timespan(1) & OUT.time<=timespan(2);
    TSERIES.time = zeros(sum(tstat),1);
    for j = 1 : J
        node = ['N',num2str(nodes(j))];
        TSERIES.(node) = ...
            zeros(sum(tstat), size(OUT.data, 2));
    end
    k = 0;
    for t = 1 : T
        if tstat(t)
            k = k + 1;
            disp(['time = ',num2str(OUT.time(t))])
            OUT = RDTSOL_LARGE(OUT,OUT.time(t));            
            TSERIES.time(k) = OUT.time(t);
            for j = 1 : J
                node = ['N',num2str(nodes(j))];
                TSERIES.(node)(k,:) = ...
                    OUT.data(nodes(j),:);
            end
        end
    end
    
elseif isfield(OUT,'datfil')
    disp(['Extracting timeseries from ',OUT.datfil])
    for n = 1 : length(OUT.names)
        disp(['data name = ',OUT.names{n}])
        name = OUT.names{n};
        T = length(OUT.(name).time);
        tstat = OUT.(name).time>=timespan(1) & OUT.(name).time<=timespan(2);
        TSERIES.(name).time = zeros(sum(tstat),1);
        for j = 1 : J
            node = ['N',num2str(nodes(j))];
            TSERIES.(name).(node) = ...
                zeros(sum(tstat), size(OUT.(name).data, 2));
        end
        k = 0;
        for t = 1 : T
            if tstat(t)
                k = k + 1;
                disp(['time = ',num2str(OUT.(name).time(t))])
                OUT = RDTDAT_LARGE(OUT,name,OUT.(name).time(t));
                TSERIES.(name).time(k) = OUT.(name).time(t);
                for j = 1 : J
                    node = ['N',num2str(nodes(j))];
                    TSERIES.(name).(node)(k,:) = ...
                        OUT.(name).data(nodes(j),:);
                end
            end
        end
    end
else, disp('Unrecognised file type')
end