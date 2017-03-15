%               //// runtimes_plot \\\\
% Plots cool figures of the cumulative and instantaneous runtime ratios of
% a TuflowFV sim
% Useful for determining the fast and slow parts of your model
%
% input:
%       log - the full path to the log file in windows format (right click
%             ultraedit and copy file name/path)
%       secplot - if anything is in here, another plot of external and
%                 internal timesteps is made
% 
% output:
%       pretty pictures
%
% TDevlin Nov2015
%


function runtimes_plot(log,varargin)
    secplot=false;
    if nargin==2
        secplot=true;
    end
    % sort out paths
    if isunix
        log = strrep(log,'L:\','/mnt/coastal/');
        log = strrep(log,'P:\','/mnt/catchments/');
        log = strrep(log,'K:\','/mnt/rivers/');
        log = strrep(log,'O:\','/mnt/oceanics/');
        log = strrep(log,'\\blaydos','');
        log = strrep(log,'Y:\','/scratch3/');
        log = strrep(log,'X:\','/scratch2/');
        log = strrep(log,'\','/');
    end

    % allocate locals
    st = cell(1000000,1); % simulation time
    rt = cell(1000000,1); % real time
    iT = cell(1000000,1);
    eT = cell(1000000,1);
    a=0;
    check = double('t = ');
    check2 = double('Entering timestep loop');
    check3  = double('END TIME == ');
    start = false;

    % Open file and step through lines
    fid = fopen(log);
    while ~feof(fid)
        fline = fgets(fid);
        if start % are we in the timesteps yet?
            if length(fline)>3 
                if double(fline(1:4))==check % is this a timestep line?
                    a=a+1;
                    rt{a} = fline(67:end-4);
                    st{a} = fline(5:23);
                    iT{a} = fline(32:38);
                    eT{a} = fline(43:47);
                end
            end
        else
            if length(fline)>21
                if double(fline(1:22))==check2 % are we in the timesteps yet?
                    start = true;
                end
            end
            if length(fline)>12
                if double(fline(1:12))==check3
                    tend = datenum(fline(13:end),'dd/mm/yyyy HH:MM:SS');
                end
            end
        end
    end

    % Turn strings into numbers and calculate Ratios
    st = datenum(st(1:a),'dd/mm/yyyy HH:MM:SS');
    rt = str2double(rt(1:a));
    iT = str2double(iT(1:a));
    eT = str2double(eT(1:a));
    inst = diff(st)*24*3600./diff(rt);
    cumu = (st-st(1))*24*3600./rt;

    % Plot
    figure
    set(gca,'position',[0.1 0.1 0.85 0.85])
    line('xdata',st(2:end),'ydata',inst,'color','b','displayname','Instantaneous Runtime')
    line('xdata',st,'ydata',cumu,'color','r','displayname','Cumulative Runtime')
    dynamicDateTicks(gca,[],'dd/mmm')
    legend(gca,'show')
    xlabel('Simulation Time')
    ylabel('Runtime Ratio')
    if exist('tend','var')
        tdif = (tend-st(end));
        tleft=datevec(tdif./cumu(end));
        if all(tleft<1)
            title('Finished');
        else
            names = {'Years','Months','Days','Hours','Minutes','Seconds'};
            ii = find(tleft>0,1,'first');
            this = [];
            for aa=ii:6
                this = [this , sprintf('%2.0f %s ',tleft(aa),names{aa})];
            end
            title([this 'Remaining']);
        end
    end
    if secplot
        figure
        set(gca,'position',[0.1 0.1 0.85 0.85])
        line('xdata',st,'ydata',iT,'color','b','displayname','Internal Loop')
        line('xdata',st,'ydata',eT,'color','r','displayname','External Loop')
        dynamicDateTicks(gca,[],'dd/mmm')
        legend(gca,'show')
        xlabel('Simulation Time')
        ylabel('Timestep')
    end
end