% /////// mywave ///////
% very similar to mytide but y axis starts at 0
% produces a small image (axes) showing the current position in the wave signal
% you would usually source the times and wave heights from a wave buoy
%
% h = mytide(ax,tim_s,tim_e,tim_dat,wav_dat,tim_now,varargin)
%
% inputs
%   tim_s = time of start pf plot (matlab convention)
%   tim_e = time of end pf plot (matlab convention)
%   tim_dat = vector of times corresponding to inputed tides (matlab convention)
%   wav_dat = tide levels corresponding to time vector
%   tim_now = current time (matlab convention)
%
% outputs
% h_wav = handle to axes containing information of where in wave tseries you are
%
% variable arguments
% tim_prior = sets axes limit to extend this far (hours) before t_now, default: 12 hours
% tim_after = sets axes limit to extend this far (hours) before t_now, default: tim_prior
% position = vector (normalized to "parent" axes) [xp yp width hight] or string 'topright' 'topleft'
%   'bottomleft' 'bottomright' (default). the width (position(3))
% backing = creates an axes behind the mytide so lables are easily visible, default: false
%
% JN November 2011


function h = mywave(ax,tim_s,tim_e,tim_dat,wav_dat,tim_now,varargin)

% have we created a mywave axes previously
f = get(ax,'parent');
h = findobj(f,'Tag','mywave');
if isempty(h)
    new = true;
elseif length(h) == 1
    new = false;
else
    error('only one mywave permitted per figure')
end

if new
    % defaults
    pos = 'topright';
    isvec = false;
    height = 0.25;
    width = 0.25;
    x_buff = 0.025;
    y_buff = x_buff;
    backing = false;
    
    % variable arguments
    if mod(length(varargin),2) > 0
        error('Expecting variable arguments as descriptor/value pairs')
    end
    
    for i = 1 : 2 : length(varargin)
        varargtyp{i} = varargin{i};
        varargval{i} = varargin{i+1};
        switch lower(varargtyp{i})
            case 'position'
                pos = varargin{i+1};
                if isnumeric(pos)
                    isvec = true;
                end
                case 'backing'
                backing = varargval{i};
                if ~islogical(backing)
                    error('expecting logical input for backing')
                end
            otherwise
                error('unexpected variable argument type')
        end
    end
    
    % checks
    if tim_s ~= round(tim_s) | tim_e ~= round(tim_e)
        error('start and end times must be whole days')
    end
    
 %   if sum(tim_dat < (tim_s - 3/24)) <= 2
 %       error('supply more tide data at start of record to correspond with your start time')
 %   end
    
 %   if sum(tim_dat > (tim_s + 3/24)) <= 2
%        error('supply more tide data at end of record to correspond with your end time')
 %   end
    
    % make parent figure current
    f = get(ax,'Parent');
    figure(f);
    
    % pos of "parent" axes within figure
    f_pos = get(gcf,'position');
    yx_ratio = f_pos(4) / f_pos(3);
    
    ax_pos = get(ax,'position');
    xlim = get(ax,'XLim');
    ylim = get(ax,'YLim');
    dx = diff(xlim);
    dy = diff(ylim);
    fx = ax_pos(3) / dx;
    fy = ax_pos(4) * yx_ratio / dy;
    if abs(fx - fy)/fx > 0.01
        display('WARNING inconsistent scaling between x & y axis')
    end
    
    % replace  pos(3) with 'lim'
    % if isvec
    %     pos(3) = lim / dx;
    % end
    
    % convert any string pos to a normalized pos vector within "parent" axis
    if ~isvec
        switch lower(pos)
            case 'bottomright'
                xp = 1 - x_buff - width;
                yp = y_buff;
            case 'topright'
                xp = 1 - x_buff - width;
                yp = 1 - y_buff - height;
            case 'topleft'
                xp = x_buff;
                yp = 1 - y_buff - height;
            case 'bottomleft'
                xp = x_buff;
                yp = y_buff;
            otherwise
                error('unexpected pos string')
        end
        pos = [xp yp width height];
    end
    
    % update mywave pos vector to reference the parent figure
    pos(1) = ax_pos(1) + ax_pos(3) * pos(1);
    pos(2) = ax_pos(2) + ax_pos(4) * pos(2);
    pos(3) = pos(3) * ax_pos(3);
    pos(4) = pos(4) * ax_pos(4);
    
    % create tide signal
    dt = mean(diff(tim_dat));
    
    its = find(tim_dat >= tim_s,1,'first');
    ite = find(tim_dat <= tim_e,1,'last');
    tim_plot = tim_dat(its:ite);
    wav_plot = wav_dat(its:ite);
    
    % where are you on plot
    wav_now = interp1(tim_plot,wav_plot,tim_now);
    
    % fiddly fiddly
    xlim = [tim_s - 6/24 tim_e + 6/24];
%     lim_tmp = ceil(max(wav_plot));
%     ylim = [0 lim_tmp];
ylim = [0 max(wav_plot) + max(wav_plot)/10];

    % generate plot
    h = axes('parent',f,'Tag','mywave');
    set(h,'position',pos)
    set(h,'Nextplot','add')
    set(h,'color',[0.9 0.9 0.9])
    % set(h,'XColor',[0.9 0.9 0.9])
    % set(h,'YColor',[0.9 0.9 0.9])
    set(h,'XLim',xlim)
    set(h,'YLim',ylim);
    set(h,'XTick',tim_s:tim_e)
    datetick(h,'x','dd','keepticks','keeplimits')
%     set(h,'YTick',ylim(1)+ 1:ylim(2)-1)
    
    plot(h,tim_plot,wav_plot,'-','color','r','LineWidth',1,'Tag','line')
    tmp = plot(h,tim_now,wav_now,'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',4,'Tag','marker');
    set(tmp,'XDataSource','tim_now','YDataSource','wav_now')
    
    % plot(h,tim_plot,wav_plot,'-','color','r','LineWidth',2)
    % plot(h,tim_now,wav_now,'o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',5)
    
    % backing
    % do this prior to setting limits because the ticklables prior todateticking them take up space
    if backing
        tit = get(h,'TightInset');
        posb(1) = pos(1) - tit(1);
        posb(2) = pos(2) - tit(2);
        %         posb(3) = pos(3) + tit(1) + tit(3);
        posb(3) = pos(3) + tit(1) + 0;
        %         posb(4) = pos(4) + tit(2) + tit(4);
        posb(4) = pos(4) + tit(2) + 0;
        hb = axes('parent',f,'position',posb,'color','w','XTick',[],'YTick',[],'XColor','w','YColor','w');
        % -- set behind mytide
        kids = get(f,'children');
        ib = kids == hb;
        i = kids == h;
        kids(ib) = h;
        kids(i) = hb;
        set(f,'children',kids)
    end
    
else
    % leave everything and update the moving marker (position in tide)
    h_lin = findobj(h,'Tag','line');
    tim_plot = get(h_lin,'XData');
    tid_plot = get(h_lin,'YData');
    % -- where are you on plot
    wav_now = interp1(tim_plot,tid_plot,tim_now);
    h_mark = findobj(h,'Tag','marker');
    refreshdata(h_mark,'caller');
end
% set(h,'FontWeight','bold')
