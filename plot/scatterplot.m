% /////// scatterplot ///////
% plots currents onto scatter plot.
%
% inputs
%   ax = axes to plot into
%   vx = x-component of velocity (m/s)
%   vy = y-component of velocity (m/s)
%
% optional inputs as discriptor value pairs
%   name / name = DisplayName ie DATA or MODEL to appear in legend
%   color / color vec  = either [0 0 1], [0 1 0] or [1 0 0];
%   alpha / transperancy 0-1 transparency of markers
%   draw / logical = true to draw radial plots, false to draw onto existing ones
%   markersize / size of triangle marker base, fraction of axis limits
%   rmax / rmax =  radial distance limit on plot
%   rticks /rticks =  # ticks indicating distance

function scatterplot(ax,vx,vy,varargin)

% defaults
name = 'SERIES 1';
color = [0 0 1];
draw = true;
rmax = 100;
rticks = 5;
markersize = 1/200;
alp = 0.1;

% variable arguments
if mod(length(varargin),2) > 0
    error('Expecting variable arguments as descriptor/value pairs')
end

for i = 1 : 2 : length(varargin)
    varargtyp{i} = varargin{i};
    varargval{i} = varargin{i+1};
    switch lower(varargtyp{i})
        case 'name'
            name = varargval{i};
        case 'color'
            col = varargval{i};
            if ~isnumeric(col)
                error('color input must be 1x3 rgb vector')
            end
        case 'alpha'
            alp = varargval{i};
            if alp <= 0 || alp > 1
                error('input "alpha" must be within (0 1]')
            end
        case 'markersize'
            markersize = varargval{i};
            if markersize ~= round(markersize)
                error('input "markersize" must be an integer')
            end
        case 'draw'
            draw = varargval{i};
            if ~islogical(draw)
                error('expecting logical input for "draw"')
            end
        case 'rmax'
            rmax = varargval{i};
        case 'rticks'
            rticks = varargval{i};
        otherwise
            error('unexpected variable argument type')
    end
end

% make inputted axes current
f = get(ax,'parent');
set(f,'CurrentAxes',ax)
set(ax,'NextPlot','add')
set(ax,'visible','off')

% h = hggroup;
% set(h,'parent',ax)

if draw
    % Create polar axis to work with
    polar(ax,rmax,rticks,false);
    axis fill
    frac = 0.6;
    dy_new = 2*rmax / frac;
    t_buff = dy_new * 0.1;
    ylim(2) = rmax + t_buff;
    ylim(1) = ylim(2) - dy_new;
    set(ax,'YLim',ylim);
    axis equal
    
    pos = get(ax,'position');
    
    % Create current magnitude legend
    %     frac = 0.30;
    frac = 0.40;
    h_frq = axes;
    polar(h_frq,rmax,rticks,true);
    axis fill;
    axis equal;
    set(gca,'visible','off')
    xp = pos(1) + (pos(3)/2-frac*pos(4))/2;
    yp = pos(2);
    width = frac * pos(4);
    height = frac * pos(4);
    
    set(h_frq,'position',[xp yp width height])
end

% plot data
if ~isempty(vx)
    if length(vx) ~= length(vy)
        error('mismatch between x and y component velocities')
    end
    
    if size(vx,2) == 1
        vx = vx';
    end
    
    if size(vy,2) == 1
        vy = vy';
    end
    
    nm = length(vx);
    
    % marker size
    xlim = get(ax,'xlim');
    mx = markersize * diff(xlim);
    my = mx;
    tri_x = repmat([-mx/2;mx/2;0],1,nm);
    tri_y = repmat([-my/2;-my/2;my/2],1,nm);
    
    tri_x = tri_x + repmat(vx,3,1);
    tri_y = tri_y + repmat(vy,3,1);
    
    p = patch(tri_x,tri_y,col);
    set(p,'Edgecolor','none')
    if alp < 1
        set(gcf,'Renderer','OpenGL')
        set(p,'AlphaDataMapping','none')
        set(p,'FaceAlpha',alp)
    end
    alpha(p,alp)
    set(p,'DisplayName',name)
end

% /////// nested functions ///////

function [hpol,rmax] = polar(ax,rmax,rticks,leg)

f = get(ax,'parent');
set(f,'CurrentAxes',ax);

h = hggroup;
set(h,'parent',ax)
ha = get(h,'Annotation');
hl = get(ha,'LegendInformation');
set(hl,'IconDisplayStyle','off')

% cax = get(ax,'parent');
% type = get(cax,'Type');
% switch type
%     case 'axes'
%         % good
%     otherwise
%         error('parent of "ax" must be an axes')
% end

tc = get(ax,'xcolor');         % line color
ls = get(ax,'gridlinestyle');  % line style

% make a radial grid
rmin = 0;

% define a circle
th = 0:pi/50:2*pi;
xunit = cos(th);
yunit = sin(th);

% now really force points on x/y axes to lie on them exactly
inds = 1:(length(th)-1)/4:length(th);
xunit(inds(2:2:4)) = zeros(2,1);
yunit(inds(1:2:5)) = zeros(3,1);

% draw radial circles
c89 = cos(89*pi/180);
s89 = sin(89*pi/180);
rinc = (rmax-rmin)/rticks;
for i=(rmin+rinc):rinc:rmax
    hhh = line(xunit*i,yunit*i,'linestyle',ls,'color',tc,'linewidth',1,'parent',h);
    %     set(cax,'Layer','Bottom')
    if leg
        text((i+rinc/20)*c89,(i+rinc/20)*s89,...
            [num2str(i) ' m/s'],'verticalalignment','bottom','parent',h,...
            'fontweight','bold','Color',[0 0 0])
    end
end
set(hhh,'linestyle','-') % Make outer circle solid

% spokes
if ~leg
    % -- plot spokes
    th = (1:4)*2*pi/8;
    %     th = (1:8)*2*pi/16;
    cst = cos(th);
    snt = sin(th);
    cs = [-cst; cst];
    sn = [-snt; snt];
    line(rmax*cs,rmax*sn,'linestyle',ls,'color',tc,'linewidth',1,'parent',h)
    txt = {'SW','S','SE','E','NE','N','NW','W'};
    %     txt = {'WSW','SW','SSW','S','SSE','SE','ESE','E','ENE','NE','NNE','N','NNW','NW','WNW','W'};
    
    % -- annotate spokes
    rt = 1.1*rmax;
    for i = 1:length(th)
        
        text(rt*cst(i),rt*snt(i),txt{i+length(th)},...
            'horizontalalignment','center','parent',h,...
            'fontweight','bold');
        
        text(-rt*cst(i),-rt*snt(i),txt{i},'horizontalalignment','center','parent',h,...
            'fontweight','bold')
    end
end
