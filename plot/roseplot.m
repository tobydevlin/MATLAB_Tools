% /////// roseplot ///////
% based on 'windroseplot.m but plots the "pieces of pie" the same width as the interval in theta_bin
% all inputs are in nautical convention, which type of nautical convention is specied by the input 'type'
%
% inputs
%   ax = axes to plot roseplot on
%   mag = magnitude component of data
%   theta = directional component of data (nautical)
%   mag_bin = bins for mag, eg mag_bin = [0:2:4] creates 3 bins [<0 & 0-2,2-4,>4]
%   theta_bin = [0:dtheta:360], binning as performed as above on [theta_bin(1):theta_bin(end-1)]
%   rmax = percentage limit on plot
%   rticks = #ticks determines percentage intervals
%   type = 'wind', 'wave' or 'current' - determines the convention
%
%   variable arguments as descriptor/value pairs:
%   'freq_table'/freq - specified frequency table with size [length(theta_bin) length(mag_bin)]
%
% calls on the nested function hist2D_JN

function freq=roseplot(ax,mag,theta,mag_bin,theta_bin,rmax,rticks,type,varargin)

freq=[];
if mod(nargin,2)>0
    error('Expecting variable arguments as descriptor/value pairs')
end
for i = 1 : 2 : nargin-8
    varargtyp{i} = varargin{i};
    varargval{i} = varargin{i+1};
    switch varargtyp{i}
        case 'freq_table'
            freq = varargval{i};
        otherwise
            error('unexpected variable argument type')
    end
end

if ~isvector(theta_bin)&&~isvector(mag_bin)
    error('theta_bin and mag_bin must be vectors')
end
I = length(theta_bin);
J = length(mag_bin);

% warnings
imag = isnan(mag);
itheta = isnan(theta);
if sum(imag) > 0
    display([num2str(sum(imag)/length(imag)*100) ' % of magnitudes are NaN'])
end
if sum(itheta) > 0
    display([num2str(sum(itheta)/length(itheta)*100) ' % of theta values are NaN'])
end

% frequency matrix
% hist = hist2d_JN(mag,theta,mag_bin,theta_bin(1:end-1));
if isempty(freq)
    hist = hist2d_JN(mag,theta,mag_bin,theta_bin);
    freq = hist / sum(hist(:))*100;
end

%if size(freq,2)~=J, error('size(freq,2)~=length(mag_bin)'), end
%if size(freq,1)~=I, error('size(freq,1)~=length(theta_bin)'), end
if ~issorted(mag_bin)&&~issorted(theta_bin)
    error('mag_bin and theta_bin must be monotonic')
end

% make parent figure current as with axes
f = get(ax,'Parent');
figure(f);
set(f,'CurrentAxes',ax)

% convert to cartesian with the difference between new bins same as
% original ie. not (90    45     0   315   270   225   180   135    90)
theta_bin = 270 - theta_bin;
theta_bin = theta_bin - 180;
cumfreq = [zeros(I,1),cumsum(freq,2)];   % cumsum so pieces of pie build outwards with the different velocity bins

% Create polar axis to work with
polar(ax,[],[],rmax,rticks,false);


axis fill
frac = 0.6;
dy_new = 2*rmax / frac;
t_buff = dy_new * 0.1;
ylim(2) = rmax + t_buff;
ylim(1) = ylim(2) - dy_new;
set(ax,'YLim',ylim);
axis equal

uns = get(ax,'units');
% set(ax,'units','normalized');
pos = get(ax,'position');
xlim = get(ax,'XLim');
% colors
c = rgb(length(mag_bin));
colormap(c);
set(ax,'CLim',[mag_bin(1) mag_bin(end) + mode(diff(mag_bin))])

% Loop through directions
for i = 1 : I-1                % 'I-1' because the bins starts with 0 (0-45) and ends with 360 where the last bin is (theta_bin(end-1)-360)
    % Loop through magnitudes
    for j = 1 : J
        thet = [linspace(theta_bin(i),theta_bin(i+1),10),...
            linspace(theta_bin(i+1),theta_bin(i),10),...
            theta_bin(i)];
        r = [repmat(cumfreq(i,j),1,10),...      % inner arc of shape
            repmat(cumfreq(i,j+1),1,10),...     % outer arc of shape
            cumfreq(i,j)];                      % this is why there is a column of zeros in cumfreq so a 'piece of pie' shape is plotted
        x = r.*cosd(thet);
        y = r.*sind(thet);
        patch(x,y,c(j,:))
    end
end

% Create frequency legend
frac = 0.30;

h_frq = axes;
set(h_frq,'units',uns);
polar(h_frq,[],[],rmax,rticks,true);
axis fill;
xp = pos(1) + (pos(3)/3-frac*pos(4))/2;
yp = pos(2);
width = frac * pos(4);
height = frac * pos(4);

set(h_frq,'position',[xp yp width height])
axis equal


% Create magnitude legend symetrical with freq legend
frac = 0.05; % width of colorbar as fraction of parent axes width

h_cbar = colorbar('peer',ax);
set(h_cbar,'units',uns);
set(h_cbar,'Location','EastOutside')
xp = pos(1) + pos(3)/2.1;
yp = pos(2) + 0.15 * pos(4) - frac * pos(4) /2;
width1 = rmax / diff(xlim) * pos(3);
% width1=pos(3);
% width2 = 0.8 * pos(4);
% width = max(width1,width2);
width = width1;
height = frac * pos(4);
pos = get(h_cbar,'Position');

set(h_cbar,'Position',[pos(1) pos(2)+pos(4)./4 pos(3)./2 pos(4)./2])
% set(h_cbar,'XAxisLocation','bottom')
set(h_cbar,'XTick',mag_bin)



switch lower(type)
    case 'wave'
        t_cbar = 'Hsig (m)';
    case {'wind','current'}
        t_cbar = 'Velocity (m/s)';
end
title(h_cbar,t_cbar,'FontWeight','bold')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hpol,rmax] = polar(cax,theta,rho,rmax,rticks,leg)

if ~isequal(size(theta),size(rho))
    %error('THETA and RHO must be the same size.');
    error('MATLAB:polar:InvalidInput', 'THETA and RHO must be the same size.');
end

line_style = 'auto';

% get hold state
cax = newplot(cax);

next = lower(get(cax,'NextPlot'));
hold_state = ishold(cax);

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   get(cax, 'FontSize'), ...
    'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
    'DefaultTextUnits','data')

% only do grids if hold is off
if ~hold_state
    
    % make a radial grid
    hold(cax,'on');
    rmin = 0;
    %     maxrho = max(abs(rho(:)));
    %     hhh=line([-maxrho -maxrho maxrho maxrho],[-maxrho maxrho maxrho -maxrho],'parent',cax);
    %     set(cax,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
    %     v = [get(cax,'xlim') get(cax,'ylim')];
    %     ticks = sum(get(cax,'ytick')>=0);
    %     delete(hhh);
    % % check radial limits and ticks
    %     rmin = 0; rmax = v(4); rticks = max(ticks-1,2);
    %     if rticks > 5   % see if we can reduce the number
    %         if rem(rticks,2) == 0
    %             rticks = rticks/2;
    %         elseif rem(rticks,3) == 0
    %             rticks = rticks/3;
    %         end
    %     end
    
    % define a circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);
    % now really force points on x/y axes to lie on them exactly
    inds = 1:(length(th)-1)/4:length(th);
    xunit(inds(2:2:4)) = zeros(2,1);
    yunit(inds(1:2:5)) = zeros(3,1);
    % plot background if necessary
    if ~ischar(get(cax,'color')),
        patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
            'edgecolor',tc,'facecolor',get(cax,'color'),...
            'handlevisibility','off','parent',cax);
    end
    
    % draw radial circles
    c89 = cos(89*pi/180);
    s89 = sin(89*pi/180);
    rinc = (rmax-rmin)/rticks;
    for i=(rmin+rinc):rinc:rmax
        hhh = line(xunit*i,yunit*i,'linestyle',ls,'color',tc,'linewidth',1,...
            'handlevisibility','off','parent',cax);
        set(cax,'Layer','Bottom')
        if leg
            % CHANGED BY JN SO FONTSIZE is DEFAULT
            %         text((i+rinc/20)*c89,(i+rinc/20)*s89, ...
            %             ['  ' num2str(i) '%'],'verticalalignment','bottom',...
            %             'handlevisibility','off','parent',cax,...
            %             'fontsize',16,'fontweight','bold')
            text((i-rinc/2)*c89,(i-rinc/2)*s89, ...
            ['  ' num2str(i) '%'],'verticalalignment','middle',...
            'horizontalalignment','center',...
            'handlevisibility','off','parent',cax,...
            'fontweight','bold','fontsize',8)
        end
    end
    set(hhh,'linestyle','-') % Make outer circle solid
    
    % plot spokes
    if ~leg
        %     th = (1:4)*2*pi/8;
        th = (1:8)*2*pi/16;
        cst = cos(th); snt = sin(th);
        cs = [-cst; cst];
        sn = [-snt; snt];
        line(rmax*cs,rmax*sn,'linestyle',ls,'color',tc,'linewidth',1,...
            'handlevisibility','off','parent',cax)
        %     txt = {'SW','S','SE','E','NE','N','NW','W'};
        txt = {'WSW','SW','SSW','S','SSE','SE','ESE','E','ENE','NE','NNE','N','NNW','NW','WNW','W'};
        
        % annotate spokes in degrees
        rt = 1.2*rmax;
        for i = 1:length(th)
            
            % CHANGED BY JN SO FONTSIZE is DEFAULT
            %         text(rt*cst(i),rt*snt(i),txt{i+length(th)},...
            %              'horizontalalignment','center',...
            %              'handlevisibility','off','parent',cax,...
            %              'fontsize',16,'fontweight','bold');
            
            text(rt*cst(i),rt*snt(i),txt{i+length(th)},...
                'horizontalalignment','center',...
                'handlevisibility','off','parent',cax,...
                'fontweight','bold');
            
            %         if i == length(th)
            %             loc = int2str(0)
            %         else
            %             loc = int2str(180+i*30)
            %         end
            %             text(-rt*cst(i),-rt*snt(i),txt{i},'horizontalalignment','center',...
            %                 'handlevisibility','off','parent',cax,...
            %                 'fontsize',16,'fontweight','bold')
            text(-rt*cst(i),-rt*snt(i),txt{i},'horizontalalignment','center',...
                'handlevisibility','off','parent',cax,...
                'fontweight','bold')
        end
    end
    
    % set view to 2-D
    view(cax,2);
    % set axis limits
    %     axis(cax,rmax*[-1 1 -1.15 1.15]);
end

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits',fUnits );

% transform data to Cartesian coordinates.
xx = rho.*cos(theta);
yy = rho.*sin(theta);

% plot data on top of grid
if strcmp(line_style,'auto')
    q = plot(xx,yy,'parent',cax);
else
    q = plot(xx,yy,line_style,'parent',cax);
end

if nargout > 0
    hpol = q;
end

if ~hold_state
    set(cax,'dataaspectratio',[1 1 1]), axis(cax,'off'); set(cax,'NextPlot',next);
end
set(get(cax,'xlabel'),'visible','on')
set(get(cax,'ylabel'),'visible','on')

if ~isempty(q) && ~isdeployed
    makemcode('RegisterHandle',cax,'IgnoreHandle',q,'FunctionName','polar');
end

% function Hout = hist2d(Dx,Dy,Xbin,Ybin)
% sum(Hout(:)) will always equal length(Dx) as value below the minimum bin
% values are grouped into the minimum bin and those above the highest bin
% value grouped into the highest bin.
%
% Calculates and returns the 2 Dimensional Histogram of D.
%
%eg v_bin = [0 2 4 6 8 10]   = 6 bins
%
function H = hist2d_JN(Dx,Dy,Xbin,Ybin)

Xn = length(Xbin);
Yn = length(Ybin);

N = length(Dx);
if length(Dy)~=N, error('length(Dx)~=length(Dy)'), end

H = zeros(Yn,Xn);

for i = 1:N
    x = find(Xbin <= Dx(i),1,'last');
    y = find(Ybin <= Dy(i),1,'last');
    if isempty(x), x=1; end
    if isempty(y), y=1; end
    H(y,x) = H(y,x)+1;
end

if sum(H(:)) ~= length(Dx), error('not all records included in histogram'), end