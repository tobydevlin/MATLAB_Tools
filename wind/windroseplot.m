% function windroseplot(mag_bin,theta_bin,freq,rmax,rticks)

function windroseplot(mag_bin,theta_bin,freq,rmax,rticks)

if ~isvector(theta_bin)&&~isvector(mag_bin)
    error('theta_bin and mag_bin must be vectors')
end
I = length(theta_bin)-1;
J = length(mag_bin);
if size(freq,2)~=J, error('size(freq,2)~=length(mag_bin)'), end
if size(freq,1)~=I, error('size(freq,1)~=length(theta_bin)'), end
if ~issorted(mag_bin)&&~issorted(theta_bin)
    error('mag_bin and theta_bin must be monotonic')
end

theta_bin = 270 - theta_bin;
theta_bin = theta_bin/180*pi;

% Create polar axis to work with
maxval = max(sum(freq,2));
figure
cax = subplot(15,1,1:13);
h = polar(cax,[0 0],[0 maxval],rmax,rticks,false);
delete(h);

cumfreq = [zeros(I,1),cumsum(freq,2)];

c = mag_bin(1:end);

% Loop through directions
for i = 1 : I
    % Loop through magnitudes
    for j = 1 : J
        f = 0.2;
        dthet = theta_bin(i+1)-theta_bin(i);
        thet = [linspace(theta_bin(i)+f*dthet,theta_bin(i+1)-f*dthet,10),...
            linspace(theta_bin(i+1)-f*dthet,theta_bin(i)+f*dthet,10),...
            theta_bin(i)+f*dthet]-pi;
        r = [repmat(cumfreq(i,j),1,10),...
            repmat(cumfreq(i,j+1),1,10),...
            cumfreq(i,j)];
        x = r.*cos(thet);
        y = r.*sin(thet);
        patch(x,y,c(j))
    end
end

set(cax,'handlevisibility','off')

% Create frequency legend
cax = subplot(3,2,5);
set(cax,'xlim',[-1 1],'ylim',[0 1],'color',get(gcf,'color'))
axis(cax,'off')
h = polar(cax,[0 0],[0 maxval],rmax,rticks,true);
delete(h);

% Create magnitude legend
cax = subplot(4,10,38);
set(cax,'xlim',[0 1],'ylim',[-1 1],'color',get(gcf,'color'))
axis(cax,'off')
y = linspace(-0.7,0.7,J+1);
x = zeros(1,4);
x(3:4) = 1;
for j = 1 : J
    patch(x,[y(j:j+1),y(j+1:-1:j)],c(j))
    text(-0.3,y(j),num2str(mag_bin(j)),...
        'fontsize',14,'fontweight','bold','horizontalalignment','center')
end
text(0.5,-0.9,'Magnitude (km/h)','fontsize',14,...
    'fontweight','bold','horizontalalignment','center')

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
        text((i+rinc/20)*c89,(i+rinc/20)*s89, ...
            ['  ' num2str(i) '%'],'verticalalignment','bottom',...
            'handlevisibility','off','parent',cax,...
            'fontsize',16,'fontweight','bold')
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
    rt = 1.1*rmax;
    for i = 1:length(th)
        text(rt*cst(i),rt*snt(i),txt{i+length(th)},...
             'horizontalalignment','center',...
             'handlevisibility','off','parent',cax,...
             'fontsize',16,'fontweight','bold');
%         if i == length(th)
%             loc = int2str(0)
%         else
%             loc = int2str(180+i*30)
%         end
        text(-rt*cst(i),-rt*snt(i),txt{i},'horizontalalignment','center',...
             'handlevisibility','off','parent',cax,...
             'fontsize',16,'fontweight','bold')
    end
end

% set view to 2-D
    view(cax,2);
% set axis limits
    axis(cax,rmax*[-1 1 -1.15 1.15]);
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