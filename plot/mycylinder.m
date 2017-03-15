% /////// mycylinder ///////
% h = mycylinder(ax,bot,top,inc,rmax,rticks)
% creates cylinder object
%
% outputs
%   h = handle to cyclinder object
%
% inputs
%   ax (handle)         = handle to axes to plot into
%   bot (numeric)       = bottom of cylinder
%   top (numeric)       = top of cylinder
%   inc (numeric +ve)   = size of intervals depicted in centre column
%   rmax (numeric)      = radial limit of cylinder
%   rticks (int)        = # radial ticks
%
% optional inputs as descriptor / value pairs
% ref_l / reference (numeric) = elevation to draw small disk on centre column and centre column markings around
% ref_r / fraction (numeric [0 1]) = fraction of rmax to plat reference disk
% ref_c / color vector eg. [0 0 1] or recognised matlab color string = color of reference disk
%
%
% Jesper Nielsen March 2012

function h = mycylinder(ax,bot,top,inc,rmax,rticks,varargin)

h = 1;  % TEMPORARY FLUFF

% defaults
refme = false;
ref_l = 0; % refernce level
ref_r = 0.05; % percentage of rmax
ref_c = [0.9 0.9 1];

% variable arguments
if mod(length(varargin),2) > 0
    error('Expecting variable arguments as descriptor/value pairs')
end

for i = 1 : 2 : length(varargin)
    varargtyp{i} = varargin{i};
    varargval{i} = varargin{i+1};
    switch lower(varargtyp{i})
        case {'ref_l'}
            ref_l = varargval{i};
            if ischar(ref_l)
                switch ref_l
                    case 'none'
                        refme = false;
                    otherwise
                        error('only allowable string input is "none"')
                end
            end
            refme = true;
        case {'ref_r'}
            ref_r = varargval{i};
        case {'ref_c'}
            ref_c = varargval{i};
        otherwise
            error('unexpected variable argument type')
    end
end

% make inputted axes current
f = get(ax,'parent');
set(f,'CurrentAxes',ax)
set(ax,'NextPlot','add')

% h = hggroup;
% set(h,'parent',ax)

lw = get(ax,'LineWidth');
lw = lw * 3;

% centre column
col(1,:) = [0 0 0];
col(2,:) = [0.8 0.8 0.8];

% -- work up from ref_l
nc = ceil((top - ref_l) / inc);
for aa = 1:nc
    x = [0 0];
    y = [0 0];
    z = [0 1];
    
    if aa < nc
        z = z * inc;
    else
        z = z * (top - (aa-1) * inc);
    end
    z = z + inc * (aa-1) + ref_l;
    if mod(aa,2) == 0
        line(x,y,z,'Color',col(1,:),'LineWidth',lw,'Parent',ax)
    else
        line(x,y,z,'Color',col(2,:),'LineWidth',lw,'Parent',ax)
    end
end

% -- work down from ref_l
nc = ceil(ref_l - bot / inc);
for aa = 1:nc
    x = [0 0];
    y = [0 0];
    z = [0 -1];
    
    if aa < nc
        z = z * inc;
    else
        z = z * (-(aa-1) * inc - bot);
    end
    z = z - inc * (aa-1) + ref_l;
    if mod(aa,2) == 1
        line(x,y,z,'Color',col(1,:),'LineWidth',lw,'Parent',ax)
    else
        line(x,y,z,'Color',col(2,:),'LineWidth',lw,'Parent',ax)
    end
end

% reference level
if refme
    p = linspace(0,2*pi,10);
    np = length(p);
    r = ref_r * rmax;
    x = r*cos(p);
    y = r*sin(p);
    z = ref_l * ones(np,1);
    patch(x,y,z,'FaceColor',ref_c,'EdgeColor','none','Parent',ax)
end

% polar axes at base
% polar(h,rmax,rticks,bot)
polar(ax,rmax,rticks,bot)

% /////// nested functions ///////

function [hpol,rmax] = polar(mum,rmax,rticks,z)
% mum is object to group entire cylinder into

cax = mum;
% cax = get(mum,'parent');
% type = get(cax,'Type');
% switch type
%     case 'axes'
%         % good
%     otherwise
%         error('parent of "mum" must be an axes')
% end

tc = 'k';
% tc = get(cax,'xcolor');         % line color
ls = get(cax,'gridlinestyle');  % line style

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

zunit = z * ones(length(xunit),1);

% draw radial circles
% c89 = cos(89*pi/180);
% s89 = sin(89*pi/180);
c89 = cos(68.5*pi/180);  % offset text so does not interfere with "N"
s89 = sin(68.5*pi/180);
rinc = (rmax-rmin)/rticks;
for i=(rmin+rinc):rinc:rmax
    hhh = line(xunit*i,yunit*i,zunit,'linestyle',ls,'color',tc,'linewidth',1,'parent',mum);
    %     set(cax,'Layer','Bottom')
    text((i+rinc/20)*c89,(i+rinc/20)*s89,z,...
        num2str(i),'verticalalignment','bottom','parent',mum,...
        'fontweight','bold','Color',[0.8 0.8 0.8])
end
set(hhh,'linestyle','-') % Make outer circle solid

% plot spokes
th = (1:4)*2*pi/8;
%     th = (1:8)*2*pi/16;
cst = cos(th);
snt = sin(th);
cs = [-cst; cst];
sn = [-snt; snt];
zs = z*ones(size(cs));
line(rmax*cs,rmax*sn,zs,'linestyle',ls,'color',tc,'linewidth',1,'parent',mum)
txt = {'SW','S','SE','E','NE','N','NW','W'};
%     txt = {'WSW','SW','SSW','S','SSE','SE','ESE','E','ENE','NE','NNE','N','NNW','NW','WNW','W'};

% annotate spokes in degrees
rt = 1.1*rmax;
for i = 1:length(th)
    
    text(rt*cst(i),rt*snt(i),z,txt{i+length(th)},...
        'horizontalalignment','center','parent',mum,...
        'fontweight','bold');
    
    text(-rt*cst(i),-rt*snt(i),z,txt{i},'horizontalalignment','center','parent',mum,...
        'fontweight','bold')
end
