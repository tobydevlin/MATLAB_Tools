% /////// myscalearrow ///////
% plots an arrow (to scale) with accompanying text
% be sure to use the same shape and scale as those used in your figures
%
% h_sarrow = myscalearrow(ax,lim,multi,shape)
%
% inputs
%    ax = axes to plot scale arrows into
%    lim = magnitude (m/s) or scale arrow
%    multi = scaling of arrows
%    shape_head = [a b 1] [width head, legth head, logical (true if frac of arrow length, false if abs values)]
%    shape_tail = [a 1]  [width stem, logical (true if frac of arrow length, false if abs value)];
%
% variable inputs as discriptor / value pairs
%    FaceColor = color of arrows
%    EdgeColor = edge color of arrows
%    start  = fraction of ax's limits defining starting point of scale arrow [x y]
%    text = additional text which appears after m/s
%
% JN Feburary 2012


function [h_sarrow,h_sarrow_text] = myscalearrow(ax,lim,multi,shape_head,shape_tail,varargin)

% preliminaries
xlim = get(ax,'XLim');
ylim = get(ax,'YLim');
f = get(ax,'Parent');
set(f,'CurrentAxes',ax)

% defaults
xp = xlim(1) + 0.025 * diff(xlim);
yp = ylim(1) + 0.95 * diff(ylim);
txt = '';
backcolor = 'w';

face_col = 'k';
edge_col = 'k';
linewidth = 0.5;

% variable arguments
if mod(length(varargin),2) > 0
    error('Expecting variable arguments as descriptor/value pairs')
end

for i = 1 : 2 : length(varargin)
    varargtyp{i} = varargin{i};
    varargval{i} = varargin{i+1};
    switch lower(varargtyp{i})
        case 'facecolor'
            face_col = varargval{i};
        case 'edgecolor'
            edge_col = varargval{i};
        case 'linewidth';
            linewidth = varargval{i};
        case 'text'
            txt = varargval{i};
        case 'start'
            xp = xlim(1) + varargval{i}(1) * diff(xlim);
            yp = ylim(1) + varargval{i}(2) * diff(ylim);
        case 'backcolor'
            backcolor = varargval{i};
        otherwise
            error('unexpected variable argument type')
    end
end

% scaleing of vectors
[fx,fy] = getscale(ax);
scale = sqrt(fx.^2 + fy.^2) * multi;

% create veticees to patch
[vx,vy] = arrow2(lim,0,xp,yp,'scale',scale,'shape_head',shape_head,'shape_tail',shape_tail);

% patch up your scale arrow
h_sarrow = patch(vx,vy,'w','FaceColor',face_col,'EdgeColor',edge_col,'linewidth',linewidth);

% label your scale arrow
h_sarrow_text = text(xp + 1.1*lim*scale,yp,[num2str(lim) ' m/s  ' txt],'color','k','fontweight','bold','BackgroundColor',backcolor);
