% /////// myscalebar ///////
% h = myscalebar(ax,lim,varargin)
% generates a scale bar onto the specified axis
% marks the start and end with a number and halfway with change from white to black
% returns a warning if the scaling of the "parent" axis is not equal;
% The axes which features in the input is not actually the parent. The
% scale bar is infact an axis (making use of the labeling) whos real parent is
% the figure parenting ax
%
% inputs
%    ax = axes' to put scale bar/s into
%    lim = length of scale bar in units as specified (default = metres)
% variable inputs
%   'position'/ vector (normalized to "parent" axis) [xp yp width hight] or string 'topright' 'topleft'
%   'bottomleft' 'bottomright' (default). the width (position(3))
%   is reduntant as the input 'lim' is used to set the width
%   'units'/ 'metres' (default), 'kilometres' (not yet supported),
%   'spherical' / true or false (default) whether plotting has been performed in spherical coordinates
%
% outputs
%    h = handel to scale bar axes
%
% tips
% if your scale bar is getting covered by other axis for instance within a loop try
%     kids = get(f,'children');
%     i = kids == h_sbar;
%     kids = cat(1,kids(i),kids(~i));
%     set(f,'children',kids)
%
% useful call back so scalebar automatically resizes
%     h_zoom = zoom(f);
%     set(h_zoom,'ActionPostCallback',('delete(h_sbar);h_sbar = myscalebar(ax(1),250,''spherical'',true);'));
%
% JN September 2011

function h = myscalebar(ax,lim,varargin)

% defaults
height = 0.025;      % 10% of "parent" axes' height
pos = 'bottomright';
isvec = false;
units = 'metres';
x_buff = 0.025;
y_buff = 4.5 * height;
spherical = false;

% variable arguments
if mod(length(varargin),2) > 0
    error('Expecting variable arguments as descriptor/value pairs')
end

for i = 1 : 2 : length(varargin)
    varargtyp{i} = varargin{i};
    varargval{i} = varargin{i+1};
    switch lower(varargtyp{i})
        case 'position'
            pos = varargval{i};
            if isnumeric(pos)
                isvec = true;
            end
        case 'units'
            error('varargin "units" not yet supported')
        case 'spherical'
            spherical = varargval{i};
            if ~islogical(spherical)
                error('expecting logical input for spherical input')
            end
        otherwise
            error('unexpected variable argument type')
    end
end

% loop through axes planting a scalebar in each
na = length(ax);
for aa = 1:na
    
    % make parent figure current
    f = get(ax(aa),'Parent');
    figure(f);
    
    % pos of "parent" axes within figure
    f_pos = get(f,'position');
    yx_ratio = f_pos(4) / f_pos(3);
    
    ax_pos = get(ax(aa),'position');
    xlim = get(ax(aa),'XLim');
    ylim = get(ax(aa),'YLim');
    
    if spherical
        [xlim(1) ylim(1) gz(1,:)] = ll2utm(xlim(1),ylim(1));
        [xlim(2) ylim(2) gz(2,:)] = ll2utm(xlim(2),ylim(2));
        if diff(gz(:,1)) ~= 0
            error('axes windows extends beyond a single UTM grid zone')
        end
    end
    
    dx = diff(xlim);
    dy = diff(ylim);
    fx = ax_pos(3) / dx;
    fy = ax_pos(4) * yx_ratio / dy;
    if abs(fx - fy)/fx > 0.01
        display('WARNING inconsistent scaling between x & y axis')
    end
    
    % replace  pos(3) with 'lim'
    if isvec
        pos(3) = lim / dx;
    end
    
    % convert any string pos to a normalized pos vector within "parent" axis
    if ~isvec
        width = lim / dx;
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
    
    % update scale bar pos vector to reference the parent figure
    pos(1) = ax_pos(1) + ax_pos(3) * pos(1);
    pos(2) = ax_pos(2) + ax_pos(4) * pos(2);
    pos(3) = pos(3) * ax_pos(3);
    pos(4) = pos(4) * ax_pos(4);
    
    % generate scale bar
    h(aa) = axes;
    set(h(aa),'position',pos)
    set(h(aa),'Xtick',[0 1]);
    set(h(aa),'XTickLabel',{'0';num2str(lim)})
    set(h(aa),'YTick',[]);
    title(h(aa),units)
    patch([0;0.5;0.5;0],[1;1;0;0],'k')
    patch([0.5;1;1;0.5],[1;1;0;0],'w')
    
%     % tag him
%     set(h(aa),'Tag','scalebar')
%     
%     % user data
%     UD.lim = lim;
%     UD.speriscal = spherical;
%     set(h(aa),'UserData',UD)
end



