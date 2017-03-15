% /////// mymouse ///////
%
% function mymouse
%
% When mymouse is called in the command window the position of the cursor (mouse) is
% plotted as markers (brother mice) on all axes open in the matlab session.
% Call mymouse a 2nd time in the command window to clear the mice.
% In the figure task bar the select, zoom, pan & rotate icons must be
% unselected for the mice to update.
%
% TIP:
% If you create a new figure while mymouse is active, new brother mice will
% appear in the new axes (children of the new figure). If you want to place
% the mouse over any of the new axes and have the brother mice update their
% positions you will need to cancle the mymouse session before starting a
% new one.
%
% Jesper Nielsen, February 2014

function mymouse

f = findobj(0,'Type','figure');

% are we setting up mymouse or ending it
nf = length(f);
for aa = 1:nf
    fun = get(f(aa),'WindowButtonMotionFcn');
    if isempty(fun)
        set(f(aa),'WindowButtonMotionFcn',@motion)
    else
        set(f(aa),'WindowButtonMotionFcn',{})
        h_tmp = findobj(f(aa),'Tag','mice');
        delete(h_tmp);
    end
end

% MOTION executes whenever the mouse moves within a figure window.
% The select, zoom, pan & rotate icons must be unselected in the figures task bar.
function motion(f,~)

% are any of the taskbar icons selected (this bit could be more robust)
fun = get(f,'WindowButtonDownFcn');
if ~isempty(fun)
    return
end

% make the axes which the cursor is over the current axes
set(f,'Units','normalized')
tmp = get(f,'CurrentPoint');
xp = tmp(1,1);
yp = tmp(1,2);
ax = findobj(f,'Type','axes','-not','Tag','legend','-not','Tag','Colorbar');
ax_in = [];
for aa = 1:length(ax);
    set(ax(aa),'Units','normalized');
    pos = get(ax(aa),'position');
    x(1) = pos(1);
    x(2) =  pos(1)+pos(3);
    x(3) = x(2);
    x(4) = x(1);
    y(1) = pos(2);
    y(2) = y(1);
    y(3) = pos(2) + pos(4);
    y(4) = y(3);
    if inpolygon(xp,yp,x,y)
        ax_in = ax(aa);
    end
end

if isempty(ax_in)
    return
end

% get coordinates of mouse within the axes its over
tmp = get(ax_in, 'CurrentPoint');
x = tmp(1,1);
y = tmp(1,2);

% brother mice
h_mice = findobj(0,'Tag','mice');
% -- draw the brother mice
if isempty(h_mice)
    ax_all = findobj(0,'Type','axes','-not','Tag','legend','-not','Tag','Colorbar');
    for aa = 1:length(ax_all)
        h_tmp = plot(ax_all(aa),x,y,'Marker','o','LineStyle','none','MarkerEdgeColor','k');
        set(h_tmp,'XDataSource','x','YDataSource','y','Tag','mice')
    end
else
    % -- update positions of brother mice
    refreshdata(h_mice(:),'caller')
end
