% /////// myline ///////
% Draw the line with a 'o' marker indicating the starting point onto the
% current axes. If no points are specified the user is prompted to
% draw a line onto an axes using mouse clicks.
%
% To start drawing left click, to end drawing the line right click.
%
% JN August 2014

function h = myline(varargin)

h = hggroup;

if length(varargin) == 2
    x = varargin{1};
    y = varargin{2};
    linObj = line('XData',x,'YData',y,'Color','k','parent',h);
    mrkObj = line('XData',x(1),'YData',y(1),'Marker','o','color','w','MarkerEdgeColor','k','LineStyle','none','parent',h);
    set(linObj,'HitTest','off')
    set(mrkObj,'HitTest','off')
    set(linObj,'parent',h)
    set(mrkObj,'parent',h)
    return
end


% turn off toolbar modes
plotedit off
zoom off
pan off
rotate3d off
datacursormode off
brush off

display('left click to define myline, right click to end it')
inloop = true;
set(gcf,'WindowButtonDownFcn',@wbdcb)
linObj = line('XData',[],'YData',[]);
mrkObj = line('XData',[],'YData',[]);
drawnow
    function wbdcb(src,~)
        if strcmp(get(src,'SelectionType'),'normal')
            set(src,'pointer','circle')
            cp = get(gca,'CurrentPoint');
            xinit = cp(1,1);
            yinit = cp(1,2);
            xdata = get(linObj,'XData');
            ydata = get(linObj,'YData');
            if isempty(xdata)
                mrkObj = line('XData',xinit,'YData',yinit,'Marker','o','color','w','MarkerEdgeColor','k','LineStyle','none');
                linObj = line('XData',xinit,'YData',yinit,'Color','k');
            end
            xdata = [xdata xinit];
            ydata = [ydata yinit];
            set(linObj,'XData',xdata,'YData',ydata)
            set(src,'WindowButtonMotionFcn',@wbmcb)
            set(src,'WindowButtonUpFcn',@wbucb)
        end
        
        function wbmcb(~,~)
            cp = get(gca,'CurrentPoint');
            xtmp = [xdata,cp(1,1)];
            ytmp = [ydata,cp(1,2)];
            set(linObj,'XData',xtmp,'YData',ytmp)
            drawnow
        end
        
        function wbucb(src,~)
            if strcmp(get(src,'SelectionType'),'alt')
                set(src,'Pointer','arrow')
                set(src,'WindowButtonMotionFcn','')
                set(src,'WindowButtonUpFcn','')
                set(src,'WindowButtonDownFcn','')
                inloop = false;
            else
                return
            end
        end
    end

while inloop
    drawnow
    % cycle around
end
set(linObj,'HitTest','off')
set(mrkObj,'HitTest','off')
set(linObj,'parent',h)
set(mrkObj,'parent',h)
end



