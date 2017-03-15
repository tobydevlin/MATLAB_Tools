% Function to plot a histogram of the probability distribution function and
% a line plot of the cumulative distribution function of a dataset.

function [pdf, cdf] = histplot(edges,y,varargin)

J = length(edges);
K = length(varargin) + 1;

pdf = zeros(J-1,K);
cdf = zeros(J,K);

for k = 1 : K
    if k == 1, yk = y;
    else yk = varargin{k-1};
    end
    if min(yk) < min(edges)
        disp('Warning y value out of edges span')
    elseif max(yk) > max(edges)
        disp('Warning y value out of edges span')
    end
    ycnt = histc(yk,edges);
    ycnt = ycnt(1:end-1,:);
    pdf(:,k) = ycnt / sum(ycnt) * 100;
    cdf(:,k) = cumsum([0;pdf(:,k)]);
end

cmap = colormap(jet);
figure
h = bar(edges,[pdf;zeros(1,size(pdf,2))],'histc');
set(h, 'FaceAlpha', 0.7)
ax1 = gca;
hold
xlimits = [edges(1) edges(end)];
set(ax1,'XLim', xlimits);
ylimits = get(ax1,'YLim');
yinc = (ylimits(2)-ylimits(1))/10;
set(ax1,'YTick',ylimits(1):yinc:ylimits(2),'Fontsize',12)
pbaspect(ax1, [1 1 1])
xlabel('X')
ylabel('% Events')
grid
for k = 1 : K
    label{k} = ['series',num2str(k)];
end
legend(label,'Fontsize',12)

ax2 = axes('Position',get(ax1,'Position'),...
    'XLim',xlimits,'XTick',[],...
    'XAxisLocation','bottom','YAxisLocation','right',...
    'Color','none','XColor','k','YColor','k','Fontsize',12);
pbaspect(ax2, [1 1 1])

for k = 1 : K
    m = size(cmap,1);
    if K > 1
        index = fix((k-1)/(K-1)*m)+1;
        if index < 1, index = 1; end
        if index > m, index = m; end
    else
        index = 1;
    end
    color = cmap(index,:);
    line(edges,cdf(:,k),'Parent',ax2,'Color',color,'Linewidth',2);
end
set(ax2,'YLim',[0 100],'YTick',[0 10 20 30 40 50 60 70 80 90 100])
ylabel('Cumulative %')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The following MATLAB functions were appended because minor changes were
%%% required.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hh = bar(varargin)
%BAR Bar graph.
%   BAR(X,Y) draws the columns of the M-by-N matrix Y as M groups of N
%   vertical bars.  The vector X must not have duplicate values.
%
%   BAR(Y) uses the default value of X=1:M.  For vector inputs, BAR(X,Y)
%   or BAR(Y) draws LENGTH(Y) bars.  The colors are set by the colormap.
%
%   BAR(X,Y,WIDTH) or BAR(Y,WIDTH) specifies the width of the bars. Values
%   of WIDTH > 1, produce overlapped bars.  The default value is WIDTH=0.8
%
%   BAR(...,'grouped') produces the default vertical grouped bar chart.
%   BAR(...,'stacked') produces a vertical stacked bar chart.
%   BAR(...,LINESPEC) uses the line color specified (one of 'rgbymckw').
%
%   BAR(AX,...) plots into AX instead of GCA.
%  
%   H = BAR(...) returns a vector of barseries handles in H.
%
%   Use SHADING FACETED to put edges on the bars.  Use SHADING FLAT to
%   turn them off.
%
%   Examples: subplot(3,1,1), bar(rand(10,5),'stacked'), colormap(cool)
%             subplot(3,1,2), bar(0:.25:1,rand(5),1)
%             subplot(3,1,3), bar(rand(2,3),.75,'grouped')
%
%   See also HIST, PLOT, BARH, BAR3, BAR3H.

%   Backwards compatibility
%   BAR('v6',...) creates patch objects instead of barseries
%   objects for compatibility with MATLAB 6.5 and earlier.

%   C.B Moler 2-06-86
%   Modified 24-Dec-88, 2-Jan-92 LS.
%   Modified 8-5-91, 9-22-94 by cmt; 8-9-95 WSun.
%   Copyright 1984-2005 The MathWorks, Inc. 
%   $Revision: 5.34.6.16 $  $Date: 2005/09/29 16:33:30 $

[v6,args] = usev6plotapi(varargin{:});
if v6 || ((length(args) > 1) && ...
          isa(args{end},'char') && ...
          (length(args{end}) > 3) && ...
          (strcmp(args{end}(1:4),'hist')))
  h = Lbarv6(args{:});
else
  [cax,args,nargs] = axescheck(args{:});
  error(nargchk(1,inf,nargs,'struct'));
  [args,pvpairs,msg] = parseargs(args);
  if ~isempty(msg), error(msg); end %#ok
  nargs = length(args);

  [msg,x,y] = xychk(args{1:nargs},'plot');
  if ~isempty(msg), error(msg); end %#ok
  hasXData = nargs ~= 1;
  if hasXData && length(x)>1
    sortedx = sort(x);
    if any(sortedx(2:end)-sortedx(1:end-1) == 0)
      error(id('DuplicateXValue'),...
            'XData cannot contain duplicate values.');
    end
  end
  if min(size(x))==1, x = x(:); end
  if min(size(y))==1, y = y(:); end
  n = size(y,2);

  % handle vectorized data sources and display names
  extrapairs = cell(n,0);
  if ~isempty(pvpairs) && (n > 1)
    [extrapairs, pvpairs] = vectorizepvpairs(pvpairs,n,...
                                            {'XDataSource','YDataSource','DisplayName'});
  end
  
  % Create plot
  cax = newplot(cax);

  h = []; 
  xdata = {};
  pkg = findpackage('specgraph');
  findclass(pkg,'barseries');
  listeners = getappdata(0,'SpecgraphBarListeners');
  seriesListeners = getappdata(0,'Graph2dSeriesListeners');
  % 2 is the index for the peer listener
  set(listeners(2),'enable','off');
  set(seriesListeners(end),'enable','off');
  for k=1:n
    % extract data from vectorizing over columns
    if hasXData
      xdata = {'XData', datachk(x(:,k))};
    end
    h = [h specgraph.barseries('YData',datachk(y(:,k)), ...
                               xdata{:}, pvpairs{:},...
			      extrapairs{k,:}, 'Parent', cax)];
  end
  set(h,'BarPeers',h);
  if ~ishold(cax)
      % Turn off edges when they start to overwhelm the colors
      if numel(y) > 150, 
          set(h,{'edgecolor'},get(h,{'facecolor'}));
      end
  end
  if n > 1
    set(h(2:end),'RefreshMode','auto');
  end
  set(listeners(2),'enable','on');
  set(seriesListeners(end),'enable','on');
  set(h(1),'RefreshMode','auto');
  plotdoneevent(cax,h);
  h = double(h);
end

if nargout>0, hh = h; end

function h = Lbarv6(varargin)
error(nargchk(1,inf,nargin,'struct'));
[cax,args] = axescheck(varargin{:});

[msg,x,y,xx,yy,linetype,plottype,barwidth,equal] = makebars(args{:}); %#ok
if ~isempty(msg), error(msg); end %#ok

% Create plot
cax = newplot(cax);
fig = ancestor(cax,'figure');

next = lower(get(cax,'NextPlot'));
hold_state = ishold(cax);
edgec = get(fig,'defaultaxesxcolor');
facec = 'flat';
h = []; 
cc = ones(size(xx,1),1);
if ~isempty(linetype), 
    facec = linetype; 
end
for i=1:size(xx,2)
  numBars = (size(xx,1)-1)/5;
  f = 1:(numBars*5);
  f(1:5:(numBars*5)) = [];
  f = reshape(f, 4, numBars);
  f = f';

  v = [xx(:,i) yy(:,i)];

  h=[h patch('faces', f, 'vertices', v, 'cdata', i*cc, ...
             'FaceColor',facec,'EdgeColor',edgec,'parent',cax)];
end
if length(h)==1, 
    set(cax,'clim',[1 2]), 
end
if ~equal, 
  hold(cax,'on'),
  plot(x(:,1),zeros(size(x,1),1),'*','parent',cax)
end
if ~hold_state, 
  % Set ticks if less than 16 integers
  if all(all(floor(x)==x)) && (size(x,1)<16),  
    set(cax,'xtick',x(:,1))
  end
  hold(cax,'off'), view(cax,2), set(cax,'NextPlot',next);
  set(cax,'Layer','Bottom','box','on')
  % Turn off edges when they start to overwhelm the colors
  if size(xx,2)*numBars > 150, 
    set(h,{'edgecolor'},get(h,{'facecolor'}));
  end
end

function [args,pvpairs,msg] = parseargs(args)
msg = '';
% separate pv-pairs from opening arguments
[args,pvpairs] = parseparams(args);
% check for LINESPEC or bar layout
done = false;
while ~isempty(pvpairs) && ~done
  arg = pvpairs{1};
  [l,c,m,tmsg]=colstyle(arg,'plot');
  if isempty(tmsg)
    pvpairs = pvpairs(2:end);
    if ~isempty(l) 
      pvpairs = {pvpairs{:},'LineStyle',l};
    end
    if ~isempty(c)
      pvpairs = {pvpairs{:},'FaceColor',c}; % note FaceColor, not Color
    end
    if ~isempty(m)
      pvpairs = {pvpairs{:},'Marker',m};
    end
  elseif any(strcmpi(arg(1:min(4,length(arg))),{'grou','stac'}))
    pvpairs = {pvpairs{2:end},'BarLayout',arg};
  else
    done = true; % stop looping
  end
end
% check for bar width
if length(args) > 1 && length(args{end}) == 1 && ...
      ~((length(args) == 2) && (length(args{1}) == 1) && (length(args{2}) == 1))
    pvpairs = {'BarWidth',args{end},pvpairs{:}};
    args(end) = [];
end
if isempty(args)
  msg.message = 'Must supply Y data or X and Y data as first argument(s).';
  msg.identifier = id('NoDataInputs');
else
  msg = checkpvpairs(pvpairs,false);
end

function str = id(str)
str = ['MATLAB:bar:' str];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [msg,x,y,xx,yy,linetype,plottype,barwidth,arg8] = makebars(varargin)
%MAKEBARS Make data for bar charts.
%   [MSG,X,Y,XX,YY,LINETYPE,PLOTTYPE,BARWIDTH,EQUAL] = MAKEBARS(X,Y) 
%   returns X and Y, the original data values, XX and YY, the data
%   values formatted to be plotted by one of the BARxx m-files (BAR,
%   BARH, BAR3, BAR3H).
%   
%   LINETYPE returns the color desired for the plot.
%   PLOTTYPE determines whether the plot will be grouped (PLOTTYPE=0),
%   stacked (PLOTTYPE=1), or detached (PLOTTYPE=2--only for 3-D plots).
%   BARWIDTH has the bar width (normalized to one).
%
%   [MSG,X,Y,XX,YY,LINETYPE,PLOTTYPE,BARWIDTH,ZZ] = MAKEBARS(X,Y)
%   does the same as above, except for the final parameter.
%   ZZ is the z-axis data for 3-D plots; used in BAR3 and BAR3H.
%
%   [...] = MAKEBARS(X,Y,WIDTH) or MAKEBARS(Y,WIDTH) returns the
%   specified width given in WIDTH.  The default is 0.8.
%   [...] = MAKEBARS(...,'grouped') returns the data in the form
%   so that the information will be plotted in groups.
%   [...] = MAKEBARS(...,'detached') {3-D only} returns the data
%   such that the information will be plotted detached.
%   [...] = MAKEBARS(...,'stacked') returns the data such that the
%   information will be plotted in stacked form.
%   [...] = MAKEBARS(...,'hist') creates centered bars that touch.
%   [...] = MAKEBARS(...,'histc') creates bars that touch edges.
%   EQUAL is true if spacing of the data is equal; false otherwise.
%   EQUAL is always true except for 'hist' and 'histc' plottypes.
%
%   See also HIST, PLOT, BAR, BARH, BAR3, BAR3H.

%   [...] = MAKEBARS(X,Y,[MIN MAX],'hist') creates bars where the
%   first bar includes MIN and the last includes MAX.

%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.22.4.7 $  $Date: 2005/10/28 15:54:37 $

% Initialize everything
x = []; y=[]; xx=[]; yy=[]; linetype=[]; plottype=[]; barwidth=[];
arg8 = [];

msg = nargchk(1,5,nargin);
if ~isempty(msg), return, end

barwidth = .8; % Normalized width of bar.
groupwidth = .8; % Normalized width of groups.
linetype = []; % Assume linetype is not specified
ishist = 0; % Assume this isn't a histogram
extremes = [];

nin = nargin;

if strcmp(varargin{nin},'3'), 
  threeD = 1; 
  nin = nin-1; 
  plottype = 2; % Detached plot default
else
  threeD = 0; 
  plottype = 0; % Grouped plot default
end

if ischar(varargin{nin}), % Try to parse this string as a color
  [ls,co,mark,msg] = colstyle(varargin{nin});         %#ok
  if isempty(msg), linetype = co; nin = nin - 1; end
end

if ischar(varargin{nin}), % Process 'grouped','stacked' or
                         % 'detached' string.
  kind = [lower(varargin{nin}) 'g'];
  if kind(1)=='g', % grouped
    plottype = 0;
  elseif kind(1)=='s', % stacked
    plottype = 1;
  elseif threeD && kind(1)=='d', % detached
    plottype = 2;
  elseif kind(1)=='h', % histogram
    if strcmpi(varargin{nin},'histc')
      ishist = -1; barwidth = 1;
    else
      ishist = 1; barwidth = 1;
      if nin == 4 && isnumeric(varargin{3}) && ...
            length(varargin{3}) == 2
        extremes = varargin{3};
        nin = nin-1;
      end
    end
  else
    msg.message = sprintf('Unrecognized option "%s".', ...
                          varargin{nin});
    msg.identifier = 'MATLAB:makebars:UnrecognizedOption';
    return
  end
  nin = nin-1;
end

% Parse input arguments.
if ((nin>1) && isscalar(varargin{nin}) && ~isscalar(varargin{nin-1})) || ...
      ((nin>2) && isscalar(varargin{nin}))
  % If last argument is a scalar and next to last isn't then last
  % argument must be the barwidth. Also if the third arg is a
  % scalar then that is the barwidth
  barwidth = varargin{nin};
  [msg,x,y] = xychk(varargin{1:nin-1});
else
  [msg,x,y] = xychk(varargin{1:nin});
end

% Make sure x is monotonic
[x,ndx] = sort(x);
if min(size(y))==1, y = y(ndx); else y = y(ndx,:); end

if ~isempty(msg), return, end

% Expand x to be the same size as y.
if min(size(y))>1, x = x(:,ones(size(y,2),1)); end

% Make sure vector arguments are columns
if min(size(y))==1, x = x(:); y = y(:); end

[n,m] = size(y);

% Remove y values to 0 where y is NaN;
k = find(isnan(y));
if ~isempty(k), y(k) = 0; end

if threeD,
  z = y; % Make variables consistent with 3-D bar graph
  y = x;

  if m==1 || plottype~=0, 
    groupwidth = 1;
  else
    groupwidth = min(groupwidth,m/(m+1.5));
  end

  nn = 6*n; mm = 4*m;
  zz = zeros(nn,mm);
  yy = zz;
  xx = zeros(nn,4);

  % Define xx
  xx(:,3:4) = 1;
  if plottype==0,
    xx = (xx-0.5)*barwidth*groupwidth/m;
  else
    xx = (xx-0.5)*barwidth*groupwidth;
  end

  % Define heights
  zz(2:6:nn,2:4:mm) = z;
  zz(2:6:nn,3:4:mm) = z;
  zz(3:6:nn,2:4:mm) = z;
  zz(3:6:nn,3:4:mm) = z;
        
  if plottype==1 && m>1, % Stacked
    z = cumsum(z.').';
    zz = zz + [zeros(nn,4) z(ones(6,1)*(1:n),ones(4,1)*(1:m-1))];
  end     
        
  if length(y)==1 || all(max(diff(y))==0),
    equal = []; % Special case
  else
    equal = max(abs(diff(diff(y)))) <= max(max(abs(y)))*sqrt(eps(class(y)));
    if isempty(equal), equal = 0; end % check for double-diff in above line
  end
        
  %       
  % Determine beginning of bars (t) and bar spacing (delta)
  %       
  if isempty(equal), % Scalar case and special case
    delta = ones(size(y));
    t = y - 0.5*delta;
    x = 1;
  elseif equal
    if plottype~=0, % Stacked or detached
      delta = ones(n,1) * (max(y) - min(y)) * groupwidth / (n-1);
      t = y - 0.5*delta;
    else % grouped
      delta = ones(n,1) * (max(y) - min(y)) * groupwidth / (n-1) / m ;
      t = y - 0.5*delta*m + (ones(n,1)*(0:m-1)).*delta;
    end 
    if plottype==2, x = 1:m; else x = ones(1,m); end
  else % Spacing is unequal.
    if ishist==1 % Width of bin is average of difference between points.
      dy = diff(y); dy = (dy([1 1:end],:)+dy([1:end end],:))/2;
      t = [(y(1:end-1,:)+y(2:end,:))/2;y(end,:)+(y(end,:)-y(end-1,:))/2]-dy/2;
    elseif ishist==-1 % Width of bin is difference between edges
      dy = [diff(y,1);zeros(1,size(y,2))];
      t = y + dy*groupwidth/2;
    else % Width of bin is minimum of difference between points.
      dy = ones(n,1)*min(diff(y));
      t = y;
    end 
    if plottype~=0, % Stacked or detached
       delta = dy * groupwidth;
       t = t - delta/2;
    else % Grouped
       delta = dy * groupwidth / m ;
       t = t - delta/2*m + (ones(n,1)*(0:m-1)).*delta;
    end
    if plottype==2, x = 1:m; else x = ones(1,m); end
  end     
  t = t(:,ones(4,1)*(1:m));
  delta = delta(:,ones(4,1)*(1:m));
  yy(1:6:nn,:) = t + (1-barwidth)/2.*delta;
  yy(2:6:nn,:) = t + (1-barwidth)/2.*delta;
  yy(3:6:nn,:) = t + (1+barwidth)/2.*delta;
  yy(4:6:nn,:) = t + (1+barwidth)/2.*delta;
  yy(5:6:nn,:) = t + (1-barwidth)/2.*delta;
        
  % Insert NaN's so distinct bars are drawn
  ii1 = [(1:6:nn) (4:6:nn) (5:6:nn) (6:6:nn)];
  ii2 = 6:6:nn;
  xx(ii1,[1 4]) = nan;
  xx(ii2,[2 3]) = nan;
  yy(ii1,[(1:4:mm) (4:4:mm)]) = nan;
  yy(ii2,[(2:4:mm) (3:4:mm)]) = nan;
  zz(ii1,[(1:4:mm) (4:4:mm)]) = nan;
  zz(ii2,[(2:4:mm) (3:4:mm)]) = nan;
  arg8 = zz;
else

  nn = 5*n;
  yy = zeros(nn+1,m);
  xx = yy;
  if plottype && (m>1)
       yc = cumsum(y,2);
       ys = [zeros(n,1),yc(:,1:end-1)];     

       yy(2:5:nn,:) = ys;
       yy(3:5:nn,:) = yc;
       yy(4:5:nn,:) = yc;
       yy(5:5:nn,:) = ys;
  else
       yy(3:5:nn,:) = y;
       yy(4:5:nn,:) = y;
  end
%%% IAT change to make bars use full group width for 'histc' plot
  if m==1 || plottype || ishist==-1, 
    groupwidth = 1;
  else
    groupwidth = min(groupwidth,m/(m+1.5));
  end
%%% End IAT change

  equal = max(abs(diff(diff(x)))) <= max(max(abs(x)))*sqrt(eps(class(x)));
  if ishist==0,
    if isempty(equal), equal = 0; end
    if length(x)==1 || all(max(diff(x))==0), equal = []; end
  else
    if isempty(equal), equal = 1; end
  end

  %
  % Determine beginning of bars (t) and bar spacing (delta)
  %
  dt = [];
  t = x;
  if isempty(equal), % Scalar case and special case
    delta = 1;
    if (ishist~=-1), dt = 0.5*delta; end
  elseif equal
    if plottype,
      delta = ones(n,1) * (max(x) - min(x)) * groupwidth / (n-1);
      if (ishist~=-1), dt = 0.5*delta; end
    else
      delta = ones(n,1) * (max(x) - min(x)) * groupwidth / (n-1) / m ;
      t = x + (ones(n,1)*(0:m-1)).*delta;
      if (ishist~=-1), dt = 0.5*delta*m; end
    end
  else % Spacing is unequal.
    if ishist==1, % Width of bin is average of difference between points.
      dx = diff(x); dx = (dx([1 1:end],:)+dx([1:end end],:))/2;
      t = [(x(1:end-1,:)+x(2:end,:))/2;x(end,:)+(x(end,:)-x(end-1,:))/2]-dx/2;
    elseif ishist==-1, % Width of bin is difference between edges
      dx = diff(x,1); dx = dx([(1:end) end],:);
%      t = x + dx*groupwidth/2;
    else % Width of bin is minimum of difference between points.
      dx = ones(n,1)*min(diff(x));
    end
    if plottype,
      delta = dx * groupwidth;
      if (ishist~=-1), dt = delta/2; end
    else
      delta = dx * groupwidth / m ;
      t = t + (ones(n,1)*(0:m-1)).*delta;
      if (ishist~=-1), dt = delta/2*m; end
    end
  end
  if (length(dt)>0), t = t - dt; end

  xx(1:5:nn,:) = t;
  xx(2:5:nn,:) = t + (1-barwidth)/2.*delta;
  xx(3:5:nn,:) = t + (1-barwidth)/2.*delta;
  xx(4:5:nn,:) = t + (1+barwidth)/2.*delta;
  xx(5:5:nn,:) = t + (1+barwidth)/2.*delta;

  % ensure the bars cover the entire hist data range
  if ~isempty(extremes) && size(xx,2) == 1
    xx(2) = min(xx(2), extremes(1));
    xx(3) = min(xx(3), extremes(1));
    xx(nn) = max(xx(nn), extremes(2));
    xx(nn-1) = max(xx(nn-1), extremes(2));
  end
  xx(nn+1,:) = xx(nn,:);

  if (ishist~=1) && (ishist~=-1), equal = 1; end
  arg8 = equal;
end


