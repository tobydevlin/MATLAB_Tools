% Box and Whisker Plots
% Takes a profile file in and plots the box and whisker plots at each
% location from the depth averaged data as specified.
% Takes an expression as the variable if required - usual expression rules
% apply
%
% Input:
%           profil:     valid TFV profile file
%           variable:   expression or variable (1 at a time)
%           ref:        depth averaging reference
%           range:      depth averaging range
%           labels:     Xaxis Labels corresponding to locations
%
% Output: 
%           f:          Handle to Figure created
%
%
% TD Sep2015



function f=boxandwhisker(profil,variable,ref,range,labels)


cbox = clr(2);
cmean = clr(1);

% Parse Inputs
if nargin<5
    labels={};
    if nargin<4
        range = [0 1];
        if nargin<3
            ref = 'sigma';
        end
    end
end       

% Check File
nci = netcdf.open(profil);
if ~strcmpi(netcdf.getAtt(nci,netcdf.getConstant('NC_GLOBAL'),'Type'),'TUFLOWFV Profile Output')
    error('Not A Profile File')
end
grids = netcdf.inqGrps(nci);
np = length(grids);
  
% Setup Figure    (people can edit their own afterwards)
f = figure('color','w','inverthardcopy','off','units','centimeters','position',[10 10 16 10],'paperunits','centimeters','paperposition',[0 0 16 10]);
ax = axes('parent',f,'units','centimeters','position',[1.1 1.1 13.9 8],'color','none','yaxislocation','left','ylim',[0 200],'box','off'); %#ok<*NASGU>


for aa = 1 : np
    
    %{
    Qmin: min
    Q1: first quartile
    Q2: median
    Q3: 3rd quartile
    Qmax: max
    Qmean: average
    %}
    
    % Get Location Name and Data
    location{aa} = netcdf.inqGrpName(grids(aa)); %#ok<*AGROW>
    data = gettseries(profil,location{aa},variable,ref,range);
    
    
    % Find stats as above
    tmpset = data(~isnan(data));
    Q = prctile(tmpset,[25 50 75]);
    Qmin = min(tmpset(:));
    Qmax = max(tmpset(:));
    Qmean = mean(tmpset(:));
    
    % Box and Whisker
    x=(aa-1)+[0.5,0.5,0.75,0.75,0.25,0.25,0.5,0.5,0.5,0.75,0.75,0.25,0.25,0.5];
    y=[Qmax,Q(3),Q(3),Q(2),Q(2),Q(1),Q(1),Qmin,Q(1),Q(1),Q(2),Q(2),Q(3),Q(3)];
    plot(x,y,'color',cbox)
    
    hold on
    
    % Extrema
    x=aa-1+[0.5,0.5];
    y=[Qmin,Qmax];
    plot(x,y,'x','color',cbox)
    
    % Average
    x=aa-1+0.5;
    y=Qmean;
    plot(x,y,'d','markerfacecolor',cmean,'markeredgecolor',cmean);
    
end

netcdf.close(nci);

if isempty(labels)
    labels=location;
end
set(gca,'XTick',0.5:1:np)
set(gca,'XTickLabel',labels)

end