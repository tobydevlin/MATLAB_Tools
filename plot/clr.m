function [rgb] = clr(varargin)
   

    % -- CLR
    % 
    % Returns the RGB of some colors when you specify 
    % a string input. With no argument, it returns 
    % an RGB matrix of all those colors, [nclrs * 3].
    % use input 'demo' to see a plot of all the colors
    %
    % TD 2015

    demo = false;
    if nargin==1
        val = varargin{1};
        if strcmpi(val,'demo')
            demo = true;
        end
    elseif nargin==0
        val = 'all';
    end
    
    switch lower(val)
        case {'orange','or',1}
            rgb = [255/255 104/255 0];
        case {'aqua',2}
            rgb = [26/255 189/255 201/255];
        case {'green','g','grn',3}
            rgb = [11/255 168/255 16/255];
        case {'blue','b',4}
            rgb = [72/255 76/255 196/255];
        case {'lime','lightgreen','lgrn',5}
            rgb = [181/255 230/255 20/255];
        case {'dkblue','darkblue','dkb',6}
            rgb = [22/255 73/255 112/255];
        case {'dkred','darkred','dkr',7}
            rgb = [169/255 16/255 11/255];
        case {'lblue','lightblue','ltblue','lb',8}
            rgb = [0/255 222/255 255/255];
        case {'yellow','y',9}
            rgb = [255/255 230/255 0/255];
        case {'jn','jespergreen','jespergrn','land',10}
        	rgb = [170/255 209/255 183/255];
        case {'stevengreen','stevegreen','land2','se',11}
            rgb=[181/255 212/255 82/255];
        case {12}
        	rgb = [0 0 0];
        case {'all','demo'}
            rgb = [181/255 230/255  20/255 ;
                   255/255 104/255  0      ;
                   26/255  189/255  201/255;
                   22/255  73/255   112/255;
                   169/255 16/255   11/255 ;
                   11/255  168/255  16/255 ;
                   72/255  76/255   196/255;
                   255/255 230/255  0/255  ;
                   0/255   222/255  255/255;
                   170/255 209/255 183/255;
                   181/255 212/255 82/255;
                   0 0 0];
        otherwise
            return
    end
    
    if demo
        for aa = 1 : length(rgb)
            plot([aa aa],'Color',rgb(aa,:),'linewidth',4)
            hold on
        end
    end