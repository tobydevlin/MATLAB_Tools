% Script to write csv of loads due to bcs between two timesteps
% Handles extrapolation the same way as TUFLOW FV.
% Requires a full path to an fvc (can be an include file) and it will read
% all Q and QC boundaries within this and summarize
% All BCS must have a bc header specified to know what to read.
% Also requires time bounds to integrate between.
%
%
% fvc =  'L:\B22036.L.iat.TamarScenarios\Modelling\TUFLOWFV\bc\catchment\inflows_v03\inflows_control.fvc';
% ts  =  '01/12/2014 00:00:00'; % dd/mm/yyyy HH:MM:SS
% te  =  '01/03/2015 00:00:00'; % dd/mm/yyyy HH:MM:SS
%
% Output units are the concentrations integrated over the volume, so no
% units changes are made ie. mMol/m3 concentrations will result in the 
% output of this in mMol. Variables such as temperature will require 
% changing in excel to make sense.
%
% TDEVLIN FEB 2017
%
% Revisions
%   1.01 : Bug fixes for Adi (TJD 08/03/2017)
%
%
function check_Q_loads(fvc,ts,te)

inblock = false;
ts = datenum(ts,'dd/mm/yyyy HH:MM:SS');
te = datenum(te,'dd/mm/yyyy HH:MM:SS');
fidout = fopen(['Flow_Load_Sum_' datestr(ts,'yyyymmddHHMMSS') '_' datestr(te,'yyyymmddHHMMSS') '.csv'],'W');
isfirs = true;
fid = fopen(fvc);
while ~feof(fid)
    
    line = strtrim(strrep(fgetl(fid),' ',''));
    tmp = strsplit(line,'!');
    line = tmp{1};
    tmp = strsplit(line,{',','=='},'CollapseDelimiters',true);
    
    
    if inblock
        if strncmpi(line,'ENDBC',5)
            pat = fileparts(fvc);
            if ~exist(fullfile(pat,fil),'file')
                error(['Could not find file: ' fullfile(pat,fil)]); % TODO, make this a warning and set data to NaN in output csv... 
            end
            fidsub = fopen(fullfile(pat,fil));
            hedline = strsplit(fgetl(fidsub),',');
            fmt = strsplit(repmat('%s,',1,length(hedline)),',');
            if isfirs
                fprintf(fidout,'%s,%s,%s,%s','FileName','X','Y','NS');
            end
            for aa = 2:length(heds)
                fmt(strncmpi(hedline,strtrim(heds{aa}),length(heds{aa}))) = {'%f'};
                if isfirs
                    if aa==2
                        nm='VOL (m3)';
                    else
                        nm = sprintf('S%d',aa-2);
                    end
                    fprintf(fidout,',%s',nm);
                end
            end
            if isfirs
                fprintf(fidout,'\r\n');
                isfirs = false;
            end
            
            fmt = [strjoin(fmt,'')];
            data = textscan(fidsub, fmt, 'delimiter', ',');
            fclose(fidsub);
            
            itime = find(strncmpi(hedline,heds{1},length(heds{1})));
            time = datenum(data{itime},'dd/mm/yyyy HH:MM:SS');
            
            tnew = [ts;time(time>ts & time<te);te];
            
            
            fprintf(fidout,'%s,%s,%s,%s',fil,x,y,ns);
            
            for aa = 2:length(heds)
                ih = find(strncmpi(hedline,strtrim(heds{aa}),length(heds{aa})));
                lch = length(ih);
                if lch>1
                    display(['Warning! There are ' num2str(lch) ' Matches to header ''' heds{aa} ''' in file ''' fil ''' Using the first one.'])
                end
                
                if isempty(ih)
                    tmpdat = ones(length(tnew),1).*defa{1}(aa-1);
                else
                    tmpdat = interp1(time,data{ih(1)},tnew,'linear',inf).*scals{1}(aa-1);
                    if isnan(tmpdat(1)); tmpdat(1)=tmpdat(2);end
                    if isnan(tmpdat(end)); tmpdat(end)=tmpdat(end-1);end
                end
                if aa==2
                    tmpdat = tmpdat.*86400;%  = m3/s *s/d = m3/d
                    flow=tmpdat; 
                else
                    tmpdat = tmpdat.*flow; %mM/m3 * m3/d
                end
                OUT=trapz(tnew,tmpdat);
                if any(isnan(OUT))
                    display('NaN Detected')
                end
                fprintf(fidout,',%f',OUT);
            end
            fprintf(fidout,'\r\n');
            inblock=false;
        elseif strncmpi(line,'BCHEADER',8)
            heds = tmp(2:end);
            defa = {NaN(1,length(heds))};
            scals = {ones(1,length(heds))};
        elseif strncmpi(line,'BCSCALE',7)
            tmp = strsplit(line,'==');
            scals = textscan(tmp{2},'%f','delimiter',',');
        elseif strncmpi(line,'BCDEFAULT',9)
            tmp = strsplit(line,'==');
            defa = textscan(tmp{2},'%f','delimiter',',');
        end
    else
        if strncmpi(line,'BC==Q',5)
            inblock = true;
            if tmp{2}=='QC'
                typ = 'QC';
                x=tmp{3}; y=tmp{4};
                ns = 'None'; fil = tmp{5};
            elseif tmp{2}=='Q'
                x='None'; y='None';
                ns = tmp{3}; fil = tmp{4};
                typ='Q';
            end
        end
    end
end

fclose(fidout);
    
    
end
    
    
    
    
% If you are reading this, please make sure you have commented any changes    
    
    
    
    
    
    
    
    