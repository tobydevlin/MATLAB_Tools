

% Start User Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnam = '474A.tpi' % constituent file (name amplitude phase)
start_date = '01/08/2009 00:00' % dd/mm/yyyy HH:MM
end_date = '01/10/2009 00:00' % dd/mm/yyyy HH:MM
interval = 0.25; % hours
latitude = 0.80; % decimal degrees
% End User Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(fnam,'r');
C = textscan(fid, '%4s %f %f');
fclose(fid);

J = length(C{1});
names = char(C{1});
if size(names,2)<4
    names = [names,repmat(' ',J,4-size(names,2))];
end
tidecon = zeros(J,4);
tidecon(:,1) = C{2};
tidecon(:,3) = C{3};

time = [datenum(start_date,'dd/mm/yyyy HH:MM'):...
    interval/24:...
    datenum(end_date,'dd/mm/yyyy HH:MM')]';
T = length(time);

const = t_getconsts(mean(time));
freq = zeros(J,1);
for j = 1 : J
    ind = strmatch(names(j,:),const.name);
    freq(j) = const.freq(ind);
end

tpred = t_predic(time,names,freq,tidecon,'latitude',latitude);