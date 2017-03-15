% Function to read mid/mif files into a structure
% The resulting structure is number of objects long, with all the
% associated fieldnames attached if the mid can be found.
% Three fields are appended:
%   coords :  for the coordinates of the object
%   type   :  for the type of object (polyline, region, point, text)
%   string :  for text objects.
%
% Input:
%           miffil : the full path to the .mif file
%           midfil : the full path to the .mid file, will assume the same
%                    name as the mif file if not specified.
%
% Output:
%           mif    : structure that is no. objects long, with all fields.
%
%
% TDevlin Jan 2017
%
function mif = read_mifmidv2(miffil, midfil)


if nargin==1
    midfil = strrep(strrep(miffil,'.mif','.mid'),'.MIF','.MID');
end

miffid = fopen(miffil);
if miffid == -1, error('ERROR:read_mifmidv2:fopen', 'Unable to open %s', miffil), end
disp(['Reading ', miffil])

fieldon = false;
k=0;
while feof(miffid)~=1
    lin = fgetl(miffid);
    if length(lin)>2
        switch lower(deblank(lin(1:3)))
            case 'pli' 
                k=k+1;
                tmp = textscan(miffid,'%15f %15f');
                coords(k) = {cell2mat(tmp)};
                type(k) = {'poly'};
                string(k) = {''};
            case 'reg'
                k=k+1;
                lin = fgetl(miffid);
                tmp = textscan(miffid,'%15f %15f');
                coords(k) = {cell2mat(tmp)};
                type(k) = {'region'};
                string(k) = {''};
            case 'poi'
                k = k+1;
                tmp = textscan(miffid,'%15f %15f');
                coords(k) = {cell2mat(tmp)};
                if isempty(coords{k})
                    tmp = textscan(lin(6:end),'%15f %15f');
                    coords(k) = {cell2mat(tmp)};
                end
                type(k) = {'point'};
                string(k) = {''};
            case 'lin'
                k = k+1;
                tmp = textscan(miffid,'%15f %15f');
                coords(k) = {cell2mat(tmp)};
                if isempty(coords{k})
                    tmp = textscan(lin(6:end),'%15f %15f');
                    coords(k) = {cell2mat(tmp)};
                end
                type(k) = {'line'};
                string(k) = {''};
            case 'tex'
                k = k+1;
                lin = fgetl(miffid);
                string(k) = strrep(deblank(lin),'"','');
                tmp = textscan(miffid,'%15f %15f');
                coords(k) = cell2mat(tmp);
                type(k) = 'text';
            case 'col'
                tmp = textscan(lin,'%s %f');
                N = tmp{2};
                tmp = textscan(miffid,'%s %s',N);
                val = tmp{2};
                fmt = tmp{2};
                fields = tmp{1};
                iin = strncmpi(tmp{2},'inte',4); % these are integers..
                ifl = strncmpi(tmp{2},'floa',4); % these are floats
                ich = strncmpi(tmp{2},'char',4); % these are char(#)
                fmt(iin|ifl) = {'%f'};
                fmt(ich) = {'%q'};
                
                fieldon = true;
      end
    end  
end
fclose(miffid);

if fieldon
    fields = [fields(:);'coords';'type';'string'];
    midfid = fopen(midfil);
    if midfid == -1, error('ERROR:read_mifmidv2:fopen', 'Unable to open %s', midfil), end
    disp(['Reading ', midfil])
    data = cell(k,length(fields));
    for aa = 1:k
        tmp = textscan(midfid,strjoin(fmt,''),1,'delimiter',',');
        coo = coords{aa};
        if isempty(coo)
            coo = {coo};
        end
        data(aa,:) = [tmp,coo,type(aa),string(aa)];
    end
    fclose(midfid);
else
    fields = {'coords';'type';'string'};
    data = cell(k,length(fields)+3);
    for aa = 1:k
        coo = coords{aa};
        if isempty(coo)
            coo = {coo};
        end
        data(aa,:) = [coo,type(aa),string(aa)];
    end
end

mif = cell2struct(data,fields',2);

