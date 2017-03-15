% /////// read_mifmid ///////
% Read a .mif and the accompanying .mid file and store the information in a
% matlab structure. The coordinates are stored in the .mif file and all
% other information stored in the .mid file. The information in the .mid
% file corresponds to the columns in the MapInfo tables.
%
%
% TO DATE ONLY POLYGONS ARE SUPPORTED
%
% Jesper Nielsen, August 2014

function MIF = read_mifmid(miffil)

% check files exist
fid1 = fopen(miffil);
if fid1 == -1
    error('Unable to open %s', miffil)
end
midfil = strrep(lower(miffil),'.mif','.mid');
fid2 = fopen(midfil);
if fid2 == -1
    error('Unable to open %s', midfil)
end

% what information is stored in the .mid file
while ~feof(fid1)
    lin = fgetl(fid1);
    if strcmp(lin(1:7),'Columns')
        C = textscan(lin,'%s%d');
        nh = C{2}; % number of headers / columns in table
        C = textscan(fid1,'%s%s',nh);
        headers = C{1};
        type = C{2}; % leave this up to myformat - either character or numeric
        break
    end
end

% load in the information stored in the .mid file
fmat = myformat(midfil,'Delimiter',',');
MID = textscan(fid2,fmat,'Delimiter',',');
fclose(fid2);

% continue through the .mif file extracting the coords whilst
% assigning the information from the .mid to the structure.
MIF = struct();

i = 1;
j = 1;
k = 1;
l = 1;
m = 1;
% read through file top to bottom
fgetl(fid1); % skip these lines
fgetl(fid1);
fgetl(fid1);
while ~feof(fid1)
    lin = fgetl(fid1);
    switch lower(deblank(lin(1:3)))
        case 'pli'
            error('only regions are supported at present')
            tmp = textscan(fid1,'%15f %15f');
            MIF.poly(i).coordinates = cell2mat(tmp);
            i = i + 1;
            m = m + 1;
        case 'reg'
            r_name = ['REG_' num2str(j)];
            C = textscan(lin,'%s%d');
            np = C{2};     % # polygons in "region"
            for aa = 1:np
                p_name = ['POLY_' num2str(aa)];
                lin = fgetl(fid1);
                C = textscan(lin,'%d'); % # of points in individual polygon
                tmp = textscan(fid1,'%f%f',C{1});
                MIF.(r_name).coordinates.(p_name) = cell2mat(tmp);
                % -- assign the .mid info
                for bb = 1:nh
                    h_name = validfieldname(headers{bb});
                    tmp = MID{bb}(m);
                    if iscell(tmp)
                        tmp = tmp{1};
                        if ischar(tmp)
                            tmp = strrep(tmp,'"','');
                        end
                    end
                    MIF.(r_name).(h_name) = tmp;
                end
                fgetl(fid1); % jump to next line from texscan
            end
            j = j + 1;
            m = m + 1;
            fgetl(fid1); % skip the Pen
            fgetl(fid1); % skip the Brush
            fgetl(fid1); % skip the Centre
        case 'poi'
            error('only regions are supported at present')
            tmp = textscan(fid1,'%15f %15f');
            MIF.point(k).coordinates = cell2mat(tmp);
            k = k + 1;
            m = m + 1;
        case 'tex'
            error('only regions are supported at present')
            lin = fgetl(fid);
            MIF.text(l).string = strrep(deblank(lin),'"','');
            tmp = textscan(fid1,'%15f %15f');
            MIF.text(l).coordinates = cell2mat(tmp);
            l = l + 1;
            m = m + 1;
    end
end
fclose(fid1);
end

function name2 = validfieldname(name1)
name2 = name1;
ind = regexp(name1,'\W'); % alphabetic, numeric, or underscore characters only as valid field names
name2(ind) = '';
if ~isempty(ind)
    display(['WARNING: ' name1 ' in .mif file converted to ' name2])
end
if isempty(regexp(name2(1),'[a-zA-Z]','once'));
    error('headers / column names must begin with a character')
end
end