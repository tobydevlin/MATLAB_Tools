% Read an SMS .2DM file

function MESH = RD2DM(fnam)

MESH = struct();

ND = [200000,4];
E3T = [200000,5];
E4Q = [200000,6];
E6T = [200000,8];
E8Q = [200000,10];
E9Q = [200000,11];

fid = fopen(fnam);
if fid == -1, error('ERROR:RD2DM:fopen', 'Unable to open %s', fnam), end
disp(['Reading ', fnam])

i = 0;
while feof(fid)~=1
    i = i + 1;
    lin = fgetl(fid);
    switch lower(deblank(lin(1:3)))
        case 'mes'
        case 'num'
        case 'nd'
            if ~isfield(MESH,'ND')
                ind = 0;
                MESH.ND = zeros(ND);
            end
            ind = ind + 1;
            MESH.ND(ind,:) = strread(lin(3:end),'%n',ND(2));
        case 'e3t'
            if ~isfield(MESH,'E3T')
                ie3 = 0;
                MESH.E3T = zeros(E3T);
            end
            ie3 = ie3 + 1;
            tmp=strread(lin(4:end),'%n',E3T(2));
            MESH.E3T(ie3,1:length(tmp)) = tmp;
        case 'e4q'
            if ~isfield(MESH,'E4Q')
                ie4 = 0;
                MESH.E4Q = zeros(E4Q);
            end
            ie4 = ie4 + 1;
            tmp=strread(lin(4:end),'%n',E4Q(2));
            MESH.E4Q(ie4,1:length(tmp)) = tmp;
        case 'e6t'
            if ~isfield(MESH,'E6T')
                ie6 = 0;
                MESH.E6T = zeros(E6T);
            end
            ie6 = ie6 + 1;
            tmp=strread(lin(4:end),'%n',E6T(2));
            MESH.E6T(ie6,1:length(tmp)) = tmp;
        case 'e8q'
            if ~isfield(MESH,'E8Q')
                ie8 = 0;
                MESH.E8Q = zeros(E8Q);
            end
            ie8 = ie8 + 1;
            tmp=strread(lin(4:end),'%n',E8Q(2));
            MESH.E8Q(ie8,1:length(tmp)) = tmp;
        case 'e9q'
            if ~isfield(MESH,'E9Q')
                ie9 = 0;
                MESH.E9Q = zeros(E9Q);
            end
            ie9 = ie9 + 1;
            tmp=strread(lin(4:end),'%n',E9Q(2));
            MESH.E9Q(ie9,1:length(tmp)) = tmp;
        case 'ns'
            if ~isfield(MESH,'NS')
                ins = 1;
                MESH.NS{ins} = [];
            end
            MESH.NS{ins} = [MESH.NS{ins};strread(lin(4:end),'%n',100)];
            if MESH.NS{ins}(end-1) < 0
                MESH.NS{ins}(end-1) = -MESH.NS{ins}(end-1);
                MESH.NS{ins}(end) = [];
                ins = ins + 1;
                MESH.NS{ins} = [];
            elseif MESH.NS{ins}(end) < 0
                MESH.NS{ins}(end) = -MESH.NS{ins}(end);
                ins = ins + 1;
                MESH.NS{ins} = [];
            end        
        case 'beg'
            break
        otherwise
            error(['ERROR:RD2DM:UnknownLine', 'Unknown line encountered,',num2str(i)])
    end  
end
fclose(fid);
if isfield(MESH,'ND'), MESH.ND = MESH.ND(1:ind,:); end
if isfield(MESH,'E3T'), MESH.E3T = MESH.E3T(1:ie3,:); end
if isfield(MESH,'E4Q'), MESH.E4Q = MESH.E4Q(1:ie4,:); end
if isfield(MESH,'E6T'), MESH.E6T = MESH.E6T(1:ie6,:); end
if isfield(MESH,'E8Q'), MESH.E8Q = MESH.E8Q(1:ie8,:); end
if isfield(MESH,'E9Q'), MESH.E9Q = MESH.EQ(1:ie9,:); end
if isfield(MESH,'NS'), MESH.NS = MESH.NS(1:ins-1); end