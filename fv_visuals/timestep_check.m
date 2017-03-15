% Mesh Scanner
% Checks for cells in a mesh that will slow things down.
% Write a file similar to the dt outputs from TUFLOWFV
%
% Inputs:
%           twodm :  path to mesh file
%           zfil  :  path to bathy file (cell centres)
%       spherical :  true|false  if mesh in spherical coordinates
%
% Outputs:
%           Writes a file called %meshname%__dt_check.csv
%           which contains two 'timesteps'. The first one scales on dx only
%           and the second scales on sqrt(cell depth).
%           They dont correspond to actual timesteps, and are just
%           indicative of relative cell timestep constraints.
%
%
% eg.
%       twodm = '../geo/Mesh_000.2dm';
%       zfil = '../geo/Mesh_000_centres.txt';
%       spherical = true;
%
% TDEVLIN March 2017

function timestep_check(twodm, zfil, spherical)

    vel = 1; % Assume current velocity for DT1

    % -- Sort Out Bathy First
    if ~exist(zd,'file')

        % Inspect bathy from mesh perhaps...
        error('Not implemented to handle bathy on mesh yet')

    else
        % Open File
        fid = fopen(zfil);

        % Get Header Line and find columns
        heds = fgetl(fid);
        heds = strsplit(heds,',');
        isz = strcmpi(heds,'Z');
        isid = strcmpi(heds,'ID');

        % Figure out format
        fmt = heds;   fmt(:) = {'%s'};
        fmt(isz|isid) = {'%f'};
        fmt = strjoin(fmt,'');

        % Get Data
        dat = textscan(fid,fmt,'delimiter',',');
        zids = dat{isid};
        zcel = dat{isz};
        % Cleanup
        fclose(fid); clearvars dat fmt heds idz isid

    end


    % -- Now Read 2DM
    fid = fopen(twodm);
    options = {'ND ','E3T','E4Q'};

    % Get Counts of Cells and Nodes
    cnts = [0,0,0];
    while ~feof(fid)
        line = fgetl(fid);
        cnts = cnts+strncmpi(line,options,3);
    end
    node = zeros(cnts(1),2);
    nc = sum(cnts(2:3));
    cell = zeros(nc,4);

    % Now Read Cell and Node info
    frewind(fid);
    while ~feof(fid)
        line = fgetl(fid);
        ix = find(strncmpi(line,options,3));
        if ~isempty(ix)
        switch ix
            case 1 % ND
                tmp = textscan(line,'ND %f %f %f %*[^\n]',...
                                    'CollectOutput',true);
                id = tmp{1}(1);
                node(id,:) = tmp{1}(2:3);
            case 2 % E3T
                tmp = textscan(line,'E3T %f %f %f %f %f %*[^\n]',...
                                    'CollectOutput',true);
                id = tmp{1}(1);
                cell(id,:) = [tmp{1}(2:4) 0];
            case 3 % E4Q
                tmp = textscan(line,'E4Q %f %f %f %f %f %f %*[^\n]',...
                                    'CollectOutput',true);
                id = tmp{1}(1);
                cell(id,:) = tmp{1}(2:5);
        end
        end
    end
    fclose(fid);

    % If spherical then convert nodes to meters (UTM)
    utnode = node;
    if spherical
        [utnode(:,1),utnode(:,2)] = ll2utm(node(:,1),node(:,2));
    end


    % -- Calculate Cell Area./Face Length to get dx and compute dt

    % preallocate
    L = zeros(4,1);
    dt = zeros(nc,2);
    x = zeros(nc,1);
    y = zeros(nc,1);

    % loop through cells
    for aa = 1 : nc

        % Loop through edges
        A = 0;
        L(:) = 0.001;
        for bb = 1 : 4

            % First Point
            n1=cell(aa,bb);
            if n1==0
                bb=3;
                break
            end

            % Second Point
            n2=cell(aa,rem(bb,4)+1);
            if n2==0
                n2=cell(aa,1);
            end

            % Get Length, Area, and Centroid
            points = utnode([n1;n2],:);
            x(aa) = x(aa)+node(n1,1);
            y(aa) = y(aa)+node(n1,2);
            L(bb) = sqrt(sum(diff(points).^2));
            A = A + det(points)./2;
            
        end
        x(aa) = x(aa)./bb;
        y(aa) = y(aa)./bb;

        % Get Dx and sqrt(h)
        dx = min(A./L);
        h = sqrt(1-zcel(zids==aa));

        % Calculate Dt
        dt(aa,:) = [dx./vel,dx./h];    
    end

    % -- Write Output
    fid = fopen(strrep(twodm,'.2dm','_dt_check.csv'),'W');
    fprintf(fid,'%s,%s,%s,%s\r\n','X','Y','DT1','DT2');
    datout = [x,y,dt];
    fprintf(fid,'%f,%f,%f,%f\r\n',datout');
    fclose(fid);


end






