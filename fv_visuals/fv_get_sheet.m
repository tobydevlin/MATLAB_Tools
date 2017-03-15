% /////// fv_get_sheet ///////
% obj = fv_get_sheet(obj)
% Processes TUFLOW-FV results (3D & 2D) from netcdf output into 2D values.
%
% options for depth averaging:
%   sigma:     [s1 s2]     - average from s1*depth above the bed up to s2*depth above the bed
%   elevation: [e1 e2]     - average from e1 metres up to e2 metres (refereced to model datum)
%   height:    [h1 h2]     - average from h1 metres above the bed up to h2 metres above the bed
%   depth:     [d1 d2]     - average from d1 metres below the surface down to d2 metres below the surface
%   top:       [t1 t2]     - average from t1 layer down to t2 layer. 1 = top layer
%   bot:       [b1 b2]     - average from b1 layer up to b2 layer. 1 = bottom layer
%
% Where a 2D cell is dry the processed results are set to NaN
% Where 3D model results do not exist between the specified depth averaging limits the processed results are also set to NaN but the stat output is unchanged
% Where 3D model results exist between a subset of the depth averaging limits then the 3D results are depth avearged over the subset
%
% inputs
%   obj = handle to fvres_sheet object generated by fvres_sheet.m
%
% outputs
%   obj.results_cell = structure with fields for each variables processed
%   obj.WORK = structure containing variables required during the processing, created on first call
%
% Customized cell centre model results can be processed with fv_get_sheet.
% When obj.feed = true the customized results within the
% structure obj.results_custom are processed as per obj.ref and obj.range
% If obj.results_custum contain 2D results then they are not processed and
% merely transfered to obj.results_cell with the dry cells set to NaN's
%
% Jesper Nielsen, Copyright (C) BMTWBM 2014

function obj = fv_get_sheet(obj)

% do the hard work once
if obj.WORK.refresh
    
    nci = obj.Nci(1); % more than 1 results files can be used to generate customized results. In which case the stat and layerface_Z are taken from the 1st resfil specified
    
    % clear away any old variables in structure
    obj.ResultsCell = struct(); % obj.results_node is emptied in fvres2node
    
    % check the means of depth averaging
    fv_check_dave(obj.Ref,obj.Range)
    
    % info info info
    names = {'NL';'idx2';'idx3'};
    TMP = netcdf_get_var(nci,'names',names,'timestep',1);
    nl = TMP.NL;                                                         % # of layers in each cell
    nlf = sum(nl) + length(nl);                                          % total # of layer faces
    idx2 = TMP.idx2;
    idx3 = TMP.idx3;
    nl_max = max(nl);
    nc2 = length(idx3);
    nc3 = length(idx2);                                                     % # cells, sum(ne*number of layers in element)
    
    % ready yourself for netcdf.getVar (even when feeding in data you want layerface_Z and stat)
    variables = cat(1,obj.Variables,{'layerface_Z';'stat'});
    nv = length(variables);
    [~,~,~,unlimdimid] = netcdf.inq(nci);
    varid = zeros(nv,1);
    i_ud = zeros(nv,1);
    is_2d = false(nv,1);
    is_3d = false(nv,1);
    is_zl = false(nv,1);
    is_bed = false(nv,1);
    nb = 0;
    for aa = 1:nv
        v_name = variables{aa};
        varid(aa) = netcdf.inqVarID(nci,v_name);
        [~, ~, dimids, ~] = netcdf.inqVar(nci,varid(aa));
        nd = length(dimids);
        START.(v_name) = zeros(nd,1);
        for bb = 1:nd
            [dimname,dimlen] = netcdf.inqDim(nci,dimids(bb));
            if dimids(bb) == unlimdimid;
                COUNT.(v_name)(bb) = 1;
                i_ud(aa) = bb;
            else
                COUNT.(v_name)(bb) = dimlen;
            end
            switch dimname
                case 'NumCells2D'
                    is_2d(aa) = true;
                case 'NumCells3D'
                    is_3d(aa) = true;
                case 'NumLayerFaces3D'
                    is_zl(aa) = true;
                case 'NumSedFrac'
                    is_bed(aa) = true;
                    nb = dimlen;
            end
        end
    end
    is_2d(is_bed) = false;
    
    % 3D variables (V_x, TSS, SAL etc.) are 2D and require no processing when the simulation was performed in 2D
    % If this is the case then the layerface_Z variable is not extracted.
    if nc2 == nc3
        is_2d(is_3d) = true;
        is_3d(:) = false;
    end
    
    % when feeding in customized results
    if obj.Feed
        variables_custom = fieldnames(obj.ResultsCustom);
        nvc = length(variables_custom);
        pad = false(nvc,1);
        is_2d = cat(1,pad,is_2d);
        is_3d = cat(1,pad,is_3d);
        is_zl = cat(1,pad,is_zl);
        is_bed = cat(1,pad,is_bed);
        varid = cat(1,zeros(nvc,1),varid);
        i_ud = cat(1,zeros(nvc,1),i_ud);
        for aa = 1:nvc
            v_name = variables_custom{aa};
            [dim1,dim2] = size(obj.ResultsCustom.(v_name));
            if dim1 == nc2 || (dim1 == 1 && dim2 == nc2) % second condition is for bed variables which have been either indexed into a sigle sediment class or have been summed across the sediment classes
                is_2d(aa) = true;
            elseif dim1 == nc3
                is_3d(aa) = true;
            else
                error('The 1st dimension of the customized results must be consistent with either NumCells2D or NumCells3D')
            end
        end
        variables = cat(1,variables_custom,variables);
    end
    
    % indexing required for processing
    if any(is_3d)
        % assign 3D cells for given 2D cell to specific column
        ir = false(nl_max,nc2,sum(is_3d));
        for aa = 1:nc2
            ir(1:nl(aa),aa,:) = true;
        end
        
        % assign 3D layers for given 2D cell to specific column
        il = cat(1,true(1,nc2,1),ir(:,:,1));
        
        % index into top layers
        itop = double(idx3) + (0:nc2-1)';
        
        % index into bottom layers
        ibot = itop(2:end) - 1;
        ibot(nc2) = nlf;
        
        % idx4 (align depths with layerfaces)
        idx4 = ones(nlf,1);
        k = 1;
        for aa = 1:nc2
            kk = k + nl(aa);
            idx4(k:kk) = aa;
            k = kk + 1;
        end
    else
        itop = [];
        ibot = [];
        ir = [];
        idx4 = [];
        il = [];
    end
    
    % preallocate and store
    obj.WORK.mod3 = zeros(nc3,1,sum(is_3d),'single');
    obj.WORK.mod2 = zeros(nc2,1,sum(is_2d)-1,'single');
    obj.WORK.bed = zeros(nb,nc2,sum(is_bed),'single');
    obj.WORK.res = zeros(nl_max,nc2,sum(is_3d),'single');
    obj.WORK.laym = zeros(nl_max+1,nc2,1,'single');
    
    % store small dicky variables
    v = {'variables';'START';'COUNT';'varid';'is_2d';'is_3d';'is_bed';'is_zl';'i_ud';'ir';'il';'idx2';'idx4';'itop';'ibot';'nl'};
    for aa = 1:length(v)
        eval(['obj.WORK.(v{aa}) = ' v{aa} ';'])
    end
    obj.WORK.refresh = false;
else
    variables = obj.WORK.variables; START = obj.WORK.START; COUNT = obj.WORK.COUNT; varid = obj.WORK.varid;
    is_2d = obj.WORK.is_2d; is_3d = obj.WORK.is_3d; is_bed = obj.WORK.is_bed; is_zl = obj.WORK.is_zl;
    i_ud = obj.WORK.i_ud; ir = obj.WORK.ir; il = obj.WORK.il;
    idx2 = obj.WORK.idx2; idx4 = obj.WORK.idx4;
    itop = obj.WORK.itop; ibot = obj.WORK.ibot; nl = obj.WORK.nl;
end
% changes in the following 3 variables outside fv_get_sheet will not effect what work has needed to be done above
ref = obj.Ref;
range = obj.Range;
bedref = obj.BedRef;
bedfrac = obj.BedFrac;
feed = obj.Feed;
nci = obj.Nci(1); % more than 1 results files can be used to generate customized results
it = obj.TimeStep(1);   % as above

% extract & process all cells for given timestep
% -- retrieve results
nv = length(variables);
i = 1;
j = 1;
k = 1;
for aa = 1 : nv
    v_name = variables{aa};
    if ~feed || ~ismember(v_name,{'custom';'custom_x';'custom_y'}) 
        START.(v_name)(i_ud(aa)) = it - 1;
    end
    if is_2d(aa)
        switch v_name
            case 'stat'
                stat_raw = netcdf.getVar(nci,varid(aa),START.(v_name),COUNT.(v_name));
                stat = stat_raw == -1;
            otherwise
                if feed
                    obj.WORK.mod2(:,1,i) = obj.ResultsCustom.(v_name);
                else
                    obj.WORK.mod2(:,1,i) = netcdf.getVar(nci,varid(aa),START.(v_name),COUNT.(v_name));
                end
                i = i+1;
        end
    elseif is_bed(aa)
        obj.WORK.bed(:,:,j) = netcdf.getVar(nci,varid(aa),START.(v_name),COUNT.(v_name));
        j = j+1;
    elseif is_3d(aa)
        if feed
            obj.WORK.mod3(:,1,k) = obj.ResultsCustom.(v_name);
        else
            obj.WORK.mod3(:,1,k) = netcdf.getVar(nci,varid(aa),START.(v_name),COUNT.(v_name));
        end
        k = k+1;
    elseif is_zl(aa) && any(is_3d)
        lay = netcdf.getVar(nci,varid(aa),START.(v_name),COUNT.(v_name));
    end
end

% -- stamping of results, when soft dry cells are still set to NaN
if ~isempty(obj.StampVec)
    switch obj.StampType
        case 'hard'
            stat = obj.StampVec;
        case 'soft'
            stat = all([stat obj.StampVec],2);
    end
end
% -- results in dry cells are set to NaN (Bed variables are not effected by stat)
if any(is_2d)
    obj.WORK.mod2(~stat,:,:) = NaN;
end
if any(is_3d)
    stat3 = stat(idx2);
    obj.WORK.mod3(~stat3,:,:) = NaN;
end

% -- process (depth average) 3D results
if any(is_3d)
    obj.WORK.res(ir) = obj.WORK.mod3;
    top = lay(itop);
    bot = lay(ibot);
    % -- limits of depth averaging, d1 is below d2
    switch ref
        case 'sigma'
            d = top - bot;
            d1 = bot + range(1) * d;
            d2 = bot + range(2) * d;
        case 'elevation'
            d1 = max(bot,range(1));
            d2 = min(top,range(2));
        case 'height'
            d1 = bot + range(1);
            d2 = min(bot+range(2),top);
        case 'depth'
            d1 = max(top-range(2),bot);
            d2 = top - range(1);
        case 'top'
            id1 = min(itop+range(2),ibot);
            id2 = itop + range(1) - 1;
            i_nolay = nl < range(1);
            id2(i_nolay) = 1;   % dummy assignment
            d1 = lay(id1);
            d2 = lay(id2);
        case 'bot'
            id1 = ibot - range(1) + 1;
            i_nolay = nl < range(1);
            id1(i_nolay) = 2;   % dummy assignment
            id2 = max(ibot-range(2),itop);
            d1 = lay(id1);
            d2 = lay(id2);
    end
    % -- engine room
    d = d1 - d2;
    lay = max(lay,d1(idx4));
    lay = min(lay,d2(idx4));
    obj.WORK.laym(il) = lay;
    frc = bsxfun(@rdivide,diff(obj.WORK.laym),d');
    out = bsxfun(@times,obj.WORK.res,frc);
    dav = sum(out,1);
    
    % -- cells beyond depth averaging limits are set to NaN, stat variable is unchanged
    switch ref
        case {'elevation','height','depth'}
            i3 = d1 > top | d2 < bot; % could be replaced by d <= 0
            dav(1,i3,:) = NaN;
        case {'top','bot'}
            dav(1,i_nolay,:) = NaN;
    end
end

% store away
i = 1;
j = 1;
k = 1;
for aa = 1 : nv
    v_name = variables{aa};
    if is_2d(aa)
        switch v_name
            case 'stat'
                obj.ResultsCell.stat = stat;
            otherwise
                obj.ResultsCell.(v_name) = obj.WORK.mod2(:,1,i);
                i = i+1;
        end
    elseif is_bed(aa)
        if isempty(bedref)
            v_name = strcat(v_name,'_TOTAL');
            obj.ResultsCell.(v_name) =sum(obj.WORK.bed(:,:,j),1); % it is possible that we are overwriting the _TOTAL variable generated by TUFLOW-FV
        else
            if bedfrac
                tol = realmin('single') * 1000;
                tot = sum(obj.WORK.bed(:,:,j),1);
                tot(tot<tol) = Inf;
            end
            for bb = 1:length(bedref)
                v_name = strcat(v_name,'_SED_',num2str(bedref(bb)));
                if bedfrac
                    obj.ResultsCell.(v_name) = obj.WORK.bed(bedref(bb),:,j) ./ tot;
                else
                    obj.ResultsCell.(v_name) = obj.WORK.bed(bedref(bb),:,j);
                end
            end
        end
        j = j + 1;
    elseif is_3d(aa)
        obj.ResultsCell.(v_name) = dav(1,:,k)';
        k = k+1;
    end
end