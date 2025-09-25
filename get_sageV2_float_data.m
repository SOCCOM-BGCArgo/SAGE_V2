function d = get_sageV2_float_data(file_path)
% This function extracts profiling float data from *Sprof.nc files and
% generates a structure including a data matrix. MAIN USE IS FOR MBARI
% SAGEV2 QC GUI
%
% INPUTS:
%   file_path - path to Sprof file or to raw ODV file (if available the
%               ODV adjusted data will also be extracted.
%
% OUTPUTS:
%   d.hdr  - cell_array
%   d.data - data matrix
%   d.info - meta info
%   pH_PO  - pH pump offset structure with fields:
%     hdr  - hdr for pH.PO.POdata matrix : {'cycle' 'CpActP' 'spline offset'
%             'spline SSR' 'resid std' 'Zish' 'lin offset' 'lin diff'
%             'lin_resid std' 'lin Zish'}
%     data - pump offset data matix
%     pH   - offset corrected pH for d.data matrix [cycle, P, raw pH, offset corrected raw pH 
%     offset_corr_type -  'linear', 'poly','mix','none'
%
% UPDATES:
% 5/15/25 TM: Added adjusted salinity to indexing and removed exclusion.
% 8/19/25 LG: - Changed kordi to kiost to follow Korean DAC naming convention
%             - Modified logic for the initial file check on lines 112 and 114 to prevent
%               isempty from being used on a logical value, which causes a matlab error.


% *************************************************
% TESTING
% file_path = 'C:\temp\4903274_Sprof.nc'; % wn1357
%file_path = 'C:\temp\5904855_Sprof.nc'; % ua12559 %cpAct P changes after cycle &
%file_path = 'C:\temp\5906028_Sprof.nc'; % ua19065 % SBE83 ONLY
%file_path = 'C:\Users\jplant\Documents\MATLAB\ARGO\DATA\SPROF\aoml\5906568_Sprof.nc'; % ua19065 % SBE83 ONLY

%file_path = 'C:\temp\4903274.TXT'; % wn1357
%fn = '5904855.TXT'; % ua12559 %cpAct P changes after cycle 7. pH bad 1-
%fn = '5906028.TXT'; % ua19065 % SBE83 ONLY
%fn = '1902459.TXT'; % ua19065 % SBE83 ONLY
% % fn  = '5905069.TXT'; %ua12558	12558	
% fn  = '2903473.TXT';	%ua21441	PH	17-	4
% fn  = '1902383.TXT';	%wn1351	PH	73-	4
% fn  = '5906526.TXT'; % ua20108	good pump offset example	
% fn  = '5906568.TXT';
 % fn  = '5906339.TXT';
 % fn  = '5906339.TXT';
% file_path = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATVIZ\',fn];
%file_path = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATVIZ\NO_WMO_un0948.TXT';
%file_path = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATVIZ\5906316.TXT';
% file_path = 'C:\Users\lgrady\Documents\3902506_Sprof.nc';

%file_path = 'C:\temp\junk.ttt'; 

% *******************************************************************
% INIT
d.data    = []; % default output if fuction fails
d.hdr     = '';
d.info    = [];
info.log = {'No issues to log'};

% pH_PO.offset_corr_type = 'poly'; % NOTE: pH pump offeset correction type is hardwired to poly for now

fp     = file_path; clear file_path;

gdac_urls = {'https://data-argo.ifremer.fr/'; ...
             'https://usgodae.org/pub/outgoing/argo/'};

%mat_dir   = '\\atlas\Chem\ARGO_PROCESSING\DATA\FLOATS'; % to get CpActP
log_ct     = 1;
tf_CpAct   = 0;
tf_Oxygen  = 0;
tf_Nitrate = 0;
tf_pH      = 0;

% THIS CAN BE ADDED TO INCLUDE ANY DESIRED PARAMETERS
keep_vars = {.... % Sprof params, FV params, output params
    'JULD'               'SDN'                'SDN'
    'CYCLE_NUMBER'       'Station'            'Cycle'          
    'LONGITUDE'          'Lon [°E]'           'Lon'
    'LATITUDE'           'Lat [°N]'           'Lat'
    'PRES'               'Pressure[dbar]'     'P'
    'PSAL'               'Salinity[pss]'      'S'
    'TEMP'               'Temperature[°C]'    'T'
    'DOXY'               'Oxygen[µmol/kg]'    'O2'
    'NITRATE'            'Nitrate[µmol/kg]'   'NO3'
    'PH_IN_SITU_TOTAL'   'pHinsitu[Total]'    'pH'
    };


% DAC REFERENCE CODES FOUND
dac_codes = { ...
    'AO|MB|SI|UW|PM|WH|NA' 'aoml';
    'BO'                   'bodc';
    'CI|ME'                'meds';
    'CS'                   'csiro';
    'GE|IF'                'coriolis';
    'HZ|NM'                'csio';
    'IN'                   'incois';
    'JA|JM'                'jma';
    'KM'                   'kma';
    'KO'                   'kiost'}; %LG note: check korean float compatibility


% *******************************************************************
% CHECK FOR FILE EXISTANCE & FILE TYPE
if ~isfile(fp)
    str         = sprintf('File not found: %s',fp);
    info.log{log_ct,1} = str; disp(str)
    log_ct      = log_ct+1;
    return
end

% CHECK FILE TYPE
[fd, fn, ext] = fileparts(fp);
if ~isempty(regexp(ext,'nc$','once')) && ~isempty(regexp(fn,'\d{7}_Sprof','once'))
    ftype = 'Sprof';
elseif ~isempty(regexp(ext,'TXT$','once')) && ~isempty(regexp(fn,'^\d{7}|^NO\_WMO','once')) %LG note: changed logic here to prevent matlab error
    ftype = 'ODV';
    % Try and create mat dir from file path dir
    mat_dir = regexprep(fd,'FLOATVIZ','FLOATS');
    if isfolder(mat_dir)
        fprintf('*.mat cycle dir found for CpActP values!\n')
    else % MBARI ONLY!
        str = '*.mat cycle dir not found for MBARI ODV files!';
        info.log{log_ct,1} = str; disp(str)
        log_ct      = log_ct+1;
        return
    end
else
    str         = sprintf('File type is not recognized: %s', [fn,ext]);
    info.log{log_ct,1} = str; disp(str)
    log_ct      = log_ct+1;
    return
end

info.keep_vars = keep_vars;
info.fn        = [fn,ext];
info.source    = fd;
info.ftype     = ftype;
info.QFbad     = 8; % MBARI ODV TXT FILE
info.QFmiss    = 1; % MBARI ODV TXT FILE
if strcmp(ftype,'Sprof')
    info.QFbad   = 4; % Argo Sprof
    info.QFmiss  = 9; % Argo Sprof  
end

% BUILD FLOAT DATA MATRIX HEADER (Add Adjused & QF cols)
meta = keep_vars(1:4,3)'; % sdn - lat
raw  = keep_vars(5:end,3)'; % raw params
tmp  = [raw, regexprep(raw,'(\w$)','$1_adj')]; % raw & adj
QF   = regexprep(tmp,'\w+','QF'); % quality flag cols
hdr  = [meta,'QF',reshape([tmp;QF],1,size(QF,2)*2)]; %QF added for position QC
rraw = size(raw,2); % use for indexing adj, QF to output
%clear meta raw tmp QF

% ************************************************************************
% ************************************************************************
% GET THE DATA
switch ftype
    case 'Sprof'
        nc_info    = ncinfo(fp);
        nc_vars    = {nc_info.Variables.Name}'; % Sprof variables in a col
        tP         = strcmp('PRES',nc_vars);    % use to get matrix size
        Sprof_rc   = nc_info.Variables(tP).Size;
        [tf,~]     = ismember(keep_vars(:,1), nc_vars); % In master list
        % float vars in master, use to check existance in loop below
        avail_vars = keep_vars(tf,1); %params on master in *.nc

        if any(strcmp(avail_vars,'DOXY')) % set sensor existance flags
            tf_Oxygen = 1;
        end
        if any(strcmp(avail_vars,'NITRATE'))
            tf_Nitrate = 1;
        end
        if any(strcmp(avail_vars,'PH_IN_SITU_TOTAL'))
            tf_pH = 1;
        end

        % GET WMO & DAC FROM DATA CENTRE CODE
        info.wmo = ncread(fp,'PLATFORM_NUMBER',[1 1],[7 1])'; % wmo as 7 char str
        DC = ncread(fp,'DATA_CENTRE',[1,1],[inf,1])'; % only need 1st listing
        tf = ~cellfun(@isempty, regexp(dac_codes(:,1),DC,'once'));
        if any(tf)
            info.dac = dac_codes{tf,2};
        else
            str = sprintf(['The data center responsible for %s could', ...
                'not be identified in the Sprof file'], info.fn);
            info.log{log_ct,1} = str; disp(str)
            log_ct      = log_ct+1;
            info.dac    = '';
        end

        % GET META FILE FROM THE GDAC & THEN GET NEEDED INFO IN META FILE
        % i.e. Cpactivation depth if it exists
        meta_fn = regexprep(info.fn,'Sprof','meta');
        meta_fps = fullfile(gdac_urls,'dac', info.dac, info.wmo, meta_fn);
        meta_fp = '';
        for ct = 1: size(meta_fps,1)
            try
                fprintf('getting meta file from: %s\n', meta_fps{ct})
                meta_fp = websave(fullfile(info.source,meta_fn),meta_fps{ct});
                break
            catch
                str = sprintf('Meta file grab failed from: %s', meta_fps{ct});
                info.log{log_ct,1} = str; disp(str)
                log_ct      = log_ct+1;
            end
        end

        if isfile(meta_fp)
            %meta_info = ncinfo(meta_fp);
            CMN = ncread(meta_fp,'CONFIG_MISSION_NUMBER');
            Config_PN = cellstr(ncread(meta_fp,'CONFIG_PARAMETER_NAME')'); % param names
            Config_PV = ncread(meta_fp,'CONFIG_PARAMETER_VALUE'); % Values
            tP        = strncmp(Config_PN,'CONFIG_CPAct',12); % find index for CpAct name
            if any(tP)
                % Cp Activation pressures for each CONFIG MISSION NUMBER
                % (cycle # not known) This would require match up with tech
                % file. "match" with Sprof at end
                info.CpActP = [CMN, Config_PV(tP,:)']; % mission #, CpAct P
                tf_CpAct    = 1;
            else
                str = sprintf('CONFIG_CPAct parameter not found in %s', meta_fn);
                info.log{log_ct,1} = str; disp(str)
                log_ct      = log_ct+1;
            end
            delete(meta_fp)
        else
            str = sprintf('Meta file not found at %s', meta_fp);
            info.log{log_ct,1} = str; disp(str)
            log_ct      = log_ct+1;
        end
        clear meta_fn meta_fps ct meta_fp meta_info CMN
        clear Config_PN Config_PV tP


        data    = ones(prod(Sprof_rc), size(hdr,2)) * NaN; % predim
        for vct = 1: size(keep_vars,1) % step through potential params on master list
            param = keep_vars{vct,1};

            if ~any(strcmp(param, avail_vars)) % skip if float doesn't have param
                continue
            
            % NEED TO GET PIVOT YEAR FOR FINAL SDN CALC
            elseif strcmp(param, 'JULD') % deal with time
                juld   = ncread(fp,'JULD');
                t1     = strcmp(nc_vars,'JULD');
                attr   = nc_info.Variables(t1).Attributes;
                t2     = strcmp({attr.Name},'units');
                dstr   = regexp(attr(t2).Value,'\d.+\d','match', 'once');
                ref_yr = datenum(dstr,'yyyy-mm-dd HH:MM:SS');
                SDN    = ones(Sprof_rc(1),1) * (ref_yr + juld)';
                clear juld t1 t2 attr dstr ref_yr
                data(:,strcmp(keep_vars{vct,3}, hdr)) = SDN(:);
                clear juld t1 t2 dstr ref_yr

            elseif ~isempty(regexp(param,'^CYC|^LAT|^LON','once'))
                X = ones(Sprof_rc(1),1) * ncread(fp,param)';
                data(:,strcmp(keep_vars{vct,3}, hdr)) = X(:);
                if strncmp(param,'LAT',3)
                    str = ncread(fp,'POSITION_QC');
                    str = regexprep(str',' ','7')'; % set ' ' to 7 tmp
                    QF  = sscanf(str,'%1f');
                    QF(QF==7) = NaN;
                    QF = ones(Sprof_rc(1),1) * QF';
                    % shift logical +1 for POS QC COL
                    data(:,circshift(strcmp(keep_vars{vct,3}, hdr),1)) = QF(:);
                    clear str QF
                end

            else
                tPARAM = strcmp(hdr, keep_vars{vct,3});
                X      = ncread(fp,param);
                XA     = ncread(fp,[param,'_ADJUSTED']);
                XQF    = ncread(fp,[param,'_QC']);
                XAQF   = ncread(fp,[param,'_ADJUSTED_QC']);
                
                % FOR PTS, if no ADJ replace with RAW for SAGE  
                if ~isempty(regexp(param,'^PRES|^PSAL|^TEMP','once')) 
                    tADJ = isnan(XA) & ~isnan(X); % logical matrix
                    if any(tADJ(:))
                        XA(tADJ) = X(tADJ); 
                    end
                end

                % ADD DATA TO OUTPUT MATRIX
                data(:, tPARAM) = X(:); % raw
                data(:, circshift(tPARAM, rraw*2)) = XA(:); % adj
                
                % DEAL WITH QC FLAGS (char arrays to num & space to NaN
                tmp = sscanf(regexprep(XQF(:)',' ','7'),'%1f'); 
                tmp(tmp==7) = NaN;
                data(:, circshift(tPARAM, 1)) = tmp; % raw QF

                tmp = sscanf(regexprep(XAQF(:)',' ','7'),'%1f'); 
                tmp(tmp==7) = NaN;
                data(:, circshift(tPARAM, rraw*2+1)) = tmp; % adj QF
            end 
            clear X XA XQF XAQF tADJ tmm
        end

    % *********************************************************************
    case 'ODV' % MBARI ONLY REALLY
        % GET WMO IF AVAILABLE
        wmo = regexp(fn,'\d{7}','match','once');
        if isempty(wmo)
            wmo = regexp(info.fn,'^NO_WMO.+(?=\.TXT)','match','once');
        end
        info.wmo = wmo;

        draw   = get_FloatViz_data(fp); % RAW
        draw.data(draw.data == -1e10) = NaN; % fill values to NaN

        adj_fp = fullfile([fd,'\QC'],[fn,'QC',ext]);
        dadj   = get_FloatViz_data(adj_fp); % ADJUSTED DATA
        if isempty(dadj)
            str = fprintf('NO adjusted data file was found at:%s',adj_fp);
            info.log{log_ct,1} = str; disp(str)
            log_ct      = log_ct+1;
        else
            dadj.data(dadj.data == -1e10) = NaN;
        end


        data    = ones(size(draw.data,1), size(hdr,2)) * NaN; % predim
        for vct = 1: size(keep_vars,1) % step through potential params on master list
            param = keep_vars{vct,2};
            tHDR  = strcmp(hdr, keep_vars{vct,3}); % param location output hdr
            tFVR  = strcmp(draw.hdr, param); % param location in raw FV hdr
            if ~isempty(dadj) 
            tFVA  = strcmp(dadj.hdr, param); % param location in adj FV hdr
            end
            

            if ~any(tFVR) % skip if float doesn't have raw param
                continue

            elseif ~isempty(regexp(param,'^SDN|^CYC|^LAT|^LON','once'))
                data(:,tHDR) = draw.data(:,tFVR);
                if strncmp(param,'Lat',3)
                    % shift logical +1 for POS QC COL
                    data(:,circshift(tHDR,1)) = draw.data(:,circshift(tFVR,1));
                    clear str QF
                end

            else % state variables
                if strcmp(param,'Oxygen[µmol/kg]') % set sensor existance flags
                    tf_Oxygen = 1;
                elseif strcmp(param,'Nitrate[µmol/kg]')
                    tf_Nitrate = 1;
                elseif strcmp(param,'pHinsitu[Total]')
                    tf_pH = 1;
                end

                data(:,tHDR) = draw.data(:,tFVR); % raw
                data(:,circshift(tHDR,1)) = draw.data(:,circshift(tFVR,1)); % raw QC

                % FOR PTS, if no ADJ replace with RAW for SAGE
                if isempty(dadj) && ~isempty(regexp(param,'^PRES|^PSAL|^TEMP','once'))
                    data(:,circshift(tHDR, rraw*2)) = draw.data(:,tFVR);
                    data(:,circshift(tHDR, rraw*2+1)) = 1; % ODV un-inspected/missing
                elseif isempty(dadj) % BGC adj params = NaN, QF = 1
                    data(:,circshift(tHDR, rraw*2+1)) = 1; % ODV un-inspected/missing
                else % adj data exist
                    data(:,circshift(tHDR, rraw*2)) = dadj.data(:,tFVA); % adj
                    data(:,circshift(tHDR,rraw*2+1)) = ...
                        dadj.data(:,circshift(tFVA,1)); % adj QC
                end
            end
        end

        % Flip cycle matrices to match Sprof Pdir (flip vs sort to preserve
        % order)
        tC = strcmp(hdr,'Cycle');
        uC = unique(data(:,tC));
        for ct = 1:size(uC,1)
            t1 = data(:,tC) == uC(ct);
            data(t1,:) = flip(data(t1,:));
        end

        % GET CpActP for MBARI MAT FILES
        mat_fd  = fullfile(mat_dir, wmo, filesep);
        tmp    = dir([mat_fd, wmo,'*mat']);
        if isempty(tmp)
            f = warndlg(['*.mat directory is empty! Please rebuild ',...
                'repository and try again'],'Warning');
            str = '*.mat cycle dir is empty!';
            info.log{log_ct,1} = str; disp(str)
            log_ct      = log_ct+1;
            return
        end




        fnames = {tmp.name}';
        CpActP = nan(size(fnames,1),2);
        fprintf('Getting CpAct P for each cycle ......\n');
        tf_CpAct  = 1; % should be there except for SOLO's
        for pct = 1:size(fnames,1)
            mat_fp = fullfile(mat_fd, fnames{pct});
            load(mat_fp,'INFO')

            if isfield(INFO,'CpActivationP')
                CpActP(pct,1) = INFO.cast;
                CpActP(pct,2) = INFO.CpActivationP;
            else
                str = sprintf('No CpActivation P found for float: %s',wmo);
                info.log{log_ct,1} = str; disp(str)
                log_ct      = log_ct+1;
                tf_CpAct    = 0; % reset to zero & stop looking
                break
            end
                
        end

        if tf_CpAct == 1
            info.CpActP = CpActP;
        end
        
        clear mat_fd tmp fnames CpActP mat_fp
        clear tC uC ct t1 tNO_NO3_pH tFVR tFVA tHDR draw dadj
end

% ************************************************************************
% ************************************************************************
% DO SOME FINAL CHECKS/TASKS BEFORE ASSIGNING OUTPUTS
% CLEAN UP DATA - FOR SAGE ONLY NEED DATA CO-LOCATED WITH NO3 & PH
% FOR ODV FLIP EACH CYCLE BY CYCLE SO PRESS DESCENDING TO MATCH
% SPROF PRES DIRECTION
% ************************************************************************
% ************************************************************************

% remove un-needed data including fill values from nc square matrix
tNO_NO3_pH = isnan(data(:, strcmp(hdr,'NO3'))) &  isnan(data(:, strcmp(hdr,'pH')));
data(tNO_NO3_pH,:) = [];
pct_removed = sum(tNO_NO3_pH)./ size(tNO_NO3_pH,1)*100;
if pct_removed == 100
    str = sprintf(['%0.1f percent of data were excluded for SAGE! ',...
        'Most likely no NO3 or pH sensors on board!'], ...
        sum(tNO_NO3_pH)./ size(tNO_NO3_pH,1)*100);
    info.log{log_ct,1} = str; disp(str)
    log_ct      = log_ct+1;

elseif pct_removed > 0
    str = sprintf(['%0.1f percent of data were excluded for SAGE! ',...
        '(non complementary NO3 or pH & *.nc fill)'], ...
        sum(tNO_NO3_pH)./ size(tNO_NO3_pH,1)*100);
    info.log{log_ct,1} = str; disp(str)
    log_ct      = log_ct+1;
end


% for now only use RAW PTS since that is what NO3 & pH are calculated from
% This can be changed in the future if need be
tf = ~cellfun(@isempty,regexp(hdr,'^P_adj|T_adj','once'));
tf = tf + circshift(tf,1); % get quality flag cols too
data = data(:,~tf);
hdr  = hdr(~tf);

% build Index structure
I.SDN = find(strcmp(hdr,'SDN') == 1);
I.STA = find(strcmp(hdr,'Cycle') == 1);
I.LON = find(strcmp(hdr,'Lon') == 1);
I.LAT = find(strcmp(hdr,'Lat') == 1);
I.P   = find(strcmp(hdr,'P') == 1);
I.S   = find(strcmp(hdr,'S') == 1);
I.SA   = find(strcmp(hdr,'S_adj') == 1);
I.T   = find(strcmp(hdr,'T') == 1);
I.O   = find(strcmp(hdr,'O2') == 1);
I.N   = find(strcmp(hdr,'NO3') == 1);
I.PH  = find(strcmp(hdr,'pH') == 1);
I.OA  = find(strcmp(hdr,'O2_adj') == 1); %adj
I.NA  = find(strcmp(hdr,'NO3_adj') == 1); %adj
I.PHA = find(strcmp(hdr,'pH_adj') == 1); %adj


% CHECK IF ALL NO3 or PH data are bad or missing from the start. If so call sensor bad
if tf_pH ==1
    tbad = isnan(data(:,I.PH)) | data(:,I.PH+1) == info.QFbad | ...
        data(:,I.PH+1) == info.QFmiss;
    if sum(~tbad) == 0 % no valid data
        tf_pH = 0;
        str = sprintf(['pH sensor exists but all data are bad so setting ',...
                       'sensor flag to bad\n']);
        info.log{log_ct,1} = str; disp(str)
        log_ct      = log_ct+1;

    end
end

if tf_Nitrate ==1
    tbad = isnan(data(:,I.N)) | data(:,I.N+1) == info.QFbad | ...
        data(:,I.N+1) == info.QFmiss;
    if sum(~tbad) == 0 % no valid data
        tf_Nitrate = 0;
        str = sprintf(['Nitrate sensor exists but all data are bad so ', ...
                       'setting sensor flag to bad\n']);
        info.log{log_ct,1} = str; disp(str)
        log_ct      = log_ct+1;

    end
end

if tf_Oxygen ==1 % Check adjusted Oxygen
    tbad = isnan(data(:,I.OA)) | data(:,I.OA+1) == info.QFbad | ...
        data(:,I.OA+1) == info.QFmiss;
    if sum(~tbad) == 0 % no valid data
        tf_Oxygen = 0;
        str = sprintf(['Oxygen sensor exists but all data are bad so ', ...
                       'setting sensor flag to bad\n']);
        info.log{log_ct,1} = str; disp(str)
        log_ct      = log_ct+1;

    end
end



% build track matrix [SDN Cycle Lon Lat ....]
[~,ia,~]        = unique(data(:,I.STA));
info.ftrack_hdr = [hdr([I.SDN,I.STA,I.LON,I.LAT]),'Pmin' 'Pmax'];
info.ftrack     = data(ia,[I.SDN,I.STA,I.LON,I.LAT]);
info.ftrack     = [info.ftrack nan(size(info.ftrack,1),2)];

% loop through cycles and get min and mix profile depths
for ct = 1:size(info.ftrack,1)
    t1 = data(:,I.STA) == info.ftrack(ct,2);
    info.ftrack(ct,5) = min(data(t1,I.P),[],1,'omitnan');
    info.ftrack(ct,6) = max(data(t1,I.P),[],1,'omitnan');
end


% Check CP Activation matrix size with float track & truncate to shortest
if tf_CpAct == 1 && size(info.ftrack,1) ~= size(info.CpActP,1)
    C = setxor(info.ftrack(:,2), info.CpActP(:,1)); % non commom cycle(s)
    str = sprintf(['WARNING: CpActP array size (%d) does not match ', ...
        'float track array size (%d) for %s. Truncating to common ', ...
        'cycles'],size(info.CpActP,1), size(info.ftrack,1), info.wmo);
    info.log{log_ct,1} = str; disp(str)
    log_ct      = log_ct+1;
    for ct = 1:numel(C) % check all 3 for ease of operation
        t1 = info.ftrack(:,2) == C(ct);
        info.ftrack(t1,:) = [];
        t2 = info.CpActP(:,1) == C(ct);
        info.CpActP(t2,:) = [];
        t3 = data(:,2) == C(ct);
        data(t3,:) = [];
    end
    %clear C str t1 t2 t3
end

% get some info for initial populating of GUI edit fields
info.CycleRange = [min(info.ftrack(:,2)), max(info.ftrack(:,2))];

ind = strcmp(info.ftrack_hdr,'Pmin');
info.PminRange  = [min(info.ftrack(:,ind), [], 1,'omitnan'), ...
                   max(info.ftrack(:,ind), [], 1,'omitnan')];

ind = strcmp(info.ftrack_hdr,'Pmax');
info.PmaxRange  = [min(info.ftrack(:,ind), [], 1,'omitnan'), ...
                   max(info.ftrack(:,ind), [], 1,'omitnan')];


% ***********************************************************************
% ALWAYS CALCULATE pH PUMP OFFSETS IF AVAILABLE - 1 VALUE PER CYCLE &
% ADJUSTED RAW pH FOR ALL FLAVORS FOR DATA SET
% THE GUI WILL DETERMINE IF THEY GET USED OR NOT

% Default values
pH_PO.tf     = 0; % No pH pump offset calculation
pH_PO.POhdr  = {}; % header for cycle offset values
pH_PO.POdata = []; % cycle offset values
pH_PO.hdr    = {}; % header for raw offset corrected values
pH_PO.data   = []; % raw offset corrected values same row count as float data

% only do so if valid pH data & CpAct values exist
if tf_pH && tf_CpAct %&& ~strcmp(pH_PO.offset_corr_type,'none')
    pH_PO.tf = 1; % offsets will get calculated

    % start with raw uncorr x2 & fill in with offset if available
    pH_PO.hdr  = [hdr([I.STA, I.P]),'pH PO none','pH PO poly',...
                  'pH PO lin','pH PO mixed'];
    corr_data  = data(:,[I.STA, I.P, I.PH, I.PH, I.PH, I.PH]); 
   
    po_data   = NaN(size(info.CpActP,1), 10);
    po_data(:,1:2) = info.CpActP; % cycle CpActP
    pH.PO.POhdr = {'cycle','CpActP'};
    tf_hdr    = 0 ; % switch to add hdr from offset calcs
    for ct = 1:size(info.CpActP,1)
        t1     = corr_data(:,1) == info.CpActP(ct,1); % cycle subset
        tmp    = data(t1,[I.P,I.PH, I.PH+1]); % PRES pH pHQC
        tmpT   = data(t1,I.T); % for pH TCOR
        CpActP = info.CpActP(ct,2);
        out    = calc_pH_pump_offset(tmp, CpActP, info.QFbad);


        % all bad data or other failure - no offset calculated
        % offset corrected float data = raw
        if isempty(out) 
            fprintf('No pH pump offset calculated for cycle %d\n', ...
                info.CpActP(ct,1));
            % po_data(ct,1:2) = info.CpActP(ct,:);
            continue
        end

        if tf_hdr == 0 % finish building header
            % out.hdr = {'poly offset' 'poly SSR' 'poly resid std' 'poly Zish', ...
            %     'lin offset' 'lin diff' 'lin_resid std' 'lin Zish'};
            pH_PO.POhdr = [pH.PO.POhdr, out.hdr]; % augmented header
            tf_hdr = 1;
        end
        po_data(ct,3:10) = out.data; % calculated offsets
        tP = tmp(:,1) <= CpActP; % pH pump offset upper cycle subset

        % NEED TO APPLY pH TCOR TO OFFSET CORRECTION TOO SINCE WE ARE IN pH
        % space but it is really a k0 issue
        tf_max = tmp(:,1) == max(tmp(tP,1));
        REF_T  = tmpT(tf_max);
        TCOR   = (REF_T + 273.15) ./ (tmpT + 273.15); % Ratio in Kelvin
        %TCOR   = ones(size(TCOR)); % TESTING




        % Sometimes the linear approach will work and poly will not & the
        % reverse. Check for NaN's & fill each corr column seperately
        if ~isnan(out.data(1,1))
            tmp1     = tmp(:,2); % raw pH
            tmp1(tP) = tmp1(tP) - TCOR(tP).*out.data(1,1); %poly PO corrected
            corr_data(t1,4) = tmp1;
        else
            fprintf(['Interp poly pH offset estimate returned NaN for ',...
                'cycle %d. Filled in with uncorrected raw\n'],info.CpActP(ct,1));
        end

        if ~isnan(out.data(1,5))
            tmp1     = tmp(:,2); % raw pH
            tmp1(tP) = tmp1(tP) - TCOR(tP).*out.data(1,5); %linear PO corrected
            corr_data(t1,5) = tmp1;
        else
            fprintf(['Linear pH offset estimate returned NaN for ',...
                'cycle %d. Filled in with uncorrected raw\n'],info.CpActP(ct,1));
        end

        if ~isnan(mean(out.data(1,[1,5]),'omitnan'))
            tmp1     = tmp(:,2); % raw pH
            tmp1(tP) = tmp1(tP) - TCOR(tP).*mean(out.data(1,[1,5]),'omitnan'); % mix - avg of both
            corr_data(t1,6) = tmp1;
        else
            fprintf(['Linear & poly pH offset estimates returned NaN for ',...
                'cycle %d. Filled in with uncorrected raw\n'],info.CpActP(ct,1));
        end
    end
    pH_PO.data   = corr_data; % sta, P,pump offset raw & corrected raw data
    pH_PO.POdata = po_data; % cycle offsets & stats
    clear tmp tmp1 tp po_data corr_data tf_hdr
end

% % Check
% ax(1) =subplot(1,2,1);
% plot(data(:,I.PH),data(:,I.P),'bo'), set(gca,'YDir','Reverse')
% 
% ax(2) =subplot(1,2,2);
% plot(corr_data, data(:,I.P),'r*'), set(gca,'YDir','Reverse')
% 
% x_lim = [min([data(:,I.PH);corr_data]), ...
%     max([data(:,I.PH);corr_data])];
% xlim(x_lim);



% **********************************************************************
% ASSIGN OUTPUTS
info.I          = I;
info.tf_Oxygen  = tf_Oxygen;
info.tf_Nitrate = tf_Nitrate;
info.tf_pH      = tf_pH;
info.tf_CpAct   = tf_CpAct;

d.hdr   = hdr;
d.data  = data;
d.info  = info;
d.pH_PO = pH_PO;

%clearvars -except d










 

