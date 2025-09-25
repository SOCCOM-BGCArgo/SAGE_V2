function d = get_CpActPRES(wmo, source, GDAC)

% THIS FUNCTION TRIES TO EXTRACT CP ACTIVATION PRESSURE FOR EACH PROFILE
% FROM GDAC FILES OR FROM MBARI *.MAT PROFILE FILES IF THE VALUES EXSIST

% INPUTS:
%  wmo    - 7 digit wmo number as a string
%  source - either 'GDAC' or 'MBARI' as a string
%  gdac   - 'US' or 'EU' as a string to choose GDAC
%     US = https://usgodae.org/pub/outgoing/argo/
%     EU = https://data-argo.ifremer.fr/
%
% OUTPUTS:
%   d.hdr  - 1 x 3 cell array {'Cycle number' 'Mission number' 'CpActPRESS'}
%   d.data - n x 3  matrix
%   d.info - meta data

% TESTING
%wmo     = '5904855'; %mbari aoml APEX
% wmo     = '7901028' % ifremer
% source  = 'GDAC';
% source = 'MBARI';
% GDAC   = 'US';

% Default output if function fails
d.hdr = {'Cycle number' 'Mission number' 'CpActPRESS'};
d.data = [];
d.info = [];

% Temporary place for *.nc files
dest_fd = fullfile(userpath,'TEMP\'); % for meta & prof files
if ~isfolder(dest_fd)
    mkdir(dest_fd);
end

info.WMO    = wmo;
info.source = source;
info.GDAC   = source;
info.tf_CpActPRES = 0; % true if data exist
info.msg    = '';


switch source

    case 'GDAC'
        US_https     = 'https://usgodae.org/pub/outgoing/argo/';
        EU_https     = 'https://data-argo.ifremer.fr/';
        if strcmp(GDAC,'US')
            https_target = US_https;
        else
            https_target = EU_https;
        end
        info.GDAC = https_target;

        % GET META FILE PATH FROM META FILE INDEX FILE
        fprintf('Finding meta file path .....\n');
        meta_url = [https_target,'ar_index_global_meta.txt'];
        str      = webread(meta_url);  %meta file list as big string
        filt_str = sprintf('\\w+/%s/%s_meta\\.nc',wmo,wmo);
        meta_fp  = regexp(str, filt_str,'match','once'); %extracted path
        if isempty(meta_fp)
            str = sprintf(['Meta file for WMO: %s not found - no data can be ',...
                'extracted!'],wmo);
            info.msg = str;
            fprintf('%s\n', str);
            return
        end

        [fd, fn, ext] = fileparts(meta_fp);

        % DOWNLOAD *META.NC & *PROF.NC FILES
        fprintf('Downloading meta file from the GDAC (%s): %s .....\n',...
            https_target, [fn,ext]);
        source_path  = [https_target,'dac/', meta_fp];
        dest_path    = [dest_fd, fn,ext];
        dest_meta_fp = websave(dest_path, source_path);

        fprintf('Downloading prof file from the GDAC: %s .....\n', ...
            [regexprep(fn, 'meta','prof'),ext]);
        dest_prof_fp = websave(regexprep(dest_path,'meta','prof'), ...
            regexprep(source_path,'meta','prof'));

        % GET CYCLE NUMBER & CONFIG MISSION NUMBER FROM PROF FILE
        prof_info   = ncinfo(dest_prof_fp); % you can look at available parameters here
        prof_cycles = ncread(dest_prof_fp ,'CYCLE_NUMBER'); % Values
        prof_cfg_MN = ncread(dest_prof_fp ,'CONFIG_MISSION_NUMBER'); % Values

        data = ones(size(prof_cycles,1),3) * NaN;
        data(:,1:2) = [prof_cycles, prof_cfg_MN];

        % GET CpAct P & CONFIG MISSION NUMBER FROM META FILE
        meta_info = ncinfo(dest_meta_fp); % you can look at available parameters here
        %meta_vars = {meta_info.Variables.Name}';
        meta_cfg_PN = cellstr(ncread(dest_meta_fp ,'CONFIG_PARAMETER_NAME')'); % param names
        meta_cfg_PV = ncread(dest_meta_fp ,'CONFIG_PARAMETER_VALUE'); % Values
        tP          = strncmp(meta_cfg_PN,'CONFIG_CPAct',12); % find index for CpAct name
        CpActP      = meta_cfg_PV(tP,:)'; % 157 x 1 Cp Activation pressures for each cycle (actual cycle # not known)
        %meta_cfg_MN = ncread(dest_meta_fp ,'CONFIG_MISSION_NUMBER'); % Values
        info.FloatType = strtrim(ncread(dest_meta_fp,'PLATFORM_TYPE')');

        % DATA EXIST? DIMMENSION CORRECT?
        if ~isempty(CpActP) && size(CpActP,1) == size(prof_cycles,1)
            data(:,3) = CpActP;
            info.tf_CpActPRES = 1;
        elseif ~isempty(CpActP)
            str = sprintf(['WARNING: Cycle number array size in %s does not ', ...
                'match  CpAct P array size in %s! Data will not be extracted'], ...
                [regexprep(fn,'meta','prof'), ext]  ,[fn, ext]);
            info.msg = str;
            fprintf('%s\n', str);
        end

        d.data = data;
        d.info = info;
        
        delete(dest_meta_fp)
        delete(dest_prof_fp)

        clearvars -except d


    case 'MBARI'
        fd = fullfile('\\atlas\Chem\ARGO_PROCESSING\DATA\FLOATS\',wmo,'\');
        if ~isfolder(fd)
            str = sprintf('Directory not found: %s\n',fd);
            info.msg = str;
            fprintf('%s\n', str);
            return
        end
        tmp     = dir([fd,wmo,'*mat']);
        fnames  = {tmp.name}'; %.mat profile file list
        rfnames = size(fnames,1);
        data    = ones(rfnames,3 ) * NaN;

        for ct = 1:rfnames
            fnp = fullfile(fd,fnames{ct});
            load(fnp,'INFO')
            data(ct,1) = INFO.cast;

            if ct == 1
                info.FloatType = INFO.float_type;
                info.mbari_id  = INFO.name;
            end

            if isfield(INFO,'CpActivationP') && ~isempty(INFO.CpActivationP)
                data(ct,3) = INFO.CpActivationP;
            end

        end
        if sum(~isnan(data(:,3))) > 0
            info.tf_CpActPRES = 1; % true if data exist
        end

        d.data = data;
        d.info = info;
        clearvars -except d
end


