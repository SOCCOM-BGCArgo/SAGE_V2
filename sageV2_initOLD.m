function d = sageV2_init()
%
% This function is triggered when SAGEv2 is first run in startupFcn(app).
% Most non float specifc, user configurable information such as
% directories, default file type, BIC error limits. The existance of needed
% folders based on those defined in dirs is tested. If at MBARI, additional
% MBARI specifc info is loaded

% ************************************************************************
%                  SET DIRECTORIES & ADD TO PATH
% ************************************************************************
%data_default = 'SPROF'; % 'ODV' or 'SPROF'
data_default  = 'ODV'; % 'ODV' or 'SPROF'
user_path     = userpath; % i.e. 'C:\Users\jplant\Documents\MATLAB'

d.dirs.sage      = fullfile(user_path,'Apps\SAGE_V2\');
d.dirs.esper     = fullfile(user_path,'ARGO_PROCESSING\MFILES\ESPER\');
d.dirs.caynb     = fullfile(user_path,'ARGO_PROCESSING\MFILES\CANYON_B\');
d.dirs.WOA       = fullfile(user_path,'ARGO_PROCESSING\DATA\WOA2023\');
d.dirs.glodap    = fullfile(user_path,'ARGO_PROCESSING\DATA\GLODAPv2\');
d.dirs.glodap_fn = 'GLODAPv2.2023_Merged_Master_File.mat';
d.dirs.QClists   = fullfile(user_path,'ARGO_PROCESSING\DATA\CAL\QC_LISTS\'); % user needs to define this
d.dirs.QClog     = fullfile(user_path,'ARGO_PROCESSING\DATA\CAL\');            % user needs to define this
d.dirs.QClog_fn  = 'FloatQCList_log.txt'; % this should be in QC lists folder I think not in cal
d.dirs.sprof     = fullfile(user_path,'ARGO\DATA\SPROF\'); % MBARI/OUTSIDE WORLD, user defined

% MBARI SPECIFC
d.dirs.FV        = fullfile(user_path,'ARGO_PROCESSING\DATA\FLOATVIZ\'); % MBARI
d.dirs.FVlocal   = fullfile(user_path,'ARGO_PROCESSING\DATA\FLOATVIZ\'); % MBARI -legacy for argo2ODV_LIAR
d.dirs.cal       = fullfile(user_path,'ARGO_PROCESSING\DATA\CAL\'); % MBARI
d.dirs.mat       = fullfile(user_path,'ARGO_PROCESSING\DATA\FLOATS\'); % mbari CpActP
d.dirs.bottle    = fullfile(user_path,'ARGO_PROCESSING\DATA\SHIPBOARD\'); % MBARI
d.dirs.bottle_fn = 'BottleData_lookup_table.txt'; % MBARI
d.dirs.BadSensor_fn = 'bad_sensor_list.txt'; % MBARI
d.dirs.BadSample_fn = 'bad_sample_list.txt'; % MBARI
d.dirs.temp        = fullfile('C:\temp\');     % MBARI

% SET DEFAULT EDIT BOX RANGES (CYCLE, PRESURE RANGES, ETC)
d.defaults.CycleMin = 1;
d.defaults.CycleMax = 500;
d.defaults.Pmin     = 1280;
d.defaults.Pmax     = 1520;
d.defaults.XOver_km = 30;
d.defaults.QCATable = [];

% SET DEFAULT FLOAT DATA PATH HERE
if strcmp(data_default,'ODV')
    d.dirs.float     = d.dirs.FV; %default = ODV files
    %d.dirs.QClists   = fullfile(d.dirs.cal,'QC_LISTS'); % CQ LIST FILES
elseif strcmp(data_default,'SPROF')
    d.dirs.float     = d.dirs.sprof; %default = sprof files
    %d.dirs.QClists = fullfile(d.dirs.sprof,'QC_LISTS'); % CQ LIST FILES
else
    fprintf(['Default data source (%s) not recognized ', ...
        '(should be ODV or SPROF)\n'], data_default);
    return
end

% TEST FOLDER EXISTANCE
PathStr = path; % long string of all file paths in the Matlab search path
flds    = fieldnames(d.dirs);
% remove any file names from the list - make file name end in "fn"
tg      = cellfun(@isempty, regexp(flds,'fn$','once'));
flds    = flds(tg);
for ct = 1:size(flds) % step through saveV2 paths
    fd = d.dirs.(flds{ct});
    if ~isfolder(fd) % folder not found
        fprintf('WARNING: folder not found: %s\n', fd); 
        answer = questdlg({'Would you like to create:', fd}, ...
            'FOLDER NOT FOUND!','No');      
        if strcmp(answer,'Yes') % Handle response
            mkdir(fd)
        else
            continue
        end
    else
        tf_path  = contains(PathStr, fd, 'IgnoreCase', ispc);
        if ~tf_path
            addpath(fd)
            fprintf('%s was added to the Matlab path\n', fd)
        end
    end
end


% ************************************************************************
% ADD ADDITIONAL CODE HERE TO SET DEFAULT VALUES CHECK/POPULATE MISSING REF DATA ETC
% d.error_limits.Nitrate = 0.5; % For BIC calculations to limimit over fits
% d.error_limits.pH      = 0.005;

d.error_limits.Nitrate = 0.3; % TESTING these are the SAGE limits
d.error_limits.pH      = 0.004;

% Set significance level for calc adjustment slope segment different from
% zero. set a low bar otherwise will force slope to zero more than desired
d.sig_level  = 0.25; %75%,  normally 0.05 for 95%          


% ************************************************************************


% If appropriate load MBARI master list & bottle data look up table
d.mlist = [];
d.bottle_lookup = {};
mbari_fp = fullfile(d.dirs.cal,'MBARI_float_list.mat');
if isfile(mbari_fp)
    mbari = load(mbari_fp);
    d.mlist = mbari.d;

    % Next try & load bottle lookup table
    if isfile(fullfile(d.dirs.bottle, d.dirs.bottle_fn))
        fid = fopen(fullfile(d.dirs.bottle, d.dirs.bottle_fn));
        tmp   = textscan(fid, '%s %s %s %s %f %f %s',...
            'Delimiter', '\t','HeaderLines',1);
        fclose(fid);
        d.bottle_lookup = ...
            [tmp{1,1}, tmp{1,2}, tmp{1,3}, tmp{1,4}, num2cell(tmp{1,5}), ...
            num2cell(tmp{1,6}), tmp{1,7}];
       % Check for floats without asoc bottle files & remove
        tf = ~cellfun(@isempty, d.bottle_lookup(:,7));
        if any(tf)
            fprintf(['%d floats were removed from the bottle lookup ', ...
                'table list because no bottle files exist for them\n'], ...
                sum(~tf));
            d.bottle_lookup = d.bottle_lookup(tf,:);
        end
        fprintf('Bottle lookup table loaded.\n')
    else
        fprintf('Bottle lookup file not found.\n')
    end
else
    fprintf('MBARI float list not found.\n')
end

clearvars -except d





