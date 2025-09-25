function d = get_WOA23_local(data_source,track, depth_bnds, ocean_var)
%
% UPDATED FROM get_WOA_local.m TO APPLY TO WOA2023 DATASET.
%
% PURPOSE: 
%   Extract a subset from the monthly WOA2023 climatology for a given
%   variable and then interpolate data along the provided track. Access
%   different data sets by adjusting paths and entries in the WOA_info
%   cell array. My use is for profiling float data. Must Build 4D matrix
%   from 3D monthly variable files first.
%   [month x depth x lat x lon] => [time x depth x lat x lon]
%   [1 x m x n x p]             => [12 x m x n x p]
%
% NOTES: 
%   WOA nutrients are reported in units of umol/kg
%   Depths below monthly climatologies (800 for NO3,PO4,Si; 1500 for T,S,O
%      are filled with annual climatology values even though T & S have
%      seasonal values too
%
% USAGE:
%	data = get_WOA23_local(data_source,track, depth_bnds, ocean_var)
%
% INPUTS:
%   data_source = path to local WOA data repo
%   track      = n x 3 matrix [Matlab_SDN, Lat, Lon]
%   depth_bnds = depth bounds [min depth  max depth]
%	ocean_var  = a string, ocean parameter
%                must be one of these: T S O2 O2sat NO3 Si PO4 AOU
%
% OUTPUTS:
%	data =   a data structure.
%       data.d   = data matrix subset for WOA2018 variable [depth x time]
%       data.sd  = data matrix subset for WOA2018 variable sea [depth x time]
%       data.Z   = depth array subset
%              
% EXAMPLES:
%   data = get_WOA23_local(track, [0 1000], 'NO3')
%
% CHANGES
%    04/25/2020 JP carried along and interpolated the std dev of the
%    measurement to0
%    03/05/2024 JP updated to WOA 2023, removed nctoolbox depandancies,
%    interpolation bug fixes & enhancements

plot_it = 0; % 0 to turn off plotting

% ************************************************************************
% ************************************************************************
% TESTING 
% %data_source = fullfile(userpath, 'WOA2023\');
% data_source = fullfile(userpath, 'WOA2018\'); % TESTING OLD
% % % %fvn = '5904855.TXT'; % southern Ocean along Antarctic Just W of Drake Passage
% % % %fvn = '7901106.TXT'; %Hawaii
% %fvn = '5906531.TXT'; %NPac
% % % fvn = '4903591.TXT'; %NPac
% % fvn = '1902457.TXT';
% %fvn = '5905993.TXT'; % 12886
% fvn = '7900825.TXT'; % 21836
% fvd = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATVIZ\';
% fvp = fullfile(fvd,fvn);
% %d        = get_sageV2_float_data(fvp);
% d        = get_FloatViz_data(fvp);
% filt_str = '^SDN|^Lon|^Lat|';
% tf       = ~cellfun(@isempty, regexp(d.hdr, filt_str, 'once'));
% data     = d.data(:,tf);
% hdr      = d.hdr(tf);
% hdr      = hdr([1 3 2]); %lon, lat to lat,lon
% data     = data(:,[1 3 2]); %lon, lat to lat,lon
% [~,ia,~] = unique(data(:,1));
% track    = data(ia,:);


% depth_bnds = [0 2100];
% ocean_var  = 'NO3';
% depth_bnds = [0 2];
% ocean_var  = 'O2sat';

% data_source = fullfile(userpath, 'WOA2023\');
% fd   = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATVIZ\';
% fn   =  '1902644.TXT'; %
% %fn = '5906531.TXT'; %NPac
% fp   = fullfile(fd,fn);
% d    = get_sageV2_float_data(fp);
% I    = d.info.I; % get predefined fdata col idx struct
% idx  = [I.SDN,I.STA,I.LON,I.LAT,I.P,I.S,I.T,I.OA];
% data = d.data(:,idx);
% hdr  = d.hdr(idx);
% [uCycles,ia,~] = unique(data(:,2));
% track = data(ia,[1,4,3]);
% depth_bnds = [0 2100];
% ocean_var  = 'NO3';
% 
% *************************************************************************
% ************************************************************************

%      
% ***********************************************************************
% **********************       WOA 2023 INFO       **********************
% WOA2018 file name format:  woaYY_[DECA]_[v][tp][ft][gr].[form_end]
%    [v] = oceanographic variable (t s o n i p)
%           t = temperature
%           s = salinity
%           o = oxygen
%           O = % Oxygen saturation
%           n = nitrate
%           i = silicate
%           p = phosphate
%    [tp]= time averaging period
%           00 – annual statistics, all data used;
%           01 to 12 – monthly statistics (01 = Jan, 12 – Dec);
%           13 to 16 – seasonal statistics:
%               13 – North Hemisphere winter (January - March);
%               14 – North Hemisphere spring (April - June);
%               15 – North Hemisphere summer (July - September);
%               16 – North Hemisphere autumn (October - December);
%    [ft] = field type = an = Objectively analyzed climatological mean
%    [ft] = field type = sd = std dev of statistical mean 
%    [gr] = the grid size (01 – 1-degree grid resolution this use)
%     YY  = WOA year: 23 18 13
%
%  NOTE: Monthly nitrate, phosphate, silica only to 800m
%        Monthly oxygen only to 1500m
%        
% dimensions = time x depth x lat x lon
% lat -89.5 to 89.5
% lon -179.5 to +179.5

% ***********************************************************************
% SET SOME VARIABLES, PATHS, STRINGS, TEMPLATES
% ***********************************************************************
if regexp(data_source,'2018|2013','once')
    fn_str = 'woa18_%s_%s%02d_01.nc'; % decade, param, month for monthly TESTING
    tf2023 = 0;
elseif regexp(data_source,'2023','once')
    tf2023 = 1;
    fn_str = 'woa23_%s_%s%02d_01.nc';
end

%O2_volume  =  22.3916; % L/ mole O2 @ STP

% VARIABLES EXTRACTION TABLES
if tf2023 == 1
    WOA_info(1,:) = {'nitrate'     'n',  'all'   'NO3'};
    WOA_info(2,:) = {'oxygen'      'o'   'all'   'O2' };
    WOA_info(3,:) = {'phosphate'   'p'   'all'   'PO4'};
    WOA_info(4,:) = {'silicate'    'i'   'all'   'Si' };
    WOA_info(5,:) = {'temperature' 't'   'B5C2'  'T' };
    WOA_info(6,:) = {'salinity'    's'   'B5C2'  'S' };
    WOA_info(7,:) = {'o2sat'       'O'   'all'   'O2sat'};
    WOA_info(8,:) = {'AOU'         'A'   'all'   'AOU'};
else
    WOA_info(1,:) = {'nitrate'    , 'n', 'all', 'NO3'};
    WOA_info(2,:) = {'oxygen'     , 'o', 'all', 'O2' };
    WOA_info(3,:) = {'phosphate'  , 'p', 'all', 'PO4'};
    WOA_info(4,:) = {'silicate'   , 'i', 'all', 'Si' };
    WOA_info(5,:) = {'temperature', 't', 'A5B7', 'T' };
    WOA_info(6,:) = {'salinity'   , 's', 'A5B7', 'S' };
    WOA_info(7,:) = {'o2sat'      , 'O', 'all', 'O2sat'};
    WOA_info(8,:) = {'AOU'        , 'A', 'all', 'AOU'};
end

tf_var = strcmp(WOA_info(:,4), ocean_var);

% month_label = {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' ...
%                'Oct' 'Nov' 'Dec'};

% ***********************************************************************
% CHECK INPUTS
% ***********************************************************************
s1 = ['ocean_var must be a character string from this list: ', ...
      'T S O2 O2sat NO3 Si PO4 AOU'];  
s2 = 'Check track input: must have 3 rows or columns [SDN Lat Lon]';
s3 = 'Check depth range input: should be [min_depth max_depth]';
pos_fixes = track;

[r,c] = size(pos_fixes);
if r ~= 3 && c ~= 3
    disp(s2)
    return
elseif length(depth_bnds) ~= 2 || depth_bnds(1) > depth_bnds(2)
    disp(s3)
    return
% not a sting or variable not in the designated list
elseif ~ischar(ocean_var) || ~any(tf_var)
    disp(s1)
    return
end
if r == 3 && c ~= 3
    pos_fixes = pos_fixes'; % now  n x 3
end

% SDN LAT LON
% SET ANY MISSING VALUES (ODV FILL VALUES) IN LAT OR LON TO NaN
t1 = pos_fixes(:,2) == -1e10; 
pos_fixes(t1,2:3) = NaN;
clear d t1

% CONVERT LON TO -180 to +180 IF NECESSARY
t1 = pos_fixes(:,3) > 180; % WOA2018 -180 + 180
pos_fixes(:,3) = pos_fixes(:,3) - (t1*360);

% QUICK FIX IF pos_fixes DATA IS NOT IN WOA GRID [-179.5 to 179.5] 9/12/13
% THIS SHOULD BOUND EVERY POSITION IN LONGITUDE. NO FLOATS NEAR LAT MAXES
% BUT PLAY THE SAME GAME. IN THESE SITUATATIONS IT'S A CRUDE NEAREST
% NEIGHBOR APPROACH VS. LINEAR INTERPOLATION
pos_fixes(pos_fixes(:,3) <= -179.5, 3) = -179.49;
pos_fixes(pos_fixes(:,3) >=  179.5, 3) = 179.49;
pos_fixes(pos_fixes(:,2) <= -89.5, 2)  = -89.49;
pos_fixes(pos_fixes(:,2) >=  89.5, 2)  = 89.49; 

% TEST FOR MERIDIAN CROSSING. CAN USE TO DECREASE SIZE OF SUBSET GRAB IF
% NEEDED FOR PROCESSING SPEED. NOT USED FOR NOW
cross180 = 0;
tPOS = any(sign(pos_fixes(:,3)) >= 0);
tNEG = any(sign(pos_fixes(:,3)) < 0);
if tPOS && tNEG
    fprintf('+180 / -180 crossing detected!\n');
    cross180 = 1;
end

% float track bounds
lon_bnds = [min(pos_fixes(:,3)) max(pos_fixes(:,3))]; % lon
lat_bnds = [min(pos_fixes(:,2)) max(pos_fixes(:,2))]; % lat 

clear s1 s2 r c t1

% ***********************************************************************
%   EXTRACT WOA2023 MONTHLY CLIMATOLOGY SUBSETS & BUILD ANNUAL MATRIX (4D)
%       FROM INDIVIDUAL MONTHS: [LON x LAT x DEPTH x TIME]
%             ncread appears to squeeze singlton time dimension
% ***********************************************************************

% LOCAL PATH
data_path = fullfile(data_source, WOA_info{tf_var,1}, filesep);

var_name = sprintf('%s_an', WOA_info{tf_var,2});
sea_name  = regexprep(var_name,'_an','_sea'); % standard error of analysis

annual_flag = 0;
month_error_str = 'The requested %s range is not in the grid domain [%0.1f %0.1f]\n';
% *************************************************************************
% GET NetCDF INDEXING BEFORE EXTRACTING MONTHLY & ANNUAL CLIMATOLOGIES

% **** MONTHLY INDICES ****
fname = sprintf(fn_str, WOA_info{tf_var,3}, WOA_info{tf_var,2}, 1);
fp    = fullfile(data_path, fname);

lon = ncread(fp,'lon');
lat = ncread(fp,'lat');
z   = ncread(fp,'depth');

lon_range = [min(lon) max(lon)]; % WOA bounds
lat_range = [min(lat) max(lat)];
z_range   = [min(z) max(z)];

tZ = z >= depth_bnds(1) & z <= depth_bnds(2); % float subset index
if ~any(tZ)
    fprintf(month_error_str, var_name, z_range);
    return
end

% float subset index
tLON = lon >= lon_bnds(1)-1 & lon <= lon_bnds(2)+1; % add +/- 1 to insure bounding
if ~any(tLON)
    fprintf(month_error_str, var_name, lon_range)
    return
end

% float subset index
tLAT = lat >= lat_bnds(1)-1 & lat <= lat_bnds(2)+1;
if ~any(tLON)
    fprintf(month_error_str, var_name, lon_range)
    return
end

lon_inds  = [find(tLON == 1,1,'first'), find(tLON == 1,1,'last')]; % monthly
lat_inds  = [find(tLAT == 1,1,'first'), find(tLAT == 1,1,'last')];
z_inds    = [find(tZ == 1,1,'first'), find(tZ == 1,1,'last')];
start_idx = [lon_inds(1) lat_inds(1) z_inds(1) 1]; % extraction start vector
count_idx = [diff(lon_inds)+1 diff(lat_inds)+1 diff(z_inds)+1 1]; % extraction count vector

lon = lon(tLON); lat = lat(tLAT); Z   = z(tZ);

% **** DEEP ANNUAL ****
if depth_bnds(2) > max(z) % need to fill in with annual data at depth so get it
    annual_flag = 1;
    dfname    = sprintf(fn_str, WOA_info{tf_var,3}, WOA_info{tf_var,2}, 0);
    fp       = fullfile(data_path, dfname);
    zd       = ncread(fp,'depth'); % annual clim Z
    zd_range = [min(zd) max(zd)];
    tZD = zd > max(z) & zd <= depth_bnds(2);

    if any(tZ)
        zd_inds   = [find(tZD == 1,1,'first'), find(tZD == 1,1,'last')];
        dstart_idx = [lon_inds(1) lat_inds(1) zd_inds(1) 1];
        dcount_idx = [diff(lon_inds)+1 diff(lat_inds)+1 diff(zd_inds)+1 1]; %  # of indices to grab
        d_an  = ncread(fp, var_name, dstart_idx, dcount_idx);  % objectively analyzed climatology

        if tf2023 % not in WOA2018
            d_sea = ncread(fp, sea_name, dstart_idx, dcount_idx); % Standard error of the analysis
        end

        Z = [z; zd(tZD)];
    else
        fprintf(month_error_str, var_name, z_range);
        annual_flag = 0;
    end 
end

tSURF = Z <= max(z(tZ)); % monthly depth indices
tDEEP = ~tSURF;

AN  = nan([count_idx(1:2), size(Z,1), 14]); % predim: Make 4D add time dimmension +2 for bookends
if tf2023
    SEA = nan([count_idx(1:2), size(Z,1), 14]); % predim: Make 4D add time dimmension +2 for bookends
end

% STEP THROUGH MONTHLY CLIMATOLOGIES
% Add ends in time once 4D built to catch days < 15 & > 349
t = datenum(2022,1:12,15) - datenum(2022,1,0); % year day array for mid months 


for fct = 1:12 
    %fname = sprintf(fn_str, WOA_info{tf_var,3}, WOA_info{tf_var,2}, fct+3)
    fname = sprintf(fn_str, WOA_info{tf_var,3}, WOA_info{tf_var,2}, fct);
    fp    = fullfile(data_path, fname);

    % Objectively analyzed climatology lon x lat x depth x time
    % ncread squeezes singelton dimmensions
    an  = ncread(fp, var_name, start_idx, count_idx);  % objectively analyzed climatology


    AN(:,:, tSURF,fct+1) = an; % add to time dimmension, +1 to leave room for bookends

    if tf2023
        sea  = ncread(fp, sea_name, start_idx, count_idx); % Standard error of the analysis
        SEA(:,:,tSURF,fct+1) = sea;
    end

    if annual_flag == 1
        AN(:,:,tDEEP,fct+1) = d_an;

        if tf2023
            SEA(:,:,tDEEP,fct+1) = d_sea;
        end
    end
end

% CREATE BOOKENDS IN TIME SO YOU CAN INTERPOLATE PAST GRID BOUNDS
% Year days go from 15 to 349 so bound AN ends in time dimmension for
% easier interpolation. mid month, start month or end month??

% TIME
t = [0, t, 366]'; % transpose to column so  "topkrows" can be used
AN_ends = mean(AN(:,:,:,[2,13]),4); % avg of months 1 & 12
AN(:,:,:,1) = AN_ends;
AN(:,:,:,14) = AN_ends;

if tf2023
    SEA_ends = mean(SEA(:,:,:,[2,13]),4); % avg of months 1 & 12
    SEA(:,:,:,1)  = SEA_ends;
    SEA(:,:,:,14) = SEA_ends;
end

% clear AN_ends SEA_ends an sea fp fct tSURF tDEEP d_an d_sea z zd tZD tZ
% clear AN_ends SEA_ends an sea fp fct tSURF tDEEP d_an d_sea z zd tZD tZ
% clear dstart_idx dcount_idx start_idx count_idx lon_inds lat_inds z_inds
% clear dz_inds tLAT tLON tNEG tPOS zd_inds


% ***********************************************************************
% ***********************************************************************
% DATA MERGED & SUBSETTED! NOW DO THE INTERPOLATION AT POSITION FIXES
% track
dvec   = datevec(pos_fixes(:,1));
yr_day = pos_fixes(:,1) - datenum(dvec(:,1),1,0);
ANi   = nan(size(Z,1), size(pos_fixes,1));
SEAi  = ANi;

% predim array count vectors for logical tests. These are the size of
% the netcdf subset matirix
lon_cts = (1:size(lon,1))'; 
lat_cts = (1:size(lat,1))';
day_cts = (1:size(t,1))';

for pct = 1 :size(pos_fixes,1)
    % GET SUBSET INDICES FOR PROFILE(pct)
    aLON     = lon - pos_fixes(pct,3); % LON
    [~,idx]  = topkrows(abs(aLON),2,'ascend'); % bounding anomalies
    [tLON,~] = ismember(lon_cts, idx); % bounding profile logical array
    wt_lon  = (max(lon(tLON)) - pos_fixes(pct,3)) ./ diff(lon(tLON));

    aLAT     = lat - pos_fixes(pct,2); % LAT
    [~,idx]  = topkrows(abs(aLAT),2,'ascend'); % 
    [tLAT,~] = ismember(lat_cts, idx);
    wt_lat   = (max(lat(tLAT)) - pos_fixes(pct,2)) ./ diff(lat(tLAT));

    aDay    = t - yr_day(pct); % TIME
    [~,idx] = topkrows(abs(aDay),2,'ascend'); %
    [tT,~]  = ismember(day_cts, idx);
    wt_day  = (max(t(tT)) - yr_day(pct)) ./ diff(t(tT));

    % STEP DOWN THROUGH DIMENSIONS - I'M SURE THERE IS A BETTER WAY BUT
    % NEED TO ALSO DEAL WITH NaN's IN WOA PROFILE DATA DUE TO BATHYMETRY OR
    % LAND MASKS
    tmpA = AN(tLON ,tLAT,:,tT);

    if tf2023
        tmpS = SEA(tLON ,tLAT,:,tT);
    end
   
    % Matrix sum excluding NaN's yields 0 if all elements NaN - need to check
    % & reset to NaN

    % ********************************************************************
    % START WITH TIME DIMMENSION
    wt_day1 = ones(size(tmpA(:,:,:,1)))*wt_day;
    wt_day2 = 1 - wt_day1;
    tnan    = isnan(tmpA);
    wt_day1(~tnan(:,:,:,1) & tnan(:,:,:,2)) = 1; % wt back to 100% for side with value
    wt_day2(~tnan(:,:,:,2) & tnan(:,:,:,1)) = 1; % wt back to 100% for side with value
    
    tmpT  = sum(cat(4,tmpA(:,:,:,1) .* wt_day1, tmpA(:,:,:,2) .* wt_day2),4,'omitnan'); %AN
    tZero = sum(tmpA == 0,4); % check if any true zero values to begin with
    tmpT(~tZero & tmpT == 0) = NaN; % reset any "new zero's" back to NaN

    if tf2023
        tmpTS  = sum(cat(4,tmpS(:,:,:,1) .* wt_day1, tmpS(:,:,:,2) .* wt_day2),4,'omitnan'); %SEA
        tmpTS(~tZero & tmpT == 0) = NaN; % reset any "new zero's" back to NaN
    end

    % ********************************************************************
    % LONGITUDE NEXT
    wt_lon1 = ones(size(tmpT(1,:,:)))*wt_lon;
    wt_lon2 = 1 - wt_lon1;
    tnan    = isnan(tmpT);
    wt_lon1(~tnan(1,:,:) & tnan(2,:,:)) = 1; % wt back to 100% for side with value
    wt_lon2(~tnan(2,:,:) & tnan(1,:,:)) = 1; % wt back to 100% for side with value

    tmpLO = squeeze(sum(cat(1,tmpT(1,:,:) .* wt_lon1, tmpT(2,:,:) .* wt_lon2),1,'omitnan'));
    tZero = squeeze(sum(tmpT == 0,1)); % check if any true zero values to begin with
    tmpLO(~tZero & tmpLO == 0) = NaN; % reset any "new zero's" back to NaN

    if tf2023
        tmpLOS = squeeze(sum(cat(1,tmpTS(1,:,:) .* wt_lon1, tmpTS(2,:,:) .* wt_lon2),1,'omitnan'));
        tmpLOS(~tZero & tmpLO == 0) = NaN; % reset any "new zero's" back to NaN
    end

    % ********************************************************************
    % LATITUDE LAST
    % if only 1 depth level original 'an; gets squeezed so 1 x n instead of
    % n x1 so transpose
    if sum(tSURF) == 1 
        tmpLO = tmpLO';
    end
    wt_lat1 = ones(size(tmpLO(1,:)))*wt_lat;
    wt_lat2 = 1 - wt_lat1;
    tnan    = isnan(tmpLO);
    wt_lat1(~tnan(1,:) & tnan(2,:)) = 1; % wt back to 100% for side with value
    wt_lat2(~tnan(2,:) & tnan(1,:)) = 1; % wt back to 100% for side with value

    %ANi(:,pct) = squeeze(tmpLO(1,:) * wt_lat + tmpLO(2,:) *(1-wt_lat));
    tmpLA = sum(cat(1,tmpLO(1,:) .* wt_lat1, tmpLO(2,:) .* wt_lat2),1,'omitnan');
    tZero = squeeze(sum(tmpLO == 0,1)); % check if any true zero values to begin with
    tmpLA(~tZero & tmpLA == 0) = NaN; % reset any "new zero's" back to NaN
    %clear tmpA tmpT tmpLO
    ANi(:,pct) = tmpLA;

    if tf2023
        tmpLAS = sum(cat(1,tmpLOS(1,:) .* wt_lat1, tmpLOS(2,:) .* wt_lat2),1,'omitnan');
        tmpLAS(~tZero & tmpLA == 0) = NaN; % reset any "new zero's" back to NaN
        SEAi(:,pct) = tmpLAS;
    end

%    tZZ = Z > 700 & Z <850; % TESTING TO LOOK AT ODD WOA DATA
%   % tZZ = Z > 300 & Z <500; % TESTING TO LOOK AT ODD WOA DATA
%    if any(ANi(tZZ,pct) < 40)
%     %if any(isnan(ANi(:,pct)))
%         disp('hello')
%         sz = size(tmpA);
%         jp = reshape(permute(tmpA, [3,1,2,4]), sz(3), 2^3);
%         figure
%         plot(jp, Z*ones(1,size(jp,2))), set(gca,'YDir','Reverse')
%         pause
%     end


end

% FUNCTION OUTPUT
d.Z = Z;
d.d = ANi;   % objectively analyzed climatology
if tf2023
    d.sea = SEAi;  % Standard error of the analysis
end

% CHECK THE RESULTS
if plot_it == 1
    F10 = figure(10);
    F10.Units = 'inches';
    F10.Position = [2,2,10,7];

    WOA_YR_str = regexp(data_source,'WOA\d+','once','match');

    ax(1) = subplot(2,2,1);
    plot(ANi,Z*ones(1,size(ANi,2)),'o')
    title([WOA_YR_str,' - new interp function'])
    xlabel(ocean_var)
    ax(1).YDir = 'Reverse';


    % Compare to exisiting function (using WOA2023)
    data = get_WOA2013_local(data_source,track, depth_bnds, ocean_var);

    ax(2) = subplot(2,2,2);
    plot(data.d, data.Z*ones(1,size(data.d,2)),'o')
    title([WOA_YR_str,' - old interp function'])
    xlabel(ocean_var)
    ax(2).YDir = 'Reverse';
    ax(2).XLim = ax(1).XLim;

    ax(3) = subplot(2,2,3);
    %ax(3).Position = ax(3).Position + [0.15 0 -0.3 0];

    h = earthimage('watercolor',[224 224 224]/256);
    tstr = sprintf('Platform track');
    title(tstr)

    hold on
    dvec = datevec(pos_fixes(:,1));
    dec_year = dvec(:,1) + (pos_fixes(:,1) - datenum(dvec(:,1),1,0))/365;

    tt = datetime(pos_fixes(:,1),'ConvertFrom','datenum');
    scatter(pos_fixes(:,3), pos_fixes(:,2), 36, dec_year,'filled')
    colorbar
    ax(3).FontSize = 14;
    hold off
 
    ax(4) = subplot(2,2,4);
    plot(tt, ANi(1,:), 'bo', tt, data.d(1,:),'r*');
    ylabel('Shallowest sample');
    legend('New', 'Old')

end

clearvars -except d % don't need d anymore   

