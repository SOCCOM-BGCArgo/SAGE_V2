function d = get_sageV2_ref_data(data, tfQF, Refs, WOA_path, app)
% This function estimates various static & dynamic NO3 & pH climatologies
% for the given data which are used as reference data to correct the raw
% float data with
%
% INPUTS:
%    data = a matix of float data
%           [MatlabSDN, Cycle, Lon, Lat, Pres, Salt, Temp, Oxygen] must be
%           n x 8
%           Time   = GMT
%           Pres   = dbar
%           Oxygen = µmol/kg & must be corrected 1st
%                    If oxygen is all NAN only ESPER NO O2 & WOA will have
%                    data, Bad oxygen data should be set to NaN first
%
%   tfQF = logical matrix, nx4 [S, O2, NO3, pH], where n = size(data,1) and 
%          n= 0 if param is missing, NaN or bad. Used to speed up ESPER &
%          CANYON B calcs, mostly for solo or Navis
%
%   Refs = a string or cell array of strings identifying the climatology to use
%          'all' - all available climatologies (CANYONB, ESPER_NN, ESPER_LIR,
%              ESPER_MIX, WOA ESPER_MIX_NO_02)
%          'CANYONB'
%          {'CANYONB' 'ESPER_NN'}
%
%   WOA_path = [] for default path or top level WOA data directory. 
%              i.e 'C:\Users\jplant\Documents\MATLAB\WOA2023\'
%
%   app -  This can either be empty or app. Used to send staus mesages to
%          GUI - probaly a better way
%
%   NOTE: ONLY GOOD QUALITY INPUT DATA WILL YIELD GOOD OUTPUTS!
%         Quality flag checks should be performed on inputs either before
%         or after the function is invoked.
%
% OUTPUTS:
%   d = a structure with fields for each Ref computation

% **********************************************************************
%TESTING
% fd = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATVIZ\';
% %fn =  '5904855.TXT'; % Southern Ocean along antarctic coast
% %fn =  '4903591.TXT'; % NEAR OSP
% %fn =  '5904685.TXT'; %
% fn =  'NO_WMO_un0948.TXT'; %
% fp = fullfile(fd,fn);
% d        = get_sageV2_float_data(fp);
% 
% I    = d.info.I; % get predefined fdata col idx struct
% idx  = [I.SDN,I.STA,I.LON,I.LAT,I.P,I.S,I.T,I.OA];
% data = d.data(:,idx);
% hdr  = d.hdr(idx);
% WOA_path = [];
% app      = [];
% 
% % filt_str = '^SDN|^Cycle|^Lon|^Lat|^P_adj|^S_adj|^T_adj|^O2_adj';
% % tf       = ~cellfun(@isempty, regexp(d.hdr, filt_str, 'once'));
% % % hdr      = d.hdr(tf);
% % data     = d.data(:,tf);
% % clear d tf file_str
% % WOA_path = [];
% 
% 
% %Refs = 'ESPER_NN';
% %Refs = 'CANYONB';
% %Refs = {'CANYONB' 'ESPER_NN'};
% Refs = 'all';
% %Refs = 'WOA';

% **********************************************************************

% ************************************************************************
% CHECK & SET INPUTS

% function called without app structure input - no message will be returned
% to the GUI
tf_app = 1;
if isempty(app) 
    tf_app = 0;
end

d =[]; % default output
default_WOA = fullfile('C:\Users\jplant\Documents\MATLAB\WOA2023');
%default_WOA = fullfile('C:\Users\jplant\Documents\MATLAB\WOA2018');

% THE MASTER LIST CAN BE ADDED TO OVER TIME AS NEW ALGORITHMS ARISE
ref_master = {'CANYONB' 'ESPER_NN' 'ESPER_LIR' 'ESPER_MIX' , ...
              'ESPER_MIX_NoO2' 'WOA'};

ref_tmp = NaN(size(data(:,1))); %predim output ref array template

if isempty(WOA_path)
    WOA_path = default_WOA;
end

if size(Refs,1) > size(Refs,2) % if cell array make sure 1 x n
    Refs = Refs';
end

if ischar(Refs) && strcmpi(Refs,'all') % shortcut input
    ref_list = ref_master;
else
    tg       = ismember(ref_master, Refs);
    ref_list = ref_master(tg);
end

if isempty(ref_list)
    fprintf('Refs input is not a valid option: ')
    disp(Refs)
    fprintf('Valid options:\n')
    disp(ref_master);
    return
end

% **********************************************************************
%     1         2     3    4     5     6    7       8
% [MatlabSDN, Cycle, Lon, Lat, Pres, Salt, Temp Oxygen]
% GET DATA

if any(strncmp(ref_list,'ESPER',5)) % ESPER PREP
    Z       = sw_dpth(data(:,5), data(:,4)); % Pres, Lat
    dvec    = datevec(data(:,1));
    % very crude decimal year for OA 30.41*12 =~ 365
    dec_yr  = dvec(:,1) +(dvec(:,2)*30.41)/365 + dvec(:,3)/365;

    %DesireVar      = [3,5]; % pH, NO3
    OutCoords      = [data(:,3), data(:,4), Z]; % Lon, Lat, Depth
    PredictorTypes = [1 2 6]; % PSAL, TEMP, DOXY_ADJ,
    Measurements   = [data(:,6), data(:,7), data(:,8)]; % S,T, O2

    %TESTING EFFECT OF SALINITY ON ESPER
    %Measurements(:,1) = Measurements(:,1) - 0.1;

    Equations      = 7; % S, T, O2

    PredictorTypes_NoO2 = [1, 2]; % PSAL, TEMP,
    Measurements_NoO2   = [data(:,6), data(:,7)]; % S, T
    Equations_NoO2      = 8;
end


if sum(~isnan(Measurements(:,3))) == 0
    srch_str = 'CANYONB|ESPER_NN|ESPER_LIR|ESPER_MIX$';
    tbad       = ~cellfun(@isempty, regexp(ref_list, srch_str,'once'));

    if any(tbad)
        %ref_list = ref_list(tg);
        ref_list = ref_list(~tbad);
        fprintf(['No O2 so O2 dependent references have been removed ',...
            'from the list!']);
    end
    clear tbad
end

% STEP THROUGH FLAVORS
tgO2 = ~isnan(data(:,8)); % any NaN's in O2 file? (all NaN already taken care of)
if sum(~tgO2) >0
    fprintf(['WARNING: Bad O2 values detected (%d) in reference ', ...
        'input data matrix will be exlcuded where applicable\n'], sum(tgO2));
end

% GOOD DATA ARRAYS
tg_S = tfQF(:,1); % S
tg_O2 = tfQF(:,1) & tfQF(:,2); % S, O2
tg_N  = tfQF(:,1) & tfQF(:,2) & tfQF(:,3); % S, O2, NO3
tg_PH = tfQF(:,1) & tfQF(:,2) & tfQF(:,4); % S, O2, PH

for ref_ct = 1:size(ref_list,2)
    ref = ref_list{ref_ct};

    switch ref

        %gtime,lat,lon,pres,temp,psal,doxy,
        case 'CANYONB' % umol /kg
            fprintf('Generating CANYON-B NO3 & pH estimates....\n');

            d.(ref).Nitrate = ref_tmp;
            if any(tg_N)
                CB = CANYONB(data(tg_N,1), data(tg_N,4), data(tg_N,3), ...
                    data(tg_N,5), data(tg_N,7), data(tg_N,6), data(tg_N,8),...
                    {'NO3'});
                d.(ref).Nitrate(tg_N) = CB.NO3;
            end

            d.(ref).pH = ref_tmp;
            if any(tg_PH)
                CB = CANYONB(data(tg_PH,1), data(tg_PH,4), data(tg_PH,3), ...
                    data(tg_PH,5), data(tg_PH,7), data(tg_PH,6), data(tg_PH,8),...
                    {'pH'});

                % CANYON B pH is consistent with pHtotal calculated from DIC and TA
                % We want it to be consitent with pHtotal measured
                % spectrophotometrically using puriﬁed m-cresol purple
                % indicator dye.
                % Stealing linear regresssion from guts of ESPER_NN line 581
                % Solving equation 1 in the LIRv2 paper to counter 3->4 adjustment.
                %Est=(Est+0.3168)/1.0404; (ESPER pN => CANYONB pH) - We want
                %opposite

                d.(ref).pH(tg_PH) = CB.pH * 1.0404 - 0.3168; % Closer to spec pH now
            end
            clear CB

        case 'ESPER_NN' 
            fprintf('Generating ESPER-LIR nutrient estimates....\n');

            d.(ref).Nitrate = ref_tmp;
            if any(tg_N)
                [Estimates,~] = ESPER_LIR(5, OutCoords(tg_N,:), ...
                    Measurements(tg_N,:), PredictorTypes, 'Equations', ...
                    Equations, 'EstDates', dec_yr(tg_N));
                d.(ref).Nitrate(tg_N) = Estimates.nitrate;
            end

            d.(ref).pH = ref_tmp;
            if any(tg_PH)
                [Estimates,~] = ESPER_LIR(3, OutCoords(tg_PH,:), ...
                    Measurements(tg_PH,:), PredictorTypes, 'Equations', ...
                    Equations, 'EstDates', dec_yr(tg_PH));
                d.(ref).pH(tg_PH)      = Estimates.pH;
            end

        case 'ESPER_LIR'
            fprintf('Generating ESPER-NN nutrient estimates....\n');

            d.(ref).Nitrate = ref_tmp;
            if any(tg_N)
                [Estimates,~] = ESPER_NN(5, OutCoords(tg_N,:), ...
                    Measurements(tg_N,:), PredictorTypes, 'Equations', ...
                    Equations, 'EstDates', dec_yr(tg_N));
                d.(ref).Nitrate(tg_N) = Estimates.nitrate;
            end

            d.(ref).pH = ref_tmp;
            if any(tg_PH)
                [Estimates,~] = ESPER_NN(3, OutCoords(tg_PH,:), ...
                    Measurements(tg_PH,:), PredictorTypes, 'Equations', ...
                    Equations, 'EstDates', dec_yr(tg_PH));
                d.(ref).pH(tg_PH)      = Estimates.pH;
            end

        case 'ESPER_MIX'
            fprintf('Generating ESPER-MIX nutrient estimates....\n');

            d.(ref).Nitrate = ref_tmp;
            if any(tg_N)
                [Estimates,~] = ESPER_Mixed(5, OutCoords(tg_N,:), ...
                    Measurements(tg_N,:), PredictorTypes, 'Equations', ...
                    Equations, 'EstDates', dec_yr(tg_N));
                d.(ref).Nitrate(tg_N) = Estimates.nitrate;
            end

            d.(ref).pH = ref_tmp;
            if any(tg_PH)
                [Estimates,~] = ESPER_Mixed(3, OutCoords(tg_PH,:), ...
                    Measurements(tg_PH,:), PredictorTypes, 'Equations', ...
                    Equations, 'EstDates', dec_yr(tg_PH));
                d.(ref).pH(tg_PH)      = Estimates.pH;
            end

        case 'ESPER_MIX_NoO2'
            fprintf('Generating ESPER-MIX No O2 nutrient estimates....\n');

            d.(ref).Nitrate = ref_tmp;
            if any(tg_S & tfQF(:,3))
                [Estimates,~] = ESPER_Mixed(5, OutCoords(tg_S & tfQF(:,3),:), ...
                    Measurements_NoO2(tg_S & tfQF(:,3),:), PredictorTypes_NoO2, 'Equations', ...
                    Equations_NoO2, 'EstDates', dec_yr(tg_S & tfQF(:,3),:));
                d.(ref).Nitrate(tg_S & tfQF(:,3)) = Estimates.nitrate;
            end

            d.(ref).pH = ref_tmp;
            if any(tg_S & tfQF(:,4))
                [Estimates,~] = ESPER_Mixed(3, OutCoords(tg_S & tfQF(:,4),:), ...
                    Measurements_NoO2(tg_S & tfQF(:,4),:), PredictorTypes_NoO2, 'Equations', ...
                    Equations_NoO2, 'EstDates', dec_yr(tg_S & tfQF(:,4),:));
                d.(ref).pH(tg_S & tfQF(:,4))  = Estimates.pH;
            end

        case 'WOA' % NO3 ONLY
            fprintf('EXTRACTING WOA NO3 profiles float cycle locations....\n');
            % base on cycle, in profile time could change
            [uCycles,ia,~] = unique(data(:,2)); 
            WOA_data = data(ia,[1,4,3]);
            
            % INPUTS:
            %   data_source = path to local WOA2018 data repo
            %   track      = n x 3 matrix [Matlab_SDN, Lat, Lon]
            %   depth_bnds = depth bounds [min depth  max depth]
            %	ocean_var  = a string, ocean parameter
            %                must be one of these: T S O2 O2sat NO3 Si PO4 AOU

            dWOA = get_WOA23_local(WOA_path, WOA_data , [0 2100], 'NO3');
            
            % INTERPLATE TO FLOAT PROFILE PRESSURE AXES
            WOA_NO3 = nan(size(data(:,1))); % NO3, pH but pH will be empty
            for pct = 1: size(uCycles,1)
                t1 = data(:,2) == uCycles(pct);
                WOA_NO3(t1) = interp1(dWOA.Z, dWOA.d(:,pct), data(t1,5));
            end
            d.(ref).Nitrate = WOA_NO3; % pH all NaN
    end

    if tf_app == 1
        str = sprintf('Loading Ref: %s ...',ref);
        disp(str)
        app.SelectfileButton.Text = str;
        drawnow
    end
end

