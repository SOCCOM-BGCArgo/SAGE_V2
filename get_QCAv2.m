function QCA = get_QCAv2(qc_path)
% RETURN QC ADJUSTMENTS FOR A GIVEN PARAMETER
%
% INPUTS:
%   qc_path    - path to qc adjustment list file
%   float_name - mabri float name
%   data_type  - 'NO3' , 'PH', or 'O2'
%
% OUTPUT:
%   QCA       - A STRUCTURE OF ADJUSTMENTS 
%      .offset - for pH only
%      .data   - a matrix [CYLE GAIN OFFSET DRIFT]

% *************************************************************************
% TESTING
% qc_path = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\', ...
%            'CAL\QC_LISTS\5905107_FloatQCList.txt'];

% qc_path = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\', ...
%            'CAL\QC_LISTS\NO_WMO_un0948_FloatQCList.txt'];

% qc_path = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\', ...
%     'CAL\QC_LISTS\5906568_FloatQCList_JP.txt'];

% qc_path = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\', ...
%     'CAL\QC_LISTS\5906568_FloatQCList.txt'];

% qc_path = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\', ...
%     'CAL\QC_LISTS\5906568_FloatQCList.txt'];

% qc_path = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\', ...
%     'CAL\QC_LISTS\5906000_FloatQCList.txt'];

% qc_path = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\', ...
%     'CAL\QC_LISTS\5905100_FloatQCList.txt'];



% 
% *************************************************************************
% CHANGES:
% 08/02/2017 - modified code to now use float-specific FloatQCList files, located in DATA\CAL\QC_LISTS\
% 10/30/2017 - get CHL correction if it exists
% 04/03/2021 - updated for GOBGC name changes, removed 2nd input. Get WMO
%              from file name
% 09/13/24  - JP updated field names to match QC list
% 01/21/24  - SB updated parsing of past QC to do a case sensitive regexp
%             for 'PREVIOUS'. Issues occured with parsing when 'previous'
%             was used in QC comments.


% *************************************************************************
% predim ouput
QCA.Oxygen    = [];
QCA.Nitrate   = [];
QCA.pH_offset = '';
QCA.pH        = [];
QCA.QC_info   = []; % QC meta info


% GET FILE NAME & WMO FROM FILE NAME (\w get underscore too)
fn  = regexp(qc_path,'\w+\.txt','once','match');
wmo = regexp(fn,'\w+(?=_Float)','once','match');

fid = fopen(qc_path);
if fid < 0 %invalid file identifier
    fprintf('QC LIST FILE NOT FOUND : %s\n', qc_path);
    return
end

% FIND SPECIFIC SECTIONS
tline = '';

% FIND SPECIFIC DATA TYPE
while ischar(tline)
    tmp = regexp(tline,',','split');

    if regexp(tmp{1},'PREVIOUS','once') % Done with parsing current QC
        break

    elseif strcmp(tmp{1},'Oxygen')
        QCA.Oxygen  =[QCA.Oxygen; str2double(tmp(2:end))];

    elseif strcmp(tmp{1},'Nitrate') % cycle gain offset, drift
        QCA.Nitrate  =[QCA.Nitrate; str2double(tmp(2:end))];

    elseif regexp(tline,'pH.+offset','once') % string defining correction type
        QCA.pH_offset = regexp(tline,'linear|poly|mixed','once','match');
        if isempty(QCA.pH_offset)
            QCA.pH_offset = 'none';
        end

    elseif strcmp(tmp{1},'pH')  % cycle, offset, drift  or  cycle gain offset, drift
        QCA.pH =[QCA.pH; str2double(tmp(2:end))];

    % Re-make QC_info structures that can be easily suucked into sageV2
    % These can be used to "remember" last QC settings evantually
    elseif strcmp(tmp{1},'QC parameter')
        tmp_info = regexp(tmp(3:end),',','split'); %skip param ID, already grabbed
        info = reshape(tmp_info, 2, numel(tmp_info)/2)'; %8x2
        for ct = 1:size(info,1)
            if strcmp(info{ct,1},'PresRange')
                d.(info{ct,1}{1}) = sscanf(info{ct,2}{1},'[%f %f]')';
            elseif contains(info{ct,1}{1},{'User','Date','Ref'})
                d.(info{ct,1}{1}) = info{ct,2}{1};
            elseif contains(info{ct,1}{1},{'Node','BIC','RMSE'})
                d.(info{ct,1}{1}) = str2double(info{ct,2}{1});
            end
        end
        QCA.QC_info.(tmp{2}) = d;


    end
    tline = fgetl(fid);
end

fclose(fid);

% Backwards compatibility for older  pH QC lines (n x 3 with no gain col)
if isfield(QCA,'pH') && size(QCA.pH,2) == 3 % old style. no gain col so add one
    QCA.pH = [QCA.pH(:,1), ones(size(QCA.pH(:,1))), QCA.pH(:,2:3)];
end



