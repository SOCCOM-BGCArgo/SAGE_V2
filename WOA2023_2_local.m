function WOA2023_2_local(dest)
% This will dowload al WOA2023 parameter files to your local computer based
% on the info stored in the "WOA_info" cell array. Modify the array to
% change desired variables. For example in column 3 'B5C2' represents the
% 2015-2022 average for T & S.

% ***********************************************************************
% TESTING
%dest = '';
% ***********************************************************************

if isempty(dest)
    %dest = fullfile(userpath, 'ARGO_PROCESSING\DATA\WOA2023\');
    dest = fullfile(userpath, 'WOA2023\'); % store out side of Argo processing
end

fprintf(['Downloading all WOA 2023 data files will take some time ', ...
    'and requires ~ 4.8Gb of disk space\n'])

% VARIABLES EXTRACTION TABLE
WOA_info(1,:) = {'nitrate'     'n'  'all' };
WOA_info(2,:) = {'oxygen'      'o'  'all' };
WOA_info(3,:) = {'phosphate'   'p'  'all' };
WOA_info(4,:) = {'silicate'    'i'  'all' };
WOA_info(5,:) = {'temperature' 't'  'B5C2'}; % 2015-2022
WOA_info(6,:) = {'salinity'    's'  'B5C2'}; % 2015-2022
WOA_info(7,:) = {'o2sat'       'O'  'all' };
WOA_info(8,:) = {'AOU'         'A'  'all' };

catalogue_str = ['https://www.ncei.noaa.gov/thredds-ocean/fileServer/',...
    'woa23/DATA/PARAM/netcdf/XXX/1.00/'];

fn_str        = 'woa23_%s_%s%02d_01.nc';

for ct = 1:size(WOA_info,1)
    fprintf('Downloading WOA 2023 %s files ......\n',WOA_info{ct,1})

    pd = fullfile(dest, WOA_info{ct,1},'\'); % destination for param files
    if ~isfolder(pd)
        mkdir(pd)
    end

    source_fd = regexprep(catalogue_str,'PARAM', WOA_info{ct,1});
    source_fd = regexprep(source_fd,'XXX', WOA_info{ct,3});

    for fct = 1:17 % 16 files
        fn      = sprintf(fn_str, WOA_info{ct,3}, WOA_info{ct,2}, fct-1);
        fp      = fullfile(source_fd,fn);
        dest_fn = fullfile(dest,WOA_info{ct,1},fn);
        fprintf('    %s\n',fn)
        ofn     = websave(dest_fn,fp);
    end
end

% LASTLY GET WOA 2023 DOCUMENTATION:
fprintf('Downloading WOA 2023 product documentation ......\n');
doc_url = ['https://www.ncei.noaa.gov/data/oceans/woa/WOA23/',...
    'DOCUMENTATION/WOA23_Product_Documentation.pdf'];
dest_fn = fullfile(dest,'WOA23_Product_Documentation.pdf');
ofn     = websave(dest_fn, doc_url);

clearvars



