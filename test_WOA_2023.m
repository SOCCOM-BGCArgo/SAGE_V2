% test WOA 2023
WOA_info = {'nitrate' 'n'};
YR = {'23' '18' '13'}; % choos WOA , ie 2018 or 2023 or 2013
data_source = fullfile(userpath, ['WOA20',YR,'\']);
MO = 1:12;
fn_str = 'woa%s_all_%s%02d_01.nc';
% woa18_all_n00_01.nc

% lon_bnds = [-153.50000 -150.50000];
% lat_bnds = [44.5 46.5];

lon_bnds = [-153.316-3	-151.27+3];
lat_bnds = [45.067-3 45.984+3];






for ct = 1: size(MO,2) % step through months
    F1 = figure(1);
    F1.Units = 'inches';
    F1.Position = [2.1875 5.5833 12.1979 4.7917];
    
    for dt = 1:size(YR,2) % step thrugh data sets

        yr = YR{dt};
        data_source = fullfile(userpath, ['WOA20',yr,'\']);

        fn = sprintf(fn_str, yr, WOA_info{2}, MO(ct));
        fp = fullfile(data_source, WOA_info{1}, fn);
        lon = ncread(fp,'lon');
        lat = ncread(fp,'lat');
        z   = ncread(fp,'depth');

        tLON = lon >= lon_bnds(1) & lon <= lon_bnds(2);
        tLAT = lat >= lat_bnds(1) & lat <= lat_bnds(2);

        lon_inds  = [find(tLON == 1,1,'first'), find(tLON == 1,1,'last')]; % monthly
        lat_inds  = [find(tLAT == 1,1,'first'), find(tLAT == 1,1,'last')];
        start_idx = [lon_inds(1) lat_inds(1) 1 1]; % extraction start vector
        count_idx = [diff(lon_inds)+1 diff(lat_inds)+1 inf 1]; % extraction count vector
        lon = lon(tLON); lat = lat(tLAT);
        an  = ncread(fp, 'n_an', start_idx, count_idx);  % objectively analyzed climatology
        %an  = ncread(fp, 'n_mn', start_idx, count_idx);  % objectively analyzed climatology

        s_an =size(an);
        jp = reshape(permute(an,[3,1,2]),s_an(3),s_an(1)*s_an(2)); %lon lat Z => Z lon lat
        % I checked individual profiles that I am permuting, reshaping
        % correctly

        ax(dt) = subplot(1,3,dt);
        plot(jp, z*ones(1,size(jp,2)),'o-'), set(gca,'YDir','Reverse')
        tcell = { sprintf('%s',fn), ...
            sprintf('lon bounds = [%0.1f %0.1f]',lon_bnds), ...
            sprintf('lat bounds = [%0.1f %0.1f]',lat_bnds)};
        title(tcell,'interpreter', 'none')
        ylabel('Depth, m')
        xlabel('Nitrate, µM')
        ylim([0 800]);
        xlim([5 45])

        %if MO(ct) == 7, pause, end

        save_fn = sprintf('WOA_%02d.png',MO(ct));
        save_fp = fullfile('C:\temp',[save_fn,'.png']);
    end
    print(F1, save_fp,'-r160','-dpng');
    pause
end