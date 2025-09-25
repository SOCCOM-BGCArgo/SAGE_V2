
% FIND DEAD FLOATS ON 980 LIST
% pH_pumpoffset_980
% pH980 = pH_pumpoffset_980_floats; % make a smaller name
% load('C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\CAL\MBARI_float_list.mat')
% tf_dead = cell2mat(d.list(:,13)) == 1;
% tf_APEX = strcmp(d.list(:,4),'APEX');
% tf_pH   = cell2mat(d.list(:,15)) == 1;
% mlist   =  d.list(tf_dead & tf_APEX & tf_pH,:);
% Lia = ismember(str2double(mlist(:,3)), pH980);
% dlist = mlist(Lia,:);
% 
% pause


fd1 = 'C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATVIZ\QC\';
fd2 = '\\atlas\chem\ARGO_PROCESSING\DATA\FLOATVIZ\QC\';

%wmo = '5906568';
%wmo = '5905099';
%wmo = '5905991';
wmo = '5906031';


fn = sprintf('%sQC.TXT',wmo);
fp1 = fullfile(fd1, fn); % local - offset corrected
fp2 = fullfile(fd2, fn); % chem version - no corection

d1 = get_FloatViz_data(fp1);
d1.data(d1.data == -1e10) = NaN;

d2 = get_FloatViz_data(fp2);
d2.data(d2.data == -1e10) = NaN;


iSTA  = find(strcmp('Station',d1.hdr) == 1);
iSDN  = find(strcmp('SDN',d1.hdr)   == 1);
iLAT  = find(strcmp('Lat [°N]',d1.hdr)   == 1);
iLON  = find(strcmp('Lon [°E]',d1.hdr)   == 1);
iP    = find(strcmp('Pressure[dbar]',     d1.hdr)  == 1);
iT    = find(strcmp('Temperature[°C]',    d1.hdr)  == 1);
iS    = find(strcmp('Salinity[pss]',      d1.hdr)  == 1);
iO    = find(strcmp('Oxygen[µmol/kg]',    d1.hdr)  == 1);
iN    = find(strcmp('Nitrate[µmol/kg]',   d1.hdr)  == 1);
iPH   = find(strcmp('pHinsitu[Total]',    d1.hdr)  == 1);

d1.data(d1.data(:,iPH+1) == 8, iPH) = NaN;
d2.data(d2.data(:,iPH+1) == 8, iPH) = NaN;

uC1 = unique(d1.data(:,iSTA));

% BUILD PROFILE COLORMAP
cmap = colormap;
plt_cmap = interp1((1:size(cmap,1))', cmap, linspace(1, size(cmap,1), size(uC1,1))');

F10          = figure(10);
F10. Units   = 'inches';
F10.Position = [2 2 8 8];

ax(1)      = subplot(3,2,1:2);
ax(1).Position = ax(1).Position + [0 0.05 0 -0.05];
ax(1).Title.String = wmo;
ax(1).YLabel.String = 'Surface pH';
%ax(1).XLabel.String = 'Cycle';
hold(ax(1),'on');

ax(2)      = subplot(3,2,3:4);
ax(2).Position = ax(2).Position + [0 0.12 0 -0.05];
ax(2).YLabel.String = {'Surface pH' 'PO cor - 980 cor'};
ax(2).XLabel.String = 'Cycle';
hold(ax(2),'on');

ax(3)      = subplot(3,2,5);
ax(3).Position = ax(3).Position + [0 0 0 0.1];
ax(3).YDir = 'Reverse';
ax(3).YLabel.String = 'Pressure, dbar';
ax(3).XLabel.String = {'adj pH' 'pump offset corrected'};
hold(ax(3),'on');

ax(4)      = subplot(3,2,6);
ax(4).Position = ax(4).Position + [0 0 0 0.1];
ax(4).YDir = 'Reverse';
ax(4).YLabel.String = 'Pressure, dbar';
ax(4).XLabel.String = {'adj pH' 'corrected at 950m'};
hold(ax(4),'on');

% PLOT PROFILES
for ct = 1:size(uC1,1)
    t1   = d1.data(:,iSTA) == uC1(ct);
    tmp1 = d1.data(t1,:);
    tP1  = tmp1(:,iP) < 30;

    t2   = d2.data(:,iSTA) == uC1(ct);
    tmp2 = d2.data(t2,:);
    tP2  = tmp2(:,iP) < 30;

    plot(ax(1),tmp1(tP1,iSTA), tmp1(tP1,iPH),'ko',...
        'MarkerfaceColor','b', 'MarkerSize',6)
    plot(ax(1),tmp2(tP2,iSTA), tmp2(tP2,iPH),'ko',...
        'MarkerfaceColor','r', 'MarkerSize',6)

    % Get surf diff
    pH_diff = tmp1(tP1,iPH) - tmp2(tP2,iPH);
    plot(ax(2),tmp1(tP1,iSTA), pH_diff,'ko',...
        'MarkerfaceColor','g', 'MarkerSize',6)

    plot(ax(3),tmp1(:,iPH), tmp1(:,iP),'Color', plt_cmap(ct,:));
    plot(ax(4),tmp2(:,iPH), tmp2(:,iP),'Color', plt_cmap(ct,:));

end

hold(ax(1),'off');
hold(ax(2),'off');
hold(ax(3),'off');
hold(ax(4),'off');

ax(1).FontSize = 12;
ax(2).FontSize = 12;
ax(3).FontSize = 12;
ax(4).FontSize = 12;

L1 = legend(ax(1),'PO cor','980 cor','Location','SouthEast');

x_lims = sort([ax(3).XLim, ax(4).XLim]);
x_lim  = [x_lims(1), x_lims(end)]; 
ax(3).XLim = x_lim;
ax(4).XLim = x_lim;

save_fd = 'C:\Users\jplant\Documents\MATLAB\Apps\SAGE_V2\PLOTS';
save_fn = sprintf('pH_PO_corr_comp_%s.png',wmo);
save_fp = fullfile(save_fd, save_fn);
print(F10, save_fp,'-dpng');









