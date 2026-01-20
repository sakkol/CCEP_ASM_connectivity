% Figure 3 for channel based functional connectviity comparison

data_root = fullfile(cd,'neural_data');
project_name = 'CCEP_PrePost';

figure('position',[0 0 1800 750])
for s = 1:2
    if s==1, sbj_ID = 'P2';
    elseif s==2, sbj_ID = 'P4';
    end
    subplot(1,2,s)

    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    load(fullfile(Sbj_Metadata.results,'corr_matrix',[Sbj_Metadata.sbj_ID, '_corr_perm.mat']),'p_values','significant_pairs')

    % load chan info
    AllBlockInfo = readtable(fullfile(Sbj_Metadata.project_root,[Sbj_Metadata.project_name '_BlockInfo.xlsx'])); % "F:\HBML\PROJECTS_DATA\CL_Train\CL_Train_BlockInfo.xlsx"
    preoffmed = AllBlockInfo.BlockList{strcmp(AllBlockInfo.sbj_ID,sbj_ID) & strcmp(AllBlockInfo.task_type,'preoffmedrest')};
    preonmed = AllBlockInfo.BlockList{strcmp(AllBlockInfo.sbj_ID,sbj_ID) & strcmp(AllBlockInfo.task_type,'preonmedrest')};
    off_info = load(fullfile(Sbj_Metadata.iEEG_data, preoffmed, [preoffmed '_info.mat']));
    on_info = load(fullfile(Sbj_Metadata.iEEG_data, preonmed, [preonmed '_info.mat']));
    [~,off_nonartefact_idx] = get_info_nonartefactchans(off_info.info);
    [~,on_nonartefact_idx] = get_info_nonartefactchans(on_info.info);
    good_ch_idx = off_nonartefact_idx & on_nonartefact_idx;
    soz_idx = strcmp(off_info.info.channelinfo.chan_info,'SOZ') & strcmp(on_info.info.channelinfo.chan_info,'SOZ');
    EPZ_idx = strcmp(off_info.info.channelinfo.chan_info,'EPZ') & strcmp(on_info.info.channelinfo.chan_info,'EPZ');
    IZ_idx = strcmp(off_info.info.channelinfo.chan_info,'IZ') & strcmp(on_info.info.channelinfo.chan_info,'IZ');
    Healthy_idx = strcmp(off_info.info.channelinfo.chan_info,'Healthy') & strcmp(on_info.info.channelinfo.chan_info,'Healthy');
    soz_idx = soz_idx & good_ch_idx;
    EPZ_idx = EPZ_idx & good_ch_idx;
    IZ_idx = IZ_idx & good_ch_idx;
    Healthy_idx = Healthy_idx & good_ch_idx;
    chans = on_info.info.channelinfo.Label;

    onmedcorr = load(fullfile(Sbj_Metadata.iEEG_data, preonmed,'HFB_files','corr_matrix.mat'),'corrmatrix');
    offmedcorr = load(fullfile(Sbj_Metadata.iEEG_data, preoffmed,'HFB_files','corr_matrix.mat'),'corrmatrix');

    % plot and save the results
    r=onmedcorr.corrmatrix - offmedcorr.corrmatrix;
    r(~significant_pairs)=NaN; % filter based on permutation results
    % change order
    new_order = [find(soz_idx);find(EPZ_idx);find(IZ_idx);find(Healthy_idx)];
    r=r(new_order,new_order);
    isupper = logical(triu(ones(size(r)),1));
    r(isupper) = 0;
    for i=1:size(r,1),r(i,i)=NaN;end

    h = imagesc(r);
    set(gca,'color',0*[1 1 1]);
    set(h, 'AlphaData', ~isnan(r))
    % clim_lim = [-max(max(abs(r))) max(max(abs(r)))];
    % clim(clim_lim)
    % set(gca,'XTick',1:length(new_order),'XTickLabel',chans(new_order));
    % set(gca,'YTick',1:length(new_order),'YTickLabel',chans(new_order));
    set(gca,'Colormap',master_ColorMaps('bwr'))
    clim_lim = [-max(max(abs(r))) max(max(abs(r)))];
    % clim(clim_lim)
    clim([-0.18 0.18])

    % for i=1:length(new_order)
    %     xline([i-.5, i+.5], 'k-','LineWidth',.2);
    %     yline([i-.5, i+.5], 'k-','LineWidth',.2);
    % end
    col = [length(find(soz_idx)), length(find(soz_idx))+length(find(EPZ_idx)),...
        length(find(soz_idx))+length(find(EPZ_idx))+length(find(IZ_idx))];
    % length(find(soz_idx))+length(find(EPZ_idx))+length(find(IZ_idx))+length(find(Healthy_idx))];
    xline(gca, col+.5, 'r-','LineWidth',2);
    yline(gca, col+.5, 'r-','LineWidth',2);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);

end

set(gcf,'Color','w')
set(gcf, 'InvertHardcopy', 'off')
savedir = fullfile(cd,'Figures');
print(fullfile(savedir,['Figure3_' char(datetime('today','Format','uuuu-MM-dd')) '.png']),'-dpng','-r300')


%% Add colormap

figure('Units','normalized','Position',[0 0 .5 1])
% ax = axes;
colormap(master_ColorMaps('bwr'));
% set(gca,'ColorScale','log')
clim([-0.18 0.18])
c=colorbar(gca);
% c.Label.String = 'Percentages';c.Label.Rotation=270;
c.FontSize=30;
c.FontWeight='bold';
c.LineWidth = 1;
% c.Limits=[0 25];
% c.Ticks = [0,5,10,15,20,25];
c.TickDirection='both';
c.Location='eastoutside';

a=get(c); %gets properties of colorbar
a = a.Position; %gets the positon and size of the color bar
% set(c,'Units','normalized','Position',[.1 .1 .5 .03])% To change size
set(c,'Units','normalized','Position',[.7 .1 .02 .5])% To change size
set(gca,'Visible','off')
% ax.Visible = 'off';
set(gcf,'color','w')
set(gcf, 'InvertHardcopy', 'off') % to set saving white color as white

print(fullfile(savedir,['Figure3_colorbar_' char(datetime('today','Format','uuuu-MM-dd'))]),'-dpng','-r300')
