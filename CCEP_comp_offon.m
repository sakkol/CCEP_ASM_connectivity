function CCEP_comp_offon(Sbj_Metadata, offmed_block, onmed_block, elec_pair, elec_type)
% This function is to find the N1 and N2 peaks within the data

%% Load CCEP_proc_plot2023 outputs
ref='bp';
xx=strsplit(elec_pair,'-');
StimCh1=xx{1};
StimCh2=xx{2};clear xx

off_ccep = load(fullfile(Sbj_Metadata.iEEG_data,offmed_block, 'CCEPdir', [offmed_block '_ccep_' ref '_' StimCh1 '-' StimCh2 '.mat']));
on_ccep = load(fullfile(Sbj_Metadata.iEEG_data,onmed_block, 'CCEPdir', [onmed_block '_ccep_' ref '_' StimCh1 '-' StimCh2 '.mat']));
off_info = load(fullfile(Sbj_Metadata.iEEG_data, offmed_block, [offmed_block '_info.mat']));
importance = {'SOZ','EPZ','IZ','Healthy'};
off_info_chans_split = strsplit_SA(off_ccep.CCEPs.label);
[off_info_chans_classes] = assign_classes(off_info_chans_split(:,1), off_info_chans_split(:,2), off_info.info.channelinfo.Label, off_info.info.channelinfo.chan_info,importance);
% on_info = load(fullfile(Sbj_Metadata.iEEG_data, onmed_block, [onmed_block '_info.mat']));

% smooth data a bit
% off_ccep = master_smoothData(off_ccep.CCEPs,0.02);
% on_ccep = master_smoothData(on_ccep.CCEPs,0.02);
off_ccep=off_ccep.CCEPs;
on_ccep=on_ccep.CCEPs;

%% Comparing on vs off
% Prepare output info (peaks) to be saved
peak_info=table;
peak_info.label = off_ccep.label;
peak_info.off_peak_times = cell([height(peak_info.label),1]);
peak_info.off_peak_amps = cell([height(peak_info.label),1]);
peak_info.off_N1_time = cell([height(peak_info.label),1]);
peak_info.off_N1_amp = cell([height(peak_info.label),1]);
peak_info.off_N2_time = cell([height(peak_info.label),1]);
peak_info.off_N2_amp = cell([height(peak_info.label),1]);
peak_info.off_N1_time_all = cell([height(peak_info.label),1]);
peak_info.off_N1_amp_all = cell([height(peak_info.label),1]);
peak_info.off_N2_time_all = cell([height(peak_info.label),1]);
peak_info.off_N2_amp_all = cell([height(peak_info.label),1]);
peak_info.off_rms = cell([height(peak_info.label),1]);
peak_info.off_rmsN1 = cell([height(peak_info.label),1]);
peak_info.off_rmsN2 = cell([height(peak_info.label),1]);

peak_info.on_peak_times = cell([height(peak_info.label),1]);
peak_info.on_peak_amps = cell([height(peak_info.label),1]);
peak_info.on_N1_time = cell([height(peak_info.label),1]);
peak_info.on_N1_amp = cell([height(peak_info.label),1]);
peak_info.on_N2_time = cell([height(peak_info.label),1]);
peak_info.on_N2_amp = cell([height(peak_info.label),1]);
peak_info.on_N1_time_all = cell([height(peak_info.label),1]);
peak_info.on_N1_amp_all = cell([height(peak_info.label),1]);
peak_info.on_N2_time_all = cell([height(peak_info.label),1]);
peak_info.on_N2_amp_all = cell([height(peak_info.label),1]);
peak_info.on_rms = cell([height(peak_info.label),1]);
peak_info.on_rmsN1 = cell([height(peak_info.label),1]);
peak_info.on_rmsN2 = cell([height(peak_info.label),1]);

peak_info.ttestresults{1}=[];

% choosing colors
cmap=master_ColorMaps('RdBu',20);
off_col = cmap(3,:);
on_col = cmap(18,:);
alpha_c = 0.05/length(off_ccep.label);

% save directory
res_dir=fullfile(Sbj_Metadata.results,'off_on_peaks');
if ~isfolder(res_dir); mkdir(res_dir); end

% Start plotting all channels separately
for el = 1:length(off_ccep.label)

    %% Calculate averages and peaks (N1:10-30ms; N2:80-250ms)
    % for offmed:
    off_for_avg=squeeze(off_ccep.trial(:,el,:));
    off_avg = mean(off_for_avg,1);

    % for onmed:
    on_for_avg=squeeze(on_ccep.trial(:,el,:));
    on_avg = mean(on_for_avg,1);

    %% Calculate peaks (either positive or negative)
    % offmed: get all peaks, N1 and N2 peaks 
    [peak_info.off_peak_amps{el},peak_info.off_peak_times{el},...
        peak_info.off_N1_amp{el},peak_info.off_N1_time{el},...
        peak_info.off_N2_amp{el},peak_info.off_N2_time{el}] = ...
        ccep_N1N2(off_avg,off_ccep.time);
    for t = 1:size(off_for_avg,1)
    [~,~,...
        peak_info.off_N1_amp_all{el}(t),peak_info.off_N1_time_all{el}(t),...
        peak_info.off_N2_amp_all{el}(t),peak_info.off_N2_time_all{el}(t)] = ...
        ccep_N1N2(off_for_avg(t,:),off_ccep.time);
    end
    peak_info.off_rms{el} = ccep_rms(off_avg,off_ccep.time);
    peak_info.off_rmsN1{el} = ccep_rms(off_avg,off_ccep.time,'N1');
    peak_info.off_rmsN2{el} = ccep_rms(off_avg,off_ccep.time,'N2');

    % onmed: get all peaks, N1 and N2 peaks 
    [peak_info.on_peak_amps{el},peak_info.on_peak_times{el},...
        peak_info.on_N1_amp{el},peak_info.on_N1_time{el},...
        peak_info.on_N2_amp{el},peak_info.on_N2_time{el}] = ...
        ccep_N1N2(on_avg,on_ccep.time);
    for t = 1:size(on_for_avg,1)
    [~,~,...
        peak_info.on_N1_amp_all{el}(t),peak_info.on_N1_time_all{el}(t),...
        peak_info.on_N2_amp_all{el}(t),peak_info.on_N2_time_all{el}(t)] = ...
        ccep_N1N2(on_for_avg(t,:),on_ccep.time);
    end
    peak_info.on_rms{el} = ccep_rms(on_avg,on_ccep.time);
    peak_info.on_rmsN1{el} = ccep_rms(on_avg,on_ccep.time,'N1');
    peak_info.on_rmsN2{el} = ccep_rms(on_avg,on_ccep.time,'N2');

    peak_info_toplot={peak_info.off_N1_time_all{el},peak_info.on_N1_time_all{el};peak_info.off_N2_time_all{el},peak_info.on_N2_time_all{el}};
    peak_info_toplot(:,:,2)={peak_info.off_N1_amp_all{el},peak_info.on_N1_amp_all{el};peak_info.off_N2_amp_all{el},peak_info.on_N2_amp_all{el}};

    %% Start plotting
    figure('Units','normalized','Position', [0 0  1 1]);
    off_ccep = master_smoothData(off_ccep,0.005);
    on_ccep = master_smoothData(on_ccep,0.005);
    % for offmed:
    off_for_avg=squeeze(off_ccep.trial(:,el,:));
    off_avg = mean(off_for_avg,1);

    % for onmed:
    on_for_avg=squeeze(on_ccep.trial(:,el,:));
    on_avg = mean(on_for_avg,1);

    %% First the CCEP time-series themselves
    % first the shaded error bars
    subplot(2,6,[1:4,7:10])
    shadedErrorBar(off_ccep.time,off_avg,stderr(off_for_avg),'lineprops',{'Color',off_col,'LineWidth',1});
    hold on
    shadedErrorBar(on_ccep.time,on_avg,stderr(on_for_avg),'lineprops',{'Color',on_col,'LineWidth',1});

    % place N1s and N2s
    N1N2_str = {'\color{black}N1:';['\color[rgb]{' num2str(off_col) '}' 'Off-Med: '];['\color[rgb]{' num2str(on_col) '}' 'On-Med: '];...
        '\color{black}N2:';['\color[rgb]{' num2str(off_col) '}' 'Off-Med: '];['\color[rgb]{' num2str(on_col) '}' 'On-Med: ']};
    
    if ~isempty(peak_info.off_N1_time{el})
    plot([peak_info.off_N1_time{el} peak_info.off_N1_time{el}], ylim,'Color',off_col,'LineWidth',1)
    N1N2_str{2} = [N1N2_str{2} num2str(peak_info.off_N1_amp{el},2) 'uV at ' num2str(peak_info.off_N1_time{el},3)];
    end
    if ~isempty(peak_info.off_N2_time{el})
    plot([peak_info.off_N2_time{el} peak_info.off_N2_time{el}], ylim,'Color',off_col,'LineWidth',1)
    N1N2_str{5} = [N1N2_str{5} num2str(peak_info.off_N2_amp{el},2) 'uV at ' num2str(peak_info.off_N2_time{el},3)];
    end
    if ~isempty(peak_info.on_N1_time{el})
    plot([peak_info.on_N1_time{el} peak_info.on_N1_time{el}], ylim,'Color',on_col,'LineWidth',1)
    N1N2_str{3} = [N1N2_str{3} num2str(peak_info.on_N1_amp{el},2) 'uV at ' num2str(peak_info.on_N1_time{el},3)];
    end
    if ~isempty(peak_info.on_N2_time{el})
    plot([peak_info.on_N2_time{el} peak_info.on_N2_time{el}], ylim,'Color',on_col,'LineWidth',1)
    N1N2_str{6} = [N1N2_str{6} num2str(peak_info.on_N2_amp{el},2) 'uV at ' num2str(peak_info.on_N2_time{el},3)];
    end

    text(1,0,N1N2_str,'Units','normalized','FontSize',15,'VerticalAlignment','bottom','HorizontalAlignment','right')
    
    % general stuff
    plot(xlim,[0 0],'k')
    xlim([-0.1 0.4])
    set(gca, 'FontSize',12);ylims = ylim;
    plot([0 0], ylim,'k')
    set(gca, 'FontSize',12,'FontWeight','bold');
    xlabel('Time (s)');
    ylim(ylims)
    legend({'Off-Meds','On-Meds'},'FontSize',14)

    title({['Stim: ' elec_pair ' (Type: ' elec_type ') || Rec: ' off_ccep.label{el} ' (Type: ' off_info_chans_classes{el} ')']; ...
        ['CCEP (mean+SEM)  ||  Off-Meds: ' num2str(size(off_for_avg,1)) ' trials; On-Meds: ' num2str(size(on_for_avg,1)) ' trials']})

    %% Now the peaks as scatters
    % Top row for latency
    subplot(2,6,5)
    box_scatter_simple(peak_info_toplot(1,:,1),{'N1'},{'Off-Med','On-Med'},[off_col,0.5;on_col,0.5],0)
    [significant,p,~,stats] = ttest2(peak_info_toplot{1,1,1},peak_info_toplot{1,2,1},'alpha',alpha_c,'tail','both');
    peak_info.ttestresults{el}.N1_time.h=significant;
    peak_info.ttestresults{el}.N1_time.p=p;
    peak_info.ttestresults{el}.N1_time.stats=stats;
    text(0.5,1,['p=' num2str(p)],'Units','normalized','FontSize',10,'VerticalAlignment','middle','HorizontalAlignment','center')
    box off
    ylabel('Latency (s)','FontWeight','bold','FontSize',14)
    set(gca,'FontSize',12)
    subplot(2,6,6)
    box_scatter_simple(peak_info_toplot(2,:,1),{'N2'},{'Off-Med','On-Med'},[off_col,0.5;on_col,0.5],0)
    [significant,p,~,stats] = ttest2(peak_info_toplot{2,1,1},peak_info_toplot{2,2,1},'alpha',alpha_c,'tail','both');
    peak_info.ttestresults{el}.N2_time.h=significant;
    peak_info.ttestresults{el}.N2_time.p=p;
    peak_info.ttestresults{el}.N2_time.stats=stats;
    text(0.5,1,['p=' num2str(p)],'Units','normalized','FontSize',10,'VerticalAlignment','middle','HorizontalAlignment','center')
    box off
    set(gca,'FontSize',12)

    % Second row for amplitudes
    subplot(2,6,11)
    box_scatter_simple(peak_info_toplot(1,:,2),{'N1'},{'Off-Med','On-Med'},[off_col,0.5;on_col,0.5],0)
    [significant,p,~,stats] = ttest2(peak_info_toplot{1,1,2},peak_info_toplot{1,2,2},'alpha',alpha_c,'tail','both');
    peak_info.ttestresults{el}.N1_amp.h=significant;
    peak_info.ttestresults{el}.N1_amp.p=p;
    peak_info.ttestresults{el}.N1_amp.stats=stats;
    text(0.5,1,['p=' num2str(p)],'Units','normalized','FontSize',10,'VerticalAlignment','middle','HorizontalAlignment','center')
    box off
    ylabel('Amplitude (uV)','FontWeight','bold','FontSize',14)
    set(gca,'FontSize',11)
    subplot(2,6,12)
    box_scatter_simple(peak_info_toplot(2,:,2),{'N2'},{'Off-Med','On-Med'},[off_col,0.5;on_col,0.5],0)
    [significant,p,~,stats] = ttest2(peak_info_toplot{2,1,2},peak_info_toplot{2,2,2},'alpha',alpha_c,'tail','both');
    peak_info.ttestresults{el}.N2_amp.h=significant;
    peak_info.ttestresults{el}.N2_amp.p=p;
    peak_info.ttestresults{el}.N2_amp.stats=stats;
    text(0.5,1,['p=' num2str(p)],'Units','normalized','FontSize',10,'VerticalAlignment','middle','HorizontalAlignment','center')
    box off
    set(gca,'FontSize',11)

    % Saving
    fprintf('Printing %d. of %d recording channels.\n',el,length(off_ccep.label))
    print(fullfile(res_dir,['Stim_' elec_pair '_' off_ccep.label{el} '_peaks.png']),'-dpng','-r100')

    close all
end

% Save the calculated peak info for later use
save(fullfile(res_dir,['peak_info_' elec_pair '.mat']),'peak_info')

end