function CCEP_offon_collect_drectns(Sbj_Metadata, offmed_block, onmed_block)
% CCEP_offon_collect: to collect and see what changes in what percentages
% in CCEP N1 and N2 latencies and amplitudes in which direction based on
% stimulation and recording electrode character. This version includes the
% directions of how they change as well.

%% on which electrodes were the analysis was run:
% 'peak_info_LA1-LA2.mat'
fileList = dir(fullfile(Sbj_Metadata.results,'off_on_peaks','peak_info_*.mat'));
analyzed_chans = cellfun(@(x) erase(x, {'peak_info_','.mat'}), {fileList.name}', 'UniformOutput', false);

%% find good channels on both blocks
ref='bp';
% load a random CCEP to get bipolar names
off_ccep = load(fullfile(Sbj_Metadata.iEEG_data,offmed_block, 'CCEPdir', [offmed_block '_ccep_' ref '_' analyzed_chans{1} '.mat']));
% on_ccep = load(fullfile(Sbj_Metadata.iEEG_data,onmed_block, 'CCEPdir', [onmed_block '_ccep_' ref '_' StimCh1 '-' StimCh2 '.mat']));
bp_chans = off_ccep.CCEPs.label;clear off_ccep

% now find non-artefactual chans in bipolar data
off_info = load(fullfile(Sbj_Metadata.iEEG_data, offmed_block, [offmed_block '_info.mat']));
on_info = load(fullfile(Sbj_Metadata.iEEG_data, onmed_block, [onmed_block '_info.mat']));

bp_good_chans_off = get_info_goodchans_bp(off_info.info,bp_chans);
bp_good_chans_on = get_info_goodchans_bp(on_info.info,bp_chans);

bp_good_chans = bp_good_chans_off(ismember(bp_good_chans_off,bp_good_chans_on));
clear bp_good_chans_off bp_good_chans_on

%% get channel characters
% % find characters of good channels
importance = {'SOZ','EPZ','IZ','Healthy'};
bp_good_chans_split = strsplit_SA(bp_good_chans);
[bp_good_chans_classes] = assign_classes(bp_good_chans_split(:,1), bp_good_chans_split(:,2), off_info.info.channelinfo.Label, off_info.info.channelinfo.chan_info,importance);

% now find character of each pair of analyzed electrodes
events_chans = strcat(off_info.info.events.StimCh1,repmat({'-'},[length(off_info.info.events.StimCh1),1]),off_info.info.events.StimCh2);
analyzed_chans_character = cell([length(analyzed_chans),1]);
for c = 1:length(analyzed_chans_character)
    xxx=find(ismember(events_chans,analyzed_chans{c}));
    analyzed_chans_character(c) = off_info.info.events.elec_type(xxx(1));
end

%% is there any analyzed pair that is actually bad: remove them
bad_analyzed_chans_idx = ~ismember(analyzed_chans,bp_good_chans);
bad_analyzed_chans = analyzed_chans(bad_analyzed_chans_idx);
if ~isempty(bad_analyzed_chans)
warning('%s: these channels are pairs as bad, but analyses were run nevertheless, removing from here.',strjoin(bad_analyzed_chans,','))
end
analyzed_chans = analyzed_chans(~bad_analyzed_chans_idx);
analyzed_chans_character = analyzed_chans_character(~bad_analyzed_chans_idx);
analyzed_chans_split = strsplit_SA(analyzed_chans);

%%

% create criss-cross
for_crisscross = {'SOZ','EPZ','IZ','Healthy'};
for x1=1:4
    for x2=1:4
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N1_time = [];'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp = [];'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N2_time = [];'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp = [];'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N1_time_dir = [];'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp_dir = [];'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N2_time_dir = [];'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp_dir = [];'])
    end
end

% loop stim elec based on their seizure character
for c = 1:length(analyzed_chans)

    % load CCEP_comp_offon stat output
    load(fullfile(Sbj_Metadata.results,'off_on_peaks',['peak_info_' analyzed_chans{c} '.mat']),'peak_info')
    % first remove bad channels
    peak_info = peak_info(ismember(peak_info.label,bp_good_chans),:);
    % remove stimulation channels
    peak_info_spl = strsplit_SA(peak_info.label);
    peak_info = peak_info(~any(ismember(peak_info_spl,analyzed_chans_split(c,:)),2),:);
    peak_info_spl = strsplit_SA(peak_info.label);

    % find classes for recording channels
    [peakinfo_chan_classes] = assign_classes(peak_info_spl(:,1), peak_info_spl(:,2), off_info.info.channelinfo.Label, off_info.info.channelinfo.chan_info,importance);
    
    % assign direction of change for each
    peak_info.N1_time_dir = zeros([height(peak_info),1]);peak_info.N1_amp_dir = zeros([height(peak_info),1]);
    peak_info.N2_time_dir = zeros([height(peak_info),1]);peak_info.N2_amp_dir = zeros([height(peak_info),1]);

    peak_info.N1_time_dir(nanmean(cell2mat(peak_info.off_N1_time_all),2) < nanmean(cell2mat(peak_info.on_N1_time_all),2)) = 1; % on>off = 1
    peak_info.N1_time_dir(nanmean(cell2mat(peak_info.off_N1_time_all),2) > nanmean(cell2mat(peak_info.on_N1_time_all),2)) = 2; % on<off = 2
    peak_info.N2_time_dir(nanmean(cell2mat(peak_info.off_N2_time_all),2) < nanmean(cell2mat(peak_info.on_N2_time_all),2)) = 1;
    peak_info.N2_time_dir(nanmean(cell2mat(peak_info.off_N2_time_all),2) > nanmean(cell2mat(peak_info.on_N2_time_all),2)) = 2;
    peak_info.N1_amp_dir((nanmean(cell2mat(peak_info.off_N1_amp_all),2) - nanmean(cell2mat(peak_info.on_N1_amp_all),2))>0) = 1;
    peak_info.N1_amp_dir((nanmean(cell2mat(peak_info.off_N1_amp_all),2) - nanmean(cell2mat(peak_info.on_N1_amp_all),2))<0) = 2;
    peak_info.N2_amp_dir((nanmean(cell2mat(peak_info.off_N2_amp_all),2) - nanmean(cell2mat(peak_info.on_N2_amp_all),2))>0) = 1;
    peak_info.N2_amp_dir((nanmean(cell2mat(peak_info.off_N2_amp_all),2) - nanmean(cell2mat(peak_info.on_N2_amp_all),2))<0) = 2;

    for p = 1:height(peak_info)
        % first gather the significance
        eval([analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N1_time = [' ...
            analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N1_time; peak_info.ttestresults{p}.N1_time.h];'] )
        eval([analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N1_amp = [' ...
            analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N1_amp; peak_info.ttestresults{p}.N1_amp.h];'] )
        eval([analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N2_time = [' ...
            analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N2_time; peak_info.ttestresults{p}.N2_time.h];'] )
        eval([analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N2_amp = [' ...
            analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N2_amp; peak_info.ttestresults{p}.N2_amp.h];'] )

        eval([analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N1_time_dir = [' ...
            analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N1_time_dir; peak_info.N1_time_dir(p)];'] )
        eval([analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N1_amp_dir = [' ...
            analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N1_amp_dir; peak_info.N1_amp_dir(p)];'] )
        eval([analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N2_time_dir = [' ...
            analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N2_time_dir; peak_info.N2_time_dir(p)];'] )
        eval([analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N2_amp_dir = [' ...
            analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N2_amp_dir; peak_info.N2_amp_dir(p)];'] )
        
        % if strcmp([analyzed_chans_character{c} '2' peakinfo_chan_classes{p}], 'EPZ2EPZ')
        %     fprintf('%s to %s\n',analyzed_chans{c},peak_info.label{p})
        % end
        % if strcmp([analyzed_chans_character{c} '2' peakinfo_chan_classes{p}], 'SOZ2SOZ')
        %     fprintf('%s to %s\n',analyzed_chans{c},peak_info.label{p})
        % end
    end
end

% "Inc" is for off-to-on increase
N1_time_percs=[];N1_time_dirInc=[];N1_time_dirDec=[];
N1_amp_percs=[];N1_amp_dirInc=[];N1_amp_dirDec=[];
N2_time_percs=[];N2_time_dirInc=[];N2_time_dirDec=[];
N2_amp_percs=[];N2_amp_dirInc=[];N2_amp_dirDec=[];

for x1=1:4
    for x2=1:4
        % calculate the overall percentages here
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N1_time_perc = '...
            'nansum(' for_crisscross{x1} '2' for_crisscross{x2} '_N1_time) / length(~isnan(' for_crisscross{x1} '2' for_crisscross{x2} '_N1_time));'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp_perc = '...
            'nansum(' for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp) / length(~isnan(' for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp));'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N2_time_perc = '...
            'nansum(' for_crisscross{x1} '2' for_crisscross{x2} '_N2_time) / length(~isnan(' for_crisscross{x1} '2' for_crisscross{x2} '_N2_time));'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp_perc = '...
            'nansum(' for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp) / length(~isnan(' for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp));'])

        eval(['N1_time_percs(x1,x2) = ' for_crisscross{x1} '2' for_crisscross{x2} '_N1_time_perc;'])
        eval(['N1_amp_percs(x1,x2) = ' for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp_perc;'])
        eval(['N2_time_percs(x1,x2) = ' for_crisscross{x1} '2' for_crisscross{x2} '_N2_time_perc;'])
        eval(['N2_amp_percs(x1,x2) = ' for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp_perc;'])

        % process the directions
        % multiply tstat with h value
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N1_time_dir = ' ...
            for_crisscross{x1} '2' for_crisscross{x2} '_N1_time_dir .* ' for_crisscross{x1} '2' for_crisscross{x2} '_N1_time;'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp_dir = ' ...
            for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp_dir .* ' for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp;'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N2_time_dir = ' ...
            for_crisscross{x1} '2' for_crisscross{x2} '_N2_time_dir .* ' for_crisscross{x1} '2' for_crisscross{x2} '_N2_time;'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp_dir = ' ...
            for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp_dir .* ' for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp;'])
        % now calculate the numbers
        eval(['N1_time_dirInc(x1,x2) = sum(' for_crisscross{x1} '2' for_crisscross{x2} '_N1_time_dir == 1)/length('  for_crisscross{x1} '2' for_crisscross{x2} '_N1_time_dir);'])
        eval(['N1_amp_dirInc(x1,x2) = sum(' for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp_dir == 2)/length('  for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp_dir);'])
        eval(['N2_time_dirInc(x1,x2) = sum(' for_crisscross{x1} '2' for_crisscross{x2} '_N2_time_dir == 1)/length('  for_crisscross{x1} '2' for_crisscross{x2} '_N2_time_dir);'])
        eval(['N2_amp_dirInc(x1,x2) = sum(' for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp_dir == 2)/length('  for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp_dir);'])
        eval(['N1_time_dirDec(x1,x2) = sum(' for_crisscross{x1} '2' for_crisscross{x2} '_N1_time_dir == 2)/length('  for_crisscross{x1} '2' for_crisscross{x2} '_N1_time_dir);'])
        eval(['N1_amp_dirDec(x1,x2) = sum(' for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp_dir == 1)/length('  for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp_dir);'])
        eval(['N2_time_dirDec(x1,x2) = sum(' for_crisscross{x1} '2' for_crisscross{x2} '_N2_time_dir == 2)/length('  for_crisscross{x1} '2' for_crisscross{x2} '_N2_time_dir);'])
        eval(['N2_amp_dirDec(x1,x2) = sum(' for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp_dir == 1)/length('  for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp_dir);'])
    end
end

%% 4 Figures for each N1-N2 amp and latency
% save directory
res_dir=fullfile(Sbj_Metadata.results,'offon_collect_drectns');
if ~isfolder(res_dir); mkdir(res_dir); end

%% N1 latency
figure('position',[0 0 1800 400])
subplot(1,3,1)
image(100*N1_time_percs','CDataMapping','scaled')
colormap(master_ColorMaps('hawaii'));
colorbar
xticks(1:4)
xticklabels(for_crisscross)
xlabel('Stimulation Channels')
yticks(1:4)
yticklabels(for_crisscross)
ylabel('Recording Channels')
for x1=1:4
    for x2=1:4
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N1_time_percs(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('All changed combined')

subplot(1,3,2)
image(100*N1_time_dirInc','CDataMapping','scaled')
colormap(master_ColorMaps('hawaii'));
colorbar
clim([0 100*max([N1_time_dirInc,N1_time_dirDec],[],'all')])
xticks(1:4)
xticklabels(for_crisscross)
xlabel('Stimulation Channels')
yticks(1:4)
yticklabels(for_crisscross)
ylabel('Recording Channels')
for x1=1:4
    for x2=1:4
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N1_time_dirInc(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('Increased')

subplot(1,3,3)
image(100*N1_time_dirDec','CDataMapping','scaled')
colormap(master_ColorMaps('hawaii'));
colorbar
clim([0 100*max([N1_time_dirInc,N1_time_dirDec],[],'all')])
xticks(1:4)
xticklabels(for_crisscross)
xlabel('Stimulation Channels')
yticks(1:4)
yticklabels(for_crisscross)
ylabel('Recording Channels')
for x1=1:4
    for x2=1:4
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N1_time_dirDec(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('Decreased')
sgtitle('Change from off to on meds N1-latency')

print(fullfile(res_dir,[Sbj_Metadata.sbj_ID '_offon_dir_N1-latency.png']),'-dpng','-r300')

%% N1 amplitude

figure('position',[0 0 1800 400])
subplot(1,3,1)
image(100*N1_amp_percs','CDataMapping','scaled')
colormap(master_ColorMaps('hawaii'));
colorbar
xticks(1:4)
xticklabels(for_crisscross)
xlabel('Stimulation Channels')
yticks(1:4)
yticklabels(for_crisscross)
ylabel('Recording Channels')
for x1=1:4
    for x2=1:4
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N1_amp_percs(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('All changed combined')

subplot(1,3,2)
image(100*N1_amp_dirInc','CDataMapping','scaled')
colormap(master_ColorMaps('hawaii'));
colorbar
clim([0 100*max([N1_amp_dirInc,N1_amp_dirDec],[],'all')])
xticks(1:4)
xticklabels(for_crisscross)
xlabel('Stimulation Channels')
yticks(1:4)
yticklabels(for_crisscross)
ylabel('Recording Channels')
for x1=1:4
    for x2=1:4
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N1_amp_dirInc(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('Increased')

subplot(1,3,3)
image(100*N1_amp_dirDec','CDataMapping','scaled')
colormap(master_ColorMaps('hawaii'));
colorbar
clim([0 100*max([N1_amp_dirInc,N1_amp_dirDec],[],'all')])
xticks(1:4)
xticklabels(for_crisscross)
xlabel('Stimulation Channels')
yticks(1:4)
yticklabels(for_crisscross)
ylabel('Recording Channels')
for x1=1:4
    for x2=1:4
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N1_amp_dirDec(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('Decreased')
sgtitle('Change from off to on meds N1-Amplitude')

print(fullfile(res_dir,[Sbj_Metadata.sbj_ID '_offon_dir_N1-amplitude.png']),'-dpng','-r300')


%% N2 latency

figure('position',[0 0 1800 400])
subplot(1,3,1)
image(100*N2_time_percs','CDataMapping','scaled')
colormap(master_ColorMaps('hawaii'));
colorbar
xticks(1:4)
xticklabels(for_crisscross)
xlabel('Stimulation Channels')
yticks(1:4)
yticklabels(for_crisscross)
ylabel('Recording Channels')
for x1=1:4
    for x2=1:4
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N2_time_percs(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('All changed combined')

subplot(1,3,2)
image(100*N2_time_dirInc','CDataMapping','scaled')
colormap(master_ColorMaps('hawaii'));
colorbar
clim([0 100*max([N2_time_dirInc,N2_time_dirDec],[],'all')])
xticks(1:4)
xticklabels(for_crisscross)
xlabel('Stimulation Channels')
yticks(1:4)
yticklabels(for_crisscross)
ylabel('Recording Channels')
for x1=1:4
    for x2=1:4
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N2_time_dirInc(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('Increased')

subplot(1,3,3)
image(100*N2_time_dirDec','CDataMapping','scaled')
colormap(master_ColorMaps('hawaii'));
colorbar
clim([0 100*max([N2_time_dirInc,N2_time_dirDec],[],'all')])
xticks(1:4)
xticklabels(for_crisscross)
xlabel('Stimulation Channels')
yticks(1:4)
yticklabels(for_crisscross)
ylabel('Recording Channels')
for x1=1:4
    for x2=1:4
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N2_time_dirDec(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('Decreased')
sgtitle('Change from off to on meds N2-latency')

print(fullfile(res_dir,[Sbj_Metadata.sbj_ID '_offon_dir_N2-latency.png']),'-dpng','-r300')

%% N2 amplitude

figure('position',[0 0 1800 400])
subplot(1,3,1)
image(100*N2_amp_percs','CDataMapping','scaled')
colormap(master_ColorMaps('hawaii'));
colorbar
xticks(1:4)
xticklabels(for_crisscross)
xlabel('Stimulation Channels')
yticks(1:4)
yticklabels(for_crisscross)
ylabel('Recording Channels')
for x1=1:4
    for x2=1:4
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N2_amp_percs(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('All changed combined')

subplot(1,3,2)
image(100*N2_amp_dirInc','CDataMapping','scaled')
colormap(master_ColorMaps('hawaii'));
colorbar
clim([0 100*max([N2_amp_dirInc,N2_amp_dirDec],[],'all')])
xticks(1:4)
xticklabels(for_crisscross)
xlabel('Stimulation Channels')
yticks(1:4)
yticklabels(for_crisscross)
ylabel('Recording Channels')
for x1=1:4
    for x2=1:4
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N2_amp_dirInc(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('Increased')

subplot(1,3,3)
image(100*N2_amp_dirDec','CDataMapping','scaled')
colormap(master_ColorMaps('hawaii'));
colorbar
clim([0 100*max([N2_amp_dirInc,N2_amp_dirDec],[],'all')])
xticks(1:4)
xticklabels(for_crisscross)
xlabel('Stimulation Channels')
yticks(1:4)
yticklabels(for_crisscross)
ylabel('Recording Channels')
for x1=1:4
    for x2=1:4
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N2_amp_dirDec(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('Decreased')
sgtitle('Change from off to on meds N2-Amplitude')

print(fullfile(res_dir,[Sbj_Metadata.sbj_ID '_offon_dir_N2-amplitude.png']),'-dpng','-r300')



%% save the data as well
collected=[];
for x1=1:4
    for x2=1:4
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N1_time.perc = ' for_crisscross{x1} '2' for_crisscross{x2} '_N1_time_perc;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N1_time.all = ' for_crisscross{x1} '2' for_crisscross{x2} '_N1_time;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N1_amp.perc = ' for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp_perc;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N1_amp.all = ' for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N2_time.perc = ' for_crisscross{x1} '2' for_crisscross{x2} '_N2_time_perc;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N2_time.all = ' for_crisscross{x1} '2' for_crisscross{x2} '_N2_time;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N2_amp.perc = ' for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp_perc;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N2_amp.all = ' for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp;'])

        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N1_time_dir.all = ' for_crisscross{x1} '2' for_crisscross{x2} '_N1_time_dir;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N1_time_dir.Inc =  N1_time_dirInc;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N1_time_dir.Dec =  N1_time_dirDec;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N1_amp_dir.all = ' for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp_dir;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N1_amp_dir.Inc =  N1_amp_dirInc;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N1_amp_dir.Dec =  N1_amp_dirDec;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N2_time_dir.all = ' for_crisscross{x1} '2' for_crisscross{x2} '_N2_time_dir;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N2_time_dir.Inc =  N2_time_dirInc;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N2_time_dir.Dec =  N2_time_dirDec;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N2_amp_dir.all = ' for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp_dir;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N2_amp_dir.Inc =  N2_amp_dirInc;'])
        eval(['collected.' for_crisscross{x1} '2' for_crisscross{x2} '.N2_amp_dir.Dec =  N2_amp_dirDec;'])
    end
end

save(fullfile(res_dir,[Sbj_Metadata.sbj_ID '_offon_allDirections.mat']),'collected','bp_good_chans_classes','bp_good_chans','analyzed_chans_character','analyzed_chans')

end

