function CCEP_offon_collect(Sbj_Metadata, offmed_block, onmed_block)
% CCEP_offon_collect: to collect and see what changes in what percentages
% in CCEP N1 and N2 latencies and amplitudes in which direction based on
% stimulation and recording electrode character.

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

    for p = 1:height(peak_info)
        eval([analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N1_time = [' ...
            analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N1_time; peak_info.ttestresults{p}.N1_time.h];'] )
        eval([analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N1_amp = [' ...
            analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N1_amp; peak_info.ttestresults{p}.N1_amp.h];'] )
        eval([analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N2_time = [' ...
            analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N2_time; peak_info.ttestresults{p}.N2_time.h];'] )
        eval([analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N2_amp = [' ...
            analyzed_chans_character{c} '2' peakinfo_chan_classes{p} '_N2_amp; peak_info.ttestresults{p}.N2_amp.h];'] )
    end

end

N1_time_percs=[];
N1_amp_percs=[];
N2_time_percs=[];
N2_amp_percs=[];

for x1=1:4
    for x2=1:4
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
    end
end

figure('position',[0 0 1200 1000])
subplot(2,2,1)
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
        eval(['text(' num2str(x2) ',' num2str(x1) ',[''%'' num2str(100*N1_time_percs(x2,x1),4)],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('...N1-latency')


subplot(2,2,2)
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
        eval(['text(' num2str(x2) ',' num2str(x1) ',[''%'' num2str(100*N1_amp_percs(x2,x1),4)],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('...N1-amplitude')

subplot(2,2,3)
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
        eval(['text(' num2str(x2) ',' num2str(x1) ',[''%'' num2str(100*N2_time_percs(x2,x1),4)],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('...N2-latency')

subplot(2,2,4)
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
        eval(['text(' num2str(x2) ',' num2str(x1) ',[''%'' num2str(100*N2_amp_percs(x2,x1),4)],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('...N2-amplitude')

sgtitle('Percentage of electrodes with significant change in...')

% save directory
res_dir=fullfile(Sbj_Metadata.results,'offon_collect');
if ~isfolder(res_dir); mkdir(res_dir); end

print(fullfile(res_dir,[Sbj_Metadata.sbj_ID '_offon_all.png']),'-dpng','-r300')

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
    end
end

save(fullfile(res_dir,[Sbj_Metadata.sbj_ID '_offon_all.mat']),'collected','bp_good_chans_classes','bp_good_chans','analyzed_chans_character','analyzed_chans')

% print out info about the channels
prt=['Number of total recording channels: ',num2str(height(off_info.info.channelinfo)) '\n'];
prt=[prt,'Number of total recording pairs: ',num2str(length(bp_chans)) '\n'];
prt=[prt,'Number of total recording non-artefactual pairs: ',num2str(length(bp_good_chans)) '\n'];
prt=[prt,'Number of stimulation pairs included in analysis: ',num2str(length(analyzed_chans)) '\n'];
prt=[prt,'Number of stimulation pairs that were removed due to being artefactual: ',num2str(length(bad_analyzed_chans)) '\n'];
prt=[prt,'Number of stimulation pairs that were EPZ: ',num2str(sum(strcmp(bp_good_chans_classes,'EPZ'))) '\n'];
prt=[prt,'Number of stimulation pairs that were SOZ: ',num2str(sum(strcmp(bp_good_chans_classes,'SOZ'))) '\n'];
prt=[prt,'Number of stimulation pairs that were IZ: ',num2str(sum(strcmp(bp_good_chans_classes,'IZ'))) '\n'];
prt=[prt,'Number of stimulation pairs that were Healthy: ',num2str(sum(strcmp(bp_good_chans_classes,'Healthy'))) '\n'];
prt=[prt,'Number of analyzed stimulation pairs that were EPZ: ',num2str(sum(strcmp(analyzed_chans_character,'EPZ'))) '\n'];
prt=[prt,'Number of analyzed stimulation pairs that were SOZ: ',num2str(sum(strcmp(analyzed_chans_character,'SOZ'))) '\n'];
prt=[prt,'Number of analyzed stimulation pairs that were IZ: ',num2str(sum(strcmp(analyzed_chans_character,'IZ'))) '\n'];
prt=[prt,'Number of analyzed stimulation pairs that were Healthy: ',num2str(sum(strcmp(analyzed_chans_character,'Healthy'))) '\n'];

% save some info as text
fileID = fopen(fullfile(res_dir,[Sbj_Metadata.sbj_ID '_offon_all.txt']),'w');
fprintf(fileID,prt);
fclose(fileID);
type(fullfile(res_dir,[Sbj_Metadata.sbj_ID '_offon_all.txt']))

end




% if strcmp(analyzed_chans_character{c},'EPZ')
%     % EPZ 2 EPZ
%     EPZ2EPZ_off_N1_amp = [EPZ2EPZ_off_N1_amp;peak_info.ttestresults(strcmp(peakinfo_chan_classes,'EPZ'))];
%     EPZ2EPZ_off_N1_time = [EPZ2EPZ_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'EPZ'))];
%     EPZ2EPZ_off_N2_amp = [EPZ2EPZ_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'EPZ'))];
%     EPZ2EPZ_off_N2_time = [EPZ2EPZ_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'EPZ'))];
%     % EPZ 2 IZ
%     EPZ2IZ_off_N1_amp = [EPZ2IZ_off_N1_amp;peak_info.off_N1_amp(strcmp(peakinfo_chan_classes,'IZ'))];
%     EPZ2IZ_off_N1_time = [EPZ2IZ_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'IZ'))];
%     EPZ2IZ_off_N2_amp = [EPZ2IZ_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'IZ'))];
%     EPZ2IZ_off_N2_time = [EPZ2IZ_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'IZ'))];
%     % EPZ 2 SOZ
%     EPZ2SOZ_off_N1_amp = [EPZ2SOZ_off_N1_amp;peak_info.off_N1_amp(strcmp(peakinfo_chan_classes,'SOZ'))];
%     EPZ2SOZ_off_N1_time = [EPZ2SOZ_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'SOZ'))];
%     EPZ2SOZ_off_N2_amp = [EPZ2SOZ_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'SOZ'))];
%     EPZ2SOZ_off_N2_time = [EPZ2SOZ_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'SOZ'))];
%     % EPZ 2 healthy
%     EPZ2Heathy_off_N1_amp = [EPZ2Heathy_off_N1_amp;peak_info.off_N1_amp(strcmp(peakinfo_chan_classes,'Healthy'))];
%     EPZ2Healthy_off_N1_time = [EPZ2Healthy_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'Healthy'))];
%     EPZ2Healthy_off_N2_amp = [EPZ2Healthy_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'Healthy'))];
%     EPZ2Healthy_off_N2_time = [EPZ2Healthy_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'Healthy'))];
% elseif strcmp(analyzed_chans_character{c},'IZ')
%     % IZ 2 EPZ
%     IZ2EPZ_off_N1_amp = [IZ2EPZ_off_N1_amp;peak_info.off_N1_amp(strcmp(peakinfo_chan_classes,'EPZ'))];
%     IZ2EPZ_off_N1_time = [IZ2EPZ_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'EPZ'))];
%     IZ2EPZ_off_N2_amp = [IZ2EPZ_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'EPZ'))];
%     IZ2EPZ_off_N2_time = [IZ2EPZ_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'EPZ'))];
%     % IZ 2 IZ
%     IZ2IZ_off_N1_amp = [IZ2IZ_off_N1_amp;peak_info.off_N1_amp(strcmp(peakinfo_chan_classes,'IZ'))];
%     IZ2IZ_off_N1_time = [IZ2IZ_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'IZ'))];
%     IZ2IZ_off_N2_amp = [IZ2IZ_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'IZ'))];
%     IZ2IZ_off_N2_time = [IZ2IZ_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'IZ'))];
%     % IZ 2 SOZ
%     IZ2SOZ_off_N1_amp = [IZ2SOZ_off_N1_amp;peak_info.off_N1_amp(strcmp(peakinfo_chan_classes,'SOZ'))];
%     IZ2SOZ_off_N1_time = [IZ2SOZ_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'SOZ'))];
%     IZ2SOZ_off_N2_amp = [IZ2SOZ_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'SOZ'))];
%     IZ2SOZ_off_N2_time = [IZ2SOZ_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'SOZ'))];
%     % IZ 2 healthy
%     IZ2Heathy_off_N1_amp = [IZ2Heathy_off_N1_amp;peak_info.off_N1_amp(strcmp(peakinfo_chan_classes,'Healthy'))];
%     IZ2Healthy_off_N1_time = [IZ2Healthy_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'Healthy'))];
%     IZ2Healthy_off_N2_amp = [IZ2Healthy_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'Healthy'))];
%     IZ2Healthy_off_N2_time = [IZ2Healthy_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'Healthy'))];
% elseif strcmp(analyzed_chans_character{c},'SOZ')
%     % SOZ 2 EPZ
%     SOZ2EPZ_off_N1_amp = [SOZ2EPZ_off_N1_amp;peak_info.off_N1_amp(strcmp(peakinfo_chan_classes,'EPZ'))];
%     SOZ2EPZ_off_N1_time = [SOZ2EPZ_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'EPZ'))];
%     SOZ2EPZ_off_N2_amp = [SOZ2EPZ_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'EPZ'))];
%     SOZ2EPZ_off_N2_time = [SOZ2EPZ_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'EPZ'))];
%     % SOZ 2 IZ
%     SOZ2IZ_off_N1_amp = [SOZ2IZ_off_N1_amp;peak_info.off_N1_amp(strcmp(peakinfo_chan_classes,'IZ'))];
%     SOZ2IZ_off_N1_time = [SOZ2IZ_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'IZ'))];
%     SOZ2IZ_off_N2_amp = [SOZ2IZ_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'IZ'))];
%     SOZ2IZ_off_N2_time = [SOZ2IZ_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'IZ'))];
%     % SOZ 2 SOZ
%     SOZ2SOZ_off_N1_amp = [SOZ2SOZ_off_N1_amp;peak_info.off_N1_amp(strcmp(peakinfo_chan_classes,'SOZ'))];
%     SOZ2SOZ_off_N1_time = [SOZ2SOZ_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'SOZ'))];
%     SOZ2SOZ_off_N2_amp = [SOZ2SOZ_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'SOZ'))];
%     SOZ2SOZ_off_N2_time = [SOZ2SOZ_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'SOZ'))];
%     % SOZ 2 healthy
%     SOZ2Heathy_off_N1_amp = [SOZ2Heathy_off_N1_amp;peak_info.off_N1_amp(strcmp(peakinfo_chan_classes,'Healthy'))];
%     SOZ2Healthy_off_N1_time = [SOZ2Healthy_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'Healthy'))];
%     SOZ2Healthy_off_N2_amp = [SOZ2Healthy_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'Healthy'))];
%     SOZ2Healthy_off_N2_time = [SOZ2Healthy_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'Healthy'))];
% elseif strcmp(analyzed_chans_character{c},'Healthy')
%     % Healthy 2 EPZ
%     Healthy2EPZ_off_N1_amp = [Healthy2EPZ_off_N1_amp;peak_info.off_N1_amp(strcmp(peakinfo_chan_classes,'EPZ'))];
%     Healthy2EPZ_off_N1_time = [Healthy2EPZ_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'EPZ'))];
%     Healthy2EPZ_off_N2_amp = [Healthy2EPZ_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'EPZ'))];
%     Healthy2EPZ_off_N2_time = [Healthy2EPZ_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'EPZ'))];
%     % Healthy 2 IZ
%     Healthy2IZ_off_N1_amp = [Healthy2IZ_off_N1_amp;peak_info.off_N1_amp(strcmp(peakinfo_chan_classes,'IZ'))];
%     Healthy2IZ_off_N1_time = [Healthy2IZ_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'IZ'))];
%     Healthy2IZ_off_N2_amp = [Healthy2IZ_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'IZ'))];
%     Healthy2IZ_off_N2_time = [Healthy2IZ_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'IZ'))];
%     % Healthy 2 SOZ
%     Healthy2SOZ_off_N1_amp = [Healthy2SOZ_off_N1_amp;peak_info.off_N1_amp(strcmp(peakinfo_chan_classes,'SOZ'))];
%     Healthy2SOZ_off_N1_time = [Healthy2SOZ_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'SOZ'))];
%     Healthy2SOZ_off_N2_amp = [Healthy2SOZ_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'SOZ'))];
%     Healthy2SOZ_off_N2_time = [Healthy2SOZ_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'SOZ'))];
%     % Healthy 2 healthy
%     Healthy2Heathy_off_N1_amp = [Healthy2Heathy_off_N1_amp;peak_info.off_N1_amp(strcmp(peakinfo_chan_classes,'Healthy'))];
%     Healthy2Healthy_off_N1_time = [Healthy2Healthy_off_N1_time;peak_info.off_N1_time(strcmp(peakinfo_chan_classes,'Healthy'))];
%     Healthy2Healthy_off_N2_amp = [Healthy2Healthy_off_N2_amp;peak_info.off_N2_amp(strcmp(peakinfo_chan_classes,'Healthy'))];
%     Healthy2Healthy_off_N2_time = [Healthy2Healthy_off_N2_time;peak_info.off_N2_time(strcmp(peakinfo_chan_classes,'Healthy'))];
% end