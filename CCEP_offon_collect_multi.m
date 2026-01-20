function CCEP_offon_collect_multi(Sbj_Metadatas)
% Based on CCEP_offon_collect; to collect and see what changes in what percentages
% in CCEP N1 and N2 latencies and amplitudes in which direction based on
% stimulation and recording electrode character.

comb_sbjIDs = '';
chan_count = 0;
bp_good_chans_classes_all=[];
bad_analyzed_chans_all=[];
bp_chans_all=[];
bp_good_chans_all=[];
analyzed_chans_character_all=[];
analyzed_chans_all=[];

for_crisscross = {'SOZ','EPZ','IZ','Healthy'};
for x1=1:4
    for x2=1:4
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N1_time = [];'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N1_amp = [];'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N2_time = [];'])
        eval([for_crisscross{x1} '2' for_crisscross{x2} '_N2_amp = [];'])
    end
end


%% Loop all subjects to gather analyzed channels
for s = 1:length(Sbj_Metadatas)

Sbj_Metadata = Sbj_Metadatas{s};
sbj_ID = Sbj_Metadata.sbj_ID;
comb_sbjIDs = [comb_sbjIDs '_' sbj_ID];
AllBlockInfo = readtable(fullfile(Sbj_Metadata.project_root,[Sbj_Metadata.project_name '_BlockInfo.xlsx'])); % "F:\HBML\PROJECTS_DATA\CL_Train\CL_Train_BlockInfo.xlsx"

offblocks = find(strcmp(AllBlockInfo.sbj_ID,sbj_ID) & strcmp(AllBlockInfo.task_type,'offmed'));
onblocks = find(strcmp(AllBlockInfo.sbj_ID,sbj_ID) & strcmp(AllBlockInfo.task_type,'onmed'));

offmed_block = AllBlockInfo.BlockList{offblocks(1)};
onmed_block = AllBlockInfo.BlockList{onblocks(1)};

%% on which electrodes were the analysis was run:
% 'peak_info_LA1-LA2.mat'
fileList = dir(fullfile(Sbj_Metadata.results,'off_on_peaks','peak_info_*.mat'));
analyzed_chans = cellfun(@(x) erase(x, {'peak_info_','.mat'}), {fileList.name}', 'UniformOutput', false);

%% find good channels on both blocks
ref='bp';
% load a random CCEP to get bipolar names
off_ccep = load(fullfile(Sbj_Metadata.iEEG_data,offmed_block, 'CCEPdir', [offmed_block '_ccep_' ref '_' analyzed_chans{1} '.mat']));
bp_chans = off_ccep.CCEPs.label;clear off_ccep
bp_chans_all = [bp_chans_all;bp_chans];

% now find non-artefactual chans in bipolar data
off_info = load(fullfile(Sbj_Metadata.iEEG_data, offmed_block, [offmed_block '_info.mat']));
on_info = load(fullfile(Sbj_Metadata.iEEG_data, onmed_block, [onmed_block '_info.mat']));

% maybe change here to include all channels or line 87, 89
bp_good_chans_off = get_info_goodchans_bp(off_info.info,bp_chans);
bp_good_chans_on = get_info_goodchans_bp(on_info.info,bp_chans);

bp_good_chans = bp_good_chans_off(ismember(bp_good_chans_off,bp_good_chans_on));
chan_count = chan_count + height(off_info.info.channelinfo);
bp_good_chans_all = [bp_good_chans_all;bp_good_chans];
clear bp_good_chans_off bp_good_chans_on

%% get channel characters
% % find characters of good channels
importance = {'SOZ','EPZ','IZ','Healthy'};
bp_good_chans_split = strsplit_SA(bp_good_chans);
bp_good_chans_classes_all = [bp_good_chans_classes_all ;...
                            assign_classes(bp_good_chans_split(:,1), bp_good_chans_split(:,2), off_info.info.channelinfo.Label, off_info.info.channelinfo.chan_info,importance)];

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
bad_analyzed_chans_all = [bad_analyzed_chans_all;bad_analyzed_chans];
% analyzed_chans = analyzed_chans(~bad_analyzed_chans_idx);
analyzed_chans_all = [analyzed_chans_all;analyzed_chans];
% analyzed_chans_character = analyzed_chans_character(~bad_analyzed_chans_idx);
analyzed_chans_character_all = [analyzed_chans_character_all;analyzed_chans_character];
analyzed_chans_split = strsplit_SA(analyzed_chans);

%% create criss-cross
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
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N1_time_percs(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
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
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N1_amp_percs(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
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
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N2_time_percs(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
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
        eval(['text(' num2str(x2) ',' num2str(x1) ',[num2str(100*N2_amp_percs(x2,x1),4) ''%''],''Color'',''k'',''HorizontalAlignment'',''center'');'])
    end
end
title('...N2-amplitude')

sgtitle('Percentage of electrodes with significant change in...')

% save directory
res_dir=fullfile(Sbj_Metadata.project_root,'COMBINED_RESULTS','offon_collect_multi');
if ~isfolder(res_dir); mkdir(res_dir); end

print(fullfile(res_dir,[comb_sbjIDs '_offon_all.png']),'-dpng','-r300')

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

save(fullfile(res_dir,[comb_sbjIDs '_offon_all.mat']),'collected','bp_good_chans_classes_all','bp_good_chans_all','analyzed_chans_character_all','analyzed_chans_all')

% print out info about the channels
prt=['Number of total recording channels: ',num2str(chan_count) '\n'];
prt=[prt,'Number of total recording pairs: ',num2str(length(bp_chans_all)) '\n'];
prt=[prt,'Number of total recording non-artefactual pairs: ',num2str(length(bp_good_chans_all)) '\n'];
prt=[prt,'Number of recording pairs included in analysis: ',num2str(length(analyzed_chans_all)) '\n'];
prt=[prt,'Number of recording pairs that were SOZ: ',num2str(sum(strcmp(bp_good_chans_classes_all,'SOZ'))) '\n'];
prt=[prt,'Number of recording pairs that were EPZ: ',num2str(sum(strcmp(bp_good_chans_classes_all,'EPZ'))) '\n'];
prt=[prt,'Number of recording pairs that were IZ: ',num2str(sum(strcmp(bp_good_chans_classes_all,'IZ'))) '\n'];
prt=[prt,'Number of recording pairs that were Healthy: ',num2str(sum(strcmp(bp_good_chans_classes_all,'Healthy'))) '\n'];
prt=[prt,'Number of stimulation pairs that were removed due to being artefactual: ',num2str(length(bad_analyzed_chans_all)) '\n'];
prt=[prt,'Number of analyzed stimulation pairs that were SOZ: ',num2str(sum(strcmp(analyzed_chans_character_all,'SOZ'))) '\n'];
prt=[prt,'Number of analyzed stimulation pairs that were EPZ: ',num2str(sum(strcmp(analyzed_chans_character_all,'EPZ'))) '\n'];
prt=[prt,'Number of analyzed stimulation pairs that were IZ: ',num2str(sum(strcmp(analyzed_chans_character_all,'IZ'))) '\n'];
prt=[prt,'Number of analyzed stimulation pairs that were Healthy: ',num2str(sum(strcmp(analyzed_chans_character_all,'Healthy'))) '\n'];

% save some info as text
fileID = fopen(fullfile(res_dir,[comb_sbjIDs '_offon_all.txt']),'w');
fprintf(fileID,prt);
fclose(fileID);
type(fullfile(res_dir,[comb_sbjIDs '_offon_all.txt']))

end
