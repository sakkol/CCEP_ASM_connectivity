function CCEP_offon_RMS_multi(Sbj_Metadatas)
% CCEP_offon_RMS: to collect and see what changes in what percentages
% in CCEP RMS for overall (0.01-0.3), around N1 (0.01-0.3) and around N2
% (0.085-0.25) based on stimulation and recording electrode character.

comb_sbjIDs = '';
chan_count = 0;
bp_good_chans_classes_all=[];
bad_analyzed_chans_all=[];
bp_chans_all=[];
bp_good_chans_all=[];
analyzed_chans_character_all=[];
analyzed_chans_all=[];

% create criss-cross
for_crisscross = {'SOZ','EPZ','IZ','Healthy'};
foronoff={'on','off'};
collected=struct;
for i1 = 1:4
    for i2 = 1:4
        for ii = 1:2
            collected.([foronoff{ii} '_' for_crisscross{i1} '2' for_crisscross{i2} 'rms']) =[];
            collected.([foronoff{ii} '_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN1']) =[];
            collected.([foronoff{ii} '_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN2']) =[];
        end
    end
end
ttestresults=struct;
ttestresultsN1=struct;
ttestresultsN2=struct;
alpha_c = 0.05/(3*4*4*565); % there are total of 565 bp channels across patients

%% start the loop for sbj metadata
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
% on_ccep = load(fullfile(Sbj_Metadata.iEEG_data,onmed_block, 'CCEPdir', [onmed_block '_ccep_' ref '_' StimCh1 '-' StimCh2 '.mat']));
bp_chans = off_ccep.CCEPs.label;clear off_ccep
bp_chans_all = [bp_chans_all;bp_chans];

% now find non-artefactual chans in bipolar data
off_info = load(fullfile(Sbj_Metadata.iEEG_data, offmed_block, [offmed_block '_info.mat']));
on_info = load(fullfile(Sbj_Metadata.iEEG_data, onmed_block, [onmed_block '_info.mat']));

bp_good_chans_off = get_info_goodchans_bp(off_info.info,bp_chans);
bp_good_chans_on = get_info_goodchans_bp(on_info.info,bp_chans);

bp_good_chans = bp_good_chans_off(ismember(bp_good_chans_off,bp_good_chans_on));
bp_good_chans_all = [bp_good_chans_all;bp_good_chans];
clear bp_good_chans_off bp_good_chans_on

%% get channel characters
importance = {'SOZ','EPZ','IZ','Healthy'};
% % find characters of good channels
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

%% Calculate RMS for each pair of stim-rec

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

    % gather RMS values
        for i2 = 1:4
            for ii = 1:2
                collected.([foronoff{ii} '_' analyzed_chans_character{c} '2' for_crisscross{i2} 'rms']) = [collected.([foronoff{ii} '_' analyzed_chans_character{c} '2' for_crisscross{i2} 'rms']);...
                    peak_info.([foronoff{ii} '_rms'])(strcmp(peakinfo_chan_classes, for_crisscross{i2}))];
                collected.([foronoff{ii} '_' analyzed_chans_character{c} '2' for_crisscross{i2} 'rmsN1']) = [collected.([foronoff{ii} '_' analyzed_chans_character{c} '2' for_crisscross{i2} 'rmsN1']);...
                    peak_info.([foronoff{ii} '_rmsN1'])(strcmp(peakinfo_chan_classes, for_crisscross{i2}))];
                collected.([foronoff{ii} '_' analyzed_chans_character{c} '2' for_crisscross{i2} 'rmsN2']) = [collected.([foronoff{ii} '_' analyzed_chans_character{c} '2' for_crisscross{i2} 'rmsN2']);...
                    peak_info.([foronoff{ii} '_rmsN2'])(strcmp(peakinfo_chan_classes, for_crisscross{i2}))];
            end
        end
end

end

% just make it more presentable :D
for i1 = 1:4
    for i2 = 1:4
        for ii = 1:2
            collected.([foronoff{ii} '_' for_crisscross{i1} '2' for_crisscross{i2} 'rms']) = cell2mat(collected.([foronoff{ii} '_' for_crisscross{i1} '2' for_crisscross{i2} 'rms']));
            collected.([foronoff{ii} '_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN1']) = cell2mat(collected.([foronoff{ii} '_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN1']));
            collected.([foronoff{ii} '_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN2']) = cell2mat(collected.([foronoff{ii} '_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN2']));
        end
    end
end
p_rms=zeros([4,4]);
p_rmsN1=zeros([4,4]);
p_rmsN2=zeros([4,4]);
h_rms=zeros([4,4]);
h_rmsN1=zeros([4,4]);
h_rmsN2=zeros([4,4]);

for i1 = 1:4
    for i2 = 1:4
        % for RMS
        [significant,p,~,stats] = ttest2(collected.(['off_' for_crisscross{i1} '2' for_crisscross{i2} 'rms']),...
            collected.(['on_' for_crisscross{i1} '2' for_crisscross{i2} 'rms']),'alpha',alpha_c,'tail','both');
        ttestresults.([for_crisscross{i1} '2' for_crisscross{i2}]).h=significant;
        ttestresults.([for_crisscross{i1} '2' for_crisscross{i2}]).p=p;
        ttestresults.([for_crisscross{i1} '2' for_crisscross{i2}]).stats=stats;
        ttestresults.([for_crisscross{i1} '2' for_crisscross{i2}]).dir=nanmean(collected.(['on_' for_crisscross{i1} '2' for_crisscross{i2} 'rms'])-...
            collected.(['off_' for_crisscross{i1} '2' for_crisscross{i2} 'rms']))/abs(nanmean(collected.(['on_' for_crisscross{i1} '2' for_crisscross{i2} 'rms'])-...
            collected.(['off_' for_crisscross{i1} '2' for_crisscross{i2} 'rms']))); % if 1, on is more; if -1, off is more
        if ~isempty(p)
            p_rms(i1,i2)=p;
            if significant
                h_rms(i1,i2)=ttestresults.([for_crisscross{i1} '2' for_crisscross{i2}]).dir;
            else
                h_rms(i1,i2)=0;
            end
        else
            p_rms(i1,i2)=NaN;
            h_rms(i1,i2)=NaN;
        end

        % for RMS-N1
        [significant,p,~,stats] = ttest2(collected.(['off_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN1']),...
            collected.(['on_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN1']),'alpha',alpha_c,'tail','both');
        ttestresultsN1.([for_crisscross{i1} '2' for_crisscross{i2}]).h=significant;
        ttestresultsN1.([for_crisscross{i1} '2' for_crisscross{i2}]).p=p;
        ttestresultsN1.([for_crisscross{i1} '2' for_crisscross{i2}]).stats=stats;
        ttestresultsN1.([for_crisscross{i1} '2' for_crisscross{i2}]).dir=nanmean(collected.(['on_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN1'])-...
            collected.(['off_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN1']))/abs(nanmean(collected.(['on_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN1'])-...
            collected.(['off_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN1'])));
        if ~isempty(p)
            p_rmsN1(i1,i2)=p;
            if significant
                h_rmsN1(i1,i2)=ttestresultsN1.([for_crisscross{i1} '2' for_crisscross{i2}]).dir;
            else
                h_rmsN1(i1,i2)=0;
            end
        else
            p_rmsN1(i1,i2)=NaN;
            h_rmsN1(i1,i2)=NaN;
        end

        % for RMS-N2
        [significant,p,~,stats] = ttest2(collected.(['off_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN2']),...
            collected.(['on_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN2']),'alpha',alpha_c,'tail','both');
        ttestresultsN2.([for_crisscross{i1} '2' for_crisscross{i2}]).h=significant;
        ttestresultsN2.([for_crisscross{i1} '2' for_crisscross{i2}]).p=p;
        ttestresultsN2.([for_crisscross{i1} '2' for_crisscross{i2}]).stats=stats;
        ttestresultsN2.([for_crisscross{i1} '2' for_crisscross{i2}]).dir=nanmean(collected.(['on_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN2'])-...
            collected.(['off_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN2']))/abs(nanmean(collected.(['on_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN2'])-...
            collected.(['off_' for_crisscross{i1} '2' for_crisscross{i2} 'rmsN2'])));
        if ~isempty(p)
            p_rmsN2(i1,i2)=p;
            if significant
                h_rmsN2(i1,i2)=ttestresultsN2.([for_crisscross{i1} '2' for_crisscross{i2}]).dir;
            else
                h_rmsN2(i1,i2)=0;
            end
        else
            p_rmsN2(i1,i2)=NaN;
            h_rmsN2(i1,i2)=NaN;
        end
    end
end


%% now plot the p-values
figure('position',[0 0 1800 400])

subplot(1,3,1)
image(p_rms','CDataMapping','scaled')
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
        text(x2,x1,num2str(p_rms(x2,x1),4),'Color','k','HorizontalAlignment','center');

        if h_rms(x2,x1) == 1
            rectangle('Position', [x2-.5, x1-.5, 1, 1], 'EdgeColor', 'r', 'LineWidth', 2);
        elseif h_rms(x2,x1) == -1
            rectangle('Position', [x2-.5, x1-.5, 1, 1], 'EdgeColor', 'b', 'LineWidth', 2);
        end
    end
end
title('RMS (10-300ms)')

subplot(1,3,2)
image(p_rmsN1','CDataMapping','scaled')
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
        text(x2,x1,num2str(p_rmsN1(x2,x1),4),'Color','k','HorizontalAlignment','center');
        if h_rmsN1(x2,x1) == 1
            rectangle('Position', [x2-.5, x1-.5, 1, 1], 'EdgeColor', 'r', 'LineWidth', 2);
        elseif h_rmsN1(x2,x1) == -1
            rectangle('Position', [x2-.5, x1-.5, 1, 1], 'EdgeColor', 'b', 'LineWidth', 2);
        end
    end
end
title('RMS (10-30ms)')

subplot(1,3,3)
image(p_rmsN2','CDataMapping','scaled')
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
        text(x2,x1,num2str(p_rmsN2(x2,x1),4),'Color','k','HorizontalAlignment','center');
        if h_rmsN2(x2,x1) == 1
            rectangle('Position', [x2-.5, x1-.5, 1, 1], 'EdgeColor', 'r', 'LineWidth', 2);
        elseif h_rmsN2(x2,x1) == -1
            rectangle('Position', [x2-.5, x1-.5, 1, 1], 'EdgeColor', 'b', 'LineWidth', 2);
        end
    end
end
title('RMS (85-250ms)')

sgtitle({'p-values of off vs on comparison for RMS/RMS-N1/RMS-N2 changes';...
        'Direction of change from off to on: blue=decreased, red=increased '})

% save directory
res_dir=fullfile(Sbj_Metadata.project_root,'COMBINED_RESULTS','offon_RMS_multi');
if ~isfolder(res_dir); mkdir(res_dir); end

print(fullfile(res_dir,[comb_sbjIDs '_offon_RMS.png']),'-dpng','-r300')

%% save the data as well

save(fullfile(res_dir,[comb_sbjIDs '_offon_RMS.mat']),'collected','ttestresults','ttestresultsN1','ttestresultsN2','p_rms','p_rmsN2','p_rmsN1','h_rms','h_rmsN2','h_rmsN1')


end