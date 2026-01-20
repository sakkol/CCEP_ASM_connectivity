function [shared_pairs_offblocks, shared_pairs_onblocks, shared_pairs,sp_types] = CCEP_findSharedPairs(Sbj_Metadata)

AllBlockInfo = readtable(fullfile(Sbj_Metadata.project_root,[Sbj_Metadata.project_name '_BlockInfo.xlsx'])); % "F:\HBML\PROJECTS_DATA\CL_Train\CL_Train_BlockInfo.xlsx"

off_uniq_pairs = [];off_block_names=[];
on_uniq_pairs = [];on_block_names=[];
for bl = 1:length(Sbj_Metadata.BlockLists)
    block_type = AllBlockInfo.task_type(strcmp(AllBlockInfo.sbj_ID,sbj_ID) & strcmp(AllBlockInfo.BlockList,Sbj_Metadata.BlockLists{bl}));

    blockname = Sbj_Metadata.BlockLists{bl};
    block_info_file = load(fullfile(Sbj_Metadata.iEEG_data, blockname, [blockname '_info.mat']));

    all_pairs=cellfun(@(x,y) strcat(x, '-', y), block_info_file.info.events.StimCh1, block_info_file.info.events.StimCh2, 'UniformOutput', false);
    if strcmp(block_type,'offmed')
        off_uniq_pairs = [off_uniq_pairs;unique(all_pairs)];
        off_block_names = [off_block_names;repmat({blockname}, size(unique(all_pairs),1), 1)];
    elseif strcmp(block_type,'onmed')
        on_uniq_pairs = [on_uniq_pairs;unique(all_pairs)];
        on_block_names = [on_block_names;repmat({blockname}, size(unique(all_pairs),1), 1)];
    end
end

% offmed_block = Sbj_Metadata.BlockLists{1};
% onmed_block = Sbj_Metadata.BlockLists{2};
% 
% off_info = load(fullfile(Sbj_Metadata.iEEG_data, offmed_block, [offmed_block '_info.mat']));
% on_info = load(fullfile(Sbj_Metadata.iEEG_data, onmed_block, [onmed_block '_info.mat']));
% 
% all_pairs=cellfun(@(x,y) strcat(x, '-', y), off_info.info.events.StimCh1, off_info.info.events.StimCh2, 'UniformOutput', false);
% off_uniq_pairs = unique(all_pairs);
% all_pairs=cellfun(@(x,y) strcat(x, '-', y), on_info.info.events.StimCh1, on_info.info.events.StimCh2, 'UniformOutput', false);
% on_uniq_pairs = unique(all_pairs);

shared_pairs = off_uniq_pairs(ismember(off_uniq_pairs,on_uniq_pairs));
shared_pairs_offblocks=[];shared_pairs_onblocks=[];
for s = 1:length(shared_pairs)
    shared_pairs_offblocks{s,1} = off_block_names(strcmp(off_uniq_pairs,shared_pairs{s}));
    shared_pairs_onblocks{s,1} = on_block_names(strcmp(on_uniq_pairs,shared_pairs{s}));
end

sp_split = strsplit_SA(shared_pairs);
importance = {'SOZ','EPZ','IZ','Healthy'};
sp_types = assign_classes(sp_split(:,1), sp_split(:,2), block_info_file.info.channelinfo.Label, block_info_file.info.channelinfo.chan_info,importance);

end