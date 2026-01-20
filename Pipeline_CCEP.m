%% Analyses_pipeline
% This file contains the collection of functions that need to be run prior to creating figures.
% Run this script where the "neural_data" folder is. 

%% Main loop to create folder structure and folder metadata
data_root = fullfile(cd,'neural_data');
project_name = 'CCEP_PrePost';
sbj_IDs = {'P1', 'P2', 'P3','P4'};
Sbj_Metadatas = cell(size(sbj_IDs));
for s = 1:length(sbj_IDs)
    sbj_ID = sbj_IDs{s};
    Sbj_Metadatas{s} = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
end
AllBlockInfo = readtable(fullfile(Sbj_Metadata.project_root,[Sbj_Metadata.project_name '_BlockInfo.xlsx'])); % "F:\HBML\PROJECTS_DATA\CL_Train\CL_Train_BlockInfo.xlsx"

%% Loop subjects to run several analyses together on individual basis

for s = 1:length(sbj_IDs)
    Sbj_Metadata = Sbj_Metadatas{s};

    % extracts time periods around CCEPs from BP referenced iEEG data
    CCEP_processing(Sbj_Metadata)

    % compare ASM-ON vs ASM-OFF in individual channels
    [shared_pairs_offblocks, shared_pairs_onblocks, shared_pairs,sp_types] = CCEP_findSharedPairs(Sbj_Metadata);
    for st = 1:length(shared_pairs)
        CCEP_comp_offon(Sbj_Metadata, shared_pairs_offblocks{st}, shared_pairs_onblocks{st}, shared_pairs{st},sp_types{st})
    end

    % choosing the blocks to analyze
    off_block_names=[];on_block_names=[];
    for bl = 1:length(Sbj_Metadata.BlockLists)
        block_type = AllBlockInfo.task_type(strcmp(AllBlockInfo.sbj_ID,sbj_ID) & strcmp(AllBlockInfo.BlockList,Sbj_Metadata.BlockLists{bl}));
        blockname = Sbj_Metadata.BlockLists{bl};
        if strcmp(block_type,'offmed')
            off_block_names = [off_block_names;{blockname}];
        elseif strcmp(block_type,'onmed')
            on_block_names = [on_block_names;{blockname}];
        end
    end
    offmed_block = off_block_names{1};
    onmed_block = on_block_names{1};

    % comparing N1/N2 amplitude and latency (combine across channels)
    CCEP_offon_collect(Sbj_Metadata, offmed_block, onmed_block)

    % off vs on RMS comparison
    CCEP_offon_RMS(Sbj_Metadata, offmed_block, onmed_block)

    % to compare the directions of N1/N2 related changes
    CCEP_offon_collect_drectns(Sbj_Metadata, offmed_block, onmed_block)

    if strcmp(sbj_ID,'P3') % there are two blocks of on-med, thus the second run here
        offmed_block = off_block_names{1}; onmed_block = on_block_names{2};
        CCEP_offon_collect(Sbj_Metadata, offmed_block, onmed_block)
        CCEP_offon_RMS(Sbj_Metadata, offmed_block, onmed_block)
        CCEP_offon_collect_drectns(Sbj_Metadata, offmed_block, onmed_block)
    end

end

%% REST OF THIS TO COMBINE THE RESULTS ACROSS PATIENTS

%% CCEP_offon_collect_multi
% to run multiple patients' data together for comparing N1/N2 amplitude and latency.
CCEP_offon_collect_multi(Sbj_Metadatas)

%% CCEP_offon_collect_drectns_multi 
% to run multiple patients' data together to compare the directions of N1/N2 related changes.
CCEP_offon_collect_drectns_multi(Sbj_Metadatas)

%% CCEP_offon_RMS_multi
% to run multiple patients' data together for RMS comparison.
CCEP_offon_RMS_multi(Sbj_Metadatas)
