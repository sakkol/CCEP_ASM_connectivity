function CCEP_offon_RMS_multi_v2(Sbj_Metadatas)
% CCEP_offon_RMS_multi_v2: Linear regression to compare RMS between ASM-ON vs ASM-OFF
% analyzing connectivity between stimulated and recorded regions
%
% Model per stim-rec pair: RMS ~ ASM_State + PatientID (main effects only)
% Tests ASM effect on connectivity for each stim-rec combination
% Family-wise Bonferroni correction across all 48 comparisons (16 cells × 3 metrics)
%
% Output: 4×4 heatmaps showing p-values (columns=stim, rows=rec)
% v2: updated to analyze the data with linear regression

% Initialize data collection structures
for_crisscross = {'SOZ','EPZ','IZ','Healthy'};

% Cell arrays to store data for each stim-rec pair
% Dimensions: [4 stim types, 4 rec types]
rms_data_stim_rec = cell(4, 4);
rms_n1_data_stim_rec = cell(4, 4);
rms_n2_data_stim_rec = cell(4, 4);
patient_id_stim_rec = cell(4, 4);

% Initialize cells
for i1 = 1:4
    for i2 = 1:4
        rms_data_stim_rec{i1, i2} = [];
        rms_n1_data_stim_rec{i1, i2} = [];
        rms_n2_data_stim_rec{i1, i2} = [];
        patient_id_stim_rec{i1, i2} = {};
    end
end

%% Start the loop for subject metadata
for s = 1:length(Sbj_Metadatas)
    
    Sbj_Metadata = Sbj_Metadatas{s};
    sbj_ID = Sbj_Metadata.sbj_ID;
    AllBlockInfo = readtable(fullfile(Sbj_Metadata.project_root,[Sbj_Metadata.project_name '_BlockInfo.xlsx']));
    
    offblocks = find(strcmp(AllBlockInfo.sbj_ID,sbj_ID) & strcmp(AllBlockInfo.task_type,'offmed'));
    onblocks = find(strcmp(AllBlockInfo.sbj_ID,sbj_ID) & strcmp(AllBlockInfo.task_type,'onmed'));
    
    offmed_block = AllBlockInfo.BlockList{offblocks(1)};
    onmed_block = AllBlockInfo.BlockList{onblocks(1)};
    
    %% Find analyzed channel pairs
    fileList = dir(fullfile(Sbj_Metadata.results,'off_on_peaks','peak_info_*.mat'));
    analyzed_chans = cellfun(@(x) erase(x, {'peak_info_','.mat'}), {fileList.name}', 'UniformOutput', false);
    
    %% Get good bipolar channels
    ref = 'bp';
    off_ccep = load(fullfile(Sbj_Metadata.iEEG_data,offmed_block, 'CCEPdir', [offmed_block '_ccep_' ref '_' analyzed_chans{1} '.mat']));
    bp_chans = off_ccep.CCEPs.label;
    clear off_ccep
    
    % Find non-artefactual channels in bipolar data
    off_info = load(fullfile(Sbj_Metadata.iEEG_data, offmed_block, [offmed_block '_info.mat']));
    on_info = load(fullfile(Sbj_Metadata.iEEG_data, onmed_block, [onmed_block '_info.mat']));
    
    bp_good_chans_off = get_info_goodchans_bp(off_info.info, bp_chans);
    bp_good_chans_on = get_info_goodchans_bp(on_info.info, bp_chans);
    
    bp_good_chans = bp_good_chans_off(ismember(bp_good_chans_off, bp_good_chans_on));
    
    %% Get channel characters for stimulation pairs
    importance = {'SOZ','EPZ','IZ','Healthy'};
    
    % Find character of each pair of analyzed (stimulated) electrodes
    events_chans = strcat(off_info.info.events.StimCh1,repmat({'-'},[length(off_info.info.events.StimCh1),1]),off_info.info.events.StimCh2);
    analyzed_chans_character = cell([length(analyzed_chans),1]);
    for c = 1:length(analyzed_chans_character)
        xxx=find(ismember(events_chans,analyzed_chans{c}));
        analyzed_chans_character(c) = off_info.info.events.elec_type(xxx(1));
    end
    
    analyzed_chans_split = strsplit_SA(analyzed_chans);
    
    %% Loop through each analyzed stim channel pair
    corrSheet = readtable(Sbj_Metadata.labelfile);
    for c = 1:length(analyzed_chans)
        
        % Load peak info for this stim pair
        load(fullfile(Sbj_Metadata.results,'off_on_peaks',['peak_info_' analyzed_chans{c} '.mat']),'peak_info')
        
        % Remove bad channels
        peak_info = peak_info(ismember(peak_info.label, bp_good_chans), :);
        
        % Remove stimulation channels themselves
        peak_info_spl = strsplit_SA(peak_info.label);
        peak_info = peak_info(~any(ismember(peak_info_spl, analyzed_chans_split(c,:)), 2), :);
        peak_info_spl = strsplit_SA(peak_info.label);
        % remove non-gray channels
        peakinfo_chan_WMvsGM = assign_classes(peak_info_spl(:,1), peak_info_spl(:,2), corrSheet.Label, lower(corrSheet.WMvsGM),{'gray','white','csf','skull'});
        peak_info = peak_info(~strcmp(peakinfo_chan_WMvsGM,'gray'),:);
        peak_info_spl = strsplit_SA(peak_info.label);

        % Get recording channel classes
        peakinfo_chan_classes = assign_classes(peak_info_spl(:,1), peak_info_spl(:,2), ...
            off_info.info.channelinfo.Label, off_info.info.channelinfo.chan_info, importance);
        
        % Extract RMS values (off and on states)
        rms_off = peak_info.off_rms;
        rms_on = peak_info.on_rms;
        rms_n1_off = peak_info.off_rmsN1;
        rms_n1_on = peak_info.on_rmsN1;
        rms_n2_off = peak_info.off_rmsN2;
        rms_n2_on = peak_info.on_rmsN2;
        
        % Get stimulation channel class index
        stim_chan_class_idx = find(strcmp(for_crisscross, analyzed_chans_character{c}));
        
        % For each recording channel class, accumulate data
        for rec_class_idx = 1:4
            rec_class = for_crisscross{rec_class_idx};
            
            % Find indices of recording channels in this class
            rec_idx = strcmp(peakinfo_chan_classes, rec_class);
            
            if sum(rec_idx) > 0
                % Combine OFF and ON data
                rms_combined = [rms_off(rec_idx); rms_on(rec_idx)];
                rms_n1_combined = [rms_n1_off(rec_idx); rms_n1_on(rec_idx)];
                rms_n2_combined = [rms_n2_off(rec_idx); rms_n2_on(rec_idx)];
                
                % ASM state indicator (0=off, 1=on)
                n_off = sum(rec_idx);
                n_on = sum(rec_idx);
                asm_state = [zeros(n_off, 1); ones(n_on, 1)];
                
                % Patient ID
                patient_ids = [repmat({sbj_ID}, n_off, 1); repmat({sbj_ID}, n_on, 1)];
                
                % Store in cell array
                rms_data_stim_rec{stim_chan_class_idx, rec_class_idx} = [rms_data_stim_rec{stim_chan_class_idx, rec_class_idx}; rms_combined];
                rms_n1_data_stim_rec{stim_chan_class_idx, rec_class_idx} = [rms_n1_data_stim_rec{stim_chan_class_idx, rec_class_idx}; rms_n1_combined];
                rms_n2_data_stim_rec{stim_chan_class_idx, rec_class_idx} = [rms_n2_data_stim_rec{stim_chan_class_idx, rec_class_idx}; rms_n2_combined];
                patient_id_stim_rec{stim_chan_class_idx, rec_class_idx} = [patient_id_stim_rec{stim_chan_class_idx, rec_class_idx}; patient_ids];
            end
        end
    end
    
end

%% Fit linear regression for each stim-rec pair and extract p-values
fprintf('\n========== LINEAR REGRESSION RESULTS PER STIM-REC PAIR ==========\n');
fprintf('Model: RMS ~ ASM_State + PatientID (main effects only)\n');
fprintf('Family-wise Bonferroni correction: α = 0.05 / 48 = %.6f\n\n', 0.05/48);

% Initialize p-value matrices for each metric
pval_rms_matrix = nan(4, 4);
pval_rms_n1_matrix = nan(4, 4);
pval_rms_n2_matrix = nan(4, 4);

% Store model results
models_rms = cell(4, 4);
models_rms_n1 = cell(4, 4);
models_rms_n2 = cell(4, 4);

% Alpha for family-wise correction
alpha_fam = 0.05 / 48;

for stim_idx = 1:4
    for rec_idx = 1:4
        
        % Get data for this stim-rec pair
        rms_data = cell2mat(rms_data_stim_rec{stim_idx, rec_idx});rms_data_stim_rec{stim_idx, rec_idx}=rms_data;
        rms_n1_data = cell2mat(rms_n1_data_stim_rec{stim_idx, rec_idx});rms_n1_data_stim_rec{stim_idx, rec_idx}=rms_n1_data;
        rms_n2_data = cell2mat(rms_n2_data_stim_rec{stim_idx, rec_idx});rms_n2_data_stim_rec{stim_idx, rec_idx}=rms_n2_data;
        patient_ids = patient_id_stim_rec{stim_idx, rec_idx};
        
        % Only fit model if we have data
        if ~isempty(rms_data)
            % Create ASM state vector (0=off, 1=on)
            n_data = length(rms_data);
            asm_state = [zeros(n_data/2, 1); ones(n_data/2, 1)];
            
            % Create tables for regression
            tbl_rms = table(categorical(asm_state), categorical(patient_ids), rms_data, ...
                           'VariableNames', {'ASM_State', 'PatientID', 'RMS'});
            tbl_rms_n1 = table(categorical(asm_state), categorical(patient_ids), rms_n1_data, ...
                              'VariableNames', {'ASM_State', 'PatientID', 'RMS_N1'});
            tbl_rms_n2 = table(categorical(asm_state), categorical(patient_ids), rms_n2_data, ...
                              'VariableNames', {'ASM_State', 'PatientID', 'RMS_N2'});
            
            % Fit models
            try
                mdl_rms = fitlm(tbl_rms, 'RMS ~ ASM_State + PatientID');
                mdl_rms_n1 = fitlm(tbl_rms_n1, 'RMS_N1 ~ ASM_State + PatientID');
                mdl_rms_n2 = fitlm(tbl_rms_n2, 'RMS_N2 ~ ASM_State + PatientID');
                
                % Store models
                models_rms{stim_idx, rec_idx} = mdl_rms;
                models_rms_n1{stim_idx, rec_idx} = mdl_rms_n1;
                models_rms_n2{stim_idx, rec_idx} = mdl_rms_n2;
                
                % Extract p-value for ASM_State effect (second row in coefficients)
                pval_rms_matrix(stim_idx, rec_idx) = mdl_rms.Coefficients.pValue(2);
                pval_rms_n1_matrix(stim_idx, rec_idx) = mdl_rms_n1.Coefficients.pValue(2);
                pval_rms_n2_matrix(stim_idx, rec_idx) = mdl_rms_n2.Coefficients.pValue(2);
                
            catch ME
                fprintf('Warning: Model fit failed for stim=%s, rec=%s\n', for_crisscross{stim_idx}, for_crisscross{rec_idx});
                fprintf('  Error: %s\n', ME.message);
            end
        end
    end
end

%% Create visualization (4x4 heatmaps with p-values)
figure('position',[0 0 1800 500])
h_rms = zeros(size(pval_rms_matrix));
h_rmsN1 = zeros(size(pval_rms_matrix));
h_rmsN2 = zeros(size(pval_rms_matrix));

% Plot 1: RMS (10-300ms)
subplot(1, 3, 1)
imagesc(pval_rms_matrix', [0.0001 0.05])
colormap(gca, master_ColorMaps('hawaii'))
cbar1 = colorbar;
cbar1.Label.String = 'p-value';
set(gca, 'XTick', 1:4, 'XTickLabel', for_crisscross,'FontSize', 12);
set(gca, 'YTick', 1:4, 'YTickLabel', for_crisscross);
xlabel('Stimulation Channel');
ylabel('Recording Channel');
title('RMS (10-300ms)');
grid off

% Add p-values and significance boxes
for stim_idx = 1:4
    for rec_idx = 1:4
        pval = pval_rms_matrix(stim_idx, rec_idx);
        if ~isnan(pval)
            if pval < 0.0001
                text(stim_idx, rec_idx, 'p<0.0001', 'Color', 'k', ...
                    'HorizontalAlignment', 'center', 'FontSize', 12);
            else
                text(stim_idx, rec_idx, num2str(pval, '%.4f'), 'Color', 'k', ...
                    'HorizontalAlignment', 'center', 'FontSize', 12);
            end

            % Draw colored box if Bonferroni-corrected significant
            if pval < alpha_fam
                % Determine direction of effect: mean(ON) - mean(OFF)
                n_each = length(rms_data_stim_rec{stim_idx, rec_idx}) / 2;
                if n_each > 0
                    rms_vals = rms_data_stim_rec{stim_idx, rec_idx};
                    mean_off = mean(rms_vals(1:n_each));
                    mean_on = mean(rms_vals(n_each+1:end));
                    if mean_on < mean_off
                        box_color = 'b'; % Blue: decrease from OFF to ON
                        h_rms(stim_idx, rec_idx) = -1;
                    else
                        box_color = 'r'; % Red: increase from OFF to ON
                        h_rms(stim_idx, rec_idx) = 1;
                    end
                else
                    box_color = 'k';
                end
                rectangle('Position', [stim_idx-0.5, rec_idx-0.5, 1, 1], ...
                    'EdgeColor', box_color, 'LineWidth', 2.5);
            end
        end
    end
end

% Plot 2: RMS-N1 (10-30ms)
subplot(1, 3, 2)
imagesc(pval_rms_n1_matrix', [0.0001 0.05])
colormap(gca, master_ColorMaps('hawaii'))
cbar2 = colorbar;
cbar2.Label.String = 'p-value';
set(gca, 'XTick', 1:4, 'XTickLabel', for_crisscross,'FontSize', 12);
set(gca, 'YTick', 1:4, 'YTickLabel', for_crisscross);
xlabel('Stimulation Channel');
ylabel('Recording Channel');
title('RMS-N1 (10-30ms)');
grid off

% Add p-values and significance boxes
for stim_idx = 1:4
    for rec_idx = 1:4
        pval = pval_rms_n1_matrix(stim_idx, rec_idx);
        if ~isnan(pval)
            if pval < 0.0001
                text(stim_idx, rec_idx, 'p<0.0001', 'Color', 'k', ...
                    'HorizontalAlignment', 'center', 'FontSize', 12);
            else
                text(stim_idx, rec_idx, num2str(pval, '%.4f'), 'Color', 'k', ...
                    'HorizontalAlignment', 'center', 'FontSize', 12);
            end

            % Draw colored box if Bonferroni-corrected significant
            if pval < alpha_fam
                % Determine direction of effect: mean(ON) - mean(OFF)
                n_each = length(rms_n1_data_stim_rec{stim_idx, rec_idx}) / 2;
                if n_each > 0
                    rms_vals = rms_n1_data_stim_rec{stim_idx, rec_idx};
                    mean_off = mean(rms_vals(1:n_each));
                    mean_on = mean(rms_vals(n_each+1:end));
                    if mean_on < mean_off
                        box_color = 'b'; % Blue: decrease from OFF to ON
                        h_rmsN1(stim_idx, rec_idx) = -1;
                    else
                        box_color = 'r'; % Red: increase from OFF to ON
                        h_rmsN1(stim_idx, rec_idx) = 1;
                    end
                else
                    box_color = 'k';
                end
                rectangle('Position', [stim_idx-0.5, rec_idx-0.5, 1, 1], ...
                    'EdgeColor', box_color, 'LineWidth', 2.5);
            end
        end
    end
end

% Plot 3: RMS-N2 (85-250ms)
subplot(1, 3, 3)
imagesc(pval_rms_n2_matrix', [0.0001 0.05])
colormap(gca, master_ColorMaps('hawaii'))
cbar3 = colorbar;
cbar3.Label.String = 'p-value';
set(gca, 'XTick', 1:4, 'XTickLabel', for_crisscross,'FontSize', 12);
set(gca, 'YTick', 1:4, 'YTickLabel', for_crisscross);
xlabel('Stimulation Channel');
ylabel('Recording Channel');
title('RMS-N2 (85-250ms)');
grid off

% Add p-values and significance boxes
for stim_idx = 1:4
    for rec_idx = 1:4
        pval = pval_rms_n2_matrix(stim_idx, rec_idx);
        if ~isnan(pval)
            if pval < 0.0001
                text(stim_idx, rec_idx, 'p<0.0001', 'Color', 'k', ...
                    'HorizontalAlignment', 'center', 'FontSize', 12);
            else
                text(stim_idx, rec_idx, num2str(pval, '%.4f'), 'Color', 'k', ...
                    'HorizontalAlignment', 'center', 'FontSize', 12);
            end


            % Draw colored box if Bonferroni-corrected significant
            if pval < alpha_fam
                % Determine direction of effect: mean(ON) - mean(OFF)
                n_each = length(rms_n2_data_stim_rec{stim_idx, rec_idx}) / 2;
                if n_each > 0
                    rms_vals = rms_n2_data_stim_rec{stim_idx, rec_idx};
                    mean_off = mean(rms_vals(1:n_each));
                    mean_on = mean(rms_vals(n_each+1:end));
                    if mean_on < mean_off
                        box_color = 'b'; % Blue: decrease from OFF to ON
                        h_rmsN2(stim_idx, rec_idx) = -1;
                    else
                        box_color = 'r'; % Red: increase from OFF to ON
                        h_rmsN2(stim_idx, rec_idx) = 1;
                    end
                else
                    box_color = 'k';
                end
                rectangle('Position', [stim_idx-0.5, rec_idx-0.5, 1, 1], ...
                    'EdgeColor', box_color, 'LineWidth', 2.5);
            end
        end
    end
end

sgtitle({'Linear Regression: ASM Effect on Connectivity (RMS ~ ASM_State + PatientID)'; ...
         sprintf('Columns = Stimulated Region | Rows = Recorded Region | Red box = p < %.6f (Bonferroni corrected)', alpha_fam)});

%% Save figure and results
res_dir = fullfile(Sbj_Metadatas{1}.project_root, 'COMBINED_RESULTS', 'offon_RMS_multi_v2');
if ~isfolder(res_dir)
    mkdir(res_dir);
end

print(fullfile(res_dir, 'RMS_regression_connectivity_pvalues.png'), '-dpng', '-r300');

%% Save results
save(fullfile(res_dir, 'RMS_regression_connectivity_results.mat'), ...
     'models_rms', 'models_rms_n1', 'models_rms_n2', ...
     'pval_rms_matrix', 'pval_rms_n1_matrix', 'pval_rms_n2_matrix', ...
     'rms_data_stim_rec', 'rms_n1_data_stim_rec', 'rms_n2_data_stim_rec', ...
     'patient_id_stim_rec', 'alpha_fam', 'for_crisscross',...
     'h_rms','h_rmsN1','h_rmsN2');

fprintf('\n\nResults saved to: %s\n', res_dir);

end
