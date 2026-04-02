function CCEP_offon_collect_drectns_perm(Sbj_Metadatas)
% CCEP_offon_collect_drectns_perm: Collect permutation p-values across all
% subjects and stim->recording pairs, apply FDR correction (Benjamini-Hochberg),
% and report % significant / direction of change per epileptogenic class pair.
%
% Direction of change: ON -> OFF
%   mean_diff > 0  => measure INCREASED from ON to OFF  (dir = 1)
%   mean_diff < 0  => measure DECREASED from ON to OFF  (dir = 2)
%
% FDR is applied separately per measure (N1_time, N1_amp, N2_time, N2_amp)
% across ALL pairs pooled from all subjects.
%
% All raw p-values, effect sizes (mean_diff), and trial-level data are
% saved in the output .mat for downstream use.

importance       = {'SOZ','EPZ','IZ','Healthy'};
measures         = {'N1_time','N1_amp','N2_time','N2_amp'};
FDR_ALPHA        = 0.05;

%% Initialise accumulators
comb_sbjIDs               = '';
chan_count                 = 0;
bp_good_chans_classes_all  = {};
bad_analyzed_chans_all     = {};
bp_chans_all               = {};
bp_good_chans_all          = {};
analyzed_chans_character_all = {};
analyzed_chans_all         = {};

% criss-cross containers: each cell holds one row per pair
% fields: stim_chan, rec_chan, stim_class, rec_class, sbj_ID,
%         raw p-values and mean_diffs for each measure,
%         trial-level data (off/on) for each measure
all_pairs = struct(...
    'stim_chan',    {{}}, ...
    'rec_chan',     {{}}, ...
    'stim_class',  {{}}, ...
    'rec_class',   {{}}, ...
    'sbj_ID',      {{}});
for m = 1:4
    all_pairs.(['p_raw_'    measures{m}]) = [];
    all_pairs.(['mean_diff_' measures{m}]) = [];  % OFF - ON
    all_pairs.(['off_trials_' measures{m}]) = {};
    all_pairs.(['on_trials_'  measures{m}]) = {};
end

% 4x4 criss-cross data structs (initialise empty)
for_crisscross = {'SOZ','EPZ','IZ','Healthy'};
cc = struct();
for x1 = 1:4
    for x2 = 1:4
        tag = [for_crisscross{x1} '2' for_crisscross{x2}];
        for m = 1:4
            cc.(tag).(measures{m}).p_raw    = [];
            cc.(tag).(measures{m}).mean_diff = [];
            cc.(tag).(measures{m}).off_trials = {};
            cc.(tag).(measures{m}).on_trials  = {};
        end
    end
end

%% Loop subjects
for s = 1:length(Sbj_Metadatas)

    Sbj_Metadata = Sbj_Metadatas{s};
    sbj_ID       = Sbj_Metadata.sbj_ID;
    comb_sbjIDs  = [comb_sbjIDs '_' sbj_ID];

    AllBlockInfo = readtable(fullfile(Sbj_Metadata.project_root, ...
        [Sbj_Metadata.project_name '_BlockInfo.xlsx']));

    offblocks    = find(strcmp(AllBlockInfo.sbj_ID, sbj_ID) & strcmp(AllBlockInfo.task_type,'offmed'));
    onblocks     = find(strcmp(AllBlockInfo.sbj_ID, sbj_ID) & strcmp(AllBlockInfo.task_type,'onmed'));
    offmed_block = AllBlockInfo.BlockList{offblocks(1)};
    onmed_block  = AllBlockInfo.BlockList{onblocks(1)};

    %% Analyzed stim channels
    fileList       = dir(fullfile(Sbj_Metadata.results,'off_on_peaks_perm','peak_info_*.mat'));
    analyzed_chans = cellfun(@(x) erase(x,{'peak_info_','.mat'}), {fileList.name}','UniformOutput',false);

    %% Good channels
    ref = 'bp';
    off_ccep_tmp = load(fullfile(Sbj_Metadata.iEEG_data, offmed_block, 'CCEPdir', ...
        [offmed_block '_ccep_' ref '_' analyzed_chans{1} '.mat']));
    bp_chans     = off_ccep_tmp.CCEPs.label; clear off_ccep_tmp
    bp_chans_all = [bp_chans_all; bp_chans];

    off_info = load(fullfile(Sbj_Metadata.iEEG_data, offmed_block, [offmed_block '_info.mat']));
    on_info  = load(fullfile(Sbj_Metadata.iEEG_data, onmed_block,  [onmed_block  '_info.mat']));

    bp_good_chans_off = get_info_goodchans_bp(off_info.info, bp_chans);
    bp_good_chans_on  = get_info_goodchans_bp(on_info.info,  bp_chans);
    bp_good_chans     = bp_good_chans_off(ismember(bp_good_chans_off, bp_good_chans_on));
    bp_good_chans_all = [bp_good_chans_all; bp_good_chans];
    chan_count        = chan_count + height(off_info.info.channelinfo);

    %% Channel classes
    bp_good_chans_split   = strsplit_SA(bp_good_chans);
    bp_good_chans_classes = assign_classes(bp_good_chans_split(:,1), bp_good_chans_split(:,2), ...
        off_info.info.channelinfo.Label, off_info.info.channelinfo.chan_info, importance);
    bp_good_chans_classes_all = [bp_good_chans_classes_all; bp_good_chans_classes];

    %% Stim channel classes
    events_chans          = strcat(off_info.info.events.StimCh1, repmat({'-'},[length(off_info.info.events.StimCh1),1]), off_info.info.events.StimCh2);
    analyzed_chans_split  = strsplit_SA(analyzed_chans);
    analyzed_chans_character = cell(length(analyzed_chans),1);
    for c_idx = 1:length(analyzed_chans_character)
        xxx = find(ismember(events_chans, analyzed_chans{c_idx}));
        analyzed_chans_character(c_idx) = off_info.info.events.elec_type(xxx(1));
    end

    %% Remove bad stim channels
    bad_idx            = ~ismember(analyzed_chans, bp_good_chans);
    bad_analyzed_chans = analyzed_chans(bad_idx);
    if ~isempty(bad_analyzed_chans)
        warning('%s: bad stim channels removed: %s', sbj_ID, strjoin(bad_analyzed_chans,','))
    end
    analyzed_chans           = analyzed_chans(~bad_idx);
    analyzed_chans_character = analyzed_chans_character(~bad_idx);
    bad_analyzed_chans_all   = [bad_analyzed_chans_all; bad_analyzed_chans];
    analyzed_chans_all       = [analyzed_chans_all;     analyzed_chans];
    analyzed_chans_character_all = [analyzed_chans_character_all; analyzed_chans_character];

    corrSheet = readtable(Sbj_Metadata.labelfile);

    %% Loop stim channels
    for c = 1:length(analyzed_chans)

        load(fullfile(Sbj_Metadata.results,'off_on_peaks_perm',['peak_info_' analyzed_chans{c} '.mat']),'peak_info')

        % Remove bad and stim channels
        peak_info      = peak_info(ismember(peak_info.label, bp_good_chans),:);
        peak_info_spl  = strsplit_SA(peak_info.label);
        peak_info      = peak_info(~any(ismember(peak_info_spl, analyzed_chans_split(c,:)),2),:);
        peak_info_spl  = strsplit_SA(peak_info.label);

        % Remove non-gray matter channels
        peakinfo_chan_WMvsGM = assign_classes(peak_info_spl(:,1), peak_info_spl(:,2), ...
            corrSheet.Label, lower(corrSheet.WMvsGM), {'gray','white','csf','skull'});
        peak_info     = peak_info(strcmp(peakinfo_chan_WMvsGM,'gray'),:);
        peak_info_spl = strsplit_SA(peak_info.label);

        % Recording channel classes
        peakinfo_chan_classes = assign_classes(peak_info_spl(:,1), peak_info_spl(:,2), ...
            off_info.info.channelinfo.Label, off_info.info.channelinfo.chan_info, importance);

        stim_class = analyzed_chans_character{c};

        %% Loop recording channels
        for p = 1:height(peak_info)

            rec_class  = peakinfo_chan_classes{p};
            perm_res   = peak_info.permresults{p};
            tag        = [stim_class '2' rec_class];

            % Accumulate into criss-cross struct and flat all_pairs
            n_ap = length(all_pairs.stim_chan) + 1;
            all_pairs.stim_chan{n_ap}  = analyzed_chans{c};
            all_pairs.rec_chan{n_ap}   = peak_info.label{p};
            all_pairs.stim_class{n_ap} = stim_class;
            all_pairs.rec_class{n_ap}  = rec_class;
            all_pairs.sbj_ID{n_ap}     = sbj_ID;

            for m = 1:4
                mn  = measures{m};
                p_r = perm_res.(mn).p;
                md  = perm_res.(mn).mean_diff;  % OFF - ON

                all_pairs.(['p_raw_'    mn])(n_ap) = p_r;
                all_pairs.(['mean_diff_' mn])(n_ap) = md;
                all_pairs.(['off_trials_' mn]){n_ap} = peak_info.(['off_' mn '_all']){p};
                all_pairs.(['on_trials_'  mn]){n_ap} = peak_info.(['on_'  mn '_all']){p};

                cc.(tag).(mn).p_raw(end+1)    = p_r;
                cc.(tag).(mn).mean_diff(end+1) = md;
                cc.(tag).(mn).off_trials{end+1} = peak_info.(['off_' mn '_all']){p};
                cc.(tag).(mn).on_trials{end+1}  = peak_info.(['on_'  mn '_all']){p};
            end
        end
    end

% end of subjects loop
end

%% FDR correction across ALL pairs, per measure
fprintf('\n=== FDR correction (Benjamini-Hochberg, q=%.2f) ===\n', FDR_ALPHA)
for m = 1:4
    mn      = measures{m};
    p_raw   = all_pairs.(['p_raw_' mn]);
    h_fdr   = fdr_bh(p_raw, FDR_ALPHA);  % logical vector, 1 = significant
    all_pairs.(['h_fdr_' mn]) = h_fdr;
    fprintf('  %s: %d / %d pairs significant after FDR\n', mn, sum(h_fdr), length(h_fdr))
end

%% Push FDR results back into criss-cross struct + compute summaries
% We need to re-match the flat all_pairs index to the cc struct entries
% We track per-(stim_class, rec_class, measure) which flat indices belong

summary = struct();  % summary.(tag).(mn): perc_sig, perc_inc, perc_dec, n_total

for x1 = 1:4
    for x2 = 1:4
        stim_cls = for_crisscross{x1};
        rec_cls  = for_crisscross{x2};
        tag      = [stim_cls '2' rec_cls];

        % Find flat indices for this class pair
        idx_tag = find(strcmp(all_pairs.stim_class, stim_cls) & strcmp(all_pairs.rec_class, rec_cls));

        for m = 1:4
            mn = measures{m};

            if isempty(idx_tag)
                cc.(tag).(mn).h_fdr      = [];
                summary.(tag).(mn).n_total  = 0;
                summary.(tag).(mn).perc_sig = NaN;
                summary.(tag).(mn).perc_inc = NaN;
                summary.(tag).(mn).perc_dec = NaN;
                continue
            end

            h_fdr_sub  = all_pairs.(['h_fdr_' mn])(idx_tag);
            mean_diffs  = all_pairs.(['mean_diff_' mn])(idx_tag);  % mean_diff = [OFF] - [ON] (line 165 in CCEP_comp_offon_perm.m)

            cc.(tag).(mn).h_fdr     = h_fdr_sub;
            cc.(tag).(mn).p_raw     = all_pairs.(['p_raw_' mn])(idx_tag);
            cc.(tag).(mn).mean_diff = mean_diffs;

            n_total = length(h_fdr_sub);
            n_sig   = sum(h_fdr_sub == 1 & ~isnan(h_fdr_sub));

            % Direction: among FDR-significant pairs
            sig_mask = logical(h_fdr_sub);
            n_inc    = sum(mean_diffs(sig_mask) > 0);   % OFF > ON  = increased going to OFF
            n_dec    = sum(mean_diffs(sig_mask) < 0);   % OFF < ON  = decreased going to OFF

            summary.(tag).(mn).n_total   = n_total;
            summary.(tag).(mn).n_sig     = n_sig;
            summary.(tag).(mn).perc_sig  = n_sig / n_total;
            summary.(tag).(mn).perc_inc  = n_inc / n_total;   % % of ALL pairs that sig. increased
            summary.(tag).(mn).perc_dec  = n_dec / n_total;   % % of ALL pairs that sig. decreased
            summary.(tag).(mn).n_inc     = n_inc;
            summary.(tag).(mn).n_dec     = n_dec;
        end
    end
end

%% Build 4x4 matrices for plotting
measure_labels = {'N1 Latency','N1 Amplitude','N2 Latency','N2 Amplitude'};
measure_units  = {'s','uV','s','uV'};

% Save directory
res_dir = fullfile(Sbj_Metadatas{1}.project_root,'COMBINED_RESULTS','offon_collect_perm');
if ~isfolder(res_dir); mkdir(res_dir); end

for m = 1:4
    mn = measures{m};
    perc_sig = nan(4,4); perc_inc = nan(4,4); perc_dec = nan(4,4);
    for x1 = 1:4
        for x2 = 1:4
            tag = [for_crisscross{x1} '2' for_crisscross{x2}];
            perc_sig(x1,x2) = summary.(tag).(mn).perc_sig;
            perc_inc(x1,x2) = summary.(tag).(mn).perc_inc;
            perc_dec(x1,x2) = summary.(tag).(mn).perc_dec;
        end
    end

    figure('Position',[0 0 1800 400])

    subplot(1,3,1)
    plot_ccep_heatmap(100*perc_sig', for_crisscross, ...
        ['% Significant (FDR q<' num2str(FDR_ALPHA) ') — ' measure_labels{m}])

    subplot(1,3,2)
    clim_val = [0, 100*max([perc_inc(:); perc_dec(:)],[],'all','omitnan')];
    plot_ccep_heatmap(100*perc_inc', for_crisscross, ...
        ['% Sig. Increased (OFF>ON) — ' measure_labels{m}], clim_val)

    subplot(1,3,3)
    plot_ccep_heatmap(100*perc_dec', for_crisscross, ...
        ['% Sig. Decreased (OFF<ON) — ' measure_labels{m}], clim_val)

    sgtitle(['ASM ON \rightarrow OFF: ' measure_labels{m} ' (perm. test, FDR corrected)'])
    print(fullfile(res_dir, [comb_sbjIDs '_offon_perm_' mn '.png']), '-dpng', '-r300')
    close
end

%% Save everything
save(fullfile(res_dir, [comb_sbjIDs '_offon_perm_results.mat']), ...
    'all_pairs', 'cc', 'summary', 'FDR_ALPHA', ...
    'bp_good_chans_classes_all','bp_good_chans_all', ...
    'analyzed_chans_character_all','analyzed_chans_all', '-v7.3')

%% Print summary text
print_channel_summary(res_dir, comb_sbjIDs, chan_count, ...
    bp_chans_all, bp_good_chans_all, analyzed_chans_all, ...
    bad_analyzed_chans_all, analyzed_chans_character_all, ...
    bp_good_chans_classes_all, measures, summary, for_crisscross, FDR_ALPHA)

end

%% =========================================================================
%  LOCAL HELPERS
%% =========================================================================

function plot_ccep_heatmap(data_matrix, class_labels, ttl, clim_val)
% data_matrix: [4 x 4], rows = rec, cols = stim (already transposed by caller)
    image(data_matrix, 'CDataMapping','scaled')
    colormap(master_ColorMaps('hawaii'));
    colorbar
    if nargin == 4 && ~isempty(clim_val)
        clim(clim_val)
    end
    xticks(1:4); xticklabels(class_labels); xlabel('Stimulation Channels')
    yticks(1:4); yticklabels(class_labels); ylabel('Recording Channels')
    for x1 = 1:4
        for x2 = 1:4
            text(x2, x1, [num2str(data_matrix(x1,x2),'%.1f') '%'], ...
                'Color','k','HorizontalAlignment','center','FontSize',9)
        end
    end
    title(ttl)
end

function h = fdr_bh(p_vals, q)
% Benjamini-Hochberg FDR correction.
% p_vals: vector of p-values (NaN treated as non-significant)
% q: FDR threshold (e.g. 0.05)
% Returns logical vector h (1 = significant)
    p_vals  = p_vals(:);
    n       = length(p_vals);
    h       = false(n,1);
    [p_sorted, sort_idx] = sort(p_vals);
    thresholds = (1:n)' * q / n;
    last_sig   = find(p_sorted <= thresholds, 1, 'last');
    if ~isempty(last_sig)
        h(sort_idx(1:last_sig)) = true;
    end
    % NaN p-values remain false
    h(isnan(p_vals)) = false;
end

function print_channel_summary(res_dir, comb_sbjIDs, chan_count, ...
    bp_chans_all, bp_good_chans_all, analyzed_chans_all, ...
    bad_analyzed_chans_all, analyzed_chans_character_all, ...
    bp_good_chans_classes_all, measures, summary, for_crisscross, FDR_ALPHA)

    prt  = sprintf('Number of total recording channels: %d\n', chan_count);
    prt  = [prt sprintf('Number of total recording pairs: %d\n',              length(bp_chans_all))];
    prt  = [prt sprintf('Number of non-artefactual recording pairs: %d\n',    length(bp_good_chans_all))];
    prt  = [prt sprintf('Number of stimulation pairs in analysis: %d\n',      length(analyzed_chans_all))];
    prt  = [prt sprintf('  Removed (artefactual): %d\n',                      length(bad_analyzed_chans_all))];
    for cls = {'SOZ','EPZ','IZ','Healthy'}
        prt = [prt sprintf('  Rec pairs that were %s: %d\n', cls{1}, sum(strcmp(bp_good_chans_classes_all,cls{1})))];
    end
    for cls = {'SOZ','EPZ','IZ','Healthy'}
        prt = [prt sprintf('  Stim pairs that were %s: %d\n', cls{1}, sum(strcmp(analyzed_chans_character_all,cls{1})))];
    end

    prt = [prt sprintf('\n=== Significance Summary (FDR q=%.2f, direction: ON->OFF) ===\n', FDR_ALPHA)];
    for m = 1:length(measures)
        mn  = measures{m};
        prt = [prt sprintf('\n%s:\n', mn)];
        prt = [prt sprintf('  %6s -> %8s : %%sig  %%inc  %%dec  (n)\n','Stim','Rec')];
        for x1 = 1:4
            for x2 = 1:4
                tag = [for_crisscross{x1} '2' for_crisscross{x2}];
                sm  = summary.(tag).(mn);
                if sm.n_total > 0
                    prt = [prt sprintf('  %6s -> %8s : %5.1f%% %5.1f%% %5.1f%%  (%d)\n', ...
                        for_crisscross{x1}, for_crisscross{x2}, ...
                        100*sm.perc_sig, 100*sm.perc_inc, 100*sm.perc_dec, sm.n_total)];
                end
            end
        end
    end

    fileID = fopen(fullfile(res_dir,[comb_sbjIDs '_offon_perm_summary.txt']),'w');
    fprintf(fileID,'%s',prt);
    fclose(fileID);
    fprintf('%s', prt)
end
