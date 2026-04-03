% Figure 4: to plot the connectivity fingerprints on vs off ASM during rest
% sessions
data_root = fullfile(cd,'data');
project_name = 'CCEP_PrePost';
sbj_IDs = {'P2', 'P4', 'P5', 'P6', 'P7'};
savedir = fullfile(cd,'Figures');
if ~isfolder(savedir), mkdir(savedir); end

for s = 1:length(sbj_IDs)
    sbj_ID = sbj_IDs{s};

    % load patient info
    Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
    AllBlockInfo = readtable(fullfile(Sbj_Metadata.project_root,[Sbj_Metadata.project_name '_BlockInfo.xlsx'])); % "F:\HBML\PROJECTS_DATA\CL_Train\CL_Train_BlockInfo.xlsx"
    preoffmed = AllBlockInfo.BlockList{strcmp(AllBlockInfo.sbj_ID,sbj_ID) & strcmp(AllBlockInfo.task_type,'preoffmedrest')};
    preonmed = AllBlockInfo.BlockList{strcmp(AllBlockInfo.sbj_ID,sbj_ID) & strcmp(AllBlockInfo.task_type,'preonmedrest')};
    on_info = load(fullfile(Sbj_Metadata.iEEG_data, preonmed, [preonmed '_info.mat']));
    off_info = load(fullfile(Sbj_Metadata.iEEG_data, preoffmed, [preoffmed '_info.mat']));
    [~,off_nonartefact_idx] = get_info_nonartefactchans(off_info.info);
    [~,on_nonartefact_idx] = get_info_nonartefactchans(on_info.info);
    good_ch_idx = off_nonartefact_idx & on_nonartefact_idx;
    chans = on_info.info.channelinfo.Label(good_ch_idx);
    channel_labels = on_info.info.channelinfo.chan_info(good_ch_idx);channel_labels = strrep(channel_labels,'Healthy','NIZ');
    onmedcorr = load(fullfile(Sbj_Metadata.iEEG_data, preonmed,'HFB_files','corr_matrix.mat'),'corrmatrix','corrmatrix_pval');
    offmedcorr = load(fullfile(Sbj_Metadata.iEEG_data, preoffmed,'HFB_files','corr_matrix.mat'),'corrmatrix','corrmatrix_pval');
    onmed_corr = onmedcorr.corrmatrix; onmed_corr = onmed_corr(good_ch_idx,good_ch_idx);
    offmed_corr = offmedcorr.corrmatrix; offmed_corr = offmed_corr(good_ch_idx,good_ch_idx);
    onmed_corr(logical(eye(size(onmed_corr)))) = 0;
    offmed_corr(logical(eye(size(offmed_corr)))) = 0;

    on_med_fc = onmed_corr;
    off_med_fc = offmed_corr;
    region_labels = {'SOZ', 'EPZ', 'IZ', 'NIZ'};
    n_regions = length(region_labels);

    % Create figure
    fig = figure('Position', [100, 100, 400, 400]);

    % Initialize data for boxplots
    region_pair_labels = {};
    on_med_data = {};
    off_med_data = {};
    p_values = [];

    % Calculate data for each region pair
    for i = 1:n_regions
        for j = 1:n_regions
            % Get channels for each region
            region_i_channels = find(ismember(channel_labels, region_labels{i}));
            region_j_channels = find(ismember(channel_labels, region_labels{j}));

            if ~isempty(region_i_channels) && ~isempty(region_j_channels)
                % Extract connectivity values
                if i == j
                    % Within-region connectivity (exclude self-connections)
                    [ri, ci] = meshgrid(region_i_channels, region_i_channels);
                    mask = ri ~= ci;
                    if sum(mask(:)) > 0
                        on_values = on_med_fc(region_i_channels, region_i_channels);
                        off_values = off_med_fc(region_i_channels, region_i_channels);

                        on_values = on_values(mask);
                        off_values = off_values(mask);
                    else
                        on_values = [];
                        off_values = [];
                    end
                else
                    % Between-region connectivity
                    on_values = on_med_fc(region_i_channels, region_j_channels);
                    off_values = off_med_fc(region_i_channels, region_j_channels);
                    on_values = on_values(:);
                    off_values = off_values(:);
                end

                if ~isempty(on_values)
                    % Store data for boxplot
                    region_pair_labels{end+1} = sprintf('%s-%s', region_labels{i}, region_labels{j});
                    on_med_data{end+1} = on_values;
                    off_med_data{end+1} = off_values;

                    % Perform statistical test (signrank for significance)
                    p_values(end+1) = signrank(on_values, off_values);
                end
            end
        end
    end

    % Create summary panel
    create_summary_visual(on_med_fc, off_med_fc, channel_labels, region_labels, p_values, region_pair_labels);

    % Make figure more compact
    set(fig, 'Color', 'white');
    % sgtitle(sprintf('Connectivity Analysis - Patient %s', patient_id), 'FontSize', 16);
    print(fullfile(savedir,[Sbj_Metadata.sbj_ID, '_fingerprint_boxplot.png']),'-dpng','-r300')
end

%% Functions
function create_summary_visual(on_med_fc, off_med_fc, channel_labels, region_labels, p_values, region_pair_labels)
% Create a circular graph visualization showing significant connections

n_regions = length(region_labels);
increased_col = [1 0 0];
decreased_col = [0.1529    0.3059    0.9098];

% Calculate region centers
theta = linspace(0, 2*pi, n_regions+1);
theta = theta(1:end-1);
radius = 0.8;
x_pos = radius * cos(theta);
y_pos = radius * sin(theta);

% Get significant changes
sig_increases = {};
sig_decreases = {};
hold on;

for i = 1:n_regions
    for j = 1:n_regions
        pair_label = sprintf('%s-%s', region_labels{i}, region_labels{j});
        idx = find(strcmp(region_pair_labels, pair_label));

        if ~isempty(idx) && p_values(idx) < 0.05
            source_channels = find(ismember(channel_labels, region_labels{i}));
            target_channels = find(ismember(channel_labels, region_labels{j}));

            % Calculate mean connectivity
            if i == j && length(source_channels) > 1
                % Within-region connectivity
                [ri, ci] = meshgrid(source_channels, source_channels);
                mask = ri ~= ci;

                if sum(mask(:)) > 0
                    on_values = on_med_fc(source_channels, source_channels);
                    off_values = off_med_fc(source_channels, source_channels);

                    on_mean = mean(on_values(mask));
                    off_mean = mean(off_values(mask));
                else
                    continue;
                end
            else
                on_mean = mean(mean(on_med_fc(source_channels, target_channels)));
                off_mean = mean(mean(off_med_fc(source_channels, target_channels)));
            end

            % Check direction of change
            if on_mean > off_mean
                sig_decreases{end+1} = [i, j, p_values(idx)];
            else
                sig_increases{end+1} = [i, j, p_values(idx)];
            end
        end
    end
end

% Draw connections
for i = 1:length(sig_increases)
    src = sig_increases{i}(1);
    tgt = sig_increases{i}(2);
    p_val = sig_increases{i}(3);

    % Calculate line thickness based on significance
    if p_val < 0.001
        lw = 3;
    elseif p_val < 0.01
        lw = 2;
    else
        lw = 1;
    end

    % Draw arrow
    arrow_x = [x_pos(src), x_pos(tgt)];
    arrow_y = [y_pos(src), y_pos(tgt)];

    % Add curve to the line
    if src ~= tgt
        % For different nodes, draw a curve
        mid_x = mean(arrow_x) + 0.2*(-y_pos(tgt) + y_pos(src));
        mid_y = mean(arrow_y) + 0.2*(x_pos(tgt) - x_pos(src));

        curve_x = [arrow_x(1), mid_x, arrow_x(2)];
        curve_y = [arrow_y(1), mid_y, arrow_y(2)];

        plot(curve_x, curve_y, 'LineWidth', lw, 'Color', increased_col);

        % Add arrowhead
        arrow_angle = atan2(curve_y(3)-curve_y(2), curve_x(3)-curve_x(2));
        arrow_length = 0.1;

        % Arrow endpoint adjustments to hit the circle edge
        end_x = x_pos(tgt) - 0.22*cos(arrow_angle);
        end_y = y_pos(tgt) - 0.22*sin(arrow_angle);

        % Draw arrow head
        ah_x = end_x - arrow_length * cos(arrow_angle - pi/6);
        ah_y = end_y - arrow_length * sin(arrow_angle - pi/6);
        plot([end_x, ah_x], [end_y, ah_y], 'LineWidth', lw, 'Color', increased_col);

        ah_x = end_x - arrow_length * cos(arrow_angle + pi/6);
        ah_y = end_y - arrow_length * sin(arrow_angle + pi/6);
        plot([end_x, ah_x], [end_y, ah_y], 'LineWidth', lw, 'Color', increased_col);
    else
        % For self-connections, draw a loop
        loop_radius = 0.2;
        loop_center_x = x_pos(src) + 0.22*cos(theta(src)+pi/4);
        loop_center_y = y_pos(src) + 0.22*sin(theta(src)+pi/4);

        loop_angles = linspace(0, 2*pi, 50);
        loop_x = loop_center_x + loop_radius * cos(loop_angles);
        loop_y = loop_center_y + loop_radius * sin(loop_angles);

        plot(loop_x, loop_y, 'LineWidth', 3, 'Color', increased_col);
    end
end

% Draw decreases (same as above but with different color)
for i = 1:length(sig_decreases)
    src = sig_decreases{i}(1);
    tgt = sig_decreases{i}(2);
    p_val = sig_decreases{i}(3);

    % Calculate line thickness based on significance
    if p_val < 0.001
        lw = 3;
    elseif p_val < 0.01
        lw = 2;
    else
        lw = 1;
    end

    % Draw arrow
    arrow_x = [x_pos(src), x_pos(tgt)];
    arrow_y = [y_pos(src), y_pos(tgt)];

    % Add curve to the line
    if src ~= tgt
        % For different nodes, draw a curve
        mid_x = mean(arrow_x) + 0.2*(-y_pos(tgt) + y_pos(src));
        mid_y = mean(arrow_y) + 0.2*(x_pos(tgt) - x_pos(src));

        curve_x = [arrow_x(1), mid_x, arrow_x(2)];
        curve_y = [arrow_y(1), mid_y, arrow_y(2)];

        plot(curve_x, curve_y, 'LineWidth', lw, 'Color', decreased_col);

        % Add arrowhead
        arrow_angle = atan2(curve_y(3)-curve_y(2), curve_x(3)-curve_x(2));
        arrow_length = 0.1;

        % Arrow endpoint adjustments to hit the circle edge
        end_x = x_pos(tgt) - 0.22*cos(arrow_angle);
        end_y = y_pos(tgt) - 0.22*sin(arrow_angle);

        % Draw arrow head
        ah_x = end_x - arrow_length * cos(arrow_angle - pi/6);
        ah_y = end_y - arrow_length * sin(arrow_angle - pi/6);
        plot([end_x, ah_x], [end_y, ah_y], 'LineWidth', lw, 'Color', decreased_col);

        ah_x = end_x - arrow_length * cos(arrow_angle + pi/6);
        ah_y = end_y - arrow_length * sin(arrow_angle + pi/6);
        plot([end_x, ah_x], [end_y, ah_y], 'LineWidth', lw, 'Color', decreased_col);
    else
        % For self-connections, draw a loop
        loop_radius = 0.2;
        loop_center_x = x_pos(src) + 0.22*cos(theta(src)-pi/4);
        loop_center_y = y_pos(src) + 0.22*sin(theta(src)-pi/4);

        loop_angles = linspace(0, 2*pi, 50);
        loop_x = loop_center_x + loop_radius * cos(loop_angles);
        loop_y = loop_center_y + loop_radius * sin(loop_angles);

        plot(loop_x, loop_y, 'LineWidth', 3, 'Color', decreased_col);
    end
end

% Create node positions
plot(x_pos, y_pos, 'o', 'MarkerSize', 45, 'LineWidth', 2.5, 'MarkerFaceColor', 'white','MarkerEdgeColor','k');

% Add region labels
for i = 1:n_regions
    text(x_pos(i), y_pos(i), region_labels{i}, 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
end

% Add legend
line_x = [-.9, -.7];
line_y = [-1.3, -1.3];
% plot(line_x, line_y, 'LineWidth', 3, 'Color', decreased_col);
% text(-.65, -1.3, 'Decreased OFF medication', 'FontSize', 14);

line_y = [-1.5, -1.5];
% plot(line_x, line_y, 'LineWidth', 3, 'Color', increased_col);
% text(-.65, -1.5, 'Increased OFF medication', 'FontSize', 14);

% Format plot
axis equal;
ylim([-1.2 1.2])
xlim([-1.2 1.2])
axis off;
% title('Significant Connectivity Changes', 'FontSize', 18);
end
