% Figure 5: visualization of network measures for multiple participants
%
% This creates a 2x4 grid of subplots comparing ON vs OFF states for
% multiple network measures. Each subplot contains individual boxplots for each participant.

data_root = fullfile(cd,'neural_data');
project_name = 'CCEP_PrePost';

% visualize across participants
results_array={[],[]};
sbj_ID = 'P2';
Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
xx=load(fullfile(Sbj_Metadata.results,'corr_matrix',[Sbj_Metadata.sbj_ID, '_graphMetrics.mat']),'results');
results_array{1}=xx.results;
sbj_ID = 'P4';
Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_ID); % 'SAkkol_Stanford'
xx=load(fullfile(Sbj_Metadata.results,'corr_matrix',[Sbj_Metadata.sbj_ID, '_graphMetrics.mat']),'results');
results_array{2}=xx.results;

patient_ids = {'P2','P4'};

% Extract data from all participants
n_participants = length(results_array);

% Define the measures to plot (8 in total for a 2x4 grid)
measures = {'degree', 'strength', 'clustering', 'local_efficiency', ...
    'path_length', 'modularity', 'global_efficiency', 'small_worldness'};
measure_titles = {'Degree', 'Strength', 'Clustering Coefficient', 'Local Efficiency',...
    'Path Length', 'Modularity', 'Global Efficiency', 'Small-Worldness'};

% Initialize data arrays for individual patients and measures
patient_data = cell(n_participants, length(measures), 2); % [patients, measures, state(ON/OFF)]
p_values = zeros(n_participants, length(measures));

% Extract data for each measure and participant
for i = 1:n_participants
    results = results_array{i};
    
    for m = 1:length(measures)
        measure = measures{m};
        
        % Check if the measure exists in the results
        if ~isfield(results, measure)
            warning('Measure %s not found in results for participant %d', measure, i);
            continue;
        end
        
        % Store ON and OFF data for this patient and measure
        patient_data{i, m, 1} = results.(measure).on;  % ON state
        patient_data{i, m, 2} = results.(measure).off; % OFF state
        
        % Store p-value if available
        if isfield(results.(measure), 'p') && ~isnan(results.(measure).p)
            p_values(i, m) = results.(measure).p;
        else
            p_values(i, m) = NaN;
        end
    end
end

%% Create the figure
fig = figure('Units','normalized','Position', [0 0.05  1 .9]);
% sgtitle('Network Measures Comparison (ASM-ON vs ASM-OFF)', 'FontWeight', 'bold');

% Create color map: Orange for on-med; Purple for off-med
offmed_color = [0.5128 0.0191 0.6551]; % purple
onmed_color = [0.9666 0.5638 0.2655]; % orange

% Iterate through measures and create subplots
for m = 1:length(measures)
    subplot(2, 4, m);
    hold on;
    
    % Variables to track plot positions and max values for significance lines
    boxPositions = [];
    allValues = [];
    
    % Plot each patient's boxplot pair
    boxWidth = 0.3; % Width of each box
    spaceBetweenPatients = 3; % Space between patient pairs
    
    for i = 1:n_participants
        % Get ON and OFF data for this patient and measure
        on_data = patient_data{i, m, 1};
        off_data = patient_data{i, m, 2};
        
        % Skip if data is empty
        if isempty(on_data) || isempty(off_data)
            continue;
        end
        
        % For scalar measures, convert to vector for boxplot
        if isscalar(on_data)
            on_data = [on_data; on_data]; % Duplicate to create a minimal boxplot
            off_data = [off_data; off_data];
        end
        
        % Track all values to set axis limits later
        allValues = [allValues; on_data; off_data];
        
        % Calculate positions for this patient's boxplots
        basePos = (i-1) * spaceBetweenPatients;
        onPos = basePos + 1;
        offPos = basePos + 2;
        
        boxPositions = [boxPositions, onPos, offPos];
        
        % Create group labels for boxplot
        data_to_plot = [on_data; off_data];
        group_labels = [repmat({sprintf('%s-ON', patient_ids{i})}, length(on_data), 1); 
                        repmat({sprintf('%s-OFF', patient_ids{i})}, length(off_data), 1)];
        
        % Create boxplot for this patient with custom colors
        boxplot(data_to_plot, group_labels, 'Positions', [onPos, offPos], ...
                'Colors', [onmed_color; offmed_color], 'symbol','')%,'MedianStyle','target','Width', boxWidth);
        % boxplot(data_to_plot, group_labels, 'Positions', [onPos, offPos], 'Width', boxWidth, ...
        %         'Colors', [onmed_color; 0.8, 0.4, 0.6], 'PlotStyle', 'compact');

        % Add significance lines if p-value is available and significant
        if ~isnan(p_values(i, m)) && p_values(i, m) < 0.05
            % Calculate height for significance line
            y_range = max(data_to_plot) - min(data_to_plot);
            if isempty(y_range) || y_range == 0
                y_range = 0.1 * max(data_to_plot);
            end
            
            % Get patient-specific boxplot height
            boxTop = max(data_to_plot) + 0.05 * y_range;
            
            % Draw significance line
            plot([onPos, offPos], [boxTop, boxTop], 'k-', 'LineWidth', 2);
            
            % Add asterisks based on significance level
            if p_values(i, m) < 0.001
                text((onPos + offPos)/2, boxTop + 0.01 * y_range, '***', 'HorizontalAlignment', 'center','FontSize',44);
            elseif p_values(i, m) < 0.01
                text((onPos + offPos)/2, boxTop + 0.01 * y_range, '**', 'HorizontalAlignment', 'center','FontSize',44);
            elseif p_values(i, m) < 0.05
                text((onPos + offPos)/2, boxTop + 0.01 * y_range, '*', 'HorizontalAlignment', 'center','FontSize',44);
            end
        end
    end

    % redo boxplot to get lines
    xxx=findobj('Tag','Median');
    % h = findobj('Tag','Box');
    tttt={};
    for ii=1:length(xxx)
        % h(ii).Color = 'k';
        xxx(ii).LineWidth=1;
        xxx(ii).LineStyle='-';
        if ii>4,continue,end
        plot(mean(get(xxx(ii),'XData')),mean(get(xxx(ii),'YData')),'diamond','Color',get(xxx(ii),'Color'),'LineWidth',6)
    end

    % Add title and customize appearance
    title(measure_titles{m}, 'FontWeight', 'bold');
    % ylabel('Value');
    set(gca,'FontSize',30)
    
    % Set the y-axis limits to ensure significance markers are visible
    if ~isempty(allValues)
        y_range = max(allValues) - min(allValues);
        if y_range == 0
            y_range = 0.1 * max(allValues);
        end
        ylim([min(allValues) - 0.1*y_range, max(allValues) + 0.2*y_range]);
    end
    
    % Adjust x-axis limits and ticks
    xlim([0, (n_participants * spaceBetweenPatients) + 1]);
    
    % Create custom x-tick labels showing only patient IDs
    xtick_positions = zeros(1, n_participants);
    xtick_labels = cell(1, n_participants);
    
    for i = 1:n_participants
        xtick_positions(i) = (i-1) * spaceBetweenPatients + 1.5; % Middle of ON-OFF pair
        xtick_labels{i} = patient_ids{i};
    end
    
    xticks(xtick_positions);
    xticklabels(xtick_labels);
    % xtickangle(45);
    xlim([0.5 5.5])
    
    % Set a common height limit across all subplots for consistent appearance
    box on;
    grid on;
    
    % Add a lighter background to differentiate subplots
    set(gca, 'Color', [0.97, 0.97, 0.97]);
    hold off;
end

% Add overall legend for ON/OFF states
h = axes('Position', [0.5, 0, 0.1, 0], 'Visible', 'off');
plot(h, [0,1], [0,0], 'Color', onmed_color, 'LineWidth', 10);
hold on;
plot(h, [0,1], [1,1], 'Color', offmed_color, 'LineWidth', 10);
legend(h, {'ASM-ON', 'ASM-OFF'}, 'Location', 'northoutside', 'Orientation', 'vertical','FontSize', 20);
xlim([-11550 -11549])
set(gca,'Visible','off')

% Add significance legend
h2 = axes('Position', [0.78, 0.02, 0.1, 0.03], 'Visible', 'off');
text(h2, 0, 0, '* p<0.05, ** p<0.01, *** p<0.001', 'FontSize', 22);

% Adjust layout
set(fig, 'Color', 'white');
tight_spacing = 0.1;
tight_subplot_dimensions = [tight_spacing, tight_spacing, tight_spacing, tight_spacing];
set(fig, 'DefaultAxesLooseInset', tight_subplot_dimensions);

% save figure
savedir = fullfile(cd,'Figures');
print(fullfile(savedir,['Figure5_GraphMeasures_' char(datetime('today','Format','uuuu-MM-dd')) '.png']),'-dpng','-r300')
