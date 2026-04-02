% CCEP_fig1_v2 (CCEP parameters inc vs dec)
% v4: using permutation testing results

sbj_IDs = {'P1', 'P2', 'P3', 'P4', 'P7'};
comb_sbjIDs = '';
for s = 1:length(sbj_IDs)
    sbj_ID = sbj_IDs{s};
    comb_sbjIDs = [comb_sbjIDs '_' sbj_ID];
end

savedir = fullfile(cd,'Figures');if ~isfolder(savedir); mkdir(savedir); end
Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_IDs{1}); % example Sbj_Metadata
res_dir=fullfile(Sbj_Metadata.project_root,'COMBINED_RESULTS','offon_collect_perm');
load(fullfile(res_dir,[comb_sbjIDs '_offon_allDirections.mat']),'collected')

% load everything
load(fullfile(res_dir, [comb_sbjIDs '_offon_perm_results.mat']),'summary')

% arrange save directory
savedir = fullfile(cd,'Figures');if ~isfolder(savedir); mkdir(savedir); end
if ~isfolder(savedir), mkdir(savedir); end

measures = {'N1_time','N1_amp','N2_time','N2_amp'};
for_crisscross = {'SOZ','EPZ','IZ','Healthy'};

for m = 1:4
    mn = measures{m};
    % perc_sig = nan(4,4); 
    perc_inc = nan(4,4); perc_dec = nan(4,4);
    for x1 = 1:4
        for x2 = 1:4
            tag = [for_crisscross{x1} '2' for_crisscross{x2}];
            % perc_sig(x1,x2) = summary.(tag).(mn).perc_sig;
            perc_inc(x1,x2) = summary.(tag).(mn).perc_inc;
            perc_dec(x1,x2) = summary.(tag).(mn).perc_dec;
        end
    end

    figure('position',[0 0 800 300])

    subplot(1,2,1)
    clim_val = [0, 12];
    plot_ccep_heatmap(100*perc_dec', for_crisscross, ...
        '', clim_val)

    subplot(1,2,2)
    plot_ccep_heatmap(100*perc_inc', for_crisscross, ...
        '', clim_val)

    print(fullfile(savedir,[char(datetime('today','Format','uuuu-MM-dd')) '_' mn '_viridis.png']),'-dpng','-r300')

end


%% colormap - horizontal
figure('Units','normalized','Position',[0 0 .5 1])
% ax = axes;
colormap(master_ColorMaps('-viridis'));
clim([0 12]);
c=colorbar(gca);
% c.Label.String = 'Percentages';c.Label.Rotation=270;
c.FontSize=30;
c.FontWeight='bold';
c.LineWidth = 1;
c.Limits=[0 12];
c.Ticks = [0,3,6,9,12];
c.TickLabels = {'0%','3%','6%','9%','12%'};
c.TickDirection='both';
c.Location='southoutside';

a=get(c); %gets properties of colorbar
a = a.Position; %gets the positon and size of the color bar
set(c,'Units','normalized','Position',[.1 .1 .5 .03])% To change size
% set(c,'Units','normalized','Position',[.5 .1 .05 .5])% To change size
set(gca,'Visible','off')
% ax.Visible = 'off';
set(gcf,'color','w')
set(gcf, 'InvertHardcopy', 'off') % to set saving white color as white

print(fullfile(savedir,[char(datetime('today','Format','uuuu-MM-dd')) '_colorbar_horz.png']),'-dpng','-r300')

%% colormap - vertical
% figure('Units','normalized','Position',[0 0 .5 1])
% % ax = axes;
% colormap(master_ColorMaps('-viridis'));
% clim([0 50]);
% c=colorbar(gca);
% % c.Label.String = 'Percentages';c.Label.Rotation=270;
% c.FontSize=30;
% c.FontWeight='bold';
% c.LineWidth = 1;
% c.Limits=[0 50];
% c.Ticks = [0,10,20,30,40,50];
% c.TickDirection='both';
% % c.Location='southoutside';
%
% a=get(c); %gets properties of colorbar
% a = a.Position; %gets the positon and size of the color bar
% % set(c,'Units','normalized','Position',[.1 .1 .5 .03])% To change size
% set(c,'Units','normalized','Position',[.5 .1 .05 .5])% To change size
% set(gca,'Visible','off')
% % ax.Visible = 'off';
% set(gcf,'color','w')
% set(gcf, 'InvertHardcopy', 'off') % to set saving white color as white
%
% print(fullfile(savedir,[char(datetime('today','Format','uuuu-MM-dd')) '_colorbar_vert.png']),'-dpng','-r300')
%

%% =========================================================================
%  LOCAL HELPERS
%% =========================================================================

function plot_ccep_heatmap(data_matrix, class_labels, ttl, clim_val)
% data_matrix: [4 x 4], rows = rec, cols = stim (already transposed by caller)
    image(data_matrix, 'CDataMapping','scaled')
    colormap(master_ColorMaps('-viridis'));
    % colorbar
    if nargin == 4 && ~isempty(clim_val)
        clim(clim_val)
    end
    set(gca, 'xtick', [])
    set(gca, 'YTick', [])
    % xticks(1:4); xticklabels(class_labels);xlabel('Stimulation Channels')
    % yticks(1:4); yticklabels(class_labels);ylabel('Recording Channels')
    for x1 = 1:4
        for x2 = 1:4
            if data_matrix(x1,x2) < 5
                text(x2, x1, [num2str(data_matrix(x1,x2),'%.1f') '%'], ...
                    'Color','k','HorizontalAlignment','center','FontSize',12,'FontWeight','bold')
            else
                text(x2, x1, [num2str(data_matrix(x1,x2),'%.1f') '%'], ...
                    'Color','w','HorizontalAlignment','center','FontSize',12,'FontWeight','bold')
            end
        end
    end
    % title(ttl)
end