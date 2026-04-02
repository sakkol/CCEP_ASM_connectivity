% Figure 2: comparing CCEP based RMS between on vs off ASM
sbj_IDs = {'P1', 'P2', 'P3', 'P4', 'P7'};
comb_sbjIDs = '';
for s = 1:length(sbj_IDs)
    sbj_ID = sbj_IDs{s};
    comb_sbjIDs = [comb_sbjIDs '_' sbj_ID];
end

savedir = fullfile(cd,'Figures');
if ~isfolder(savedir), mkdir(savedir); end
res_dir=fullfile(Sbj_Metadata.project_root,'COMBINED_RESULTS','offon_RMS_multi_v2');
load(fullfile(res_dir, 'RMS_regression_connectivity_results.mat'), ...
     'pval_rms_matrix', 'pval_rms_n1_matrix', 'pval_rms_n2_matrix', ...
     'h_rms','h_rmsN1','h_rmsN2');

%% Create figure
figure('position',[100 100 2100 600])

% loop
for sp = 1:3
    subplot(1,3,sp)
    if sp==1
        toimage=pval_rms_matrix;
        todir=h_rms;
    elseif sp ==2
        toimage=pval_rms_n1_matrix;
        todir=h_rmsN1;
    elseif sp==3
        toimage=pval_rms_n2_matrix;
        todir=h_rmsN2;
    end
    image(toimage','CDataMapping','scaled')
    colormap(master_ColorMaps('-viridis'));
    clim([0.0001 0.05])
    % set(gca,'ColorScale','log')
    set(gca, 'xtick', [])
    set(gca, 'YTick', [])
    for x1=1:4
        for x2=1:4
            if toimage(x2,x1) < 0.0001
                text(x2,x1,'p<0.0001','Color','k','HorizontalAlignment','center','FontSize',18);
            else
                if toimage(x2,x1) > 0.05
                    text(x2,x1,num2str(toimage(x2,x1),'%.4f'),'Color','w','HorizontalAlignment','center','FontSize',18);
                else
                text(x2,x1,num2str(toimage(x2,x1),'%.4f'),'Color','k','HorizontalAlignment','center','FontSize',18);
                end
            end
            % if 1, on is more; if -1, off is more. Thus, when doing on to
            % off, 1 means decrease/blue and -1 means increase/red
            if todir(x2,x1) == 1
                rectangle('Position', [x2-.5, x1-.5, 1, 1], 'EdgeColor', 'b', 'LineWidth', 4);
            elseif todir(x2,x1) == -1
                rectangle('Position', [x2-.5, x1-.5, 1, 1], 'EdgeColor', 'r', 'LineWidth', 4);
            end
        end
    end
    box off

    % if sp==3
    %     colorbar
    % end
end

print(fullfile(savedir,[char(datetime('today','Format','uuuu-MM-dd')) '.png']),'-dpng','-r300')

%% create colorbar
figure('Units','normalized','Position',[0 0 .5 1])
% ax = axes;
colormap(master_ColorMaps('-viridis'));
% set(gca,'ColorScale','log')
c=colorbar(gca);
% c.Label.String = 'p-values';c.Label.Rotation=270;
c.FontSize=30;
c.FontWeight='bold';
c.LineWidth = 1;
clim([0.0001 0.05]);
% c.Limits=[0.00000001 1];
% c.Ticks = [0.0001,0.001,0.01,0.05,1];
c.Ticks = [0.0001,0.01,0.02,0.03,0.04,0.05];
c.TickLabels = {'<0.0001','0.01','0.02','0.03','0.04','≥0.05'};

a=get(c); %gets properties of colorbar
a = a.Position; %gets the positon and size of the color bar
set(c,'Units','normalized','Position',[.5 .1 .05 .5])% To change size
set(gca,'Visible','off')
% ax.Visible = 'off';
set(gcf,'color','w')
set(gcf, 'InvertHardcopy', 'off') % to set saving white color as white

print(fullfile(savedir,[char(datetime('today','Format','uuuu-MM-dd')) '_colorbar_v2.png']),'-dpng','-r300')


