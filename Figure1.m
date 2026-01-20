% CCEP_fig1_v1: inc vs dec

sbj_IDs = {'P1', 'P2', 'P3','P4'};
comb_sbjIDs = '';
for s = 1:length(sbj_IDs)
    sbj_ID = sbj_IDs{s};
    comb_sbjIDs = [comb_sbjIDs '_' sbj_ID];
end

savedir = fullfile(cd,'Figures');if ~isfolder(savedir); mkdir(savedir); end
Sbj_Metadata = makeSbj_Metadata(data_root, project_name, sbj_IDs{1}); % example Sbj_Metadata
res_dir=fullfile(Sbj_Metadata.project_root,'COMBINED_RESULTS','offon_collect_drectns_multi');
load(fullfile(res_dir,[comb_sbjIDs '_offon_allDirections.mat']),'collected')

for_crisscross = {'SOZ','EPZ','IZ','Healthy'};
N1N2 = {'N1','N2'};
time_amp = {'time','amp'};
time_amp_title = {'Latency','Amplitude'};
inc_dec = {'Inc','Dec'};
inc_dec_title = {'Increased','Decreased'}; % this direction is off to on > will reverse this during printing

%% Figure combined
figure('position',[0 0 700 1400])
for n=1:2 %N1 vs N2
    for t=1:2 % time vs amp
        for i = 1:2 % inc vs dec

subplot(4,2,subplotno(2,subplotno(2,n,t),i))
toplot=collected.SOZ2EPZ.([N1N2{n} '_' time_amp{t} '_dir']).(inc_dec{i});
image(100*toplot','CDataMapping','scaled')
colormap(master_ColorMaps('-viridis'));clim([0 25]);
% set(gca,'ColorScale','log','Zdir','reverse')
% c=colorbar;c.Label.String = 'Percent';c.Label.Rotation=270;c.Ticks=0:5:30;
% xlabel('Stimulation Channels')
set(gca, 'xtick', [])
set(gca, 'YTick', [])
% xticks(1:4)
% yticks(1:4)
% if i==1
% xticklabels({'','','',''})
% xticklabels(for_crisscross)
% yticklabels({'','','',''});
% else
% xticklabels(for_crisscross)
% yticklabels({'','','',''});
% % yticklabels(for_crisscross)
% end
% if i==1
% ylabel([N1N2{n} ' - ' time_amp_title{t}])
% end
set(gca, 'FontSize',10,'FontWeight','bold');
if n==1 && t==1
    title(inc_dec_title{i}, 'FontSize',16,'FontWeight','bold')
end
for x1=1:4
    for x2=1:4
        if 100*toplot(x2,x1) < 11
        text(x2,x1,[num2str(100*toplot(x2,x1),'%.3g') '%'],'Color','k','HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
        else
        text(x2,x1,[num2str(100*toplot(x2,x1),'%.3g') '%'],'Color','w','HorizontalAlignment','center','FontSize',12,'FontWeight','bold');    
        end
    end
end
% title([N1N2{n} ' - ' time_amp_title{t}])
        
        
        end
    end
end

print(fullfile(savedir,['Figure1_' char(datetime('today','Format','uuuu-MM-dd')) '.png']),'-dpng','-r300')

%% colormap - horizontal
figure('Units','normalized','Position',[0 0 .5 1])
% ax = axes;
colormap(master_ColorMaps('-viridis'));
clim([0 25]);
c=colorbar(gca);
% c.Label.String = 'Percentages';c.Label.Rotation=270;
c.FontSize=30;
c.FontWeight='bold';
c.LineWidth = 1;
c.Limits=[0 25];
c.Ticks = [0,5,10,15,20,25];
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

print(fullfile(savedir,['Figure1_colorbar_horizontal_' char(datetime('today','Format','uuuu-MM-dd')) '.png']),'-dpng','-r300')
