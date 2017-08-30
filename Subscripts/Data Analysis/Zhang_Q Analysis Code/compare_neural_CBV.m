function compare_neural_CBV(x,y,xtype)
% compare_neural_CBV: generate scatter plot to compare the change of
%                     variations of CBV, neural activity at rest
% INPUT:
%       x: x axis data, should be nerual data
%       y: y axis data, should be CBV data
%       xtype: x data type: Gamma or MUA
% OUTPUT:
%       None except figures
% 
% MODIFIED BY ATW TO MATCH OTHER PLOTS

% figure('color','w','position',[100 100 300 300]); % use absolute size
figure;
% plot(x, y, 'marker','o',...
%            'MarkerFaceColor',[0.5 0.5 0.5],...
%            'MarkerEdgeColor','k',...
%            'linestyle','none',...
%            'MarkerSize',10);
plot(x, y, 'marker','o',...
           'MarkerFaceColor',[0.7 0.7 0.7],...
           'MarkerEdgeColor','k',...
           'linestyle','none');
% set(gca, 'box','off',...
%          'LineWidth',2,...
%          'xlim',[0 1.4], 'ylim', [0 1.4],...
%          'xtick',0:0.2:1.4, 'ytick',0:0.2:1.4,...
%          'xticklabel',{'0','', '', '', '', '1','',''},...
%          'yticklabel',{'0','', '', '', '', '1','',''});

hold on;

scatter(1.55,mean(y),'MarkerEdgeColor','k','MarkerFaceColor','k');
errorbar(1.55,mean(y),std(y),'k');
scatter(mean(x),1.55,'MarkerEdgeColor','k','MarkerFaceColor','k');
herrorbar(mean(x),1.55,std(x),std(x),'ko');
plot(0:2,0:2,'k--');
ylim([0 1.6])
xlim([0 1.6])
title('Pharmacology effects compared to aCSF: Musc.+Adr. Blockers')
ylabel('\DeltaR/R \sigma_{musc.+praz.+yoh.+prop} / \DeltaR/R \sigma_{aCSF}')
xlabel([xtype ' \sigma_{musc.+praz.+yoh.+prop} / ' xtype ' \sigma_{aCSF}'])
axis square;
hold off;

% line([1.4,1.4], [nanmean(y)-nanstd(y), nanmean(y)+nanstd(y)], 'color','k','linewidth',2);
% line([nanmean(x)-nanstd(x), nanmean(x)+nanstd(x)],[1.4, 1.4], 'color','k','linewidth',2);
% 
% plot([nanmean(x), 1.4],[1.4, nanmean(y)],...
%     'marker','o',...
%     'MarkerFaceColor','w',...
%     'MarkerEdgeColor','k',...
%     'linestyle','none',...
%     'MarkerSize',10);
% plot(0:0.7:1.4, 0:0.7:1.4,'linestyle','--', 'color',[0.5 0.5 0.5]);
% 
% if strcmp(xtype,'MUA')
%     xlabel('\sigma MUA Power_{treatment/vehicle}','fontsize',14,'fontname','arial');
%     ylabel('\sigma \DeltaR/R0_{treatment/vehicle}','fontsize',14,'fontname','arial');
% end
% if strcmp(xtype,'Gamma')
%     xlabel('\sigma Gamma Power_{treatment/vehicle}','fontsize',14,'fontname','arial');
%     ylabel('\sigma \DeltaR/R0_{treatment/vehicle}','fontsize',14,'fontname','arial');
% end
% end