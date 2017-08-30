function [] = PlotHRFAttributeComparison(HRFatt,testfields,const_field)
%   function [] =  ()
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: 
%   
%_______________________________________________________________
%   PARAMETERS:             
%                               
%_______________________________________________________________
%   RETURN:                     
%                               
%_______________________________________________________________

if not(iscell(const_field))
    error('Const_field must be a cell array')
end
if length(const_field)>1
    error('Input only on set of comparisons at a time, len const_field should = 1');
end

rowinds = ismember(HRFatt.dataTypes,[testfields const_field]);
colinds = ismember(HRFatt.behaviors,[testfields const_field]);
Amps = squeeze(HRFatt.Amps(rowinds,colinds,:));
TTPs = squeeze(HRFatt.TTPs(rowinds,colinds,:));
FWHMs = squeeze(HRFatt.FWHMs(rowinds,colinds,:));

NAmps = Amps./(ones(size(Amps,1),1)*Amps(1,:));
NTTPs = TTPs./(ones(size(TTPs,1),1)*TTPs(1,:));
NFWHMs = FWHMs./(ones(size(FWHMs,1),1)*FWHMs(1,:));

figure;
subplot(131)
boxplot(NAmps',1:length(testfields));
hold on;
Att_inds = repmat(1:length(testfields),size(NAmps,2),1)'; 
colors = get(gca,'ColorOrder');
scatter(Att_inds(:),NAmps(:),'MarkerEdgeColor',colors(1,:));
ylabel('HRF Amplitude / Sens. Ev. HRF Amplitude');
ylim([0 10]);
xlim([1.5 3.5])
title(['Amplitude: ' const_field ' HRF'])
ax1 = gca;
set(ax1,'XTick',min(Att_inds(:)):1:max(Att_inds(:)),...
    'XTickLabel',{'Sens.Ev.','Whisk','Rest'});
set(gca,'yscale','log')

subplot(132)
boxplot(NTTPs',1:length(testfields));
hold on;
Att_inds = repmat(1:length(testfields),size(NTTPs,2),1)'; 
colors = get(gca,'ColorOrder');
scatter(Att_inds(:),NTTPs(:),'MarkerEdgeColor',colors(1,:));
ylabel('HRF Time to Peak / Sens. Ev. Time to Peak')
ylim([0 2])
xlim([1.5 3.5])
title(['Time to Peak: ' const_field ' HRF'])
ax2 = gca;
set(ax2,'XTick',min(Att_inds(:)):1:max(Att_inds(:)),...
    'XTickLabel',{'Sens.Ev.','Whisk','Rest'});

subplot(133)
boxplot(NFWHMs',1:length(testfields));
hold on;
Att_inds = repmat(1:length(testfields),size(NFWHMs,2),1)'; 
colors = get(gca,'ColorOrder');
scatter(Att_inds(:),NFWHMs(:),'MarkerEdgeColor',colors(1,:));
ylabel('HRF FWHM / Sens. Ev. HRF FWHM')
ylim([0 2]);
xlim([1.5 3.5])
title(['Full Width Half Max: ' const_field ' HRF'])
ax3 = gca;
set(ax3,'XTick',min(Att_inds(:)):1:max(Att_inds(:)),...
    'XTickLabel',{'Sens.Ev.','Whisk','Rest'});