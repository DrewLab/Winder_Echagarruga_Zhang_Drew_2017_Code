function PlotHRFAttributeComparison_LayerDependence(HRFatt,testfields,const_field,layer)
%   function [fig_handles] = PlotHRFAttributeComparison_LayerDependence(HRFatt,testfields,const_field,layer)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Plots the HRF attributes with the layer of the implanted
%   electrode designated by color
%   
%_______________________________________________________________
%   PARAMETERS:             
%               HRFatt - [struct] contains the HRF attibutes
%
%               testfields - [cell array] contains the behavioral
%               categories to be tested
%
%               const_field - [cell array] contains the type of neural
%               activity to be used
%
%               layer - [array] numerical designators of the cortical layer
%               of the electrode for each HRF attribute
%                   1- Supragranular
%                   2- Granular
%                   3- Infragranular
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

% Setup colors
Supra = layer==1;
Gran = layer==2;
Infra = layer==3;

figure;
subplot(131)
plot(1:2,NAmps(2:end,Supra),'c','Linewidth',2,'Marker','o',...
    'MarkerEdgeColor','k','MarkerFaceColor','c');
hold on;
plot(1:2,NAmps(2:end,Gran),'b','Linewidth',2,'Marker','o',...
    'MarkerEdgeColor','k','MarkerFaceColor','b');
plot(1:2,NAmps(2:end,Infra),'Color',[0.7 0.7 0.7],'Linewidth',2,'Marker','o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7]);
ylabel('HRF Amplitude / Sens. Ev. HRF Amplitude');
ylim([0 10]);
xlim([0.5 2.5])
title(['Amplitude: ' const_field ' HRF'])
ax1 = gca;
set(ax1,'XTick',1:2,...
    'XTickLabel',{'Whisk','Rest'});
set(gca,'yscale','log')

subplot(132)
h1 = plot(1:2,NTTPs(2:end,Supra),'c','Linewidth',2,'Marker','o',...
    'MarkerEdgeColor','k','MarkerFaceColor','c');
hold on;
h2 = plot(1:2,NTTPs(2:end,Gran),'b','Linewidth',2,'Marker','o',...
    'MarkerEdgeColor','k','MarkerFaceColor','b');
h3 = plot(1:2,NTTPs(2:end,Infra),'Color',[0.7 0.7 0.7],'Linewidth',2,'Marker','o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7]);
ylabel('HRF Time to Peak / Sens. Ev. HRF Time to Peak')
ylim([0 2])
xlim([0.5 2.5])
title(['Time to Peak: ' const_field ' HRF'])
ax2 = gca;
set(ax2,'XTick',1:2,...
    'XTickLabel',{'Whisk','Rest'});
legend([h1(1),h2(1),h3(1)],{'Supragranular','Granular','Infragranular'})

subplot(133)
plot(1:2,NFWHMs(2:end,Supra),'c','Linewidth',2,'Marker','o',...
    'MarkerEdgeColor','k','MarkerFaceColor','c');
hold on;
plot(1:2,NFWHMs(2:end,Gran),'b','Linewidth',2,'Marker','o',...
    'MarkerEdgeColor','k','MarkerFaceColor','b');
plot(1:2,NFWHMs(2:end,Infra),'Color',[0.7 0.7 0.7],'Linewidth',2,'Marker','o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7]);
ylabel('HRF FWHM / Sens. Ev. HRF FWHM')
ylim([0 2]);
xlim([0.5 2.5])
title(['Full Width Half Max: ' const_field ' HRF'])
ax3 = gca;
set(ax3,'XTick',1:2,...
    'XTickLabel',{'Whisk','Rest'});