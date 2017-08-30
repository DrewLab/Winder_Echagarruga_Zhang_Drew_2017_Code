function [] = Pharmacology_PlotNeuroVsCBV(PharmData,animals,NeurType,CBVType,InfusionType)
%   function [] = Pharmacology_PlotNeuroVsCBV(PharmData,animals,NeurType,CBVType,InfusionType)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Plots a comparison of the amplitude  of the resting neural
%   and CBV fluctuations
%_______________________________________________________________
%   PARAMETERS:
%               PharmData - [struct] contains the amplitude measures and is
%               obtained from the scripts:
%                   SingleAnimalAnalysisMaster_InfusionTrials_NormByPreInfusion.m
%                   Infusion_PopulationAnalysis.m
%
%               animals - [cell array] animal IDs
%
%               NeurType - [string] the type of neural activity for
%               comparions to CBV
%
%               CBVType - [string] designates the CBV ROI to be used
%
%               InfusionType - [string] the type of infusion
%_______________________________________________________________
%   RETURN:
%
%_______________________________________________________________
figure;
AnimalNeur = [];
AnimalCBV = [];
ii=1;
for a = 1:length(animals)
    animal = animals{a};
    if not(isfield(PharmData.(NeurType),animal))
        continue;
    end
    if isempty(PharmData.(NeurType).(animal))
        continue;
    elseif isempty(PharmData.(CBVType).(animal))
        continue;
    end
    scatter(mean(PharmData.(NeurType).(animal)),mean(PharmData.(CBVType).(animal)),...
        'MarkerEdgeColor','k','MarkerFaceColor',[0.7 0.7 0.7]);
    hold on;
    AnimalNeur(ii) = mean(PharmData.(NeurType).(animal));
    AnimalCBV(ii) = mean(PharmData.(CBVType).(animal));
    ii = ii+1;
end
hold on;
scatter(1.55,mean(AnimalCBV),'MarkerEdgeColor','k','MarkerFaceColor','k');
errorbar(1.55,mean(AnimalCBV),std(AnimalCBV),'k');
scatter(mean(AnimalNeur),1.55,'MarkerEdgeColor','k','MarkerFaceColor','k');
herrorbar(mean(AnimalNeur),1.55,std(AnimalNeur),std(AnimalNeur),'ko');
plot(0:2,0:2,'k--');
ylim([0 1.6])
xlim([0 1.6])
title(['Pharmacology effects compared to aCSF: ' InfusionType])
xlabel([InfusionType ' ' NeurType '/' InfusionType ' aCSF'])
ylabel([InfusionType ' CBV/' InfusionType ' aCSF'])
axis square;
hold off;