function [CBV_fluc_ratio,MUAPower_fluc_ratio,GammaPower_fluc_ratio] = ...
    CBV_neural_fluctuations_plot(MUData,GamData, CBVData)
% CBV_neural_fluctions_plot: plot oscillations of cerebral blood volume 
%                            (CBV) signal vs. neural activity signal after 
%                            infusion of aCSF vs. muscimol + noradrenergic 
%                            antagonists (prazosin, atipamezole and propronolol,
%                            all at 1 mM) cocktail. ATW edits: added inputs
%                            and outputs to integrate data for plotting and
%                            output the results for statistical testing.
% Input:
%       MUData - MUA results from ATW experiments
%
%       GamData - Gamma power results from ATW experiments
%
%       CBVData - CBV results from ATW experiments
%
% Output:
%       CBV_fluc_ratio - ratio of the pharmacological CBV changes compared
%       to aCSF
%
%       MUAPower_fluc_ration- ratio of the pharmacological MUA changes
%       compared to aCSF
%
%       GammaPower_fluc_ratio - ratio of the pharmacological Gamma power
%       changes compared to aCSF
%
% See also:
%       compute_resting_variance.m
%       compare_neural_CBV.m

clc;

%% Load data
data = load('data.mat');

%% Calculate resting variance for CBV and neural signal
results = compute_resting_variance(data);

%% Plot results
animalList = fieldnames(data);
drugList = fieldnames(data.(animalList{1}));

% organize results for plot
for animalIdx = 1:numel(animalList) % all animals
    animal = animalList{animalIdx};
    for drugIdx = 1:numel(drugList) % all drugs
        drug = drugList{drugIdx};
        GammaPower_ave(animalIdx,drugIdx) = results.(animal).(drug).post.GammaPower.ave;
        MUAPower_ave(animalIdx,drugIdx) = results.(animal).(drug).post.MUAPower.ave;
        
        CBV_sd(animalIdx,drugIdx) = results.(animal).(drug).post.CBV.sd;
        GammaPower_sd(animalIdx,drugIdx) = results.(animal).(drug).post.GammaPower.sd;        
        MUAPower_sd(animalIdx,drugIdx) = results.(animal).(drug).post.MUAPower.sd;
    end
end     

CBV_fluc_ratio = CBV_sd(:,2)./CBV_sd(:,1);
CBV_fluc_ratio(end+1) = nonzeros(CBVData.Means);
MUAPower_fluc_ratio = MUAPower_sd(:,2)./MUAPower_sd(:,1).*(MUAPower_ave(:,2)./MUAPower_ave(:,1));
MUAPower_fluc_ratio(end+1) = nonzeros(MUData.Means);
compare_neural_CBV(MUAPower_fluc_ratio,CBV_fluc_ratio,'MUA');
set(gcf,'name','Figure 5i','numbertitle','off');

CBV_fluc_ratio = CBV_sd(:,2)./CBV_sd(:,1);
GammaPower_fluc_ratio = GammaPower_sd(:,2)./GammaPower_sd(:,1)...
                        .*(GammaPower_ave(:,2)./GammaPower_ave(:,1));
GammaPower_fluc_ratio(end+1) = nonzeros(GamData.Means);
CBV_fluc_ratio(end+1) = nonzeros(CBVData.Means);
compare_neural_CBV(GammaPower_fluc_ratio,CBV_fluc_ratio,'Gamma');
set(gcf,'name','Figure 5l','numbertitle','off');
end