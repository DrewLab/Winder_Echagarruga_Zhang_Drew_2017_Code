function [Stats] = MUAHRFPredictionStats_CorrelationCoefficient(AveR,IndR)
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


% Averaged Data
%-------------
% Test Assumptions - All distributions are not normal
[h1,p1] = adtest(AveR.MUA(:,1)); % Sensory Evoked, not normal
[h2,p2] = adtest(AveR.MUA(:,2)); % Extended Movement, not normal
[h3,p3] = adtest(AveR.MUA(:,3)); % Volitional Whisk, not normal
    % Omnibus test - Friedman
[Stats.AveR.Omnibus.pval,tbl,stat] =...
    friedman(AveR.MUA(:,1:3),1,'off');
Stats.AveR.Omnibus.ChiSq = tbl{2,5};
Stats.AveR.Omnibus.df = [tbl{2,3} tbl{3,3}];
    % Post Hoc - Wilcoxon Signed Rank
Bonferroni = 3;
        % Sensory Evoked vs. Extendend Movement
[Stats.AveR.SEvVsExtM.pval,~,stat] =...
    signrank(AveR.MUA(:,1),AveR.MUA(:,2),'method','approximate');
Stats.AveR.SEvVsExtM.Pcorrected = ...
    Stats.AveR.SEvVsExtM.pval*Bonferroni;
Stats.AveR.SEvVsExtM.Zval = stat.zval;
        % Sensory Evoked vs. Volitional Whisking
[Stats.AveR.SEvVsVW.pval,~,stat] =...
    signrank(AveR.MUA(:,1),AveR.MUA(:,3),'method','approximate');
Stats.AveR.SEvVsVW.Pcorrected = ...
    Stats.AveR.SEvVsVW.pval*Bonferroni;
Stats.AveR.SEvVsVW.Zval = stat.zval;
        % Extended Movement vs. Volitional Whisking
[Stats.AveR.ExtMVsVW.pval,~,stat] =...
    signrank(AveR.MUA(:,2),AveR.MUA(:,3),'method','approximate');
Stats.AveR.ExtMVsVW.Pcorrected = ...
    Stats.AveR.ExtMVsVW.pval*Bonferroni;
Stats.AveR.ExtMVsVW.Zval = stat.zval;

% Individual Data
%-------------
% Test Assumptions - All distributions are not normal
[h1,p1] = adtest(IndR.MUA(:,1)); % Sensory Evoked, not normal
[h2,p2] = adtest(IndR.MUA(:,2)); % Extended Movement, normal
[h3,p3] = adtest(IndR.MUA(:,3)); % Volitional Whisk, normal
    % Omnibus test
[Stats.IndR.Omnibus.pval,tbl,stat] =...
    friedman(IndR.MUA,1,'off');
Stats.IndR.Omnibus.ChiSq = tbl{2,5};
Stats.IndR.Omnibus.df = [tbl{2,3} tbl{3,3}];
Bonferroni = 6;
    % Post-Hoc, Sensory Evoked vs. Extended Movement
[Stats.IndR.SEvVsExtM.pval,~,stat] =...
    signrank(IndR.MUA(:,1),IndR.MUA(:,2),'method','approximate');
Stats.IndR.SEvVsExtM.Pcorrected = ...
    Stats.IndR.SEvVsExtM.pval*Bonferroni;
Stats.IndR.SEvVsExtM.Zval = stat.zval;
    % Post-Hoc, Sensory Evoked vs. Volitional Whisk
[Stats.IndR.SEvVsVW.pval,~,stat] =...
    signrank(IndR.MUA(:,1),IndR.MUA(:,3),'method','approximate');
Stats.IndR.SEvVsVW.Pcorrected = ...
    Stats.IndR.SEvVsVW.pval*Bonferroni;
Stats.IndR.SEvVsVW.Zval = stat.zval;
    % Post-Hoc, Sensory Evoked vs. Rest
[Stats.IndR.SEvVsRest.pval,~,stat] =...
    signrank(IndR.MUA(:,1),IndR.MUA(:,4),'method','approximate');
Stats.IndR.SEvVsRest.Pcorrected = ...
    Stats.IndR.SEvVsRest.pval*Bonferroni;
Stats.IndR.SEvVsRest.Zval = stat.zval;
    % Post-Hoc, Extended Movement vs. Volitional Whisking
[Stats.IndR.ExtMVsVW.pval,~,stat] =...
    signrank(IndR.MUA(:,2),IndR.MUA(:,3),'method','approximate');
Stats.IndR.ExtMVsVW.Pcorrected = ...
    Stats.IndR.ExtMVsVW.pval*Bonferroni;
Stats.IndR.ExtMVsVW.Zval = stat.zval;
    % Post-Hoc, Extended Movement vs. Rest
[Stats.IndR.ExtMVsRest.pval,~,stat] =...
    signrank(IndR.MUA(:,2),IndR.MUA(:,4),'method','approximate');
Stats.IndR.ExtMVsRest.Pcorrected = ...
    Stats.IndR.ExtMVsRest.pval*Bonferroni;
Stats.IndR.ExtMVsRest.Zval = stat.zval;
    % Post-Hoc, Volitional Whisking Vs. Rest
[Stats.IndR.VWVsRest.pval,~,stat] =...
    signrank(IndR.MUA(:,3),IndR.MUA(:,4),'method','approximate');
Stats.IndR.VWVsRest.Pcorrected = ...
    Stats.IndR.VWVsRest.pval*Bonferroni;
Stats.IndR.VWVsRest.Zval = stat.zval;