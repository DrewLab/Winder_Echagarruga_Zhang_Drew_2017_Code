function [Stats] = MUAHRFPredictionStats_Rsquared(AveR2,IndR2)
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


% Averaged data
%-------------
    % Test Assumptions
[h1,p1] = adtest(AveR2.MUA(:,1)); % Sensory Evoked, normal
[h2,p2] = adtest(AveR2.MUA(:,2)); % Extended Movement, not normal
[h3,p3] = adtest(AveR2.MUA(:,3)); % Volitional Whisk, normal
    % Omnibus test
[Stats.AveR2.Omnibus.pval,tbl,stat] =...
    friedman(AveR2.MUA(:,1:3),1,'off');
Stats.AveR2.Omnibus.ChiSq = tbl{2,5};
Stats.AveR2.Omnibus.df = [tbl{2,3} tbl{3,3}];
Bonferroni = 3;
    % Post-Hoc, Sensory Evoked vs. Extended Movement
[Stats.AveR2.SEvVsExtM.pval,~,stat] =...
    signrank(AveR2.MUA(:,1),AveR2.MUA(:,2),'method','approximate');
Stats.AveR2.SEvVsExtM.Pcorrected = ...
    Stats.AveR2.SEvVsExtM.pval*Bonferroni;
Stats.AveR2.SEvVsExtM.Zval = stat.zval;
    % Post-Hoc, Sensory Evoked vs. Volitional Whisk
[Stats.AveR2.SEvVsVW.pval,~,stat] =...
    signrank(AveR2.MUA(:,1),AveR2.MUA(:,3),'method','approximate');
Stats.AveR2.SEvVsVW.Pcorrected = ...
    Stats.AveR2.SEvVsVW.pval*Bonferroni;
Stats.AveR2.SEvVsVW.Zval = stat.zval;
    % Post-Hoc, Extended Movement vs. Volitional Whisking
[Stats.AveR2.ExtMVsVW.pval,~,stat] =...
    signrank(AveR2.MUA(:,2),AveR2.MUA(:,3),'method','approximate');
Stats.AveR2.ExtMVsVW.Pcorrected = ...
    Stats.AveR2.ExtMVsVW.pval*Bonferroni;
Stats.AveR2.ExtMVsVW.Zval = stat.zval;
  %-----------------
  % Individual Data
    % Test Assumptions
[h1,p1] = adtest(IndR2.MUA(:,1)); % Sensory Evoked, normal
[h2,p2] = adtest(IndR2.MUA(:,2)); % Extended Movement, not normal
[h3,p3] = adtest(IndR2.MUA(:,3)); % Volitional Whisk, normal
[h4,p4] = adtest(IndR2.MUA(:,4)); % Rest, normal
    % Omnibus test
[Stats.IndR2.Omnibus.pval,tbl,stat] =...
    friedman(IndR2.MUA,1,'off');
Stats.IndR2.Omnibus.ChiSq = tbl{2,5};
Stats.IndR2.Omnibus.df = [tbl{2,3} tbl{3,3}];
Bonferroni = 6;
    % Post-Hoc, Sensory Evoked vs. Extended Movement
[Stats.IndR2.SEvVsExtM.pval,~,stat] =...
    signrank(IndR2.MUA(:,1),IndR2.MUA(:,2),'method','approximate');
Stats.IndR2.SEvVsExtM.Pcorrected = ...
    Stats.IndR2.SEvVsExtM.pval*Bonferroni;
Stats.IndR2.SEvVsExtM.Zval = stat.zval;
    % Post-Hoc, Sensory Evoked vs. Volitional Whisk
[Stats.IndR2.SEvVsVW.pval,~,stat] =...
    signrank(IndR2.MUA(:,1),IndR2.MUA(:,3),'method','approximate');
Stats.IndR2.SEvVsVW.Pcorrected = ...
    Stats.IndR2.SEvVsVW.pval*Bonferroni;
Stats.IndR2.SEvVsVW.Zval = stat.zval;
    % Post-Hoc, Sensory Evoked vs. Rest
[Stats.IndR2.SEvVsRest.pval,~,stat] =...
    signrank(IndR2.MUA(:,1),IndR2.MUA(:,4),'method','approximate');
Stats.IndR2.SEvVsRest.Pcorrected = ...
    Stats.IndR2.SEvVsRest.pval*Bonferroni;
Stats.IndR2.SEvVsRest.Zval = stat.zval;
    % Post-Hoc, Extended Movement vs. Volitional Whisking
[Stats.IndR2.ExtMVsVW.pval,~,stat] =...
    signrank(IndR2.MUA(:,2),IndR2.MUA(:,3),'method','approximate');
Stats.IndR2.ExtMVsVW.Pcorrected = ...
    Stats.IndR2.ExtMVsVW.pval*Bonferroni;
Stats.IndR2.ExtMVsVW.Zval = stat.zval;
    % Post-Hoc, Extended Movement vs. Rest
[Stats.IndR2.ExtMVsRest.pval,~,stat] =...
    signrank(IndR2.MUA(:,2),IndR2.MUA(:,4),'method','approximate');
Stats.IndR2.ExtMVsRest.Pcorrected = ...
    Stats.IndR2.ExtMVsRest.pval*Bonferroni;
Stats.IndR2.ExtMVsRest.Zval = stat.zval;
    % Post-Hoc, Volitional Whisking Vs. Rest
[Stats.IndR2.VWVsRest.pval,~,stat] =...
    signrank(IndR2.MUA(:,3),IndR2.MUA(:,4),'method','approximate');
Stats.IndR2.VWVsRest.Pcorrected = ...
    Stats.IndR2.VWVsRest.pval*Bonferroni;
Stats.IndR2.VWVsRest.Zval = stat.zval;