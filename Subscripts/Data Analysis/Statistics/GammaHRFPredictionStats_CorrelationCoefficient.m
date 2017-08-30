function [Stats] = GammaHRFPredictionStats_CorrelationCoefficient(AveR,IndR)
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
[h1,p1] = adtest(AveR.Gamma(:,1)); % Sensory Evoked, not normal
[h2,p2] = adtest(AveR.Gamma(:,2)); % Extended Movement, not normal
[h3,p3] = adtest(AveR.Gamma(:,3)); % Volitional Whisk, normal
    % Omnibus test - Friedman
[Stats.AveR.Omnibus.pval,tbl,stat] =...
    friedman(AveR.Gamma(:,1:3),1,'off');
Stats.AveR.Omnibus.ChiSq = tbl{2,5};
Stats.AveR.Omnibus.df = [tbl{2,3} tbl{3,3}];
    % Post Hoc - Wilcoxon Signed Rank
Bonferroni = 3;
        % Sensory Evoked vs. Extendend Movement
[Stats.AveR.SEvVsExtM.pval,~,stat] =...
    signrank(AveR.Gamma(:,1),AveR.Gamma(:,2),'method','approximate');
Stats.AveR.SEvVsExtM.Pcorrected = ...
    Stats.AveR.SEvVsExtM.pval*Bonferroni;
Stats.AveR.SEvVsExtM.Zval = stat.zval;
        % Sensory Evoked vs. Volitional Whisking
[Stats.AveR.SEvVsVW.pval,~,stat] =...
    signrank(AveR.Gamma(:,1),AveR.Gamma(:,3),'method','approximate');
Stats.AveR.SEvVsVW.Pcorrected = ...
    Stats.AveR.SEvVsVW.pval*Bonferroni;
Stats.AveR.SEvVsVW.Zval = stat.zval;
        % Extended Movement vs. Volitional Whisking
[Stats.AveR.ExtMVsVW.pval,~,stat] =...
    signrank(AveR.Gamma(:,2),AveR.Gamma(:,3),'method','approximate');
Stats.AveR.ExtMVsVW.Pcorrected = ...
    Stats.AveR.ExtMVsVW.pval*Bonferroni;
Stats.AveR.ExtMVsVW.Zval = stat.zval;

% Individual Data
%-------------
% Test Assumptions - All distributions are normal
[h1,p1] = adtest(IndR.Gamma(:,1)); % Sensory Evoked, normal
[h2,p2] = adtest(IndR.Gamma(:,2)); % Extended Movement, normal
[h3,p3] = adtest(IndR.Gamma(:,3)); % Volitional Whisk, normal
[pvar,statvar] = vartestn(IndR.Gamma,'Display','off'); % Variances are equal
    % Omnibus test - Two-way ANOVA
[pval,tbl,stat] = anova2(IndR.Gamma,1,'off');
Stats.IndR.Omnibus.pval = pval(1);
Stats.IndR.Omnibus.Fval = tbl{2,5};
Stats.IndR.Omnibus.df = [tbl{2,3} tbl{3,3}];
    % Post-Hoc, Paired T-test
    Bonferroni = 6;
        % Sensory Evoked vs. Extended Movement
[~,Stats.IndR.SEvVsExtM.pval,~,stat] =...
    ttest(IndR.Gamma(:,1),IndR.Gamma(:,2));
Stats.IndR.SEvVsExtM.Pcorrected =...
    Stats.IndR.SEvVsExtM.pval*Bonferroni;
Stats.IndR.SEvVsExtM.tstat = stat.tstat;
Stats.IndR.SEvVsExtM.df = stat.df;
        % Sensory Evoked vs. Whisking
[~,Stats.IndR.SEvVsVW.pval,~,stat] =...
    ttest(IndR.Gamma(:,1),IndR.Gamma(:,3));
Stats.IndR.SEvVsVW.Pcorrected =...
    Stats.IndR.SEvVsVW.pval*Bonferroni;
Stats.IndR.SEvVsVW.tstat = stat.tstat;
Stats.IndR.SEvVsVW.df = stat.df;
        % Sensory Evoked vs. Rest
[~,Stats.IndR.SEvVsRest.pval,~,stat] =...
    ttest(IndR.Gamma(:,1),IndR.Gamma(:,4));
Stats.IndR.SEvVsRest.Pcorrected =...
    Stats.IndR.SEvVsRest.pval*Bonferroni;
Stats.IndR.SEvVsRest.tstat = stat.tstat;
Stats.IndR.SEvVsRest.df = stat.df;
        % Extended Movement vs. Whisking
[~,Stats.IndR.ExtMVsVW.pval,~,stat] =...
    ttest(IndR.Gamma(:,2),IndR.Gamma(:,3));
Stats.IndR.ExtMVsVW.Pcorrected =...
    Stats.IndR.ExtMVsVW.pval*Bonferroni;
Stats.IndR.ExtMVsVW.tstat = stat.tstat;
Stats.IndR.ExtMVsVW.df = stat.df;
        % Extended Movement vs. Rest
[~,Stats.IndR.ExtMVsRest.pval,~,stat] =...
    ttest(IndR.Gamma(:,2),IndR.Gamma(:,4));
Stats.IndR.ExtMVsRest.Pcorrected =...
    Stats.IndR.ExtMVsRest.pval*Bonferroni;
Stats.IndR.ExtMVsRest.tstat = stat.tstat;
Stats.IndR.ExtMVsRest.df = stat.df;
        % Volitional Whisking vs. Rest
[~,Stats.IndR.VWVsRest.pval,~,stat] =...
    ttest(IndR.Gamma(:,3),IndR.Gamma(:,4));
Stats.IndR.VWVsRest.Pcorrected =...
    Stats.IndR.VWVsRest.pval*Bonferroni;
Stats.IndR.VWVsRest.tstat = stat.tstat;
Stats.IndR.VWVsRest.df = stat.df;