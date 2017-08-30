function [Stats] = GammaHRFPredictionStats_Rsquared(AveR2,IndR2)
%   function [Stats] = GammaHRFPredictionStats_Rsquared(AveR2,IndR2)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Calculates the statistics for the gamma-based predictions
%   of averaged and individual CBV fluctuations
%
%_______________________________________________________________
%   PARAMETERS:
%                   AveR2 - [struct] contains the coefficients of
%                   determination for averaged CBV
%
%                   IndR2 - [struct] contains the coefficients of
%                   determination for the individual CBV
%_______________________________________________________________
%   RETURN:
%                   Stats - [struct] the results of the statistical tests
%_______________________________________________________________

% Averaged Data
%-------------
    % Test Assumptions - All distributions are normal
[h1,p1] = adtest(AveR2.Gamma(:,1)); % Sensory Evoked, normal
[h2,p2] = adtest(AveR2.Gamma(:,2)); % Extended Movement, normal
[h3,p3] = adtest(AveR2.Gamma(:,3)); % Volitional Whisk, normal
[pvar,statvar] = vartestn(AveR2.Gamma,'Display','off'); % Variances are equal
    % Omnibus test
[pval,tbl,stat] = anova2(AveR2.Gamma(:,1:3),1,'off');
Stats.AveR2.Omnibus.pval = pval(1);
Stats.AveR2.Omnibus.Fval = tbl{2,5};
Stats.AveR2.Omnibus.df = [tbl{2,3} tbl{3,3}];
    % Post-Hoc, Sensory Evoked vs. Extended Movement
Bonferroni = 3;
[~,Stats.AveR2.SEvVsExtM.pval,~,stat] =...
    ttest(AveR2.Gamma(:,1),AveR2.Gamma(:,2));
Stats.AveR2.SEvVsExtM.Pcorrected =...
    Stats.AveR2.SEvVsExtM.pval*Bonferroni;
Stats.AveR2.SEvVsExtM.tstat = stat.tstat;
Stats.AveR2.SEvVsExtM.df = stat.df;
    % Post-Hoc, Sensory Evoked vs. Volitional Whisk
[~,Stats.AveR2.SEvVsVW.pval,~,stat] =...
    ttest(AveR2.Gamma(:,1),AveR2.Gamma(:,3));
Stats.AveR2.SEvVsVW.Pcorrected =...
    Stats.AveR2.SEvVsVW.pval*Bonferroni;
Stats.AveR2.SEvVsVW.tstat = stat.tstat;
Stats.AveR2.SEvVsVW.df = stat.df;
    % Post-Hoc, Extended Movement vs. Volitional Whisk
[~,Stats.AveR2.ExtMVsVW.pval,~,stat] =...
    ttest(AveR2.Gamma(:,2),AveR2.Gamma(:,3));
Stats.AveR2.ExtMVsVW.Pcorrected =...
    Stats.AveR2.ExtMVsVW.pval*Bonferroni;
Stats.AveR2.ExtMVsVW.tstat = stat.tstat;
Stats.AveR2.ExtMVsVW.df = stat.df;
  %-----------------
  % Individual Data
    % Test Assumptions - Whisking R^2 distribution not normal
[h1,p1] = adtest(IndR2.Gamma(:,1)); % Sensory Evoked
[h2,p2] = adtest(IndR2.Gamma(:,2)); % Extended Movement
[h3,p3] = adtest(IndR2.Gamma(:,3)); % Volitional Whisk, not normal
[h4,p4] = adtest(IndR2.Gamma(:,4)); % Rest, normal
    % Omnibus test
[Stats.IndR2.Omnibus.pval,tbl,stat] =...
    friedman(IndR2.Gamma,1,'off');
Stats.IndR2.Omnibus.ChiSq = tbl{2,5};
Stats.IndR2.Omnibus.df = [tbl{2,3} tbl{3,3}];
Bonferroni = 6;
    % Post-Hoc, Sensory Evoked vs. Extended Movement
[Stats.IndR2.SEvVsExtM.pval,~,stat] =...
    signrank(IndR2.Gamma(:,1),IndR2.Gamma(:,2),'method','approximate');
Stats.IndR2.SEvVsExtM.Pcorrected = ...
    Stats.IndR2.SEvVsExtM.pval*Bonferroni;
Stats.IndR2.SEvVsExtM.Zval = stat.zval;
    % Post-Hoc, Sensory Evoked vs. Volitional Whisk
[Stats.IndR2.SEvVsVW.pval,~,stat] =...
    signrank(IndR2.Gamma(:,1),IndR2.Gamma(:,3),'method','approximate');
Stats.IndR2.SEvVsVW.Pcorrected = ...
    Stats.IndR2.SEvVsVW.pval*Bonferroni;
Stats.IndR2.SEvVsVW.Zval = stat.zval;
    % Post-Hoc, Sensory Evoked vs. Rest
[Stats.IndR2.SEvVsRest.pval,~,stat] =...
    signrank(IndR2.Gamma(:,1),IndR2.Gamma(:,4),'method','approximate');
Stats.IndR2.SEvVsRest.Pcorrected = ...
    Stats.IndR2.SEvVsRest.pval*Bonferroni;
Stats.IndR2.SEvVsRest.Zval = stat.zval;
    % Post-Hoc, Extended Movement vs. Volitional Whisking
[Stats.IndR2.ExtMVsVW.pval,~,stat] =...
    signrank(IndR2.Gamma(:,2),IndR2.Gamma(:,3),'method','approximate');
Stats.IndR2.ExtMVsVW.Pcorrected = ...
    Stats.IndR2.ExtMVsVW.pval*Bonferroni;
Stats.IndR2.ExtMVsVW.Zval = stat.zval;
    % Post-Hoc, Extended Movement vs. Rest
[Stats.IndR2.ExtMVsRest.pval,~,stat] =...
    signrank(IndR2.Gamma(:,2),IndR2.Gamma(:,4),'method','approximate');
Stats.IndR2.ExtMVsRest.Pcorrected = ...
    Stats.IndR2.ExtMVsRest.pval*Bonferroni;
Stats.IndR2.ExtMVsRest.Zval = stat.zval;
    % Post-Hoc, Volitional Whisking Vs. Rest
[Stats.IndR2.VWVsRest.pval,~,stat] =...
    signrank(IndR2.Gamma(:,3),IndR2.Gamma(:,4),'method','approximate');
Stats.IndR2.VWVsRest.Pcorrected = ...
    Stats.IndR2.VWVsRest.pval*Bonferroni;
Stats.IndR2.VWVsRest.Zval = stat.zval;