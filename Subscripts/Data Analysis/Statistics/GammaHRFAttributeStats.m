function [Stat] = GammaHRFAttributeStats(HRFattributes)
%   function [Stat] = GammaHRFAttributeStats(HRFattributes)
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

Bonferroni = 2;
Amplitudes = squeeze(HRFattributes.Amps(1,:,:));
[h1,p1] = adtest(Amplitudes(2,:)./Amplitudes(1,:)); % not normal
[h2,p2] = adtest(Amplitudes(3,:)./Amplitudes(1,:)); % not normal
[Stat.Gamma.Amp.Whisk.pval,~,stat] = ...
    signrank(Amplitudes(2,:)./Amplitudes(1,:),...
    Amplitudes(1,:)./Amplitudes(1,:),'method','approximate');
Stat.Gamma.Amp.Whisk.Zval = stat.zval;
Stat.Gamma.Amp.Whisk.Pcorrected = Stat.Gamma.Amp.Whisk.pval*Bonferroni;
[Stat.Gamma.Amp.Rest.pval,~,stat] = ...
    signrank(Amplitudes(3,:)./Amplitudes(1,:),...
    Amplitudes(1,:)./Amplitudes(1,:),'method','approximate');
Stat.Gamma.Amp.Rest.Zval = stat.zval;
Stat.Gamma.Amp.Rest.Pcorrected = Stat.Gamma.Amp.Rest.pval*Bonferroni;

% Stats - Gamma HRF Time to peak
Bonferroni = 2;
TTPs = squeeze(HRFattributes.TTPs(1,:,:));
[h1,p1] = adtest(TTPs(2,:)./TTPs(1,:)); % normal
[h2,p2] = adtest(TTPs(3,:)./TTPs(1,:)); % normal
[~,Stat.Gamma.TTP.Whisk.pval,~,stat] = ttest(TTPs(2,:)./TTPs(1,:),1);
Stat.Gamma.TTP.Whisk.tstat = stat.tstat;
Stat.Gamma.TTP.Whisk.df = stat.df;
Stat.Gamma.TTP.Whisk.Pcorrected = Stat.Gamma.TTP.Whisk.pval*Bonferroni;
[~,Stat.Gamma.TTP.Rest.pval,~,stat] = ttest(TTPs(3,:)./TTPs(1,:),1);
Stat.Gamma.TTP.Rest.tstat = stat.tstat;
Stat.Gamma.TTP.Rest.df = stat.df;
Stat.Gamma.TTP.Rest.Pcorrected = Stat.Gamma.TTP.Rest.pval*Bonferroni;

% Stats - Gamma HRF Full width Half Max
Bonferroni = 2;
FWHMs = squeeze(HRFattributes.FWHMs(1,:,:));
[h1,p1] = adtest(FWHMs(2,:)./FWHMs(1,:)); % normal
[h2,p2] = adtest(FWHMs(3,:)./FWHMs(1,:)); % normal
[~,Stat.Gamma.FWHM.Whisk.pval,~,stat] = ttest(FWHMs(2,:)./FWHMs(1,:),1);
Stat.Gamma.FWHM.Whisk.tstat = stat.tstat;
Stat.Gamma.FWHM.Whisk.df = stat.df;
Stat.Gamma.FWHM.Whisk.Pcorrected = Stat.Gamma.FWHM.Whisk.pval*Bonferroni;
[~,Stat.Gamma.FWHM.Rest.pval,~,stat] = ttest(FWHMs(3,:)./FWHMs(1,:),1);
Stat.Gamma.FWHM.Rest.tstat = stat.tstat;
Stat.Gamma.FWHM.Rest.df = stat.df;
Stat.Gamma.FWHM.Rest.Pcorrected = Stat.Gamma.FWHM.Rest.pval*Bonferroni;