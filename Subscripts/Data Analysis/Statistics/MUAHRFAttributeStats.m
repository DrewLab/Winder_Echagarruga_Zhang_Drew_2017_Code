function [Stat] = MUAHRFAttributeStats(HRFattributes)
%   function [Stat] = MUAHRFAttributeStats(HRFattributes)
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
Amplitudes = squeeze(HRFattributes.Amps(2,:,:));
[h1,p1] = adtest(Amplitudes(2,:)./Amplitudes(1,:)); % not normal
[h2,p2] = adtest(Amplitudes(3,:)./Amplitudes(1,:)); % not normal
[Stat.MUA.Amp.Whisk.pval,~,stat] = ...
    signrank(Amplitudes(2,:)./Amplitudes(1,:),...
    Amplitudes(1,:)./Amplitudes(1,:),'method','approximate');
Stat.MUA.Amp.Whisk.Zval = stat.zval;
Stat.MUA.Amp.Whisk.Pcorrected = Stat.MUA.Amp.Whisk.pval*Bonferroni;
[Stat.MUA.Amp.Rest.pval,~,stat] = ...
    signrank(Amplitudes(3,:)./Amplitudes(1,:),...
    Amplitudes(1,:)./Amplitudes(1,:),'method','approximate');
Stat.MUA.Amp.Rest.Zval = stat.zval;
Stat.MUA.Amp.Rest.Pcorrected = Stat.MUA.Amp.Rest.pval*Bonferroni;

% Stats - MUA HRF Time to peak
Bonferroni = 2;
TTPs = squeeze(HRFattributes.TTPs(2,:,:));
[h1,p1] = adtest(TTPs(2,:)./TTPs(1,:)); % normal
[h2,p2] = adtest(TTPs(3,:)./TTPs(1,:)); % normal
[~,Stat.MUA.TTP.Whisk.pval,~,stat] = ttest(TTPs(2,:)./TTPs(1,:),1);
Stat.MUA.TTP.Whisk.tstat = stat.tstat;
Stat.MUA.TTP.Whisk.df = stat.df;
Stat.MUA.TTP.Whisk.Pcorrected = Stat.MUA.TTP.Whisk.pval*Bonferroni;
[~,Stat.MUA.TTP.Rest.pval,~,stat] = ttest(TTPs(3,:)./TTPs(1,:),1);
Stat.MUA.TTP.Rest.tstat = stat.tstat;
Stat.MUA.TTP.Rest.df = stat.df;
Stat.MUA.TTP.Rest.Pcorrected = Stat.MUA.TTP.Rest.pval*Bonferroni;

% Stats - Gamma HRF Full width Half Max
Bonferroni = 2;
FWHMs = squeeze(HRFattributes.FWHMs(2,:,:));
[h1,p1] = adtest(FWHMs(2,:)./FWHMs(1,:)); % normal
[h2,p2] = adtest(FWHMs(3,:)./FWHMs(1,:)); % normal
[~,Stat.MUA.FWHM.Whisk.pval,~,stat] = ttest(FWHMs(2,:)./FWHMs(1,:),1);
Stat.MUA.FWHM.Whisk.tstat = stat.tstat;
Stat.MUA.FWHM.Whisk.df = stat.df;
Stat.MUA.FWHM.Whisk.Pcorrected = Stat.MUA.FWHM.Whisk.pval*Bonferroni;
[~,Stat.MUA.FWHM.Rest.pval,~,stat] = ttest(FWHMs(3,:)./FWHMs(1,:),1);
Stat.MUA.FWHM.Rest.tstat = stat.tstat;
Stat.MUA.FWHM.Rest.df = stat.df;
Stat.MUA.FWHM.Rest.Pcorrected = Stat.MUA.FWHM.Rest.pval*Bonferroni;