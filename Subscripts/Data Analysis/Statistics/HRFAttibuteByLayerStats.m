function [LayerStat] = HRFAttibuteByLayerStats(HRFattributes,Layer)
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

Bonferroni = 3;
Amplitudes = squeeze(HRFattributes.Amps(1,:,:));
NAmpMat = ones((size(Amplitudes,1)-1),1)*Amplitudes(1,:);
NAmps = Amplitudes(2:3,:)./NAmpMat;
Sup_Amps = NAmps(:,Layer==1);
Gr_Amps = NAmps(:,Layer==2);
Infr_Amps = NAmps(:,Layer==3);
[h1,p1] = adtest(Sup_Amps(:));
[h2,p2] = adtest(Gr_Amps(:));
[h3,p3] = adtest(Infr_Amps(:));
    % Supragranalar vs Granular Amplitude
[LayerStat.Gamma.Amps.SupVsGr.pval,~,stat] = ranksum(Sup_Amps(:),Gr_Amps(:),...
    'method','approximate');
LayerStat.Gamma.Amps.SupVsGr.Pcorrected = LayerStat.Gamma.Amps.SupVsGr.pval*Bonferroni;
LayerStat.Gamma.Amps.SupVsGr.Zval = stat.zval;
    % Supragranalar vs Infragranular Amplitude
[LayerStat.Gamma.Amps.SupVsInfr.pval,~,stat] = ranksum(Sup_Amps(:),Infr_Amps(:),...
    'method','approximate');
LayerStat.Gamma.Amps.SupVsInfr.Pcorrected = LayerStat.Gamma.Amps.SupVsInfr.pval*Bonferroni;
LayerStat.Gamma.Amps.SupVsInfr.Zval = stat.zval;
    % Granular vs Infragranular Amplitude
[LayerStat.Gamma.Amps.GrVsInfr.pval,~,stat] = ranksum(Gr_Amps(:),Infr_Amps(:),...
    'method','approximate');
LayerStat.Gamma.Amps.GrVsInfr.Pcorrected = LayerStat.Gamma.Amps.GrVsInfr.pval*Bonferroni;
LayerStat.Gamma.Amps.GrVsInfr.Zval = stat.zval;

% Stats - Layer dependency of TTP - Gamma-based HRF
Bonferroni = 3;
TTPs = squeeze(HRFattributes.TTPs(1,:,:));
NTTPMat = ones((size(TTPs,1)-1),1)*TTPs(1,:);
NTTPs = TTPs(2:3,:)./NTTPMat;
Sup_TTPs = NTTPs(:,Layer==1);
Gr_TTPs = NTTPs(:,Layer==2);
Infr_TTPs = NTTPs(:,Layer==3);
[h1,p1] = adtest(Sup_TTPs(:));
[h2,p2] = adtest(Gr_TTPs(:));
[h3,p3] = adtest(Infr_TTPs(:));
[p,stat] = vartestn([Sup_TTPs(:); Gr_TTPs(:); Infr_TTPs(:)],...
    [ones(size(Sup_TTPs(:))); 2*ones(size(Gr_TTPs(:))); 3*ones(size(Infr_TTPs(:)))],...
    'Display','off');
    % Supragranalar vs Granular Time to Peak
[~,LayerStat.Gamma.TTPs.SupVsGr.pval,~,stat] = ttest2(Sup_TTPs(:),Gr_TTPs(:));
LayerStat.Gamma.TTPs.SupVsGr.Pcorrected = LayerStat.Gamma.TTPs.SupVsGr.pval*Bonferroni;
LayerStat.Gamma.TTPs.SupVsGr.tstat = stat.tstat;
LayerStat.Gamma.TTPs.SupVsGr.df = stat.df;
    % Supragranalar vs Infragranular Time to Peak
[~,LayerStat.Gamma.TTPs.SupVsInfr.pval,~,stat] = ttest2(Sup_TTPs(:),Infr_TTPs(:));
LayerStat.Gamma.TTPs.SupVsInfr.Pcorrected = LayerStat.Gamma.TTPs.SupVsInfr.pval*Bonferroni;
LayerStat.Gamma.TTPs.SupVsInfr.tstat = stat.tstat;
LayerStat.Gamma.TTPs.SupVsInfr.df = stat.df;
    % Granular vs Infragranular Time to Peak
[~,LayerStat.Gamma.TTPs.GrVsInfr.pval,~,stat] = ttest2(Gr_TTPs(:),Infr_TTPs(:));
LayerStat.Gamma.TTPs.GrVsInfr.Pcorrected = LayerStat.Gamma.TTPs.GrVsInfr.pval*Bonferroni;
LayerStat.Gamma.TTPs.GrVsInfr.tstat = stat.tstat;
LayerStat.Gamma.TTPs.GrVsInfr.df = stat.df;

% Stats - Layer dependency of FWHM - Gamma-based HRF
Bonferroni = 3;
FWHMs = squeeze(HRFattributes.FWHMs(1,:,:));
NFWHMMat = ones((size(FWHMs,1)-1),1)*FWHMs(1,:);
NFWHMs = FWHMs(2:3,:)./NFWHMMat;
Sup_FWHMs = NFWHMs(:,Layer==1);
Gr_FWHMs = NFWHMs(:,Layer==2);
Infr_FWHMs = NFWHMs(:,Layer==3);
[h1,p1] = adtest(Sup_FWHMs(:));
[h2,p2] = adtest(Gr_FWHMs(:));
[h3,p3] = adtest(Infr_FWHMs(:));
    % Supragranalar vs Granular FWHM
[LayerStat.Gamma.FWHMs.SupVsGr.pval,~,stat] = ranksum(Sup_FWHMs(:),Gr_FWHMs(:),...
    'method','approximate');
LayerStat.Gamma.FWHMs.SupVsGr.Pcorrected = LayerStat.Gamma.FWHMs.SupVsGr.pval*Bonferroni;
LayerStat.Gamma.FWHMs.SupVsGr.Zval = stat.zval;
    % Supragranalar vs Infragranular FWHM
[LayerStat.Gamma.FWHMs.SupVsInfr.pval,~,stat] = ranksum(Sup_FWHMs(:),Infr_FWHMs(:),...
    'method','approximate');
LayerStat.Gamma.FWHMs.SupVsInfr.Pcorrected = LayerStat.Gamma.FWHMs.SupVsInfr.pval*Bonferroni;
LayerStat.Gamma.FWHMs.SupVsInfr.Zval = stat.zval;
    % Granular vs Infragranular FWHM
[LayerStat.Gamma.FWHMs.GrVsInfr.pval,~,stat] = ranksum(Gr_FWHMs(:),Infr_FWHMs(:),...
    'method','approximate');
LayerStat.Gamma.FWHMs.GrVsInfr.Pcorrected = LayerStat.Gamma.FWHMs.GrVsInfr.pval*Bonferroni;
LayerStat.Gamma.FWHMs.GrVsInfr.Zval = stat.zval;

% Stats - Layer dependency of amplitude - MUA-based HRF
Bonferroni = 3;
Amplitudes = squeeze(HRFattributes.Amps(2,:,:));
NAmpMat = ones((size(Amplitudes,1)-1),1)*Amplitudes(1,:);
NAmps = Amplitudes(2:3,:)./NAmpMat;
Sup_Amps = NAmps(:,Layer==1);
Gr_Amps = NAmps(:,Layer==2);
Infr_Amps = NAmps(:,Layer==3);
[h1,p1] = adtest(Sup_Amps(:));
[h2,p2] = adtest(Gr_Amps(:));
[h3,p3] = adtest(Infr_Amps(:));
% Supragranalar vs Granular Amplitude
[LayerStat.MUA.Amps.SupVsGr.pval,~,stat] = ranksum(Sup_Amps(:),Gr_Amps(:),...
    'method','approximate');
LayerStat.MUA.Amps.SupVsGr.Pcorrected = LayerStat.MUA.Amps.SupVsGr.pval*Bonferroni;
LayerStat.MUA.Amps.SupVsGr.Zval = stat.zval;
% Supragranalar vs Infragranular Amplitude
[LayerStat.MUA.Amps.SupVsInfr.pval,~,stat] = ranksum(Sup_Amps(:),Infr_Amps(:),...
    'method','approximate');
LayerStat.MUA.Amps.SupVsInfr.Pcorrected = LayerStat.MUA.Amps.SupVsInfr.pval*Bonferroni;
LayerStat.MUA.Amps.SupVsInfr.Zval = stat.zval;
% Granalar vs Infragranular Amplitude
[LayerStat.MUA.Amps.GrVsInfr.pval,~,stat] = ranksum(Gr_Amps(:),Infr_Amps(:),...
    'method','approximate');
LayerStat.MUA.Amps.GrVsInfr.Pcorrected = LayerStat.MUA.Amps.GrVsInfr.pval*Bonferroni;
LayerStat.MUA.Amps.GrVsInfr.Zval = stat.zval;

% Stats - Layer dependency of Time to peak - MUA-based HRF
Bonferroni = 3;
TTPs = squeeze(HRFattributes.TTPs(2,:,:));
NTTPMat = ones((size(TTPs,1)-1),1)*TTPs(1,:);
NTTPs = TTPs(2:3,:)./NTTPMat;
Sup_TTPs = NTTPs(:,Layer==1);
Gr_TTPs = NTTPs(:,Layer==2);
Infr_TTPs = NTTPs(:,Layer==3);
[h1,p1] = adtest(Sup_TTPs(:));
[h2,p2] = adtest(Gr_TTPs(:));
[h3,p3] = adtest(Infr_TTPs(:));
[p,stat] = vartestn([Sup_TTPs(:); Gr_TTPs(:); Infr_TTPs(:)],...
    [ones(size(Sup_TTPs(:))); 2*ones(size(Gr_TTPs(:))); 3*ones(size(Infr_TTPs(:)))],...
    'Display','off');
    % Supragranalar vs Granular Time to Peak
[~,LayerStat.MUA.TTPs.SupVsGr.pval,~,stat] = ttest2(Sup_TTPs(:),Gr_TTPs(:));
LayerStat.MUA.TTPs.SupVsGr.Pcorrected = LayerStat.MUA.TTPs.SupVsGr.pval*Bonferroni;
LayerStat.MUA.TTPs.SupVsGr.tstat = stat.tstat;
LayerStat.MUA.TTPs.SupVsGr.df = stat.df;
    % Supragranalar vs Infragranular Time to Peak
[~,LayerStat.MUA.TTPs.SupVsInfr.pval,~,stat] = ttest2(Sup_TTPs(:),Infr_TTPs(:));
LayerStat.MUA.TTPs.SupVsInfr.Pcorrected = LayerStat.MUA.TTPs.SupVsInfr.pval*Bonferroni;
LayerStat.MUA.TTPs.SupVsInfr.tstat = stat.tstat;
LayerStat.MUA.TTPs.SupVsInfr.df = stat.df;
    % Granular vs Infragranular Time to Peak
[~,LayerStat.MUA.TTPs.GrVsInfr.pval,~,stat] = ttest2(Gr_TTPs(:),Infr_TTPs(:));
LayerStat.MUA.TTPs.GrVsInfr.Pcorrected = LayerStat.MUA.TTPs.GrVsInfr.pval*Bonferroni;
LayerStat.MUA.TTPs.GrVsInfr.tstat = stat.tstat;
LayerStat.MUA.TTPs.GrVsInfr.df = stat.df;

% Stats - Layer dependency of FWHM - MUA-based HRF
Bonferroni = 3;
FWHMs = squeeze(HRFattributes.FWHMs(2,:,:));
NFWHMMat = ones((size(FWHMs,1)-1),1)*FWHMs(1,:);
NFWHMs = FWHMs(2:3,:)./NFWHMMat;
Sup_FWHMs = NFWHMs(:,Layer==1);
Gr_FWHMs = NFWHMs(:,Layer==2);
Infr_FWHMs = NFWHMs(:,Layer==3);
[h1,p1] = adtest(Sup_FWHMs(:));
[h2,p2] = adtest(Gr_FWHMs(:));
[h3,p3] = adtest(Infr_FWHMs(:));
    % Supragranalar vs Granular FWHM
[LayerStat.MUA.FWHMs.SupVsGr.pval,~,stat] = ranksum(Sup_FWHMs(:),Gr_FWHMs(:),...
    'method','approximate');
LayerStat.MUA.FWHMs.SupVsGr.Pcorrected = LayerStat.MUA.FWHMs.SupVsGr.pval*Bonferroni;
LayerStat.MUA.FWHMs.SupVsGr.Zval = stat.zval;
    % Supragranalar vs Infragranular FWHM
[LayerStat.MUA.FWHMs.SupVsInfr.pval,~,stat] = ranksum(Sup_FWHMs(:),Infr_FWHMs(:),...
    'method','approximate');
LayerStat.MUA.FWHMs.SupVsInfr.Pcorrected = LayerStat.MUA.FWHMs.SupVsInfr.pval*Bonferroni;
LayerStat.MUA.FWHMs.SupVsInfr.Zval = stat.zval;
    % Granular vs Infragranular FWHM
[LayerStat.MUA.FWHMs.GrVsInfr.pval,~,stat] = ranksum(Gr_FWHMs(:),Infr_FWHMs(:),...
    'method','approximate');
LayerStat.MUA.FWHMs.GrVsInfr.Pcorrected = LayerStat.MUA.FWHMs.GrVsInfr.pval*Bonferroni;
LayerStat.MUA.FWHMs.GrVsInfr.Zval = stat.zval;