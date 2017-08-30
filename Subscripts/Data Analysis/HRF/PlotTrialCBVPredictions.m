function [] = PlotTrialCBVPredictions(ProcData,GammaPred,MUAPred,CBVType,Baselines)
%   function [] = PlotTrialCBVPredictions(ProcData,GammaPred,MUAPred,CBVType,Baselines)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Plots the gamma-band and MUA predicted CBV vs the measured
%   CBV
%_______________________________________________________________
%   PARAMETERS:
%                   ProcData - [struct] Processed data from a single trial
%
%                   GammaPred - [array] The gamma-band predicted CBV for
%                   the trial
%
%                   MUAPred - [array] The MUA predicted CBV for the trial
%
%                   CBVType - [string] fieldname of 'ProcData'
%                   corresponding to the CBV ROI.
%
%                   Baseline - [struct] baselines for the animal
%                   corresponding to the ProcData structure.
%_______________________________________________________________
%   RETURN:
%
%_______________________________________________________________

% Process, filter, normalize the ProcData
Fs = ProcData.Fs.([CBVType '_fs']);
timevec = 1/Fs:1/Fs:ProcData.TrialDur;
Bin_wwf = ProcData.Bin_wwf>0.5;
whisks = timevec(Bin_wwf);
eLen = ProcData.TrialDur*Fs;
[z,p,k] = butter(4,[0.01,3]/(Fs/2));
[sos,g] = zp2sos(z,p,k);
nCBV = ProcData.(CBVType)(1:eLen)/mean(Baselines.(CBVType).Means)-1;
fCBV = filtfilt(sos,g,detrend(nCBV));

% Create the plot
figure;
subplot(212); plot(timevec,fCBV-mean(fCBV),'k','Linewidth',1.5);
hold on; 
plot(timevec,MUAPred-mean(MUAPred),'Color',[0 167/255 157/255],'Linewidth',1.5);
plot(timevec,GammaPred-mean(GammaPred),'Color',[0.4, 0.4, 0.4],'Linewidth',1.5);
legend({'Act','Gamma-predicted','MUA-predicted'},'location','southeast','orientation','horizontal');
R2 = CalculateRsquared(GammaPred - mean(GammaPred),...
    fCBV-mean(fCBV));
R22 = CalculateRsquared(MUAPred - mean(MUAPred),...
    fCBV-mean(fCBV));
title([' R2MU = ' num2str(R2) ' R2Gam = ' num2str(R22)]);
ylabel('\DeltaR/R')
xlabel('Time (s)')
hold off;
SR2 = SlidingRSquared(GammaPred-mean(GammaPred),...
    fCBV-mean(fCBV),8*30);
SR22 = SlidingRSquared(MUAPred-mean(MUAPred),...
    fCBV-mean(fCBV),8*30);
subplot(211); 
plot(timevec,max(SR22,-1),'Color',[0 167/255 157/255],'Linewidth',1.5);
hold on;
plot(timevec,max(SR2,-1),'Color',[0.4, 0.4, 0.4],'Linewidth',1.5);
ylim([0 1.1]);
scatter(whisks,ones(size(whisks)),'.','MarkerEdgeColor',[235, 176, 32]/255);
scatter(ProcData.Sol.Contra,1.05*ones(size(ProcData.Sol.Contra)),'v',...
    'MarkerEdgeColor',[1 1 1], 'MarkerFaceColor',...
    [13 115 184]/255);
scatter(ProcData.Sol.Ipsi,1.05*ones(size(ProcData.Sol.Ipsi)),'v',...
    'MarkerEdgeColor',[1 1 1], 'MarkerFaceColor',...
    [0.4 0.4 0.4]);
scatter([ProcData.Sol.Control ProcData.Sol.Tail],...
    1.05*ones(1,length(ProcData.Sol.Control)+length(ProcData.Sol.Tail)),...
    'v','MarkerEdgeColor',[1 1 1], 'MarkerFaceColor',...
    [0.7 0.7 0.7]);
hold off;
ylabel('R^2')
