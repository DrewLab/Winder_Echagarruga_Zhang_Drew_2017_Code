function [] = Supplementary_Figures_Create()
%   [] = Supplementary_Figures_Create()
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Generates the supplementary figures. Requires the that the
%   'Data' folder be in the current directory and that the files
%   'HRFs.mat', 'R2.mat', 'CBVPredictions.mat' have been created (using
%   FiguresCreate.mat')
%
%_______________________________________________________________
%   PARAMETERS:
%
%_______________________________________________________________
%   RETURN:
%
%_______________________________________________________________

%% FIGURE S1h - Comparison of auditory stim and whisking-evoked CBV
clear
display('Creating figure S1h')
animal = 'w323';
behaviors = {'VW','Control'};
dataType = 'CrossCorrROI';
handles.(dataType).alldays = figure;
subplot(212);
EventDataFileName = dir([animal filesep '*EVENTDATA_' dataType '.mat']);
load([animal filesep EventDataFileName.name]);

% Generate the average CBV plot
[~] = PlotEventTriggeredAve(animal,EventData,dataType,behaviors,handles);
xlim([-1 4])
ylabel('\DeltaR/R');
xlabel('Peri-Event Time (s)');
ylim([-0.015 0.015])
axis square;
title('')
legend({'Whisk','Auditory Stim.',},'location',[0.7 0.2 0.2 0.1])

% Get binarized whisking for auditory stimuli
dataType = 'Bin_wwf';
% EventDataFileName = ls([animal filesep '*EVENTDATA_' dataType '.mat']);
% load([animal filesep EventDataFileName]);
EventDataFileName = dir([animal filesep '*EVENTDATA_' dataType '.mat']);
load([animal filesep EventDataFileName.name]);
[DataStruct,FiltArray] = SelectBehavioralEvents(EventData.(dataType),'Control');
BehData = gt(DataStruct.Data(FiltArray,:),0);
timevec = (0:1/DataStruct.Fs:DataStruct.epoch.duration)-...
    DataStruct.epoch.offset;
timevecmat = repmat(timevec,size(BehData,1),1);
ymat = (1:size(BehData,1))'*ones(1,size(BehData,2));
subplot(211)
% subplot('Position',[0.2,0.1,0.4,0.4])
scatter(timevecmat(BehData),ymat(BehData),'k.');
xlim([-1 4])
ylabel(sprintf('Auditory Stimulation\nEvent Number'));
axis square;
title('')
set(gcf,'name','Figure S1h','numbertitle','off')
pause(0.001);

%% Figure S1i - Event Triggered Spatial Map
clear
display('Creating figure S1i')
animal = 'w316';
hem = 'LH';
mean_cmap_bounds = [-0.03 0.03];

% **Spatial Map structures were generated with the following code:
% Requires the raw .bin files containing CBV data, the output of the
% commented code was saved as 'EventTriggeredSpatialMaps.mat'.

% CBVType = 'CrossCorrROI';
% Beh = 'Contra';
% Day = '141104';
% [StimMap] = EventTriggeredCBVSpatialMap(animal,CBVType,Beh,Day,[],[]);
% Beh = 'VW';
% Poly.x = StimMap.xi;
% Poly.y = StimMap.yi;
% [WhiskMap] = EventTriggeredCBVSpatialMap(animal,CBVType,Beh,...
%     StimMap.Date,Poly,StimMap.FileDir);

load([animal filesep 'EventTriggeredSpatialMaps.mat'])
time_intervals = [0.5,1.5; 3,4];
StimIntervalMeans = mean(StimMap.Data,4);
WhiskIntervalMeans = mean(WhiskMap.Data,4);

% Plot
figure;
set(gcf,'name','Figure S1i','numbertitle','off')
subplot(221);
imagesc((StimIntervalMeans(:,:,1)-1).*StimMap.Mask)
imagesc((StimIntervalMeans(:,:,1)-1))
colormap('gray'); axis image;
title([num2str(time_intervals(1,1)) '-' num2str(time_intervals(1,2)) ' seconds'])
xlabel('Caudal');
if strcmp(hem,'LH')
    ylabel(sprintf('Sensory Evoked\nLateral'));
elseif strcmp(hem, 'RH')
    ylabel('Medial')
end
caxis(mean_cmap_bounds);
colbar = colorbar;
ylabel(colbar,'\DeltaR/R')
set(gca,'XTickLabel','','YTickLabel','');
xlim([min(StimMap.xi) max(StimMap.xi)])
ylim([min(StimMap.yi) max(StimMap.yi)])

subplot(222);
imagesc((StimIntervalMeans(:,:,2)-1).*StimMap.Mask)
imagesc((StimIntervalMeans(:,:,2)-1))
colormap('gray'); axis image;
title([num2str(time_intervals(2,1)) '-' num2str(time_intervals(2,2)) ' seconds'])
xlabel('Caudal');
if strcmp(hem,'LH')
    ylabel('Lateral');
elseif strcmp(hem, 'RH')
    ylabel('Medial')
end
caxis(mean_cmap_bounds);
colbar = colorbar;
ylabel(colbar,'\DeltaR/R')
set(gca,'XTickLabel','','YTickLabel','');
xlim([min(StimMap.xi) max(StimMap.xi)])
ylim([min(StimMap.yi) max(StimMap.yi)])

subplot(223);
imagesc((WhiskIntervalMeans(:,:,1)-1).*WhiskMap.Mask)
imagesc((WhiskIntervalMeans(:,:,1)-1))
colormap('gray'); axis image;
title([num2str(time_intervals(1,1)) '-' num2str(time_intervals(1,2)) ' seconds'])
xlabel('Caudal');
if strcmp(hem,'LH')
    ylabel(sprintf('Volitional Whisk\nLateral'));
elseif strcmp(hem, 'RH')
    ylabel('Medial')
end
caxis(mean_cmap_bounds);
colbar = colorbar;
ylabel(colbar,'\DeltaR/R')
set(gca,'XTickLabel','','YTickLabel','');
xlim([min(WhiskMap.xi) max(WhiskMap.xi)])
ylim([min(WhiskMap.yi) max(WhiskMap.yi)])

subplot(224);
imagesc((WhiskIntervalMeans(:,:,2)-1).*WhiskMap.Mask)
colormap('gray'); axis image;
title([num2str(time_intervals(2,1)) '-' num2str(time_intervals(2,2)) ' seconds'])
xlabel('Caudal');
if strcmp(hem,'LH')
    ylabel('Lateral');
elseif strcmp(hem, 'RH')
    ylabel('Medial')
end
caxis(mean_cmap_bounds);
colbar = colorbar;
ylabel(colbar,'\DeltaR/R')
set(gca,'XTickLabel','','YTickLabel','');
xlim([min(WhiskMap.xi) max(WhiskMap.xi)])
ylim([min(WhiskMap.yi) max(WhiskMap.yi)])
pause(0.001);

%% FIGURE S2b - R.M.S of Heart Rate at Rest
clear
display('Creating figure S2b...')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
RestingHeartRateRMS(animals)
set(gcf,'name','Figure S2b','numbertitle','off')
pause(0.001);

%% Figure S2c-e - Stimulus evoked Heart Rate and CBV comparisons
clear
display('Creating figures 2c-e...')
animals = {'w311','w312','w314','w315','w316','w321','w322','w323',...
    'w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
CompareHeartRate(animals,CBVType)
pause(0.001);

%% FIGURE S2f - HR POWER SPECTRUM AT REST
clear
display('Creating figure S2f...')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
RestingHeartRate_LombSpectrum(animals)
set(gcf,'name','Figure S2f','numbertitle','off')
pause(0.001);

%% FIGURE S2g - CROSS CORRELATION BETWEEN CBV AND HR
clear
display('Creating figure S2g...')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
Reshuf95 = NaN*ones(length(animals),2);
for a = 1:length(animals)
    prevdir = cd([animals{a} filesep]);
    [XC,lags,Reshuf95(a,:)] = RestingXCorr_CBVvsHR(CBVType);
    if a==1
        Xcorr_CBV_vs_HR = NaN*ones(length(animals),length(lags));
    end
    Xcorr_CBV_vs_HR(a,:) = XC;
    cd(prevdir)
end

figure; h1 = plot(lags,mean(Xcorr_CBV_vs_HR),'k');
hold on; h2 = plot(lags,mean(Xcorr_CBV_vs_HR)+std(Xcorr_CBV_vs_HR),'k:');
plot(lags,mean(Xcorr_CBV_vs_HR)-std(Xcorr_CBV_vs_HR),'k:');
h3 = plot(lags,(mean(Reshuf95))'*ones(1,size(lags,2)),'r--');
xlabel('Lags (s)')
ylabel(sprintf('Correlation Coefficient\n(\\DeltaR/R vs. heart rate)'))
legend([h1,h2,h3(1)],{'Population Mean','St.Dev','95% c.i. from reshuffled data'},...
    'location','northoutside','orientation','horizontal')
set(gcf,'name','Figure S2g','numbertitle','off')
pause(0.001);

%% FIGURE S2h - RESTING CBV VARIANCE VS HEART RATE
clear
display('Creating figure S2h...')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
[AnalysisStats] = RestingHeartRateVsCBVVariance(animals,CBVType);
Stats.RestingHRVsCBVVariance = AnalysisStats;
set(gcf,'name','Figure S2h','numbertitle','off')
pause(0.001);


%% Figure S3a - Transected Average CBV
clear
display('Generating Figure S3a...')
display('Gathering the stimulus triggered CBV for transected animals...')
animals = {'WL03','WL05','WL08'};
dataType = 'BarrelCBV';
for a = 1:length(animals)
    display([num2str(a) ' of ' num2str(length(animals)) ' animals...'])
    animal = animals{a};
    EventFile = dir([animal filesep '*_EVENTDATA.mat']);
    load([animal filesep EventFile.name])
    Criteria.Fieldname = {'Name'};
    Criteria.Comparison = {'equal'};
    Criteria.Value = {'Contra'};
    DataStruct = EventData.(dataType).stim;
    FiltArray = FilterEvents(DataStruct,Criteria);
    if a == 1
        WTAvgs = NaN*ones(length(animals),size(DataStruct.Data,2));
    end
    WTAvgs(a,:) = mean(DataStruct.Data(FiltArray,:));
    timevec1 = ((1:size(WTAvgs,2))/DataStruct.Fs)-DataStruct.epoch.offset;
end
%Plot
figure;
set(gcf,'name','Figure S3a','numbertitle','off')
subplot(2,2,1)
pretime = timevec1<0;
AvgDC = mean(WTAvgs(:,pretime),2)*ones(1,size(WTAvgs,2));
plot(timevec1,mean(WTAvgs-AvgDC),'r','linewidth',1.5);

display('Gathering stimulus triggered CBV for non-transected animals...')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
[CBV_Contra,~,timevec] = CompileEventTriggeredData(animals,CBVType,'Contra');
preinds = timevec<0;
ContraDC = mean(CBV_Contra(:,preinds),2)*ones(1,size(CBV_Contra,2));
mContra = CBV_Contra-ContraDC;
%Plot
hold on; plot(timevec,mean(mContra),'k','linewidth',1.5);
xlim([-2 4])
xlabel('Peri-stimulus time (s)')
ylabel('\DeltaR/R')
legend({'C.N.7 transected','Not-transected'},'location',[0.6 0.6 0.2 0.3])


% Stats
Fs = DataStruct.Fs;
T_Amps = NaN*ones(1,size(WTAvgs,1));
T_TTPs = NaN*ones(1,size(WTAvgs,1));
T_FWHMs = NaN*ones(1,size(WTAvgs,1));
for WT = 1:size(WTAvgs,1)
    indCBV = sgolayfilt(-1*(WTAvgs(WT,:)-AvgDC(WT,:)),3,Fs+1);
    [amps,ttps,fwhms,~] = findpeaks(indCBV,timevec1,...
        'WidthReference','halfheight');
    [amp,amp_ind] = max(amps);
    T_Amps(WT) = -1*amp;
    T_TTPs(WT) = ttps(amp_ind);
    T_FWHMs(WT) = fwhms(amp_ind);
end

N_Amps = NaN*ones(1,size(mContra,1));
N_TTPs = NaN*ones(1,size(mContra,1));
N_FWHMs = NaN*ones(1,size(mContra,1));
for mC = 1:size(mContra,1)
    indCBV = sgolayfilt(-1*(mContra(mC,:)),3,Fs+1);
    [amps,ttps,fwhms,~] = findpeaks(indCBV,timevec,...
        'WidthReference','halfheight');
    [amp,amp_ind] = max(amps);
    N_Amps(mC) = -1*amp;
    N_TTPs(mC) = ttps(amp_ind);
    N_FWHMs(mC) = fwhms(amp_ind);
end

subplot(2,3,4); scatter(ones(size(T_Amps)),-1*T_Amps,'ro');
hold on; scatter(2*ones(size(N_Amps)),-1*N_Amps,'ko');
hold off;
ylabel('\DeltaR/R_{max}')
xlim([0.5 2.5]);
ylim([0 0.04]);
ax = gca;
set(ax,'XTick',[1,2],'XTickLabel',{'Transect','Normal'});
[~,pAmp,~,stat] = ttest2(T_Amps,N_Amps);
title(['p=' num2str(round(pAmp*100)/100) '; t(' num2str(round(stat.df*10)/10) ')='...
    num2str(round(stat.tstat*100)/100)]);

subplot(2,3,5); scatter(ones(size(T_TTPs)),T_TTPs,'ro');
hold on; scatter(2*ones(size(N_TTPs)),N_TTPs,'ko');
hold off;
ylabel('Time to Peak (s)');
xlim([0.5 2.5]);
ylim([0 2]);
ax = gca;
set(ax,'XTick',[1,2],'XTickLabel',{'Transect','Normal'});
[~,pTTP,~,stat] = ttest2(T_TTPs,N_TTPs);
title(['p=' num2str(round(pTTP*100)/100) '; t(' num2str(round(stat.df*10)/10) ')='...
    num2str(round(stat.tstat*100)/100)]);

subplot(2,3,6); scatter(ones(size(T_FWHMs)),T_FWHMs,'ro');
hold on; scatter(2*ones(size(N_FWHMs)),N_FWHMs,'ko');
hold off;
ylabel('Full Width Half Max (s)')
xlim([0.5 2.5]);
ylim([0 2]);
ax = gca;
set(ax,'XTick',[1,2],'XTickLabel',{'Transect','Normal'});
[~,pFWHM,~,stat] = ttest2(T_FWHMs,N_FWHMs);
title(['p=' num2str(round(pFWHM*100)/100) '; t(' num2str(round(stat.df*10)/10) ')='...
    num2str(round(stat.tstat*100)/100)]);

%% Figure S3b - Transected cross correlations - Gamma power
clear
display('Generating figure S3b...')
display(sprintf('Calculating the cross correlation between Gamma-band power and CBV\nduring non-sensory periods of the trial for transected animals...'))
animals = {'WL03','WL05','WL08'};
CBVType = 'BarrelCBV';
CBVFilterParams.cutoff = 2;
CBVFilterParams.order = 4;
for a = 1:length(animals)
    display([num2str(a) ' of ' num2str(length(animals)) ' animals...'])
    prevdir = cd([animals{a} filesep]);
    ProcFiles = dir('*ProcData.mat');
    Filenames = {ProcFiles(:).name};
    [CC] = TrialCrossCorrelation_GamvsCBV(Filenames, CBVType, CBVFilterParams);
    if a==1
        AllCC_Transected = NaN*ones(length(animals),size(CC.Gampower.vals,2));
    end
    AllCC_Transected(a,:) = mean(CC.Gampower.vals);
    cd(prevdir)
end

display(sprintf('Calculating the cross correlation between Gamma-band power and CBV\nduring non-sensory periods of the trial for non-transected animals...'))
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
CBVFilterParams.cutoff = 2;
CBVFilterParams.order = 4;
for a = 1:length(animals)
    display([num2str(a) ' of ' num2str(length(animals)) ' animals...'])
    prevdir = cd([animals{a} filesep]);
    ProcFiles = dir('*ProcData.mat');
    Filenames = {ProcFiles(:).name};
    [CC] = TrialCrossCorrelation_GamvsCBV(Filenames, CBVType, CBVFilterParams);
    if a==1
        AllCC_Normal = NaN*ones(length(animals),size(CC.Gampower.vals,2));
    end
    AllCC_Normal(a,:) = mean(CC.Gampower.vals);
    cd(prevdir)
end

% Plot
figure;
set(gcf,'name','Figure S3b','numbertitle','off')
subplot(221);
plot(CC.Gampower.Lags,mean(AllCC_Transected),'r','Linewidth',1.5);
hold on; plot(CC.Gampower.Lags,mean(AllCC_Normal),'k','Linewidth',1.5);
ylabel('Correlation coefficient')
xlabel('Lags (s)')
legend({'C.N.7 transected','Not-transected'},'location',[0.6 0.6 0.2 0.3])

% Stats
timevec = CC.Gampower.Lags;
Fs = 30;
N_GamAmps = NaN*ones(1,size(AllCC_Normal,1));
N_GamTTPs = NaN*ones(1,size(AllCC_Normal,1));
N_GamFWHMs = NaN*ones(1,size(AllCC_Normal,1));
for NG = 1:size(AllCC_Normal,1)
    indCC_Gam = sgolayfilt(-1*(AllCC_Normal(NG,:)),3,Fs+1);
    [amps,ttps,fwhms,~] = findpeaks(indCC_Gam,timevec,...
        'WidthReference','halfheight');
    [amp,amp_ind] = max(amps);
    N_GamAmps(NG) = -1*amp;
    N_GamTTPs(NG) = ttps(amp_ind);
    N_GamFWHMs(NG) = fwhms(amp_ind);
end

T_GamAmps = NaN*ones(1,size(AllCC_Transected,1));
T_GamTTPs = NaN*ones(1,size(AllCC_Transected,1));
T_GamFWHMs = NaN*ones(1,size(AllCC_Transected,1));
for TG = 1:size(AllCC_Transected,1)
    indCC_Gam = sgolayfilt(-1*(AllCC_Transected(TG,:)),3,Fs+1);
    [amps,ttps,fwhms,~] = findpeaks(indCC_Gam,timevec,...
        'WidthReference','halfheight');
    [amp,amp_ind] = max(amps);
    T_GamAmps(TG) = -1*amp;
    T_GamTTPs(TG) = ttps(amp_ind);
    T_GamFWHMs(TG) = fwhms(amp_ind);
end

subplot(2,3,4);
scatter(ones(size(T_GamAmps)),-1*T_GamAmps,'ro');
hold on;
scatter(2*ones(size(N_GamAmps)),-1*N_GamAmps,'ko');
hold off;
ylabel('|Corr.Coef|_{max}')
xlim([0.5 2.5])
ylim([0 0.5])
ax = gca;
set(ax,'XTick',[1,2],'XTickLabel',{'Transect','Normal'});
[~,pGamAmp,~,stat] = ttest2(N_GamAmps,T_GamAmps);
title(['p=' num2str(round(pGamAmp*100)/100) '; t(' num2str(round(stat.df*10)/10) ')='...
    num2str(round(stat.tstat*100)/100)]);


subplot(2,3,5);
scatter(ones(size(T_GamTTPs)),T_GamTTPs,'ro'); hold on;
scatter(2*ones(size(N_GamTTPs)),N_GamTTPs,'ko');
hold off;
ylabel('Time to peak (s)')
xlim([0.5 2.5])
ylim([0 2])
ax = gca;
set(ax,'XTick',[1,2],'XTickLabel',{'Transect','Normal'});
[~,pGamTTP,~,stat] = ttest2(N_GamTTPs,T_GamTTPs);
title(['p=' num2str(round(pGamTTP*100)/100) '; t(' num2str(round(stat.df*10)/10) ')='...
    num2str(round(stat.tstat*100)/100)]);

subplot(2,3,6); 
scatter(ones(size(T_GamFWHMs)),T_GamFWHMs,'r');
hold on; 
scatter(2*ones(size(N_GamFWHMs)),N_GamFWHMs,'k');
hold off;
ylabel('Full Width Half Max (s)')
xlim([0.5 2.5])
ylim([0 2.5])
ax = gca;
set(ax,'XTick',[1,2],'XTickLabel',{'Transect','Normal'});
[~,pGamFWHM,~,stat] = ttest2(N_GamFWHMs,T_GamFWHMs);
title(['p=' num2str(round(pGamFWHM*100)/100) '; t(' num2str(round(stat.df*10)/10) ')='...
    num2str(round(stat.tstat*100)/100)]);

%% Figure 2c - Transected cross correlations - MUA
clear
display('Generating figure S3c...')
display(sprintf('Calculating the cross correlation between MUA and CBV\nduring non-sensory periods of the trial for transected animals...'))
animals = {'WL03','WL05','WL08'};
CBVType = 'BarrelCBV';
CBVFilterParams.cutoff = 2;
CBVFilterParams.order = 4;
for a = 1:length(animals)
    display([num2str(a) ' of ' num2str(length(animals)) ' animals...'])
    prevdir = cd([animals{a} filesep]);
    ProcFiles = dir('*ProcData.mat');
    Filenames = {ProcFiles(:).name};
    [CC] = TrialCrossCorrelation_MUAvsCBV(Filenames, CBVType, CBVFilterParams);
    if a==1
        AllCC_Transected = NaN*ones(length(animals),size(CC.MUpower.vals,2));
    end
    AllCC_Transected(a,:) = mean(CC.MUpower.vals);
    cd(prevdir)
end

display(sprintf('Calculating the cross correlation between MUA and CBV\nduring non-sensory periods of the trial for non-transected animals...'))
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
CBVFilterParams.cutoff = 2;
CBVFilterParams.order = 4;
for a = 1:length(animals)
    display([num2str(a) ' of ' num2str(length(animals)) ' animals...'])
    prevdir = cd([animals{a} filesep]);
    ProcFiles = dir('*ProcData.mat');
    Filenames = {ProcFiles(:).name};
    [CC] = TrialCrossCorrelation_MUAvsCBV(Filenames, CBVType, CBVFilterParams);
    if a==1
        AllCC_Normal = NaN*ones(length(animals),size(CC.MUpower.vals,2));
    end
    AllCC_Normal(a,:) = mean(CC.MUpower.vals);
    cd(prevdir)
end

% Plot
figure;
set(gcf,'name','Figure S3c','numbertitle','off')
subplot(221);
plot(CC.MUpower.Lags,mean(AllCC_Transected),'r','Linewidth',1.5);
hold on; plot(CC.MUpower.Lags,mean(AllCC_Normal),'k','Linewidth',1.5);
ylabel('Correlation coefficient')
xlabel('Lags (s)')
legend({'C.N.7 transected','Not-transected'},'location',[0.6 0.6 0.2 0.3])

% Stats
timevec = CC.MUpower.Lags;
Fs = 30;
N_MUAmps = NaN*ones(1,size(AllCC_Normal,1));
N_MUTTPs = NaN*ones(1,size(AllCC_Normal,1));
N_MUFWHMs = NaN*ones(1,size(AllCC_Normal,1));
for NG = 1:size(AllCC_Normal,1)
    indCC_MU = sgolayfilt(-1*(AllCC_Normal(NG,:)),3,Fs+1);
    [amps,ttps,fwhms,~] = findpeaks(indCC_MU,timevec,...
        'WidthReference','halfheight');
    [amp,amp_ind] = max(amps);
    N_MUAmps(NG) = -1*amp;
    N_MUTTPs(NG) = ttps(amp_ind);
    N_MUFWHMs(NG) = fwhms(amp_ind);
end

T_MUAmps = NaN*ones(1,size(AllCC_Transected,1));
T_MUTTPs = NaN*ones(1,size(AllCC_Transected,1));
T_MUFWHMs = NaN*ones(1,size(AllCC_Transected,1));
for TG = 1:size(AllCC_Transected,1)
    indCC_MU = sgolayfilt(-1*(AllCC_Transected(TG,:)),3,Fs+1);
    [amps,ttps,fwhms,~] = findpeaks(indCC_MU,timevec,...
        'WidthReference','halfheight');
    [amp,amp_ind] = max(amps);
    T_MUAmps(TG) = -1*amp;
    T_MUTTPs(TG) = ttps(amp_ind);
    T_MUFWHMs(TG) = fwhms(amp_ind);
end

subplot(2,3,4);
scatter(ones(size(T_MUAmps)),-1*T_MUAmps,'ro');
hold on;
scatter(2*ones(size(N_MUAmps)),-1*N_MUAmps,'ko');
hold off;
ylabel('|Corr.Coef|_{max}')
xlim([0.5 2.5])
ylim([0 0.5])
ax = gca;
set(ax,'XTick',[1,2],'XTickLabel',{'Transect','Normal'});
[~,pMUAmp,~,stat] = ttest2(N_MUAmps,T_MUAmps);
title(['p=' num2str(round(pMUAmp*100)/100) '; t(' num2str(round(stat.df*10)/10) ')='...
    num2str(round(stat.tstat*100)/100)]);


subplot(2,3,5);
scatter(ones(size(T_MUTTPs)),T_MUTTPs,'ro'); hold on;
scatter(2*ones(size(N_MUTTPs)),N_MUTTPs,'ko');
hold off;
ylabel('Time to peak (s)')
xlim([0.5 2.5])
ylim([0 2])
ax = gca;
set(ax,'XTick',[1,2],'XTickLabel',{'Transect','Normal'});
[~,pMUTTP,~,stat] = ttest2(N_MUTTPs,T_MUTTPs);
title(['p=' num2str(round(pMUTTP*100)/100) '; t(' num2str(round(stat.df*10)/10) ')='...
    num2str(round(stat.tstat*100)/100)]);

subplot(2,3,6); 
scatter(ones(size(T_MUFWHMs)),T_MUFWHMs,'r');
hold on; 
scatter(2*ones(size(N_MUFWHMs)),N_MUFWHMs,'k');
hold off;
ylabel('Full Width Half Max (s)')
xlim([0.5 2.5])
ylim([0 2.5])
ax = gca;
set(ax,'XTick',[1,2],'XTickLabel',{'Transect','Normal'});
[~,pMUFWHM,~,stat] = ttest2(N_MUFWHMs,T_MUFWHMs);
title(['p=' num2str(round(pMUFWHM*100)/100) '; t(' num2str(round(stat.df*10)/10) ')='...
    num2str(round(stat.tstat*100)/100)]);

%% Figure S4 - Trial Cross-Correlation
clear
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};

CBVType = 'CrossCorrROI';
CBVFilterParams.cutoff = 2;
CBVFilterParams.order = 4;
display('Generating figure S4...')
display('Calculating the cross correlation between MUA and CBV during non-sensory periods of the trial...')
for a = 1:length(animals)
    display([num2str(a) ' of ' num2str(length(animals)) ' animals...'])
    prevdir = cd([animals{a} filesep]);
    ProcFiles = dir('*ProcData.mat');
    Filenames = {ProcFiles(:).name};
    [CC] = TrialCrossCorrelation_MUAvsCBV(Filenames, CBVType, CBVFilterParams);
    if a==1
        AllCC_MUA = NaN*ones(length(animals),size(CC.MUpower.vals,2));
    end
    AllCC_MUA(a,:) = mean(CC.MUpower.vals);
    cd(prevdir)
end
figure;
set(gcf,'name','Figure S4','numbertitle','off')
subplot(211);
plot(CC.MUpower.Lags,mean(AllCC_MUA),'k','Linewidth',1.5);
hold on;
plot(CC.MUpower.Lags,mean(AllCC_MUA)+std(AllCC_MUA),'k:');
plot(CC.MUpower.Lags,mean(AllCC_MUA)-std(AllCC_MUA),'k:');
hold off;
ylim([-0.5 0.2])
ylabel(sprintf('Correlation Coefficient\n(MUA vs. \\DeltaR/R)'))
axis square;
ax = gca;
set(ax,'XTickLabel',[]);

% This commented code was used to calculate the cross correlation between
% LFP frequency bands and CBV. The analysis takes a long time so the
% results were saved to a structure called '*_CBVXCorr_TrialData.mat'

% % Setup the variables
% CBVType = 'CrossCorrROI';
% SpecParams.fpass = [0.2 150];
% SpecParams.tapers = [5 9];
% SpecParams.movingwin = [1 1/30];
% CBVFilterParams.cutoff = 2;
% CBVFilterParams.order = 4;
% 
% for a = 1:length(animals)
%     prevdir = cd([animals{a} filesep]);
%     ProcFiles = dir('*ProcData.mat');
%     Filenames = {ProcFiles(:).name};
%     [CC] = TrialCrossCorrelation_LFPvsCBV(Filenames,CBVType,SpecParams,...
%         CBVFilterParams);
%     cd(prevdir)
% end

% Load the cross correlations
display('Gathering the cross correlation between LFP frequency bands and CBV...')
for a = 1:length(animals)
    display([num2str(a) ' of ' num2str(length(animals)) ' animals...'])
    prevdir = cd([animals{a} filesep]);
    load([animals{a} '_LH_CBVXCorr_TrialData.mat'])
    if a==1
        AllCC = NaN*ones(length(CC.LFP.Freqs), length(CC.LFP.Lags),...
            length(animals));
    end
    AllCC(:,:,a) = mean(CC.LFP.vals,3);
    cd(prevdir)
end

GammaFreqs = CC.LFP.Freqs>40 & CC.LFP.Freqs<100;
poslags = CC.LFP.Lags>0;
GammaCC = squeeze(mean(AllCC(GammaFreqs,:,:),1))';
MaxCC = min(GammaCC(:,poslags),[],2);

Specax = subplot(212);
imagesc(CC.LFP.Lags, CC.LFP.Freqs, mean(AllCC,3));
axis xy;
caxis([-0.2 0.2])
xlabel('Lags (s)')
ylabel('LFP Frequency (Hz)')
SpecPosition = get(Specax,'Position');
ypos = SpecPosition(2);
yheight = SpecPosition(4);
cbar = colorbar('Position',[0.66, ypos, 0.02,...
    yheight]);
ylabel(cbar,sprintf('Correlation Coefficient\n(LFP vs. \\DeltaR/R)'))
colormap('parula')
axis square;
title(['Max Gamma CC=' num2str(round(mean(MaxCC)*100)/100) '\pm'...
    num2str(round(std(MaxCC)*100)/100)]);



%% FIGURE S5 - COMPARE HRF ATTRIBUTES
clear
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
display('Creating figure S5...')
dataTypes = {'Gamma','MUA'};
Behaviors = {'Contra','VW','Rest'};
load('HRFs.mat');
HRFattributes.Amps = zeros(length(dataTypes),length(Behaviors),...
    length(animals));
HRFattributes.TTPs = zeros(length(dataTypes),length(Behaviors),...
    length(animals));
HRFattributes.FWHMs = zeros(length(dataTypes),length(Behaviors),...
    length(animals));
for d = 1:length(dataTypes)
    dataType = dataTypes{d};
    for b = 1:length(Behaviors)
        beh = Behaviors{b};
        timevec = AllHRFs.(dataType).(beh).Timevec;
        Fs = AllHRFs.(dataType).(beh).Fs;
        for row = 1:size(AllHRFs.(dataType).(beh).HRFs,1)
            % Smooth the HRF more to help with peak Identification   
            indHRF = sgolayfilt(-1*(AllHRFs.(dataType).(beh).HRFs(row,:))...
                ,3,Fs+1);
            [amps,ttps,fwhms,~] = findpeaks(indHRF,timevec,...
                'WidthReference','halfheight');
            [amp,amp_ind] = max(amps);
            if isempty(amps)
                HRFattributes.Amps(d,b,row) = NaN;
                HRFattributes.TTPs(d,b,row) = NaN;
                HRFattributes.FWHMs(d,b,row) = NaN;
            else
                HRFattributes.Amps(d,b,row) = -1*amp;
                HRFattributes.TTPs(d,b,row) = ttps(amp_ind);
                HRFattributes.FWHMs(d,b,row) = fwhms(amp_ind);
            end
        end
    end
end
HRFattributes.animals = animals;
HRFattributes.dataTypes = dataTypes;
HRFattributes.behaviors = Behaviors;

% Stats - Gamma HRF Amplitude
[Stat] = GammaHRFAttributeStats(HRFattributes);
Stats.HRFAttributes.Gamma = Stat.Gamma;

% Stats - MUA HRF Amplitude
[Stat] = MUAHRFAttributeStats(HRFattributes);
Stats.HRFAttributes.MUA = Stat.MUA;

% Plot a comparison of the attributes
varfields = {'Contra','VW','Rest'};
constfield = {'Gamma'};
PlotHRFAttributeComparison(HRFattributes,varfields,constfield);
set(gcf,'name','Figure S5a','numbertitle','off')
constfield = {'MUA'};
PlotHRFAttributeComparison(HRFattributes,varfields,constfield);
set(gcf,'name','Figure S5b','numbertitle','off')

% Plot the layer dependency of the attributes and statistics
Layer = [1,1,2,3,1,2,2,1,3,1,1,3]; % 1: Supra, 2: Granular, 3: Infra
varfields = {'Contra','VW','Rest'};
constfield = {'Gamma'};
PlotHRFAttributeComparison_LayerDependence(HRFattributes,...
    varfields,constfield,Layer);
set(gcf,'name','Figure S5c','numbertitle','off')
constfield = {'MUA'};
PlotHRFAttributeComparison_LayerDependence(HRFattributes,...
    varfields,constfield,Layer);
set(gcf,'name','Figure S5d','numbertitle','off')

% Stats - Layer dependency of amplitude - Gamma-based HRF
[Stat] = HRFAttibuteByLayerStats(HRFattributes,Layer);
Stats.HRFAttributesLayerDependency = Stat;
pause(0.001);

%% Figure S6a - Low frequency power
clearvars -except Stats
display('Generating figure S6a...')
display('Plotting event triggered power < 5 Hz')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};

display('Calculating Event Triggered Spectrogram....this will take a few moments.')
for a = 1:length(animals)
    display([num2str(a) ' of ' num2str(length(animals)) ' animals...'])
    prevdir = cd([animals{a} filesep]);
%     EvFile = ls('*EVENTDATA_Gam.mat');
%     load(EvFile)
    EvFile = dir('*EVENTDATA_Gam.mat');
    load(EvFile.name)
    [StimS,~,~] = TriggeredSpecgram(EventData.Gam,'Contra');
    [WhiskS,timevec,Freqs] = TriggeredSpecgram(EventData.Gam,'VW');
    if a==1
        StimTrigSpecgram = NaN*ones(size(StimS,1),size(StimS,2),length(animals));
        WhiskTrigSpecgram = NaN*ones(size(WhiskS,1),size(StimS,2),length(animals));
    end
    StimTrigSpecgram(:,:,a) = StimS;
    WhiskTrigSpecgram(:,:,a) = WhiskS;
    cd(prevdir)
end
figure;
set(gcf,'name','Figure S6a','numbertitle','off')
subplot(211);
imagesc(timevec,Freqs,mean(StimTrigSpecgram,3))
axis xy;
ylim([0 5])
ylabel('Frequency (Hz)')
xlabel('Peri-stimulus time (s)')
caxis([-0.5 1])
h = colorbar;
ylabel(h,'\DeltaP/P')
title('Stimulus Evoked')

subplot(212);
imagesc(timevec,Freqs,mean(WhiskTrigSpecgram,3))
axis xy;
ylim([0 5])
ylabel('Frequency (Hz)')
xlabel('Peri-whisk time (s)')
caxis([-0.5 1])
h = colorbar;
ylabel(h,'\DeltaP/P')
title('Whisking Evoked')
pause(0.001);

%% Figure S6b - HRF for 10-30 Hz LFP
clearvars -except Stats
display('Generating figure S6b...')
load('HRFs.mat')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
Behaviors = {'Contra','VW','Rest'};
HRFLims = [0 5];
BetaHRF = figure;
ColorOrd = ['b','y','c'];
display('Calculating the HRF for neural power [10-30 Hz])')
for B = 1:length(Behaviors)
    display([num2str(B) ' of ' num2str(length(Behaviors)) ' Behaviors...'])
    Beh = Behaviors{B};
    for a = 1:length(animals)
        display(sprintf(['\t' num2str(a) ' of ' num2str(length(animals)) ' animals...']));
        prevdir = cd([animals{a} filesep]);
        [HRFs] = CalculateHRF_Deconvolution('Beta',CBVType,Beh);
        TimeLims = HRFs.timevec>=HRFLims(1) & HRFs.timevec<=HRFLims(2);
        if a==1
            AllHRFs.Beta.(Beh).HRFs = NaN*ones(length(animals),sum(TimeLims));
            AllHRFs.Beta.(Beh).GammaHRFs = NaN*ones(length(animals),sum(TimeLims));
        end
        AllHRFs.Beta.(Beh).HRFs(a,:) = HRFs.HRF(TimeLims);
        AllHRFs.Beta.(Beh).GammaHRFs(a,:) = HRFs.GammaHRF;
        cd(prevdir)
    end
    AllHRFs.Beta.(Beh).Timevec = HRFs.timevec(TimeLims);
    AllHRFs.Beta.(Beh).Fs = HRFs.Fs;
    mHRFs = mean(AllHRFs.Beta.(Beh).HRFs);
    figure(BetaHRF);
    plot(HRFs.timevec(TimeLims),mHRFs,ColorOrd(B));
    hold on;
end
xlabel('HRF Time (s)')
ylabel('HRF Amplitude (A.U.)')
title('HRF from Neural power (10-30 Hz)')
legend({'Sensory Evoked','Whisking','Rest'},'location','northeast')
set(gcf,'name','Figure S6b','numbertitle','off')
pause(0.001);

%% FigureS6c - HRF for 0.1-8 Hz
clearvars -except Stats AllHRFs
display('Generating Figure S6c...')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
Behaviors = {'Contra','VW','Rest'};
HRFLims = [0 5];
SubAlphaHRF = figure;
ColorOrd = ['b','y','c'];
display('Calculating the HRF for neural power [0.1-8 Hz])')
for B = 1:length(Behaviors)
    Beh = Behaviors{B};
    display([num2str(B) ' of ' num2str(length(Behaviors)) ' Behaviors...'])
    for a = 1:length(animals)
        display(sprintf(['\t' num2str(a) ' of ' num2str(length(animals)) ' animals...']));
        prevdir = cd([animals{a} filesep]);
        [HRFs] = CalculateHRF_Deconvolution('SubAlpha',CBVType,Beh);
        TimeLims = HRFs.timevec>=HRFLims(1) & HRFs.timevec<=HRFLims(2);
        if a==1
            AllHRFs.SubAlpha.(Beh).HRFs = NaN*ones(length(animals),sum(TimeLims));
            AllHRFs.SubAlpha.(Beh).GammaHRFs = NaN*ones(length(animals),sum(TimeLims));
        end
        AllHRFs.SubAlpha.(Beh).HRFs(a,:) = HRFs.HRF(TimeLims);
        AllHRFs.SubAlpha.(Beh).GammaHRFs(a,:) = HRFs.GammaHRF;
        cd(prevdir)
    end
    AllHRFs.SubAlpha.(Beh).Timevec = HRFs.timevec(TimeLims);
    AllHRFs.SubAlpha.(Beh).Fs = HRFs.Fs;
    mHRFs = mean(AllHRFs.SubAlpha.(Beh).HRFs);
    figure(SubAlphaHRF);
    plot(HRFs.timevec(TimeLims),mHRFs,ColorOrd(B));
    hold on;
end
xlabel('HRF Time (s)')
ylabel('HRF Amplitude (A.U.)')
title('HRF from Neural power (0.1-8 Hz)')
legend({'Sensory Evoked','Whisking','Rest'},'location','northeast')
set(gcf,'name','Figure S6c','numbertitle','off')
pause(0.001);

%% Figure S6d - Evaluate HRF for 10-30 Hz
clearvars -except Stats AllHRFs
display('Evaluating HRF for 10-30 Hz...')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
Behaviors = {'Contra','Str','VW','Rest'};
CBVPred = struct;
for B = 1:length(Behaviors)
    display([num2str(B) ' of ' num2str(length(Behaviors)) ' Behaviors...'])
    Beh = Behaviors{B};
    % Preallocate
    if B==1 
        AveR.Beta = NaN*ones(length(animals),length(Behaviors));
        IndR.Beta = NaN*ones(length(animals),length(Behaviors));
    end
    for a = 1:length(animals)
        display(sprintf(['\t' num2str(a) ' of ' num2str(length(animals)) ' animals...']));
        prevdir = cd([animals{a} filesep]);
        if strcmp(Beh,'Str')
            HRFs.HRF = AllHRFs.Beta.VW.HRFs(a,:);
            HRFs.timevec = AllHRFs.Beta.VW.Timevec;
        else
            HRFs.HRF = AllHRFs.Beta.(Beh).HRFs(a,:);
            HRFs.timevec = AllHRFs.Beta.(Beh).Timevec;
        end
        
        if not(isfield(CBVPred,animals{a}))
            CBVPred.(animals{a}).Beta = [];
        end
        [AveR.Beta(a,B),IndR.Beta(a,B),CBVPred.(animals{a}).Beta] ...
            = EvaluateCBVPredictionAccuracy('Beta',CBVType,Beh,HRFs,...
            CBVPred.(animals{a}).Beta);
        cd(prevdir)
    end
end

R_inds = repmat(1:length(Behaviors),length(animals),1);
figure; 
set(gcf,'name','Figure S6d','numbertitle','off')
colors = get(gca,'ColorOrder');
scatter(R_inds(:),AveR.Beta(:),'MarkerEdgeColor',colors(1,:));
hold on; scatter(1:length(Behaviors),median(AveR.Beta),...
    'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:));
scatter(R_inds(:)+0.25,IndR.Beta(:),'MarkerEdgeColor',colors(2,:));
scatter((1:length(Behaviors))+0.25,median(IndR.Beta),...
    'MarkerEdgeColor','k','MarkerFaceColor',colors(2,:));
title('Summary of R^2 values, HRF based on neural power [10-30Hz]')
ylabel('R^2')
xlim([0 length(Behaviors)+1]);
ax1 = gca;
set(ax1,'XTick',min(R_inds(:)):1:max(R_inds(:)),...
    'XTickLabel',{'Sens. Ev.', 'Ext. Move', 'Vol. Whisk', 'Rest'});
pause(0.001);

%% Figure S6e - Evaluate HRF for 0.1-8 Hz
clearvars -except Stats AllHRFs
display('Evaluating HRF for 0.1-8 Hz...')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
Behaviors = {'Contra','Str','VW','Rest'};
CBVPred = struct;
for B = 1:length(Behaviors)
    display([num2str(B) ' of ' num2str(length(Behaviors)) ' Behaviors...'])
    Beh = Behaviors{B};
    % Preallocate
    if B==1 
        AveR.SubAlpha = NaN*ones(length(animals),length(Behaviors));
        IndR.SubAlpha = NaN*ones(length(animals),length(Behaviors));
    end
    for a = 1:length(animals)
        display(sprintf(['\t' num2str(a) ' of ' num2str(length(animals)) ' animals...']));
        prevdir = cd([animals{a} filesep]);
        if strcmp(Beh,'Str')
            HRFs.HRF = AllHRFs.SubAlpha.VW.HRFs(a,:);
            HRFs.timevec = AllHRFs.SubAlpha.VW.Timevec;
        else
            HRFs.HRF = AllHRFs.SubAlpha.(Beh).HRFs(a,:);
            HRFs.timevec = AllHRFs.SubAlpha.(Beh).Timevec;
        end
        
        if not(isfield(CBVPred,animals{a}))
            CBVPred.(animals{a}).SubAlpha = [];
        end
        [AveR.SubAlpha(a,B),IndR.SubAlpha(a,B),CBVPred.(animals{a}).SubAlpha] ...
            = EvaluateCBVPredictionAccuracy('SubAlpha',CBVType,Beh,HRFs,...
            CBVPred.(animals{a}).SubAlpha);
        cd(prevdir)
    end
end

R_inds = repmat(1:length(Behaviors),length(animals),1);
figure; 
set(gcf,'name','Figure S6e','numbertitle','off')
colors = get(gca,'ColorOrder');
scatter(R_inds(:),AveR.SubAlpha(:),'MarkerEdgeColor',colors(1,:));
hold on; scatter(1:length(Behaviors),median(AveR.SubAlpha),...
    'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:));
scatter(R_inds(:)+0.25,IndR.SubAlpha(:),'MarkerEdgeColor',colors(2,:));
scatter((1:length(Behaviors))+0.25,median(IndR.SubAlpha),...
    'MarkerEdgeColor','k','MarkerFaceColor',colors(2,:));
title('Summary of R^2 values, HRF based on neural power [0.1-8Hz]')
ylabel('R^2')
xlim([0 length(Behaviors)+1]);
ax1 = gca;
set(ax1,'XTick',min(R_inds(:)):1:max(R_inds(:)),...
    'XTickLabel',{'Sens. Ev.', 'Ext. Move', 'Vol. Whisk', 'Rest'});
pause(0.001);

%% Figure S7a Correlation Coefficient - Gamma band power
clearvars -except Stats AllHRFs
display('Creating figure s7a...')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
Behaviors = {'Contra','Str','VW','Rest'};
CBVPred = struct;
for B = 1:length(Behaviors)
    Beh = Behaviors{B};
    % Preallocate
    if B==1 
        AveR.Gamma = NaN*ones(length(animals),length(Behaviors));
        IndR.Gamma = NaN*ones(length(animals),length(Behaviors));
    end
    for a = 1:length(animals)
        prevdir = cd([animals{a} filesep]);
        if strcmp(Beh,'Str')
            HRFs.HRF = AllHRFs.Gamma.VW.GammaHRFs(a,:);
            HRFs.timevec = AllHRFs.Gamma.VW.Timevec;
        else
            HRFs.HRF = AllHRFs.Gamma.(Beh).GammaHRFs(a,:);
            HRFs.timevec = AllHRFs.Gamma.(Beh).Timevec;
        end
        
        if not(isfield(CBVPred,animals{a}))
            CBVPred.(animals{a}).Gamma = [];
        end
        [AveR.Gamma(a,B),IndR.Gamma(a,B)] ...
            = EvaluateCBVPredictionAccuracy_Supplement('Gam',CBVType,Beh,HRFs,...
            CBVPred.(animals{a}).Gamma);
        cd(prevdir)
    end
end

R_inds = repmat(1:length(Behaviors),length(animals),1);
figure; 
set(gcf,'name','Figure S7a','numbertitle','off')
colors = get(gca,'ColorOrder');
scatter(R_inds(:),AveR.Gamma(:),'MarkerEdgeColor',colors(1,:));
hold on; scatter(1:length(Behaviors),median(AveR.Gamma),...
    'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:));
scatter(R_inds(:)+0.25,IndR.Gamma(:),'MarkerEdgeColor',colors(2,:));
scatter((1:length(Behaviors))+0.25,median(IndR.Gamma),...
    'MarkerEdgeColor','k','MarkerFaceColor',colors(2,:));
title('Summary of R values for average and individual events')
ylabel('R')
xlim([0 length(Behaviors)+1]);
ax1 = gca;
set(ax1,'XTick',min(R_inds(:)):1:max(R_inds(:)),...
    'XTickLabel',{'Sens. Ev.', 'Ext. Move', 'Vol. Whisk', 'Rest'});
% Stats
[RStats] = GammaHRFPredictionStats_CorrelationCoefficient(AveR,IndR);
Stats.RCompare.GammaHRF = RStats;
pause(0.001);

%% Figure S7b Correlation Coefficient - MUA
clearvars -except Stats AllHRFs
display('Creating figure s7b...')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
Behaviors = {'Contra','Str','VW','Rest'};
CBVPred = struct;
for B = 1:length(Behaviors)
    Beh = Behaviors{B};
    % Preallocate
    if B==1 
        AveR.MUA = NaN*ones(length(animals),length(Behaviors));
        IndR.MUA = NaN*ones(length(animals),length(Behaviors));
    end
    for a = 1:length(animals)
        prevdir = cd([animals{a} filesep]);
        if strcmp(Beh,'Str')
            HRFs.HRF = AllHRFs.MUA.VW.GammaHRFs(a,:);
            HRFs.timevec = AllHRFs.MUA.VW.Timevec;
        else
            HRFs.HRF = AllHRFs.MUA.(Beh).GammaHRFs(a,:);
            HRFs.timevec = AllHRFs.MUA.(Beh).Timevec;
        end
        
        if not(isfield(CBVPred,animals{a}))
            CBVPred.(animals{a}).MUA = [];
        end
        [AveR.MUA(a,B),IndR.MUA(a,B)] ...
            = EvaluateCBVPredictionAccuracy_Supplement('MUpower',CBVType,Beh,HRFs,...
            CBVPred.(animals{a}).MUA);
        cd(prevdir)
    end
end

R_inds = repmat(1:length(Behaviors),length(animals),1);
figure;
set(gcf,'name','Figure S7b','numbertitle','off')
colors = get(gca,'ColorOrder');
scatter(R_inds(:),AveR.MUA(:),'MarkerEdgeColor',colors(1,:));
hold on; scatter(1:length(Behaviors),median(AveR.MUA),...
    'MarkerEdgeColor','k','MarkerFaceColor',colors(1,:));
scatter(R_inds(:)+0.25,IndR.MUA(:),'MarkerEdgeColor',colors(2,:));
scatter((1:length(Behaviors))+0.25,median(IndR.MUA),...
    'MarkerEdgeColor','k','MarkerFaceColor',colors(2,:));
title('Summary of R values for average and individual events')
ylabel('R')
xlim([0 length(Behaviors)+1]);
ax1 = gca;
set(ax1,'XTick',min(R_inds(:)):1:max(R_inds(:)),...
    'XTickLabel',{'Sens. Ev.', 'Ext. Move', 'Vol. Whisk', 'Rest'});

[RStats] = MUAHRFPredictionStats_CorrelationCoefficient(AveR,IndR);
Stats.RCompare.MUAHRF = RStats;
pause(0.001);

%% FIGURE S7c - POWER SPECTRUM OF THE ACTUAL AND PREDICTED CBV
clearvars -except Stats AllHRFs
display('Creating figure S7c...')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
animal = 'w322';
NeurType = 'Gamma';
CBVType = 'CrossCorrROI';
BaseFile = dir([animal filesep '*Baselines.mat']);
load([animal filesep BaseFile.name]);
AnimalRow = strcmp(animals,animal);
HRF = AllHRFs.(NeurType).Contra.HRFs(AnimalRow,:);
FileStruct = dir([animal filesep '*_ProcData.mat']);
filenames = {FileStruct(:).name};

for f = 1:length(filenames)
    load([animal filesep filenames{f}])
    [~,~,FileDate,~] = GetFileInfo(filenames{f});
    strdate = ConvertDate(FileDate);
    CurrentBaseline.Gam = Baselines.Gam.(strdate);
    CurrentBaseline.(CBVType) = Baselines.(CBVType).(strdate);
    [GammaPred,~] = CalculatePredictedCBV(ProcData,'Gam',CBVType,HRF,CurrentBaseline);
    
    SpecData = detrend(GammaPred);
    
    params.Fs = ProcData.Fs.([CBVType '_fs']);
    params.tapers = [5 9];
    params.fpass = [0.025 3];
    [S,~] = mtspectrumc(SpecData,params);
    if f == 1
        S_Pred = NaN*ones(length(filenames),length(S));
    end
    S_Pred(f,:) = S;
    
    CBV_baseline = mean(CurrentBaseline.(CBVType).Means);
    normCBV = detrend(ProcData.(CBVType)/CBV_baseline-1);
    SpecData = detrend(normCBV);
    [S,fr] = mtspectrumc(SpecData,params);
    if f == 1
        S_Act = NaN*ones(length(filenames),length(S));
    end
    S_Act(f,:) = S;
end
figure; 
set(gcf,'name','Figure S7c','numbertitle','off')
loglog(fr,mean(S_Act),'k','linewidth',2);
hold on; loglog(fr,mean(S_Pred),'b','linewidth',2);
xlim([0.01 4])
xlabel('Frequency (Hz)')
ylabel('Power (s^{-1})');
legend({'Measured CBV','Predicted CBV'},'location','southwest')
pause(0.001);

%% FIGURE S7d - Power spectrum of the Gamma-based HRF residual
clearvars -except Stats AllHRFs
display('Creating figure S7d...')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
HRFNeuroType = 'Gam';
CBVType = 'CrossCorrROI';
load('CBVPredictions.mat')
StrtIndx = 2; % HRF predictions occurred on even numbered data
IndxIncr = 2;
ResidPow = cell(length(animals),1);
Freqs = cell(length(animals),1);
display('Calculating the residual power spectrum for periods of rest...')
for a = 1:length(animals)
    display(sprintf(['\t' num2str(a) ' of ' num2str(length(animals)) ' animals...']));
    RestPred = CBVPred.(animals{a}).(HRFNeuroType).Rest;
    RestFile = dir([animals{a} filesep '*RESTDATA_' CBVType '.mat']);
    load([animals{a} filesep RestFile.name]);
    
    % Get the measuredCBV events
    [DataStruct,FiltArray] = SelectBehavioralEvents(RestData.(CBVType),'Rest');
    NormRest = DataStruct.NormData(FiltArray,:);
    RestAct = NormRest(StrtIndx:IndxIncr:length(NormRest));
    
    [ResidPow{a},Freqs{a}] = RestingResidualPowerSpectrum(RestAct,RestPred,...
        DataStruct.Fs);
end
AllPow = cell2mat(ResidPow');
mPow = mean(AllPow,2);
StPow = std(AllPow,[],2);

% Power Spectrum during rest after muscimol infusion
clearvars -except Stats AllHRFs CBVPred Freqs mPow StPow
animals = {'CAN08','CAN10','CAN11','CAN12','CAN13','CAN15','CAN16','CAN18','CAN24'};
RestCriteria.Fieldname = {'Duration','PuffDistance'};
RestCriteria.Comparison = {'gt','gt'};
RestCriteria.Value = {14,5};
RestBuffer = 4;
params.tapers = [1 1];
params.fpass = [0 15];

for a = 1:length(animals)
    animal = animals{a};
    % Load the data
    prevdir = cd([animal filesep 'Muscimol' filesep]);
    RestFile = dir('*_RESTDATA_RadiusROI.mat');
    load(RestFile.name);
    DataField = fieldnames(RestData);
    DataField = DataField{1};
    AllRestData = RestData;
    ThreshFile = dir('*Thresholds.mat');
    load(ThreshFile.name)
    % Throw out periods of rest according to the RestCriteria
    [RestFiltArray] = FilterEvents(AllRestData.(DataField),RestCriteria);
    MinSinceInf = AllRestData.(DataField).MinutesSinceInfusion(RestFiltArray);
    FileDates = AllRestData.(DataField).FileDate(RestFiltArray);
    UniqueDates = unique(FileDates);
    TimeFilt = false(size(FileDates));
    for UD = 1:length(UniqueDates)
        strdate = ConvertDate(UniqueDates{UD});
        if isempty(Thresholds.ReducedMUpower_Start.(strdate))
            continue;
        end
        TimeThresh = ElapsedTime(Thresholds.ReducedMUpower_Start.(strdate),...
            0,Thresholds.Infusion_Start_Times.(strdate));
        TimeFilt = TimeFilt|MinSinceInf>TimeThresh;
    end
    RestData = AllRestData.(DataField).Data(RestFiltArray&TimeThresh);
    Fs = AllRestData.(DataField).Fs;
    % Concatenate the resting events
    Rest = cell(1,length(RestData));
    for RD = 1:length(RestData)
        RestBuffer_Ind = RestBuffer*Fs;
        clippedRest = RestData{RD}(RestBuffer_Ind:end);
        % Find any NaN in the HR data
        RestOffset = mean(clippedRest);
        Rest{RD} = [detrend(clippedRest/mean(clippedRest)-1) NaN];
    end
    JoinedRest = [Rest{:}];
    starts = [1 find(isnan(JoinedRest(1:end-1)))+1];
    stops = find(isnan(JoinedRest(1:end)))-1;
    params.Fs = Fs;
    params.err = [1 0.05];
    [S,f,Serr] = mtspectrumc_unequal_length_trials(JoinedRest',[5,1], params, [starts' stops']);
    if a == 1
        All_S = zeros(length(animals),length(S));
    end
    All_S(a,:) = S;
    cd(prevdir)
end

figure;
set(gcf,'name','Figure S7d','numbertitle','off')
subplot(211);
semilogx(Freqs{1},[mPow, mPow+StPow, mPow-StPow],'k');
xlim([0.1 2])
xlabel('Frequency (Hz)')
ylabel('Power (s^{-1}')
title('Gamma-based HRF prediction residual power spectrum')
axis square;

subplot(212); semilogx(f,[mean(All_S)', (mean(All_S)+std(All_S))',...
    (mean(All_S)-std(All_S))'],'c');
xlim([0.1 2])
xlabel('Frequency (Hz)');
ylabel('Power s^-^1')
title('Resting power spectrum - Muscimol Infusion')
axis square;
pause(0.001);

%% FIGURE S7e - Gamma-band prediction of global subtraction data
% Calculate the HRF
clearvars -except Stats AllHRFs CBVPred
display('Creating figure 6e...')
display('Calculating the Gamma-based HRF for the global subtracted CBV signal...')
load('R2.mat');
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'GlobalSubtract';
Behaviors = {'Contra','VW','Rest'};
HRFLims = [0 5];
for B = 1:length(Behaviors)
    display([num2str(B) ' of ' num2str(length(Behaviors)) ' behaviors...']);
    Beh = Behaviors{B};
    for a = 1:length(animals)
        display(sprintf(['\t' num2str(a) ' of ' num2str(length(animals)) ' animals...']))
        prevdir = cd([animals{a} filesep]);
        [HRFs] = CalculateHRF_Deconvolution('Gam',CBVType,Beh);
        TimeLims = HRFs.timevec>=HRFLims(1) & HRFs.timevec<=HRFLims(2);
        if a==1
            AllHRFs.Gamma_GlobalSubtract.(Beh).HRFs = NaN*ones(length(animals),sum(TimeLims));
            AllHRFs.Gamma_GlobalSubtract.(Beh).GammaHRFs = NaN*ones(length(animals),sum(TimeLims));
        end
        AllHRFs.Gamma_GlobalSubtract.(Beh).HRFs(a,:) = HRFs.HRF(TimeLims);
        AllHRFs.Gamma_GlobalSubtract.(Beh).GammaHRFs(a,:) = HRFs.GammaHRF;
        cd(prevdir)
    end
    AllHRFs.Gamma_GlobalSubtract.(Beh).Timevec = HRFs.timevec(TimeLims);
    AllHRFs.Gamma_GlobalSubtract.(Beh).Fs = HRFs.Fs;
end

% Predict the Global Subtracted CBV
Behaviors = {'Contra','Str','VW','Rest'};
CBVPred = struct;
display('Predicting the global subtracted CBV from the Gamma-based HRF...')
for B = 1:length(Behaviors)
    display([num2str(B) ' of ' num2str(length(Behaviors)) ' behaviors...']);
    Beh = Behaviors{B};
    % Preallocate
    if B==1 
        AveR2.Gamma_GlobalSubtract = NaN*ones(length(animals),length(Behaviors));
        IndR2.Gamma_GlobalSubtract = NaN*ones(length(animals),length(Behaviors));
    end
    for a = 1:length(animals)
        display(sprintf(['\t' num2str(a) ' of ' num2str(length(animals)) ' animals...']))
        prevdir = cd([animals{a} filesep]);
        if strcmp(Beh,'Str')
            HRFs.HRF = AllHRFs.Gamma_GlobalSubtract.VW.GammaHRFs(a,:);
            HRFs.timevec = AllHRFs.Gamma_GlobalSubtract.Contra.Timevec;
        else
            HRFs.HRF = AllHRFs.Gamma_GlobalSubtract.(Beh).GammaHRFs(a,:);
            HRFs.timevec = AllHRFs.Gamma_GlobalSubtract.(Beh).Timevec;
        end
        
        if not(isfield(CBVPred,animals{a}))
            CBVPred.(animals{a}).Gamma_GlobalSubtract = [];
        end
        [AveR2.Gamma_GlobalSubtract(a,B),IndR2.Gamma_GlobalSubtract(a,B),~] ...
            = EvaluateCBVPredictionAccuracy('Gam',CBVType,Beh,HRFs,...
            CBVPred.(animals{a}).Gamma_GlobalSubtract);
        cd(prevdir)
    end
end

% MUA prediction of global subtraction data
clearvars -except Stats AllHRFs CBVPred IndR2 AveR2
display('Calculating the MUA-based HRF for the global subtracted CBV signal...')
% Calculate the HRF
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'GlobalSubtract';
Behaviors = {'Contra','VW','Rest'};
HRFLims = [0 5];
for B = 1:length(Behaviors)
    display([num2str(B) ' of ' num2str(length(Behaviors)) ' behavior...']);
    Beh = Behaviors{B};
    for a = 1:length(animals)
        display(sprintf(['\t' num2str(a) ' of ' num2str(length(animals)) ' animals...']))
        prevdir = cd([animals{a} filesep]);
        [HRFs] = CalculateHRF_Deconvolution('MUpower',CBVType,Beh);
        TimeLims = HRFs.timevec>=HRFLims(1) & HRFs.timevec<=HRFLims(2);
        if a==1
            AllHRFs.MUA_GlobalSubtract.(Beh).HRFs = NaN*ones(length(animals),sum(TimeLims));
            AllHRFs.MUA_GlobalSubtract.(Beh).GammaHRFs = NaN*ones(length(animals),sum(TimeLims));
        end
        AllHRFs.MUA_GlobalSubtract.(Beh).HRFs(a,:) = HRFs.HRF(TimeLims);
        AllHRFs.MUA_GlobalSubtract.(Beh).GammaHRFs(a,:) = HRFs.GammaHRF;
        cd(prevdir)
    end
    AllHRFs.MUA_GlobalSubtract.(Beh).Timevec = HRFs.timevec(TimeLims);
    AllHRFs.MUA_GlobalSubtract.(Beh).Fs = HRFs.Fs;
end

% Predict the Global Subtracted CBV
display('Predicting the global subtracted CBV from the MUA-based HRF...')
Behaviors = {'Contra','Str','VW','Rest'};
CBVPred = struct;
for B = 1:length(Behaviors)
    display([num2str(B) ' of ' num2str(length(Behaviors)) ' behaviors...']);
    Beh = Behaviors{B};
    % Preallocate
    if B==1 
        AveR2.MUA_GlobalSubtract = NaN*ones(length(animals),length(Behaviors));
        IndR2.MUA_GlobalSubtract = NaN*ones(length(animals),length(Behaviors));
    end
    for a = 1:length(animals)
        display(sprintf(['\t' num2str(a) ' of ' num2str(length(animals)) ' animals...']))
        prevdir = cd([animals{a} filesep]);
        if strcmp(Beh,'Str')
            HRFs.HRF = AllHRFs.MUA_GlobalSubtract.VW.GammaHRFs(a,:);
            HRFs.timevec = AllHRFs.MUA_GlobalSubtract.Contra.Timevec;
        else
            HRFs.HRF = AllHRFs.MUA_GlobalSubtract.(Beh).GammaHRFs(a,:);
            HRFs.timevec = AllHRFs.MUA_GlobalSubtract.(Beh).Timevec;
        end
        
        if not(isfield(CBVPred,animals{a}))
            CBVPred.(animals{a}).MUA_GlobalSubtract = [];
        end
        [AveR2.MUA_GlobalSubtract(a,B),IndR2.MUA_GlobalSubtract(a,B),~] ...
            = EvaluateCBVPredictionAccuracy('MUpower',CBVType,Beh,HRFs,...
            CBVPred.(animals{a}).MUA_GlobalSubtract);
        cd(prevdir)
    end
end

% Plot Figure
figure;
set(gcf,'name','Figure S7e','numbertitle','off')
subplot(121);
scatter(ones(size(IndR2.Gamma(:,4))),IndR2.Gamma(:,4),'c','MarkerFaceColor','c');
hold on;
scatter(2*ones(size(IndR2.Gamma_GlobalSubtract(:,4))),...
    IndR2.Gamma_GlobalSubtract(:,4),'k','MarkerFaceColor','c');
scatter([1,2],[mean(IndR2.Gamma(:,4)),mean(IndR2.Gamma_GlobalSubtract(:,4))],...
    'k','MarkerFaceColor','k');
ax = gca;
ax.XTick = 1:2;
ax.XTickLabel = {'Normal','Global Subtraction'};
xlim([0.5 2.5])
ylim([-0.4 1])
ylabel('R^2')
title(sprintf('HRF predictions based on\nGamma-band power'))

subplot(122);
scatter(ones(size(IndR2.MUA(:,4))),IndR2.MUA(:,4),'c','MarkerFaceColor','c');
hold on;
scatter(2*ones(size(IndR2.MUA_GlobalSubtract(:,4))),...
    IndR2.MUA_GlobalSubtract(:,4),'k','MarkerFaceColor','c');
scatter([1,2],[mean(IndR2.MUA(:,4)),mean(IndR2.MUA_GlobalSubtract(:,4))],...
    'k','MarkerFaceColor','k');
ax = gca;
ax.XTick = 1:2;
ax.XTickLabel = {'Normal','Global Subtraction'};
xlim([0.5 2.5])
ylim([-0.4 1])
ylabel('R^2')
title('HRF predictions based on MUA')

% STATISTICS ON THE DIFFERENCE IN R^2 DUE TO GLOBAL SUBTRACTION
[h1,p1] = adtest(IndR2.Gamma_GlobalSubtract(:,4)); % Normal
[h2,p2] = adtest(IndR2.MUA_GlobalSubtract(:,4)); % Normal

[~,Stats.GlobalSubtract.Gamma.pval,~,stat] =...
    ttest(IndR2.Gamma_GlobalSubtract(:,4),IndR2.Gamma(:,4));
Stats.GlobalSubtract.Gamma.tstat = stat.tstat;
Stats.GlobalSubtract.Gamma.df = stat.df;
[~,Stats.GlobalSubtract.MUA.pval,~,stat] =...
    ttest(IndR2.MUA_GlobalSubtract(:,4),IndR2.MUA(:,4));
Stats.GlobalSubtract.MUA.tstat = stat.tstat;
Stats.GlobalSubtract.MUA.df = stat.df;
pause(0.001);

%% Figure S7f CBV variance vs prediction correlation - Gamma
clearvars -except Stats AllHRFs
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
NeurType = 'Gam';
HRFType = 'Gamma';
CBVType = 'CrossCorrROI';
HRFBeh = 'Rest';

display('Comparing the CBV variance to the prediction accuracy of the Gamma-based HRF...')
[GamSlopes,GamPvals] = CBVVarianceVsHRFPredictionAccuracy(animals,NeurType,...
    CBVType,HRFBeh,AllHRFs.(HRFType));
GamSlopeFilter = GamPvals<0.05 & GamSlopes > 0;
GamSigPos = GamSlopes(GamSlopeFilter);
GamOther = GamSlopes(not(GamSlopeFilter));
set(gcf,'name','Figure S7f','numbertitle','off')
pause(0.001);

%% Figure S7g CBV variance vs prediction correlation - MUA
clearvars -except Stats AllHRFs GamSigPos GamOther
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
HRFBeh = 'Rest';
NeurType = 'MUpower';
HRFType = 'MUA';
[MUASlopes,MUAPvals] = CBVVarianceVsHRFPredictionAccuracy(animals,NeurType,...
    CBVType,HRFBeh,AllHRFs.(HRFType));
MUASlopeFilter = MUAPvals<0.05 & MUASlopes > 0;
MUASigPos = MUASlopes(MUASlopeFilter);
MUAOther = MUASlopes(not(MUASlopeFilter));
set(gcf,'name','Figure S7g','numbertitle','off')
pause(0.001);

%% Figure S7h CBV variance vs Cross correlation - individual animals.
display('Generating figure S7h...')
% Get the whisk triggered MUA
[MUpower_VW,~,timevec] = CompileEventTriggeredData(animals,'MUpower','VW');

% Plot
figure;
set(gcf,'name','Figure S7h','numbertitle','off')
scatter(ones(size(GamSigPos)),GamSigPos,'ko','MarkerFaceColor','k');
hold on;
scatter(ones(size(GamOther)),GamOther,'MarkerEdgeColor',[0.7 0.7 0.7],...
    'MarkerFaceColor',[0.7 0.7 0.7]);
scatter(2*ones(size(MUASigPos)),MUASigPos,'ko','MarkerFaceColor','k');
scatter(2*ones(size(MUAOther)),MUAOther,'MarkerEdgeColor',[0.7 0.7 0.7],...
    'MarkerFaceColor',[0.7 0.7 0.7]);
set(gca,'Xtick',[1,2],'XTickLabel',{'Gamma-based HRF','MUA-based HRF'});
xlim([0.5 4]);
ylabel('Slope')

axes('Position',[0.6 0.7 0.3 0.2]); box on;
plot(timevec,mean(MUpower_VW(MUASlopeFilter,:)),'k');
hold on;
plot(timevec,mean(MUpower_VW(not(MUASlopeFilter),:)),'Color',[0.7 0.7 0.7]);
xlim([min(timevec), max(timevec)])
ylim([-0.2 0.5])
ylabel('\DeltaP_{MUA} / P_{MUA}');
xlabel('Peri-whisk time (s)');
pause(0.001);

%% Gather the averaged data for infusions
clearvars -except Stats
clc
animals = {'CAN08','CAN10','CAN11','CAN12','CAN13','CAN15','CAN16',...
    'CAN18','CAN24','CAN27','CAN28','CAN31'};
InfusionTypes = {'aCSF','Muscimol','Muscimol + CNQX + AP5'};
dataTypes = {'RadiusROI','MUpower','Gam'};
BehTypes = {'Contra','VW'};
AveStruct = [];
for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    InfusionField = strrep(InfusionType, ' + ','_');
    AveStruct.(InfusionField) = [];
    for dT = 1:length(dataTypes)
        dataType = dataTypes{dT};
        AveStruct.(InfusionField).(dataType) = [];
        AnimalInd = 1;
        for a = 1:length(animals)
            animal = animals{a};
            FileLocation = [animal filesep InfusionType filesep ...
                    animal '_LH_EVENTDATA_' dataType '.mat'];
            if exist(FileLocation,'file')==2
                load(FileLocation)
            elseif exist(FileLocation,'file')==0
                continue
            end
            
            for BT = 1:length(BehTypes)
                BehType = BehTypes{BT};
                if strcmp(BehType,'Contra')
                    fdname = 'stim';
                elseif strcmp(BehType,'VW')
                    fdname = 'whisk';
                end
                if AnimalInd==1
                    AveStruct.(InfusionField).(dataType).(BehType).PostAve = ...
                        NaN*ones(length(animals),...
                        length(EventData.(dataType).(fdname).Averages.TimeVec));
                    AveStruct.(InfusionField).(dataType).(BehType).PreAve = ...
                        NaN*ones(length(animals),...
                        length(EventData.(dataType).(fdname).Averages.TimeVec));
                    AveStruct.(InfusionField).(dataType).(BehType).AnimalID = ...
                        cell(1,length(animals));
                end
                if isfield(EventData.(dataType).(fdname),'Averages')
                    AveStruct.(InfusionField).(dataType).(BehType).PostAve(AnimalInd,:)...
                        = EventData.(dataType).(fdname).Averages.AllSessions.Post;
                    AveStruct.(InfusionField).(dataType).(BehType).PreAve(AnimalInd,:)...
                        = EventData.(dataType).(fdname).Averages.AllSessions.Pre;
                    AveStruct.(InfusionField).(dataType).(BehType).AnimalID{AnimalInd}...
                        = animal;
                    AveStruct.(InfusionField).(dataType).(BehType).timevec = ...
                        EventData.(dataType).(fdname).Averages.TimeVec;
                else
                    continue;
                end
            end
            AnimalInd = min(AnimalInd+1,length(animals));
        end
    end
end

for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    InfusionField = strrep(InfusionType, ' + ','_');
    for dT = 1:length(dataTypes)
        dataType = dataTypes{dT};
        for BT = 1:length(BehTypes)
            BehType = BehTypes{BT};
            NaNFilt = isnan(sum(AveStruct.(InfusionField).(dataType).(BehType).PreAve,2));
            AveStruct.(InfusionField).(dataType).(BehType).PostAve(NaNFilt,:) = [];
            AveStruct.(InfusionField).(dataType).(BehType).PreAve(NaNFilt,:) = [];
            AveStruct.(InfusionField).(dataType).(BehType).AnimalID = ...
                AveStruct.(InfusionField).(dataType).(BehType).AnimalID(not(NaNFilt));
        end
    end
end

%% Figure S7i
display('Generating figure S7i...')
% Muscimol
CBVType = 'RadiusROI';
InfusionType = 'Muscimol';
InfusionType = strrep(InfusionType, ' + ','_');
NeurType = 'MUpower';
BehType = 'Contra';

% Set the time periods of interest
NeurInds = AveStruct.(InfusionType).(NeurType).(BehType).timevec>-0.1 & ...
    AveStruct.(InfusionType).(NeurType).(BehType).timevec<0.5;
CBVInds = AveStruct.(InfusionType).RadiusROI.(BehType).timevec>0 & ...
    AveStruct.(InfusionType).RadiusROI.(BehType).timevec<3;
PreInds = AveStruct.(InfusionType).RadiusROI.(BehType).timevec<0;

% Normalize everything to the PreStim average of PreInfusion
% Puffs
SumPharmNeur = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
SumaCSFNeur = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
SumaCSFCBV = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
SumPharmCBV = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
for a = 1:size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1)
    animal = AveStruct.(InfusionType).(NeurType).(BehType).AnimalID{a};
    
    % Normalize the post infusion Pharm neural activity to
    % preinfusion amplitudes
    PrePharmNeur = mean(AveStruct.(InfusionType).(NeurType).(BehType).PreAve(a,PreInds));
    NormPharmPreNeur = AveStruct.(InfusionType).(NeurType).(BehType).PreAve(a,:)/PrePharmNeur;
    NormPharmPostNeur = AveStruct.(InfusionType).(NeurType).(BehType).PostAve(a,:)/PrePharmNeur;
    SumPharmPreNeur = sum(NormPharmPreNeur(NeurInds));
    SumPharmPostNeur = sum(NormPharmPostNeur(NeurInds));
    SumPharmNeur(a) = SumPharmPostNeur/SumPharmPreNeur;
    
    % Normalize the post infusion aCSF neural activity to
    % preinfusion amplitudes
    AnimalInd = strcmp(AveStruct.aCSF.(NeurType).(BehType).AnimalID,animal);
    PreaCSFNeur = mean(AveStruct.aCSF.(NeurType).(BehType).PreAve(AnimalInd,PreInds));
    NormaCSFPreNeur = AveStruct.aCSF.(NeurType).(BehType).PreAve(AnimalInd,:)/PreaCSFNeur;
    NormaCSFPostNeur = AveStruct.aCSF.(NeurType).(BehType).PostAve(AnimalInd,:)/PreaCSFNeur;
    SumaCSFPreNeur = sum(NormaCSFPreNeur(NeurInds));
    SumaCSFPostNeur = sum(NormaCSFPostNeur(NeurInds));
    SumaCSFNeur(a) = SumaCSFPostNeur/SumaCSFPreNeur;
    
    % Normalize the post infusion Pharm CBV to preinfusion
    % amplitudes
    PrePharmCBV = mean(AveStruct.(InfusionType).(CBVType).(BehType).PreAve(a,PreInds));
    NormPharmPreCBV = AveStruct.(InfusionType).(CBVType).(BehType).PreAve(a,:)/PrePharmCBV-1;
    NormPharmPostCBV = AveStruct.(InfusionType).(CBVType).(BehType).PostAve(a,:)/PrePharmCBV-1;
    SumPharmPreCBV = sum(NormPharmPreCBV(CBVInds));
    SumPharmPostCBV = sum(NormPharmPostCBV(CBVInds)...
        -mean(NormPharmPostCBV(PreInds)));
    SumPharmCBV(a) = SumPharmPostCBV/SumPharmPreCBV;
    
    PreaCSFCBV = mean(AveStruct.aCSF.(CBVType).(BehType).PreAve(a,PreInds));
    NormaCSFPreCBV = AveStruct.aCSF.(CBVType).(BehType).PreAve(a,:)/PreaCSFCBV-1;
    NormaCSFPostCBV = AveStruct.aCSF.(CBVType).(BehType).PostAve(a,:)/PreaCSFCBV-1;
    SumaCSFPreCBV = sum(NormaCSFPreCBV(CBVInds));
    SumaCSFPostCBV = sum(NormaCSFPostCBV(CBVInds)...
        -mean(NormaCSFPostCBV(PreInds)));
    SumaCSFCBV(a) = SumaCSFPostCBV/SumaCSFPreCBV;
    
end
NormNeur = SumPharmNeur./SumaCSFNeur;
NormCBV = SumPharmCBV./SumaCSFCBV;

% Plot Muscimol
figure;
set(gcf,'name','Figure S7i','numbertitle','off')
scatter(NormNeur,NormCBV,'ko','MarkerFaceColor',[0.7 0.7 0.7]);
ylim([0 2])
xlim([0 2])
title(['Pharmacology effects on ' BehType ' compared to aCSF: ' strrep(InfusionType,'_',' ')])
xlabel([NeurType ' ' strrep(InfusionType,'_',' ') '/' NeurType ' aCSF'])
ylabel(['CBV ' strrep(InfusionType,'_',' ') '/CBV aCSF']);
axis square;
hold on;
errorbar(1.95,mean(NormCBV),std(NormCBV),'k');
scatter(1.95,mean(NormCBV),'ko','MarkerFaceColor',[0.7 0.7 0.7]);
herrorbar(mean(NormNeur),1.95,std(NormNeur),std(NormNeur),'ko');
scatter(mean(NormNeur),1.95,'ko','MarkerFaceColor',[0.7 0.7 0.7]);

% Muscimol + CNQX + AP5

CBVType = 'RadiusROI';
InfusionType = 'Muscimol + CNQX + AP5';
InfusionType = strrep(InfusionType, ' + ','_');
NeurType = 'MUpower';
BehType = 'Contra';

% Set the time periods of interest
NeurInds = AveStruct.(InfusionType).(NeurType).(BehType).timevec>-0.1 & ...
    AveStruct.(InfusionType).(NeurType).(BehType).timevec<0.5;
CBVInds = AveStruct.(InfusionType).RadiusROI.(BehType).timevec>0 & ...
    AveStruct.(InfusionType).RadiusROI.(BehType).timevec<3;
PreInds = AveStruct.(InfusionType).RadiusROI.(BehType).timevec<0;

% Normalize everything to the PreStim average of PreInfusion
% Puffs
SumPharmNeur = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
SumaCSFNeur = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
SumaCSFCBV = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
SumPharmCBV = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
for a = 1:size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1)
    animal = AveStruct.(InfusionType).(NeurType).(BehType).AnimalID{a};
    
    % Normalize the post infusion Pharm neural activity to
    % preinfusion amplitudes
    PrePharmNeur = mean(AveStruct.(InfusionType).(NeurType).(BehType).PreAve(a,PreInds));
    NormPharmPreNeur = AveStruct.(InfusionType).(NeurType).(BehType).PreAve(a,:)/PrePharmNeur;
    NormPharmPostNeur = AveStruct.(InfusionType).(NeurType).(BehType).PostAve(a,:)/PrePharmNeur;
    SumPharmPreNeur = sum(NormPharmPreNeur(NeurInds));
    SumPharmPostNeur = sum(NormPharmPostNeur(NeurInds));
    SumPharmNeur(a) = SumPharmPostNeur/SumPharmPreNeur;
    
    % Normalize the post infusion aCSF neural activity to
    % preinfusion amplitudes
    AnimalInd = strcmp(AveStruct.aCSF.(NeurType).(BehType).AnimalID,animal);
    PreaCSFNeur = mean(AveStruct.aCSF.(NeurType).(BehType).PreAve(AnimalInd,PreInds));
    NormaCSFPreNeur = AveStruct.aCSF.(NeurType).(BehType).PreAve(AnimalInd,:)/PreaCSFNeur;
    NormaCSFPostNeur = AveStruct.aCSF.(NeurType).(BehType).PostAve(AnimalInd,:)/PreaCSFNeur;
    SumaCSFPreNeur = sum(NormaCSFPreNeur(NeurInds));
    SumaCSFPostNeur = sum(NormaCSFPostNeur(NeurInds));
    SumaCSFNeur(a) = SumaCSFPostNeur/SumaCSFPreNeur;
    
    % Normalize the post infusion Pharm CBV to preinfusion
    % amplitudes
    PrePharmCBV = mean(AveStruct.(InfusionType).(CBVType).(BehType).PreAve(a,PreInds));
    NormPharmPreCBV = AveStruct.(InfusionType).(CBVType).(BehType).PreAve(a,:)/PrePharmCBV-1;
    NormPharmPostCBV = AveStruct.(InfusionType).(CBVType).(BehType).PostAve(a,:)/PrePharmCBV-1;
    SumPharmPreCBV = sum(NormPharmPreCBV(CBVInds));
    SumPharmPostCBV = sum(NormPharmPostCBV(CBVInds)...
        -mean(NormPharmPostCBV(PreInds)));
    SumPharmCBV(a) = SumPharmPostCBV/SumPharmPreCBV;
    
    PreaCSFCBV = mean(AveStruct.aCSF.(CBVType).(BehType).PreAve(a,PreInds));
    NormaCSFPreCBV = AveStruct.aCSF.(CBVType).(BehType).PreAve(a,:)/PreaCSFCBV-1;
    NormaCSFPostCBV = AveStruct.aCSF.(CBVType).(BehType).PostAve(a,:)/PreaCSFCBV-1;
    SumaCSFPreCBV = sum(NormaCSFPreCBV(CBVInds));
    SumaCSFPostCBV = sum(NormaCSFPostCBV(CBVInds)...
        -mean(NormaCSFPostCBV(PreInds)));
    SumaCSFCBV(a) = SumaCSFPostCBV/SumaCSFPreCBV;
    
end
NormNeur = SumPharmNeur./SumaCSFNeur;
NormCBV = SumPharmCBV./SumaCSFCBV;

% Plot Muscimol
scatter(NormNeur,NormCBV,'ko','MarkerFaceColor',[1 1 1]);
title(['Pharmacology effects on ' BehType ' compared to aCSF: ' strrep(InfusionType,'_',' ')])
xlabel([NeurType ' ' strrep(InfusionType,'_',' ') '/' NeurType ' aCSF'])
ylabel(['CBV ' strrep(InfusionType,'_',' ') '/CBV aCSF']);
axis square;
errorbar(1.9,mean(NormCBV),std(NormCBV),'k');
scatter(1.9,mean(NormCBV),'ko','MarkerFaceColor',[1 1 1]);
herrorbar(mean(NormNeur),1.9,std(NormNeur),std(NormNeur),'ko');
scatter(mean(NormNeur),1.9,'ko','MarkerFaceColor',[1 1 1]);
plot(0:2,0:2,'k--');
ylim([0 2])
xlim([0 2])
hold off;
pause(0.001);

%% Figure S7j 
display('Generating figure S7j...')
% Muscimol
CBVType = 'RadiusROI';
InfusionType = 'Muscimol';
InfusionType = strrep(InfusionType, ' + ','_');
NeurType = 'Gam';
BehType = 'Contra';

% Set the time periods of interest
NeurInds = AveStruct.(InfusionType).(NeurType).(BehType).timevec>-0.1 & ...
    AveStruct.(InfusionType).(NeurType).(BehType).timevec<0.5;
CBVInds = AveStruct.(InfusionType).RadiusROI.(BehType).timevec>0 & ...
    AveStruct.(InfusionType).RadiusROI.(BehType).timevec<3;
PreInds = AveStruct.(InfusionType).RadiusROI.(BehType).timevec<0;

% Normalize everything to the PreStim average of PreInfusion
% Puffs
SumPharmNeur = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
SumaCSFNeur = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
SumaCSFCBV = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
SumPharmCBV = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
for a = 1:size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1)
    animal = AveStruct.(InfusionType).(NeurType).(BehType).AnimalID{a};
    
    % Normalize the post infusion Pharm neural activity to
    % preinfusion amplitudes
    PrePharmNeur = mean(AveStruct.(InfusionType).(NeurType).(BehType).PreAve(a,PreInds));
    NormPharmPreNeur = AveStruct.(InfusionType).(NeurType).(BehType).PreAve(a,:)/PrePharmNeur;
    NormPharmPostNeur = AveStruct.(InfusionType).(NeurType).(BehType).PostAve(a,:)/PrePharmNeur;
    SumPharmPreNeur = sum(NormPharmPreNeur(NeurInds));
    SumPharmPostNeur = sum(NormPharmPostNeur(NeurInds));
    SumPharmNeur(a) = SumPharmPostNeur/SumPharmPreNeur;
    
    % Normalize the post infusion aCSF neural activity to
    % preinfusion amplitudes
    AnimalInd = strcmp(AveStruct.aCSF.(NeurType).(BehType).AnimalID,animal);
    PreaCSFNeur = mean(AveStruct.aCSF.(NeurType).(BehType).PreAve(AnimalInd,PreInds));
    NormaCSFPreNeur = AveStruct.aCSF.(NeurType).(BehType).PreAve(AnimalInd,:)/PreaCSFNeur;
    NormaCSFPostNeur = AveStruct.aCSF.(NeurType).(BehType).PostAve(AnimalInd,:)/PreaCSFNeur;
    SumaCSFPreNeur = sum(NormaCSFPreNeur(NeurInds));
    SumaCSFPostNeur = sum(NormaCSFPostNeur(NeurInds));
    SumaCSFNeur(a) = SumaCSFPostNeur/SumaCSFPreNeur;
    
    % Normalize the post infusion Pharm CBV to preinfusion
    % amplitudes
    PrePharmCBV = mean(AveStruct.(InfusionType).(CBVType).(BehType).PreAve(a,PreInds));
    NormPharmPreCBV = AveStruct.(InfusionType).(CBVType).(BehType).PreAve(a,:)/PrePharmCBV-1;
    NormPharmPostCBV = AveStruct.(InfusionType).(CBVType).(BehType).PostAve(a,:)/PrePharmCBV-1;
    SumPharmPreCBV = sum(NormPharmPreCBV(CBVInds));
    SumPharmPostCBV = sum(NormPharmPostCBV(CBVInds)...
        -mean(NormPharmPostCBV(PreInds)));
    SumPharmCBV(a) = SumPharmPostCBV/SumPharmPreCBV;
    
    PreaCSFCBV = mean(AveStruct.aCSF.(CBVType).(BehType).PreAve(a,PreInds));
    NormaCSFPreCBV = AveStruct.aCSF.(CBVType).(BehType).PreAve(a,:)/PreaCSFCBV-1;
    NormaCSFPostCBV = AveStruct.aCSF.(CBVType).(BehType).PostAve(a,:)/PreaCSFCBV-1;
    SumaCSFPreCBV = sum(NormaCSFPreCBV(CBVInds));
    SumaCSFPostCBV = sum(NormaCSFPostCBV(CBVInds)...
        -mean(NormaCSFPostCBV(PreInds)));
    SumaCSFCBV(a) = SumaCSFPostCBV/SumaCSFPreCBV;
    
end
NormNeur = SumPharmNeur./SumaCSFNeur;
NormCBV = SumPharmCBV./SumaCSFCBV;

% Plot Muscimol
figure;
set(gcf,'name','Figure S7j','numbertitle','off')
scatter(NormNeur,NormCBV,'ko','MarkerFaceColor',[0.7 0.7 0.7]);
ylim([0 2])
xlim([0 2])
title(['Pharmacology effects on ' BehType ' compared to aCSF: ' strrep(InfusionType,'_',' ')])
xlabel([NeurType ' ' strrep(InfusionType,'_',' ') '/' NeurType ' aCSF'])
ylabel(['CBV ' strrep(InfusionType,'_',' ') '/CBV aCSF']);
axis square;
hold on;
errorbar(1.95,mean(NormCBV),std(NormCBV),'k');
scatter(1.95,mean(NormCBV),'ko','MarkerFaceColor',[0.7 0.7 0.7]);
herrorbar(mean(NormNeur),1.95,std(NormNeur),std(NormNeur),'ko');
scatter(mean(NormNeur),1.95,'ko','MarkerFaceColor',[0.7 0.7 0.7]);

% Muscimol + CNQX + AP5

CBVType = 'RadiusROI';
InfusionType = 'Muscimol + CNQX + AP5';
InfusionType = strrep(InfusionType, ' + ','_');
NeurType = 'Gam';
BehType = 'Contra';

% Set the time periods of interest
NeurInds = AveStruct.(InfusionType).(NeurType).(BehType).timevec>-0.1 & ...
    AveStruct.(InfusionType).(NeurType).(BehType).timevec<0.5;
CBVInds = AveStruct.(InfusionType).RadiusROI.(BehType).timevec>0 & ...
    AveStruct.(InfusionType).RadiusROI.(BehType).timevec<3;
PreInds = AveStruct.(InfusionType).RadiusROI.(BehType).timevec<0;

% Normalize everything to the PreStim average of PreInfusion
% Puffs
SumPharmNeur = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
SumaCSFNeur = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
SumaCSFCBV = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
SumPharmCBV = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
for a = 1:size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1)
    animal = AveStruct.(InfusionType).(NeurType).(BehType).AnimalID{a};
    
    % Normalize the post infusion Pharm neural activity to
    % preinfusion amplitudes
    PrePharmNeur = mean(AveStruct.(InfusionType).(NeurType).(BehType).PreAve(a,PreInds));
    NormPharmPreNeur = AveStruct.(InfusionType).(NeurType).(BehType).PreAve(a,:)/PrePharmNeur;
    NormPharmPostNeur = AveStruct.(InfusionType).(NeurType).(BehType).PostAve(a,:)/PrePharmNeur;
    SumPharmPreNeur = sum(NormPharmPreNeur(NeurInds));
    SumPharmPostNeur = sum(NormPharmPostNeur(NeurInds));
    SumPharmNeur(a) = SumPharmPostNeur/SumPharmPreNeur;
    
    % Normalize the post infusion aCSF neural activity to
    % preinfusion amplitudes
    AnimalInd = strcmp(AveStruct.aCSF.(NeurType).(BehType).AnimalID,animal);
    PreaCSFNeur = mean(AveStruct.aCSF.(NeurType).(BehType).PreAve(AnimalInd,PreInds));
    NormaCSFPreNeur = AveStruct.aCSF.(NeurType).(BehType).PreAve(AnimalInd,:)/PreaCSFNeur;
    NormaCSFPostNeur = AveStruct.aCSF.(NeurType).(BehType).PostAve(AnimalInd,:)/PreaCSFNeur;
    SumaCSFPreNeur = sum(NormaCSFPreNeur(NeurInds));
    SumaCSFPostNeur = sum(NormaCSFPostNeur(NeurInds));
    SumaCSFNeur(a) = SumaCSFPostNeur/SumaCSFPreNeur;
    
    % Normalize the post infusion Pharm CBV to preinfusion
    % amplitudes
    PrePharmCBV = mean(AveStruct.(InfusionType).(CBVType).(BehType).PreAve(a,PreInds));
    NormPharmPreCBV = AveStruct.(InfusionType).(CBVType).(BehType).PreAve(a,:)/PrePharmCBV-1;
    NormPharmPostCBV = AveStruct.(InfusionType).(CBVType).(BehType).PostAve(a,:)/PrePharmCBV-1;
    SumPharmPreCBV = sum(NormPharmPreCBV(CBVInds));
    SumPharmPostCBV = sum(NormPharmPostCBV(CBVInds)...
        -mean(NormPharmPostCBV(PreInds)));
    SumPharmCBV(a) = SumPharmPostCBV/SumPharmPreCBV;
    
    PreaCSFCBV = mean(AveStruct.aCSF.(CBVType).(BehType).PreAve(a,PreInds));
    NormaCSFPreCBV = AveStruct.aCSF.(CBVType).(BehType).PreAve(a,:)/PreaCSFCBV-1;
    NormaCSFPostCBV = AveStruct.aCSF.(CBVType).(BehType).PostAve(a,:)/PreaCSFCBV-1;
    SumaCSFPreCBV = sum(NormaCSFPreCBV(CBVInds));
    SumaCSFPostCBV = sum(NormaCSFPostCBV(CBVInds)...
        -mean(NormaCSFPostCBV(PreInds)));
    SumaCSFCBV(a) = SumaCSFPostCBV/SumaCSFPreCBV;
    
end
NormNeur = SumPharmNeur./SumaCSFNeur;
NormCBV = SumPharmCBV./SumaCSFCBV;

% Plot Muscimol
scatter(NormNeur,NormCBV,'ko','MarkerFaceColor',[1 1 1]);
title(['Pharmacology effects on ' BehType ' compared to aCSF: ' strrep(InfusionType,'_',' ')])
xlabel([NeurType ' ' strrep(InfusionType,'_',' ') '/' NeurType ' aCSF'])
ylabel(['CBV ' strrep(InfusionType,'_',' ') '/CBV aCSF']);
axis square;
errorbar(1.9,mean(NormCBV),std(NormCBV),'k');
scatter(1.9,mean(NormCBV),'ko','MarkerFaceColor',[1 1 1]);
herrorbar(mean(NormNeur),1.9,std(NormNeur),std(NormNeur),'ko');
scatter(mean(NormNeur),1.9,'ko','MarkerFaceColor',[1 1 1]);
plot(0:2,0:2,'k--');
ylim([0 2])
xlim([0 2])
hold off;
pause(0.001);

%% Figure S9a
display('Generating Figure S9a...')
CBVType = 'RadiusROI';
animal = 'CAN12';
filename = 'CAN12_LH_160412_11_25_3712_ProcData.mat';
prevdir = cd([animal filesep 'aCSF' filesep]);
figure;
set(gcf,'name','Figure S9a (left)','numbertitle','off')
PlotSingleTrial_Infusion(filename,CBVType)
subplot(411)
title('aCSF Infusion')
cd(prevdir)

filename = 'CAN12_LH_160407_11_32_2407_ProcData.mat';
prevdir = cd([animal filesep 'Muscimol' filesep]);
figure;
set(gcf,'name','Figure S9a (right)','numbertitle','off')
PlotSingleTrial_Infusion(filename,CBVType)
subplot(411)
title('Muscimol Infusion')
cd(prevdir)
pause(0.001);

%% Figure 9b
display('Generating Figure S9b...')
CBVType = 'RadiusROI';
animal = 'CAN28';
filename = 'CAN28_LH_161222_12_01_2422_ProcData.mat';
prevdir = cd([animal filesep 'aCSF' filesep]);
figure;
set(gcf,'name','Figure S9b (left)','numbertitle','off')
PlotSingleTrial_Infusion(filename,CBVType)
subplot(411)
title('aCSF Infusion')
cd(prevdir)

filename = 'CAN28_LH_161220_11_33_5720_ProcData.mat';
prevdir = cd([animal filesep 'Muscimol + CNQX + AP5' filesep]);
figure;
set(gcf,'name','Figure S9b (right)','numbertitle','off')
PlotSingleTrial_Infusion(filename,CBVType)
subplot(411)
title('Muscimol Infusion')
cd(prevdir)
pause(0.001);

%% Figure 9c
display('Generating Figure S9c...')
CBVType = 'RadiusROI';
animal = 'CAN31';
filename = 'CAN31_LH_170201_11_19_3501_ProcData.mat';
prevdir = cd([animal filesep 'aCSF' filesep]);
figure;
set(gcf,'name','Figure S9c (left)','numbertitle','off')
PlotSingleTrial_Infusion(filename,CBVType)
subplot(411)
title('aCSF Infusion')
cd(prevdir)

filename = 'CAN31_LH_170206_17_11_0906_ProcData.mat';
prevdir = cd([animal filesep 'Muscimol + Adrenergic Blockers' filesep]);
figure;
set(gcf,'name','Figure S9c (right)','numbertitle','off')
PlotSingleTrial_Infusion(filename,CBVType)
subplot(411)
title('Muscimol Infusion')
cd(prevdir)
pause(0.001);

%% Figure S10
% The cross correlation data were generated by the script: 
% Infusion_CrossCorrelation_InfusionRadiusvsRestofWindow.m and the results
% were saved to a .mat file.

display('Generating Figure S10...');
clear
animals = {'CAN10','CAN11','CAN12','CAN13','CAN16','CAN18'};
for a = 1:length(animals)
    load([animals{a} filesep 'aCSF' filesep animals{a} '_CBVROIXCorr_Rest.mat'])
    if a==1
        aCSFXCorr = NaN*ones(length(animals),length(CrossROI));
    end
    aCSFXCorr(a,:) = CrossROI;
    
    load([animals{a} filesep 'Muscimol' filesep animals{a} '_CBVROIXCorr_Rest.mat'])
    if a==1
        MuscXCorr = NaN*ones(length(animals),length(CrossROI));
    end
    MuscXCorr(a,:) = CrossROI;
end

figure; 
set(gcf,'name','Figure S10','numbertitle','off')
subplot(121); 
L1 = plot(lags/30,mean(aCSFXCorr),'k');
hold on;
plot(lags/30,mean(aCSFXCorr)+std(aCSFXCorr),'k:');
plot(lags/30,mean(aCSFXCorr)-std(aCSFXCorr),'k:');
L2 = plot(lags/30,mean(MuscXCorr),'c');
plot(lags/30,mean(MuscXCorr)+std(MuscXCorr),'c:');
plot(lags/30,mean(MuscXCorr)-std(MuscXCorr),'c:');
xlim([-4 4])
legend([L1,L2],{'aCSF','Musicmol'},'location','Southeast');

[aCSF_Amps,max_ind_aCSF] = max(aCSFXCorr,[],2);
[Musc_Amps,max_ind_Musc] = max(MuscXCorr,[],2);

subplot(143);
scatter(ones(size(aCSF_Amps)),aCSF_Amps,'ko');
hold on; scatter(2*ones(size(Musc_Amps)), Musc_Amps,'co');
ylim([0,1])
xlim([0.5 2.5]);
ax = gca;
set(ax, 'XTick',[],'XTickLabels',[]);
[~,p_amp,~,stat] = ttest(aCSF_Amps,Musc_Amps);
title(['p=' num2str(round(p_amp*100)/100) '; t(' num2str(stat.df) ')=' ...
    num2str(round(stat.tstat*100)/100)])

subplot(144);
scatter(ones(size(max_ind_aCSF)),lags(max_ind_aCSF)/30,'ko');
hold on;
scatter(2*ones(size(max_ind_Musc)),lags(max_ind_Musc)/30,'co');
xlim([0.5 2.5]);
ax = gca;
set(ax, 'XTick',[],'XTickLabels',[]);
ylim([-1 1]);
[~,p_lag,~,stat] = ttest(lags(max_ind_aCSF)/30,lags(max_ind_Musc)/30);
title(['p=' num2str(round(p_lag*100)/100) '; t(' num2str(stat.df) ')=' ...
    num2str(round(stat.tstat*100)/100)])