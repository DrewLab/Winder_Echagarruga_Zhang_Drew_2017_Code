function [] = FiguresCreate()
%   function [] = FiguresCreate()
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: This script plots the main figures. The data required for
%   this script are in a pre-processed form. All scripts for processing the
%   data are included with these scripts. 
%
%   The raw data were imported into MATLAB using the
%   CreateRawDataStructure_parallel.m and AddROIIntensity2RawdataFile.m
%   scripts.
%
%   Solenoid firing times, neural power, and binarization were performed
%   using ProcessRawDataFile.m
%
%   The processed data were categorized and grouped into behaviors using
%   the SingleAnimalProcessingMaster.m
%
%_______________________________________________________________
%   PARAMETERS:
%
%_______________________________________________________________
%   RETURN:
%
%_______________________________________________________________

%% Figure 1b - Total behavioral times
clear
clc
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
display('Generating figure 1b...')
TotalBehavioralTimes(animals,CBVType)
pause(0.001);

%% Figure 1d - Pixel-wise cross correlation example
clear
clc
animal = 'w316';
display('Generating figure 1d...')
PlotSpatialXCorrExample(animal)
pause(0.001);

%% Figure 1e - Behavioral CBV Variance
clear
clc
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
Behaviors = {'Str','Contra','Ipsi','Control','VW','Rest'};
display('Generating figure 1e...')
[VarianceStats] = CompareBehavioralVariance(animals,Behaviors,CBVType);
Stats.CompareCBVVariance = VarianceStats;
pause(0.001);

%% Figure 1f - Example Trial
clearvars -except Stats
clc
animal = 'w312';
FileID = 'w312_LH_140701_19_54_4001_ProcData.mat';
CBVType = 'CrossCorrROI';
prevdir = cd([animal filesep]);
display('Generating figure 1f...')
PlotSingleTrial(FileID,CBVType)
cd(prevdir)
pause(0.001);

%% Figure 2a - Stimulus triggered averages
clearvars -except Stats
clc
display('Generating figure 2a...')
display('Gather triggered MUA...')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
% Average MUA
yscale = [-1 6];
[MUpower_Contra,MUpower_ContraZ,timevec] = CompileEventTriggeredData(animals,'MUpower','Contra');
preinds = timevec<0;
ContraDC = mean(MUpower_Contra(:,preinds),2)*ones(1,size(MUpower_Contra,2));
mContra = MUpower_Contra-ContraDC;
ContraDCZ = mean(MUpower_ContraZ(:,preinds),2)*ones(1,size(MUpower_ContraZ,2));
mContraZ = MUpower_ContraZ-ContraDCZ;
figure; 
set(gcf,'name','Figure 2a','numbertitle','off')
subplot(311);
[AX,H1,H2] = plotyy([timevec', timevec', timevec'],...
    [mean(mContra)', (mean(mContra)+std(mContra))', (mean(mContra)-std(mContra))'],...
    timevec,mean(mContraZ));
scalefactor = median(mean(mContra)./mean(mContraZ));
set(AX(1),'box','off',...
    'YLim',[yscale(1)*scalefactor, yscale(2)*scalefactor],...
    'YTick',round(yscale(1)*scalefactor):1:round(yscale(2)*scalefactor))
ylabel(AX(1),'\DeltaP/P')
set(AX(2),'YLim',yscale,...
    'YTick',yscale(1):1:yscale(2),...
    'YColor',[0 0 0]);
ylabel(AX(2),'z-score')
set(H2,'visible','off');
set(H1,'Color','k')
set(AX,'xlim',[-4 6],...
    'XTickLabel',[])
axis(AX,'square')

% Average Spectrogram
display(sprintf('Calculating and assembling triggered spectrograms.\nThis will take a few moments...'))
for a = 1:length(animals)
    display([num2str(a) ' of ' num2str(length(animals)) ' animals...'])
    prevdir = cd([animals{a} filesep]);
%     EvFile = ls('*EVENTDATA_Gam.mat');
%     load(EvFile)
    EvFile = dir('*EVENTDATA_Gam.mat');
    load(EvFile.name)
    [S,timevec,Freqs] = TriggeredSpecgram(EventData.Gam,'Contra');
    if a==1
        TrigSpecgram = NaN*ones(size(S,1),size(S,2),length(animals));
    end
    TrigSpecgram(:,:,a) = S;
    cd(prevdir)
end
Specax = subplot(312);
imagesc(timevec,Freqs,mean(TrigSpecgram,3))
set(gca,'XTickLabel',[])
axis xy;
axis square;
SpecPosition = get(Specax,'Position');
ypos = SpecPosition(2);
yheight = SpecPosition(4);
cbar = colorbar('Position',[0.61, ypos, 0.02,...
    yheight]);
colormap('parula')
ylabel(cbar,'\DeltaP/P')
ylabel(sprintf('Frequency\n(Hz)'))
caxis([-0.5 1])

% Average CBV
display('Gathering stimulus triggered CBV...')
[CBV_Contra,CBV_ContraZ,timevec] = CompileEventTriggeredData(animals,CBVType,'Contra');
preinds = timevec<0;
ContraDC = mean(CBV_Contra(:,preinds),2)*ones(1,size(CBV_Contra,2));
mContra = CBV_Contra-ContraDC;
ContraDCZ = mean(CBV_ContraZ(:,preinds),2)*ones(1,size(CBV_ContraZ,2));
mContraZ = CBV_ContraZ-ContraDCZ;
subplot(313);
[AX,H1,H2] = plotyy([timevec', timevec', timevec'],...
    [mean(mContra)', (mean(mContra)+std(mContra))', (mean(mContra)-std(mContra))'],...
    timevec,mean(mContraZ));
scalefactor = median(mean(mContra)./mean(mContraZ));
yscale = [-4,1.5];
set(AX(1),'YLim',[yscale(1)*scalefactor, yscale(2)*scalefactor],...
    'YTick',round(yscale(1)*scalefactor*100)/100:0.01:(yscale(2)*scalefactor*100)/100,...
    'box','off');
ylabel(AX(1),'\DeltaR/R')
set(AX(2),'YLim',yscale,...
    'YTick',yscale(1):1:yscale(2),...
    'YColor','k');
ylabel(AX(2),'z-score')
set(AX,'XLim',[-4 6]);
set(H2,'Visible','off');
set(H1,'Color','k')
axis(AX,'square')
xlabel('Peri-stimulus time (s)')
pause(0.001);

%% Figure 2b - Whisking triggered averages
clearvars -except Stats
clc
display('Generating figure 2b...')
display('Gathering whisk-triggered MUA...')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
% Average MUA
yscale = [-1 6];
[MUpower_VW,MUpower_VWZ,timevec] = CompileEventTriggeredData(animals,'MUpower','VW');
preinds = timevec<0;
ContraDC = mean(MUpower_VW(:,preinds),2)*ones(1,size(MUpower_VW,2));
mVW = MUpower_VW-ContraDC;
ContraDCZ = mean(MUpower_VWZ(:,preinds),2)*ones(1,size(MUpower_VWZ,2));
mVWZ = MUpower_VWZ-ContraDCZ;

figure;
set(gcf,'name','Figure 2b','numbertitle','off')
subplot(311);
[AX,H1,H2] = plotyy([timevec', timevec', timevec'],[mean(mVW)', (mean(mVW)+std(mVW))', ...
    (mean(mVW)-std(mVW))'],timevec,mean(mVWZ));
scalefactor = median(mean(mVW)./mean(mVWZ));
set(AX(1),'box','off',...
    'YLim',[yscale(1)*scalefactor, yscale(2)*scalefactor],...
    'YTick',round(yscale(1)*scalefactor):1:round(yscale(2)*scalefactor))
ylabel(AX(1),'\DeltaP/P')
set(AX(2),'YLim',yscale,...
    'YTick',yscale(1):1:yscale(2),...
    'YColor',[0 0 0]);
ylabel(AX(2),'z-score')
set(H2,'visible','off');
set(H1,'Color','k')
set(AX,'xlim',[-4 6],...
    'XTickLabel',[])
axis(AX,'square')

% Average Spectrogram
display(sprintf('Calculating and assembling triggered spectrograms.\nThis will take a few moments...'))
for a = 1:length(animals)
    display([num2str(a) ' of ' num2str(length(animals)) ' animals...'])
    prevdir = cd([animals{a} filesep]);
%     EvFile = ls('*EVENTDATA_Gam.mat');
%     load(EvFile)
    EvFile = dir('*EVENTDATA_Gam.mat');
    load(EvFile.name)
    [TrigSpecgram(:,:,a),timevec,Freqs] = TriggeredSpecgram(EventData.Gam,'VW');
    cd(prevdir)
end
Specax = subplot(312);
imagesc(timevec,Freqs,mean(TrigSpecgram,3))
axis xy;
axis square;
caxis([-0.5 1])
set(gca,'XTickLabel',[]);
SpecPosition = get(Specax,'Position');
ypos = SpecPosition(2);
yheight = SpecPosition(4);
cbar = colorbar('Position',[0.61, ypos, 0.02,...
    yheight]);
colormap('parula')
ylabel(cbar,'\DeltaP/P')
ylabel(sprintf('Frequency\n(Hz)'))

% Average CBV
display('Gathering whisk triggered CBV...')
[CBV_VW,CBV_VWZ,timevec] = CompileEventTriggeredData(animals,CBVType,'VW');
preinds = timevec<0;
VWDC = mean(CBV_VW(:,preinds),2)*ones(1,size(CBV_VW,2));
mVW = CBV_VW-VWDC;
VWDCZ = mean(CBV_VWZ(:,preinds),2)*ones(1,size(CBV_VWZ,2));
mVWZ = CBV_VWZ-VWDCZ;
subplot(313);
[AX,H1,H2] = plotyy([timevec', timevec', timevec'],[mean(mVW)', (mean(mVW)+std(mVW))', ...
    (mean(mVW)-std(mVW))'],timevec,mean(mVWZ));
scalefactor = median(mean(mVW)./mean(mVWZ));
yscale = [-4,1.5];
set(AX(1),'YLim',[yscale(1)*scalefactor, yscale(2)*scalefactor],...
    'YTick',round(yscale(1)*scalefactor*100)/100:0.01:round(yscale(2)*scalefactor*100)/100,...
    'box','off');
ylabel(AX(1),'\DeltaR/R')
set(AX(2),'YLim',yscale,...
    'YTick',yscale(1):1:yscale(2),...
    'YColor','k');
ylabel(AX(2),'z-score')
set(AX,'XLim',[-4 6]);
set(H2,'Visible','off');
set(H1,'Color','k')
axis(AX,'square')
xlabel('Peri-stimulus time (s)')
pause(0.001);

%% Figure 2c - Frequency-wise cross correlation
clearvars -except Stats
clc
display('Generating figure 2c...')
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';

% Calculate cross correlation between MUA and CBV
display('Calculating cross correlation between MUA and CBV for periods of rest...')
[XC, ~, lags] = RestingXCorr(animals,'MUpower',CBVType);

% Get the lags
lagfilt = lags>0;
[~,MinInd] = min(XC(:,lagfilt),[],2);
PostLags = lags(lagfilt);
delays = PostLags(MinInd);

figure; subplot(211);
set(gcf,'name','Figure 2c','numbertitle','off')
plot(lags, [mean(XC)' (mean(XC)+std(XC))' (mean(XC)-std(XC))'],'k');
xlim([min(lags), max(lags)])
ylabel('Correlation Coefficient')
title(['Delay: ' num2str(round(mean(delays)*100)/100) '\pm'...
    num2str(round(std(delays)*100)/100)]);
axis square;

% Gather cross correlation data between LFP frequency bands and CBV. Cross
% correlations were calculated by the function 'BehavioralCrossCorrelation_LFPvsCBV.m'
display('Gathering cross correlation between LFP frequencies and CBV...')
for a = 1:length(animals)
    display([num2str(a) ' of ' num2str(length(animals)) ' animals...'])
%     XCorrFile = ls([animals{a} filesep '*_CBVXCorr.mat']);
%     load([animals{a} filesep XCorrFile])
    XCorrFile = dir([animals{a} filesep '*_CBVXCorr.mat']);
    load([animals{a} filesep XCorrFile.name])    
    
    % PreAllocate
    if a == 1
        AllXCorr = NaN*ones(size(CC.LFP.Rest.vals,1),size(CC.LFP.Rest.vals,2),...
            length(animals));
    end
    
    AllXCorr(:,:,a) = mean(CC.LFP.Rest.vals,3);
end

Specax = subplot(212);
imagesc(CC.LFP.Rest.Lags,CC.LFP.Rest.Freqs,mean(AllXCorr,3)); axis xy;
caxis([-0.25 0.25]);
ylim([0 150])
xlabel('Lags (s)')
ylabel('Frequency')
axis square;
SpecPosition = get(Specax,'Position');
ypos = SpecPosition(2);
yheight = SpecPosition(4);
cbar = colorbar('Position',[0.66, ypos, 0.02,...
    yheight]);
ylabel(cbar,'Correlation Coefficient')
colormap('parula')
pause(0.001);

%% Figure 2d - Gamma-band HRF
clearvars -except Stats
clc
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
CBVType = 'CrossCorrROI';
Behaviors = {'Contra','VW','Rest'};
HRFLims = [0 5];
GamHRF = figure;
set(gcf,'name','Figure 2d','numbertitle','off')
ColorOrd = ['b','y','c'];
display('Generating figure 2d...')
display('Calculating HRFs based on Gamma-band power...')
for B = 1:length(Behaviors)
    Beh = Behaviors{B};
    display([num2str(B) ' of ' num2str(length(Behaviors)) ' behaviors...'])
    for a = 1:length(animals)
        display(sprintf(['\t' num2str(a) ' of ' num2str(length(animals)) ' animals...']));
        prevdir = cd([animals{a} filesep]);
        [HRFs] = CalculateHRF_Deconvolution('Gam',CBVType,Beh);
        TimeLims = HRFs.timevec>=HRFLims(1) & HRFs.timevec<=HRFLims(2);
        if a==1
            AllHRFs.Gamma.(Beh).HRFs = NaN*ones(length(animals),sum(TimeLims));
            AllHRFs.Gamma.(Beh).GammaHRFs = NaN*ones(length(animals),sum(TimeLims));
        end
        AllHRFs.Gamma.(Beh).HRFs(a,:) = HRFs.HRF(TimeLims);
        AllHRFs.Gamma.(Beh).GammaHRFs(a,:) = HRFs.GammaHRF;
        cd(prevdir)
    end
    AllHRFs.Gamma.(Beh).Timevec = HRFs.timevec(TimeLims);
    AllHRFs.Gamma.(Beh).Fs = HRFs.Fs;
    mHRFs = mean(AllHRFs.Gamma.(Beh).HRFs);
    figure(GamHRF);
    plot(HRFs.timevec(TimeLims),mHRFs,ColorOrd(B));
    hold on;
end
ylabel('HRF Amplitude (A.U.)');
xlabel('HRF Time (s)');
legend({'Sensory Evoked','Whisking','Rest'},'location','southeast')
title('Gamma-based HRFs')
save('HRFs.mat','AllHRFs'); % save the result
pause(0.001);

%% Figure 2e - MUA-band HRF
clearvars -except Stats
clc
animals = {'w311','w312','w314','w315','w316','w321','w322',...
    'w323','w324','w325','w326','w327'};
load('HRFs.mat'); % Load Gamma HRF
CBVType = 'CrossCorrROI';
Behaviors = {'Contra','VW','Rest'};
HRFLims = [0 5];
MUAHRF = figure;
set(gcf,'name','Figure 2e','numbertitle','off')
ColorOrd = ['b','y','c'];
clc
display('Generating figure 2e...')
display('Calculating HRFs based on MUA...')
for B = 1:length(Behaviors)
    Beh = Behaviors{B};
    display([num2str(B) ' of ' num2str(length(Behaviors)) ' behaviors...'])
    for a = 1:length(animals)
        display(sprintf(['\t' num2str(a) ' of ' num2str(length(animals)) ' animals...']));
        prevdir = cd([animals{a} filesep]);
        [HRFs] = CalculateHRF_Deconvolution('MUpower',CBVType,Beh);
        TimeLims = HRFs.timevec>=HRFLims(1) & HRFs.timevec<=HRFLims(2);
        if a==1
            AllHRFs.MUA.(Beh).HRFs = NaN*ones(length(animals),sum(TimeLims));
            AllHRFs.MUA.(Beh).GammaHRFs = NaN*ones(length(animals),sum(TimeLims));
        end
        AllHRFs.MUA.(Beh).HRFs(a,:) = HRFs.HRF(TimeLims);
        AllHRFs.MUA.(Beh).GammaHRFs(a,:) = HRFs.GammaHRF;
        cd(prevdir)
    end
    AllHRFs.MUA.(Beh).Timevec = HRFs.timevec(TimeLims);
    AllHRFs.MUA.(Beh).Fs = HRFs.Fs;
    mHRFs = mean(AllHRFs.MUA.(Beh).HRFs);
    figure(MUAHRF);
    plot(HRFs.timevec(TimeLims),mHRFs,ColorOrd(B));
    hold on;
end
ylabel('HRF Amplitude (A.U.)');
xlabel('HRF Time (s)');
legend({'Sensory Evoked','Whisking','Rest'},'location','southeast')
title('MUA-based HRFs')
save('HRFs.mat','AllHRFs'); % save the result
pause(0.001);

%% Figure 3b - Ongoing CBV prediction example
clearvars -except Stats AllHRFs animals
clc
display('Generating figure 3b...')
animal = 'w322';
CBVType = 'CrossCorrROI';
filename = 'w322_LH_150724_16_28_1624_ProcData.mat';
load([animal filesep filename]);
[~,~,datename,~] = GetFileInfo(filename);
strdate = ConvertDate(datename);

% BaseFile = ls([animal filesep '*Baselines.mat']);
% load([animal filesep BaseFile]);
BaseFile = dir([animal filesep '*Baselines.mat']);
load([animal filesep BaseFile.name]);

NeurType = 'Gamma';
AnimalRow = strcmp(animals,animal);
HRF = AllHRFs.(NeurType).Contra.GammaHRFs(AnimalRow,:);
CurrentBaseline.Gam = Baselines.Gam.(strdate);
CurrentBaseline.(CBVType) = Baselines.(CBVType).(strdate);
[GammaPred,~] = CalculatePredictedCBV(ProcData,'Gam',CBVType,HRF,CurrentBaseline);

NeurType = 'MUA';
AnimalRow = strcmp(animals,animal);
HRF = AllHRFs.(NeurType).Contra.GammaHRFs(AnimalRow,:);
CurrentBaseline.MUpower = Baselines.MUpower.(strdate);
CurrentBaseline.(CBVType) = Baselines.(CBVType).(strdate);
[MUAPred,~] = CalculatePredictedCBV(ProcData,'MUpower',CBVType,HRF,CurrentBaseline);

PlotTrialCBVPredictions(ProcData,GammaPred,MUAPred,CBVType,CurrentBaseline)
set(gcf,'name','Figure 3b','numbertitle','off')
pause(0.001);

%% FIGURE 3C - SUMMARY OF THE ACCURACY OF THE GAMMA-BASED HRF PREDICTION
clearvars -except Stats AllHRFs animals
clc
display('Generating figure 3c...')
if exist('AllHRFs','var')==0
    load('HRFs.mat')
end
CBVType = 'CrossCorrROI';
Behaviors = {'Contra','Str','VW','Rest'};

% COMPILE THE R^2 VALUES FOR ALL ANIMALS AND BEHAVIORS
CBVPred = struct;
display('Calculating R^2 values for the gamma-based HRF...')
for B = 1:length(Behaviors)
    display([num2str(B) ' of ' num2str(length(Behaviors)) ' behaviors...'])
    Beh = Behaviors{B};
    % Preallocate
    if B==1 
        AveR2.Gamma = NaN*ones(length(animals),length(Behaviors));
        IndR2.Gamma = NaN*ones(length(animals),length(Behaviors));
    end
    for a = 1:length(animals)
        display(sprintf(['\t' num2str(a) ' of ' num2str(length(animals)) ' animals...']));
        prevdir = cd([animals{a} filesep]);
        if strcmp(Beh,'Str')
            HRFs.HRF = AllHRFs.Gamma.VW.GammaHRFs(a,:);
            HRFs.timevec = AllHRFs.Gamma.Contra.Timevec;
        else
            HRFs.HRF = AllHRFs.Gamma.(Beh).GammaHRFs(a,:);
            HRFs.timevec = AllHRFs.Gamma.(Beh).Timevec;
        end
        
        if not(isfield(CBVPred,animals{a}))
            CBVPred.(animals{a}).Gam = [];
        end
        [AveR2.Gamma(a,B),IndR2.Gamma(a,B),CBVPred.(animals{a}).Gam] ...
            = EvaluateCBVPredictionAccuracy('Gam',CBVType,Beh,HRFs,...
            CBVPred.(animals{a}).Gam);
        cd(prevdir)
    end
end

% GENERATE THE PLOT
R2_inds = repmat(1:length(Behaviors),length(animals),1);
figure; 
set(gcf,'name','Figure 3c','numbertitle','off')
scatter(R2_inds(:),AveR2.Gamma(:),'MarkerEdgeColor','k');
hold on; scatter(1:length(Behaviors),median(AveR2.Gamma),...
    'MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(R2_inds(:)+0.25,IndR2.Gamma(:),'MarkerEdgeColor','k');
scatter((1:length(Behaviors))+0.25,median(IndR2.Gamma),...
    'MarkerEdgeColor','k','MarkerFaceColor','k');
title('Summary of R^2 values for average and individual events')
ylabel('R^2')
xlim([0 length(Behaviors)+1]);
ylim([-0.2 1])
ax1 = gca;
set(ax1,'XTick', min(R2_inds(:)):1:max(R2_inds(:)),...
    'XTickLabel',{'Sens. Ev.','Ext. Move.','Whisk','Rest'});

[R2Stats] = GammaHRFPredictionStats_Rsquared(AveR2,IndR2);
Stats.R2Compare.GammaHRF = R2Stats;
pause(0.001);

%% Figure 3d - MUA-band prediction summary
clearvars -except Stats AllHRFs animals CBVPred AveR2 IndR2
clc
display('Generating figure 3d...')
if exist('AllHRFs','var')==0
    load('HRFs.mat')
end
CBVType = 'CrossCorrROI';
Behaviors = {'Contra','Str','VW','Rest'};
display('Calculating R^2 values for the MUA-based HRF...')
for B = 1:length(Behaviors)
    display([num2str(B) ' of ' num2str(length(Behaviors)) ' behaviors...'])
    Beh = Behaviors{B};
    % Preallocate
    if B==1 
        AveR2.MUA = NaN*ones(length(animals),length(Behaviors));
        IndR2.MUA = NaN*ones(length(animals),length(Behaviors));
    end
    for a = 1:length(animals)
        display(sprintf(['\t' num2str(a) ' of ' num2str(length(animals)) ' animals...']));
        prevdir = cd([animals{a} filesep]);
        if strcmp(Beh,'Str')
%             HRFs.HRF = AllHRFs.MUA.VW.HRFs(a,:);
            HRFs.HRF = AllHRFs.MUA.VW.GammaHRFs(a,:);
            HRFs.timevec = AllHRFs.MUA.VW.Timevec;
        else
%             HRFs.HRF = AllHRFs.MUA.(Beh).HRFs(a,:);
            HRFs.HRF = AllHRFs.MUA.(Beh).GammaHRFs(a,:);
            HRFs.timevec = AllHRFs.MUA.(Beh).Timevec;
        end
        if not(isfield(CBVPred.(animals{a}),'MUpower'))
            CBVPred.(animals{a}).MUpower = [];
        end
        [AveR2.MUA(a,B),IndR2.MUA(a,B),CBVPred.(animals{a}).MUpower] ...
            = EvaluateCBVPredictionAccuracy('MUpower',CBVType,Beh,HRFs,...
            CBVPred.(animals{a}).MUpower);
        cd(prevdir)
    end
end

R2_inds = repmat(1:length(Behaviors),length(animals),1);
figure;
set(gcf,'name','Figure 3d','numbertitle','off');
scatter(R2_inds(:),AveR2.MUA(:),'MarkerEdgeColor','k');
hold on; scatter(1:length(Behaviors),median(AveR2.MUA),...
    'MarkerEdgeColor','k','MarkerFaceColor','k');
scatter(R2_inds(:)+0.25,IndR2.MUA(:),'MarkerEdgeColor','k');
scatter((1:length(Behaviors))+0.25,median(IndR2.MUA),...
    'MarkerEdgeColor','k','MarkerFaceColor','k');
title('Summary of R^2 values for average and individual events')
ylabel('R^2')
xlim([0 length(Behaviors)+1]);
ax1 = gca;
set(ax1,'XTick',min(R2_inds(:)):1:max(R2_inds(:)),...
    'XTickLabel',{'Sens. Ev.','Ext. Move.','Whisk','Rest'});

[R2Stats] = MUAHRFPredictionStats_Rsquared(AveR2,IndR2);
Stats.R2Compare.MUAHRF = R2Stats;

save('R2.mat','AveR2','IndR2');
save('CBVPredictions.mat','CBVPred');
pause(0.001);

%% Figure 3a - Average CBV prediction example
clearvars -except Stats AllHRFs animals CBVPred AveR2 IndR2
clc
display('Generating figure 3a...')
CBVType = 'CrossCorrROI';
animal = 'w327';

% Calculate the CBV predictions
% ActFile = ls([animal filesep '*EVENTDATA_' CBVType '.mat']);
% load([animal filesep ActFile]);
ActFile = dir([animal filesep '*EVENTDATA_' CBVType '.mat']);
load([animal filesep ActFile.name]);
[DataStruct,FiltArray] = SelectBehavioralEvents(EventData.(CBVType),'Contra');
NormData = DataStruct.NormData(FiltArray,:);
test_inds = 2:2:size(NormData,1);
strt = round(DataStruct.epoch.offset*DataStruct.Fs);
ActDC = NormData(test_inds,strt)*ones(1,size(NormData,2)-strt+1);
Act = mean(NormData(test_inds,strt:end)-ActDC);
timevec = (1:size(Act,2))/DataStruct.Fs;
GamPred = mean(CBVPred.(animal).Gam.Contra(:,strt:end));
MUAPred = mean(CBVPred.(animal).MUpower.Contra(:,strt:end));
R2MUA = CalculateRsquared(MUAPred-mean(MUAPred),Act-mean(Act));
R2Gamma = CalculateRsquared(GamPred-mean(GamPred),Act-mean(Act));

% Plot
figure; 
set(gcf,'name','Figure 3a','numbertitle','off');
h = plot(timevec,[(Act-mean(Act)+mean(GamPred))', GamPred', MUAPred']);
title(['Gamma band: R^2 = ' num2str(R2Gamma) '; MUA: R^2 = ' num2str(R2MUA)])
xlabel('Peri-stimulus time (s)')
ylabel('\DeltaR/R')
legend({'Measured','Gamma-predicted','MUA-Predicted'},'location','southeast')
set(h(2),'Color',[0 167 157]/255,'Linewidth',1.5)
set(h(1),'Color','k','Linewidth',1.5);
set(h(3),'Color',[0.4 0.4 0.4],'Linewidth',1.5)
xlim([0 4])
pause(0.001);

%% Figure 3e - Comparison of resting predictions
clearvars -except Stats AllHRFs animals CBVPred AveR2 IndR2
clc
display('Generating figure 3e...')
figure; 
set(gcf,'name','Figure 3e','numbertitle','off');
scatter(ones(1,length(animals)), IndR2.Gamma(:,4),'co','MarkerFaceColor','c');
hold on; scatter(2*ones(1,length(animals)), IndR2.MUA(:,4),'co','MarkerFaceColor','c');
scatter([1,2],[median(IndR2.Gamma(:,4)), median(IndR2.MUA(:,4))],'ko',...
    'MarkerFaceColor','k');
xlim([0 3]);
ax1 = gca;
set(ax1,'XTick',[1,2],'XTickLabel',{'Gamma-band Power','MUA'});
ylim([-0.2 1])
ylabel('R^2')

% Statistics - Resting R^2 between neural measures
[~,Stats.R2Compare.RestingGammaVsMUA.pval,~,stat] = ttest(IndR2.Gamma(:,4),IndR2.MUA(:,4));
Stats.R2Compare.RestingGammaVsMUA.tstat = stat.tstat;
Stats.R2Compare.RestingGammaVsMUA.df = stat.df;
pause(0.001);

%% Figure 3f - Behavioral example trial
clearvars -except Stats AllHRFs animals CBVPred
clc
display('Generating figure 3f...')
CBVType = 'CrossCorrROI';
animal = 'w327';
filename = 'w327_LH_150804_13_37_4104_ProcData.mat';
load([animal filesep filename]);
[~,~,datename,~] = GetFileInfo(filename);
strdate = ConvertDate(datename);

% BaseFile = ls([animal filesep '*Baselines.mat']);
% load([animal filesep BaseFile]);
BaseFile = dir([animal filesep '*Baselines.mat']);
load([animal filesep BaseFile.name]);

NeurType = 'Gamma';
AnimalRow = strcmp(animals,animal);
HRF = AllHRFs.(NeurType).Contra.GammaHRFs(AnimalRow,:);
CurrentBaseline.Gam = Baselines.Gam.(strdate);
CurrentBaseline.(CBVType) = Baselines.(CBVType).(strdate);
[GammaPred,~] = CalculatePredictedCBV(ProcData,'Gam',CBVType,HRF,CurrentBaseline);

NeurType = 'MUA';
AnimalRow = strcmp(animals,animal);
HRF = AllHRFs.(NeurType).Contra.GammaHRFs(AnimalRow,:);
CurrentBaseline.MUpower = Baselines.MUpower.(strdate);
CurrentBaseline.(CBVType) = Baselines.(CBVType).(strdate);
[MUAPred,~] = CalculatePredictedCBV(ProcData,'MUpower',CBVType,HRF,CurrentBaseline);

PlotTrialCBVPredictions(ProcData,GammaPred,MUAPred,CBVType,CurrentBaseline)
set(gcf,'name','Figure 3f','numbertitle','off');
pause(0.001);

%% Figure 4b - Average behavioral CBV variance comparison
clearvars -except Stats animals CBVPred
clc
display('Generating figure 4b...')
CBVType = 'CrossCorrROI';
[ITV_Stats] = CompareInterTrialVariance(animals,CBVType);
set(gcf,'name','Figure 4b','numbertitle','off');
Stats.CompareInterTrialVariance = ITV_Stats;
pause(0.001);

%% Figure 4c - Gamma-band prediction residual variance comparison
clearvars -except Stats animals CBVPred
clc
display('Generating figure 4c...')
CBVType = 'CrossCorrROI';
[GamResidStats] = CompareResidualRMSE(animals,'Gam',CBVType,CBVPred);
ylabel('Gamma-band prediction RMSE / Resting residual RMSE')
set(gcf,'name','Figure 4c','numbertitle','off');
Stats.CompareResidualRMSE.GammaHRF = GamResidStats;
pause(0.001);

%% Figure 4d - MUA prediction residual variance comparison
clearvars -except Stats animals CBVPred
clc
display('Generating figure 4d...')
CBVType = 'CrossCorrROI';
[MUAResidStats] = CompareResidualRMSE(animals,'MUpower',CBVType,CBVPred);
ylabel('MUA prediction RMSE / Resting residual RMSE')
set(gcf,'name','Figure 4d','numbertitle','off');
Stats.CompareResidualRMSE.MUAHRF = MUAResidStats;
pause(0.001);

%% Figure 5b - aCSF/muscimol average response comparison
clearvars -except Stats
clc
display('Generating figure 5b...')
animal = 'CAN12';
InfusionTypes = {'aCSF','Muscimol'};
colorcodes = {'k','c'};
figure;
set(gcf,'name','Figure 5b','numbertitle','off');
for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    prevdir = cd([animal filesep InfusionType filesep]);
%     ThreshFile = ls('*Thresholds.mat');
%     load(ThreshFile)
    ThreshFile = dir('*Thresholds.mat');
    load(ThreshFile.name)
    if strcmp(InfusionType,'aCSF')
        Threshfields = fieldnames(Thresholds.Infusion_Start_Times);
        for Tf = 1:length(Threshfields)
            AverageThresh.(Threshfields{Tf}) = 20;
        end
    else
        AverageThresh = [];
        Threshfields = fieldnames(Thresholds.ReducedMUpower_Start);
        for Tf = 1:length(Threshfields)
            if isempty(Thresholds.ReducedMUpower_Start.(Threshfields{Tf}))
                continue;
            end
            AverageThresh.(Threshfields{Tf}) = ElapsedTime(Thresholds.ReducedMUpower_Start.(Threshfields{Tf}),...
                0,Thresholds.Infusion_Start_Times.(Threshfields{Tf}));
        end
    end
    EventDataFileName = dir('*EVENTDATA_MUpower.mat');
    load(EventDataFileName.name);
    subplot(211);
    PlotEventTriggeredAve_Infusion(EventData,'MUpower',...
        'Contra',AverageThresh,colorcodes{IT});
    ylabel('MUA (V^2/s)')
    ylim([0 1e-9])
    hold on;
    
    EventDataFileName = dir('*EVENTDATA_RadiusROI.mat');
    load(EventDataFileName.name);
    subplot(212);
    PlotEventTriggeredAve_Infusion(EventData,'RadiusROI',...
        'Contra',AverageThresh,colorcodes{IT});
    ylabel('\DeltaR/R')
    ylim([-0.02 0.01])
    hold on;
    cd(prevdir)
end
subplot(211); legend(InfusionTypes,'location','northeast')
pause(0.001);

%% Figure 5c - Spatial map of muscimol/aCSF CBV example
clearvars -except Stats
clc
display('Generating figure 5c...')
animal = 'CAN12';
InfusionType = 'Muscimol';
load([animal filesep 'PharmSpatialMap.mat']);
% Register the images
[~,RA,RB] = RegisterCBVImage(SampleFrames_aCSF(:,:,1),SampleFrames_Pharm(:,:,1));

% Combine the images to get the global map
[StimChangeMap,WorldPharmStim,WorldaCSFStim] = MapSpatialCBVChange(-1*PharmStimMap,RB,-1*aCSFStimMap,RA);
aCSFStdMap = sqrt(aCSFVarMap);
PharmStdMap = sqrt(PharmVarMap);
[RestChangeMap,WorldPharmRest,WorldaCSFRest] = MapSpatialCBVChange(PharmStdMap,RB,aCSFStdMap,RA);

[rows,cols] = find(mask);

figure; 
set(gcf,'name','Figure 5c','numbertitle','off');
subplot(231); imagesc(-1*WorldaCSFStim.*mask); axis image;
caxis([-0.02 0.02]); title('aCSF')
set(gca,'XTickLabel','','YTickLabel','')
% xlim([min(cols) max(cols)]);
% ylim([min(rows) max(rows)]);
cbar = colorbar('southoutside');
% ylabel(cbar,'\DeltaR/R')
xlabel(cbar,'\DeltaR/R')
colormap('parula')
subplot(232); imagesc(-1*WorldPharmStim.*mask); axis image;
caxis([-0.02 0.02]); title(InfusionType);
set(gca,'XTickLabel','','YTickLabel','')
% xlim([min(cols) max(cols)]);
% ylim([min(rows) max(rows)]);
cbar = colorbar('southoutside');
xlabel(cbar,'\DeltaR/R')
colormap('parula')
subplot(233); imagesc(StimChangeMap.*mask); axis image;
caxis([-0.01 0.02]);
set(gca,'XTickLabel','','YTickLabel','')
% xlim([min(cols) max(cols)]);
% ylim([min(rows) max(rows)]);
cbar = colorbar('southoutside');
xlabel(cbar,'\DeltaR/R_{Muscimol}-\DeltaR/R_{aCSF}')
colormap('parula')

subplot(234); imagesc(WorldaCSFRest.*mask); axis image;
caxis([0 0.02]); title('aCSF')
set(gca,'XTickLabel','','YTickLabel','')
% xlim([min(cols) max(cols)]);
% ylim([min(rows) max(rows)]);
cbar = colorbar('southoutside');
xlabel(cbar,'\sigma_{\DeltaR/R}')
colormap('parula')
subplot(235); imagesc(WorldPharmRest.*mask); axis image;
caxis([0 0.02]); title(strrep(InfusionType,'_',' '));
set(gca,'XTickLabel','','YTickLabel','')
% xlim([min(cols) max(cols)]);
% ylim([min(rows) max(rows)]);
cbar = colorbar('southoutside');
xlabel(cbar,'\sigma_{\DeltaR/R}')
colormap('parula')
subplot(236); imagesc(RestChangeMap.*mask); axis image;
caxis([-0.01 0.02]);
set(gca,'XTickLabel','','YTickLabel','')
% xlim([min(cols) max(cols)]);
% ylim([min(rows) max(rows)]);
cbar = colorbar('southoutside');
% ylabel(cbar,'\sigma_{aCSF}-\sigma_{Muscimol}')
xlabel(cbar,'\sigma_{aCSF}-\sigma_{Muscimol}')
colormap('parula')
pause(0.001);

%% Figure 5d,top - Example aCSF and muscimol resting oscillations
clearvars -except Stats
clc
display('Generating figure 5d...')
% Plot the CBV examples
animal = 'CAN12';
CBVType = 'RadiusROI';
aCSFTrial = [animal '_LH_160412_11_15_0312_ProcData.mat'];
aCSFStart = 240; % seconds
MuscTrial = [animal '_LH_160407_11_16_3407_ProcData.mat'];
MuscStart = 230; % seconds
figure; 
set(gcf,'name','Figure 5d','numbertitle','off');

% Load aCSF example trial
load([animal filesep 'aCSF' filesep aCSFTrial])
[~,~,FileDate,~] = GetFileInfo(aCSFTrial);
strdate = ConvertDate(FileDate);
% Load the baselines
% basefile = ls([animal filesep 'aCSF' filesep '*_Baselines.mat']);
% load([animal filesep 'aCSF' filesep basefile])
basefile = dir([animal filesep 'aCSF' filesep '*_Baselines.mat']);
load([animal filesep 'aCSF' filesep basefile.name])
baseline = Baselines.(CBVType).(strdate).PreInfusionMeans;
Fs = ProcData.Fs.(CBVType);
NormaCSF = ProcData.CBV.(CBVType)/baseline;

[z,p,k] = butter(4,2/(Fs/2),'low');
[sos,g] = zp2sos(z,p,k);
filtaCSF = filtfilt(sos,g,detrend(NormaCSF));
subplot(211);
plot(0:1/Fs:20,detrend(filtaCSF((aCSFStart*Fs):((aCSFStart+20)*Fs))),'k',...
    'Linewidth',1.5);

% Load aCSF example trial
load([animal filesep 'Muscimol' filesep MuscTrial])
[~,~,FileDate,~] = GetFileInfo(MuscTrial);
strdate = ConvertDate(FileDate);
% Load the baselines
% basefile = ls([animal filesep 'Muscimol' filesep '*_Baselines.mat']);
% load([animal filesep 'Muscimol' filesep basefile])
basefile = dir([animal filesep 'Muscimol' filesep '*_Baselines.mat']);
load([animal filesep 'Muscimol' filesep basefile.name])
baseline = Baselines.(CBVType).(strdate).PreInfusionMeans;
Fs = ProcData.Fs.(CBVType);
NormaCSF = ProcData.CBV.(CBVType)/baseline;

[z,p,k] = butter(4,2/(Fs/2),'low');
[sos,g] = zp2sos(z,p,k);
filtaCSF = filtfilt(sos,g,detrend(NormaCSF));
hold on;
plot(0:1/Fs:20,detrend(filtaCSF((MuscStart*Fs):((MuscStart+20)*Fs))),'c',...
    'Linewidth',1.5);
ylim([-0.03 0.03])
legend({'aCSF','Muscimol'},'location','southeast','orientation','horizontal')
ylabel('\DeltaR/R')
xlabel('Time (s)')

%% Figures 5d(bottom),e - 2PLSM muscimol comparison
clearvars -except Stats
display('Generating figure 5e...')
prevdir = cd(['CEData' filesep 'Muscimol Data for Paper_CE']);
subplot(212);
Muscimol_graph
set(gcf,'name','Figure 5e','numbertitle','off');
cd(prevdir)
pause(0.001);

%% Figure 5f - 2PLSM muscimol + CNQX + AP5 comparison
clearvars -except Stats
clc
display('Generating figure 5f...')
prevdir = cd(['CEData' filesep 'Musc_APV_CNQX_For_Paper_CE']);
Muscimol_APV_CNQX_graph
set(gcf,'name','Figure 5f','numbertitle','off');
cd(prevdir)
pause(0.001);

%% Gather the Muscimol data for the pharmacology comparisons
clearvars -except Stats
clc
animals = {'CAN08','CAN10','CAN11','CAN12','CAN13','CAN15','CAN16',...
    'CAN18','CAN24','CAN27','CAN28','CAN31'};
CBVType = 'RadiusROI';
NeurTypes = {'Gam','MUpower'};
InfusionTypes = {'Muscimol'};
PharmData = [];
for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    for a = 1:length(animals) 
        animal = animals{a};
        if not(isdir([animal filesep InfusionType filesep]))
            continue;
        end
        prevdir = cd([animal filesep InfusionType filesep]);
        % Load the CBV data
%         RestFile = ls(['*_RESTDATA_' CBVType '.mat']);
%         if isempty(RestFile)
%             error(['No rest data found for ' CBVType])
%         else
%             load(RestFile)
%         end
        RestFile = dir(['*_RESTDATA_' CBVType '.mat']);
        if isempty(RestFile)
            error(['No rest data found for ' CBVType])
        else
            load(RestFile.name)
        end
        if isempty(RestData.(CBVType).AmpRatio)
            cd(prevdir)
            continue;
        end
        PharmData.(CBVType).(animal) =...
            RestData.(CBVType).AmpRatio_SD;
        PharmData.(CBVType).Means(a) = RestData.(CBVType).AmpRatio;
        
        for NT = 1:length(NeurTypes)
            NeurType = NeurTypes{NT};
           
            % Load the neural rest data
%             RestFile = ls(['*_RESTDATA_' NeurType '.mat']);
%             if isempty(RestFile)
%                 error(['No rest data found for ' NeurType])
%             else
%                 load(RestFile)
%             end
            RestFile = dir(['*_RESTDATA_' NeurType '.mat']);
            if isempty(RestFile)
                error(['No rest data found for ' NeurType])
            else
                load(RestFile.name)
            end
            PharmData.(NeurType).(animal) = RestData.(NeurType).AmpRatio_SD;
            PharmData.(NeurType).Means(a) = RestData.(NeurType).AmpRatio;
            
        end
        cd(prevdir)
    end
end

%% Figure 5g - Resting CBV vs MUA - muscimol
display('Generating figure 5g...')
clc
CBVType = 'RadiusROI';
NeurType = 'MUpower';
animals = {'CAN08','CAN10','CAN11','CAN12','CAN13','CAN15','CAN16',...
    'CAN18','CAN24','CAN27','CAN28','CAN31'};
InfusionType = 'Muscimol';

Pharmacology_PlotNeuroVsCBV(PharmData,animals,NeurType,CBVType,InfusionType);
set(gcf,'name','Figure 5g','numbertitle','off');
ylabel('\DeltaR/R \sigma_{muscimol} / \DeltaR/R \sigma_{aCSF}')
xlabel('MUA \sigma_{muscimol} / MUA \sigma_{aCSF}')
pause(0.001);

%% Figure 5j - Resting CBV vs Gamma-power - muscimol
display('Generating figure 5j...')
clc
CBVType = 'RadiusROI';
NeurType = 'Gam';
animals = {'CAN08','CAN10','CAN11','CAN12','CAN13','CAN15','CAN16',...
    'CAN18','CAN24','CAN27','CAN28','CAN31'};
InfusionType = 'Muscimol';

Pharmacology_PlotNeuroVsCBV(PharmData,animals,NeurType,CBVType,InfusionType);
set(gcf,'name','Figure 5j','numbertitle','off');
ylabel('\DeltaR/R \sigma_{muscimol} / \DeltaR/R \sigma_{aCSF}')
xlabel('Gamma Power \sigma_{muscimol} / Gamma Power \sigma_{aCSF}')

InfusionMeans.Muscimol.MUA = nonzeros(PharmData.MUpower.Means);
InfusionMeans.Muscimol.(CBVType) = nonzeros(PharmData.(CBVType).Means);
InfusionMeans.Muscimol.Gamma = nonzeros(PharmData.Gam.Means);
pause(0.001);

%% Gather the Muscimol + CNQX + AP5 data for the pharmacology comparisons
clearvars -except Stats InfusionMeans
clc
animals = {'CAN08','CAN10','CAN11','CAN12','CAN13','CAN15','CAN16',...
    'CAN18','CAN24','CAN27','CAN28','CAN31'};
CBVType = 'RadiusROI';
NeurTypes = {'Gam','MUpower'};
InfusionTypes = {'Muscimol + CNQX + AP5'};
PharmData = [];
for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    for a = 1:length(animals) 
        animal = animals{a};
        if not(isdir([animal filesep InfusionType filesep]))
            continue;
        end
        prevdir = cd([animal filesep InfusionType filesep]);
        % Load the CBV data
%         RestFile = ls(['*_RESTDATA_' CBVType '.mat']);
%         if isempty(RestFile)
%             error(['No rest data found for ' CBVType])
%         else
%             load(RestFile)
%         end
        RestFile = dir(['*_RESTDATA_' CBVType '.mat']);
        if isempty(RestFile)
            error(['No rest data found for ' CBVType])
        else
            load(RestFile.name)
        end
        if isempty(RestData.(CBVType).AmpRatio)
            cd(prevdir)
            continue;
        end
        PharmData.(CBVType).(animal) =...
            RestData.(CBVType).AmpRatio_SD;
        PharmData.(CBVType).Means(a) = RestData.(CBVType).AmpRatio;
        
        for NT = 1:length(NeurTypes)
            NeurType = NeurTypes{NT};
           
            % Load the neural rest data
%             RestFile = ls(['*_RESTDATA_' NeurType '.mat']);
%             if isempty(RestFile)
%                 error(['No rest data found for ' NeurType])
%             else
%                 load(RestFile)
%             end
            RestFile = dir(['*_RESTDATA_' NeurType '.mat']);
            if isempty(RestFile)
                error(['No rest data found for ' NeurType])
            else
                load(RestFile.name)
            end
            PharmData.(NeurType).(animal) = RestData.(NeurType).AmpRatio_SD;
            PharmData.(NeurType).Means(a) = RestData.(NeurType).AmpRatio;
            
        end
        cd(prevdir)
    end
end

%% Figure 5h - Resting CBV vs Neuro - muscimol + CNQX + AP5
display('Generating figure 5h...')
clc
CBVType = 'RadiusROI';
NeurType = 'MUpower';
animals = {'CAN08','CAN10','CAN11','CAN12','CAN13','CAN15','CAN16',...
    'CAN18','CAN24','CAN27','CAN28','CAN31'};
InfusionType = 'Muscimol + CNQX + AP5';

Pharmacology_PlotNeuroVsCBV(PharmData,animals,NeurType,CBVType,InfusionType);
set(gcf,'name','Figure 5h','numbertitle','off');
ylabel('\DeltaR/R \sigma_{musc.+CNQX+AP5} / \DeltaR/R \sigma_{aCSF}')
xlabel('MUA \sigma_{musc.+CNQX+AP5} / MUA \sigma_{aCSF}')
pause(0.001);

%% Figure 5k - Sens. Ev. CBV vs Neuro - muscimol + CNQX + AP5
display('Generating figure 5k...')
CBVType = 'RadiusROI';
NeurType = 'Gam';
animals = {'CAN08','CAN10','CAN11','CAN12','CAN13','CAN15','CAN16',...
    'CAN18','CAN24','CAN27','CAN28','CAN31'};
InfusionType = 'Muscimol + CNQX + AP5';

Pharmacology_PlotNeuroVsCBV(PharmData,animals,NeurType,CBVType,InfusionType);
set(gcf,'name','Figure 5k','numbertitle','off');
ylabel('\DeltaR/R \sigma_{musc.+CNQX+AP5} / \DeltaR/R \sigma_{aCSF}')
xlabel('Gamma Power \sigma_{musc.+CNQX+AP5} / Gamma Power \sigma_{aCSF}')

InfusionMeans.MCA.MUA = nonzeros(PharmData.MUpower.Means);
InfusionMeans.MCA.(CBVType) = nonzeros(PharmData.(CBVType).Means);
InfusionMeans.MCA.Gamma = nonzeros(PharmData.Gam.Means);

% Stats on the reduction due to cocktail infusion compared to muscimol-only
% infusion
    % MUA
[h,p] = adtest(InfusionMeans.MCA.MUA); % normal
[h,p] = adtest(InfusionMeans.Muscimol.MUA); % normal
[h,p] = vartest2(InfusionMeans.Muscimol.MUA,InfusionMeans.MCA.MUA); % similar variance
[~,Stats.PharmaCompare.MCA.MUA.pval,~,stat] =...
    ttest2(InfusionMeans.Muscimol.MUA,InfusionMeans.MCA.MUA);
Stats.PharmaCompare.MCA.MUA.tstat = stat.tstat;
Stats.PharmaCompare.MCA.MUA.df = stat.df;
    % CBV - one sided t-test, does cocktail reduce CBV further than
    % muscimol
[h,p] = adtest(InfusionMeans.MCA.(CBVType)); % approximately normal
[h,p] = adtest(InfusionMeans.Muscimol.(CBVType)); % normal
[h,p] = vartest2(InfusionMeans.Muscimol.(CBVType),InfusionMeans.MCA.(CBVType)); % similar variance   
[~,Stats.PharmaCompare.MCA.(CBVType).pval,~,stat] =...
    ttest2(InfusionMeans.Muscimol.(CBVType),InfusionMeans.MCA.(CBVType),...
    'Tail','right');
Stats.PharmaCompare.MCA.(CBVType).tstat = stat.tstat;
Stats.PharmaCompare.MCA.(CBVType).df = stat.df;
    % Gamma-band power - one-sided, does cocktail reduce Gamma power
    % further than muscimol
[h,p] = adtest(InfusionMeans.MCA.Gamma); % normal
[h,p] = adtest(InfusionMeans.Muscimol.Gamma); % normal
[h,p] = vartest2(InfusionMeans.Muscimol.Gamma,InfusionMeans.MCA.Gamma); % similar variance
[~,Stats.PharmaCompare.MCA.Gamma.pval,~,stat] =...
    ttest2(InfusionMeans.Muscimol.Gamma,InfusionMeans.MCA.Gamma,...
    'Tail','right');
Stats.PharmaCompare.MCA.Gamma.tstat = stat.tstat;
Stats.PharmaCompare.MCA.Gamma.df = stat.df;
pause(0.001);

%% Gather the Muscimol + Adr. data for the pharmacology comparisons
clearvars -except Stats InfusionMeans
display('Generating figures 5i,l...')
animals = {'CAN08','CAN10','CAN11','CAN12','CAN13','CAN15','CAN16',...
    'CAN18','CAN24','CAN27','CAN28','CAN31'};
CBVType = 'RadiusROI';
NeurTypes = {'Gam','MUpower'};
InfusionTypes = {'Muscimol + Adrenergic Blockers'};
PharmData = [];
for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    for a = 1:length(animals) 
        animal = animals{a};
        if not(isdir([animal filesep InfusionType filesep]))
            continue;
        end
        prevdir = cd([animal filesep InfusionType filesep]);
        % Load the CBV data
%         RestFile = ls(['*_RESTDATA_' CBVType '.mat']);
%         if isempty(RestFile)
%             error(['No rest data found for ' CBVType])
%         else
%             load(RestFile)
%         end
        RestFile = dir(['*_RESTDATA_' CBVType '.mat']);
        if isempty(RestFile)
            error(['No rest data found for ' CBVType])
        else
            load(RestFile.name)
        end
        if isempty(RestData.(CBVType).AmpRatio)
            cd(prevdir)
            continue;
        end
        PharmData.(CBVType).(animal) =...
            RestData.(CBVType).AmpRatio_SD;
        PharmData.(CBVType).Means(a) = RestData.(CBVType).AmpRatio;
        
        for NT = 1:length(NeurTypes)
            NeurType = NeurTypes{NT};
           
            % Load the neural rest data
%             RestFile = ls(['*_RESTDATA_' NeurType '.mat']);
%             if isempty(RestFile)
%                 error(['No rest data found for ' NeurType])
%             else
%                 load(RestFile)
%             end
            RestFile = dir(['*_RESTDATA_' NeurType '.mat']);
            if isempty(RestFile)
                error(['No rest data found for ' NeurType])
            else
                load(RestFile.name)
            end
            PharmData.(NeurType).(animal) = RestData.(NeurType).AmpRatio_SD;
            PharmData.(NeurType).Means(a) = RestData.(NeurType).AmpRatio;
            
        end
        cd(prevdir)
    end
end

%% Figure 5i,l - Resting CBV vs Neuro - Muscimol + Adr.
CBVType = 'RadiusROI';
clc
prevdir = cd(['QZData' filesep]);
[CBVData, MUData, GamData] = CBV_neural_fluctuations_plot(PharmData.MUpower, PharmData.Gam, PharmData.(CBVType));
cd(prevdir)

InfusionMeans.Adr.MUA = nonzeros(MUData);
InfusionMeans.Adr.(CBVType) = nonzeros(CBVData);
InfusionMeans.Adr.Gamma = nonzeros(GamData);

% Stats on the reduction due to cocktail infusion compared to muscimol-only
% infusion
    % MUA - Welch's t-test
[h,p] = adtest(InfusionMeans.Adr.MUA); % normal
[h,p] = adtest(InfusionMeans.Muscimol.MUA); % normal
[h,p] = vartest2(InfusionMeans.Muscimol.MUA,InfusionMeans.Adr.MUA); % unequal variance
[~,Stats.PharmaCompare.Adr.MUA.pval,~,stat] =...
    ttest2(InfusionMeans.Muscimol.MUA,InfusionMeans.Adr.MUA,...
    'VarType','unequal');
Stats.PharmaCompare.Adr.MUA.tstat = stat.tstat;
Stats.PharmaCompare.Adr.MUA.df = stat.df;
    % CBV - one sided t-test, does cocktail reduce CBV further than
    % muscimol
[h,p] = adtest(InfusionMeans.Adr.(CBVType)); % normal
[h,p] = adtest(InfusionMeans.Muscimol.(CBVType)); % normal
[h,p] = vartest2(InfusionMeans.Muscimol.(CBVType),InfusionMeans.Adr.(CBVType)); % similar variance   
[~,Stats.PharmaCompare.Adr.(CBVType).pval,~,stat] =...
    ttest2(InfusionMeans.Muscimol.(CBVType),InfusionMeans.Adr.(CBVType),...
    'Tail','right');
Stats.PharmaCompare.Adr.(CBVType).tstat = stat.tstat;
Stats.PharmaCompare.Adr.(CBVType).df = stat.df;
    % Gamma-band power - one-sided, does cocktail reduce Gamma power
    % further than muscimol
[h,p] = adtest(InfusionMeans.Adr.Gamma); % normal
[h,p] = adtest(InfusionMeans.Muscimol.Gamma); % normal
[h,p] = vartest2(InfusionMeans.Muscimol.Gamma,InfusionMeans.Adr.Gamma); % similar variance
[~,Stats.PharmaCompare.Adr.Gamma.pval,~,stat] =...
    ttest2(InfusionMeans.Muscimol.Gamma,InfusionMeans.Adr.Gamma);
Stats.PharmaCompare.Adr.Gamma.tstat = stat.tstat;
Stats.PharmaCompare.Adr.Gamma.df = stat.df;

display('Done generating main figures.')