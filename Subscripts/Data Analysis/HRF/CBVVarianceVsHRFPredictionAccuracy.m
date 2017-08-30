function [slope,pvals] = CBVVarianceVsHRFPredictionAccuracy(animals,NeurType,...
    CBVType,HRFBeh,HRFStruct)
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
%               animals - [cell array] IDs for all animals
%
%               NeurType - [string] fieldname of predicted CBV structures
%               (_CBVPred_"NeurType".mat)
%
%               CBVType - [string] fieldname designating the CBV ROI to be
%               used.
%
%               HRFBeh - [string] behavior designating which HRF should be
%               used for the prediction.
%_______________________________________________________________
%   RETURN:
%
%_______________________________________________________________

HRFType = [NeurType HRFBeh];
SegLength = 20;
PuffBuffer = 5;
MedR = zeros(1,length(animals));

% Parameters for filtering the resting data
RestCriteria.Fieldname = {'Duration','PuffDistance'};
RestCriteria.Comparison = {'gt','gt'};
RestCriteria.Value = {8,5};
RestBuffer = 2;


animalVars = cell(1,length(animals));
animalRs = cell(1,length(animals));
animalIDArray = cell(1,length(animals));
slope = zeros(1,length(animals));
intercept = zeros(1,length(animals));
pvals = zeros(1,length(animals));

for a = 1:length(animals)
    display([num2str(a) ' of ' num2str(length(animals)) ' animals...']);
    animal = animals{a};
    prevdir = cd([animal filesep]);
    
    % Calculate Resting Variance 
%     RestFile = ls(['*_RESTDATA_' CBVType '.mat']);
%     load(RestFile);
    RestFile = dir(['*_RESTDATA_' CBVType '.mat']);
    load(RestFile.name);
    RestFs = RestData.(CBVType).Fs;
    [RestFiltArray] = FilterEvents(RestData.(CBVType),RestCriteria);
    RestData = RestData.(CBVType).NormData(RestFiltArray);
    BaselineVars = zeros(1,length(RestData));
    for RD = 1:length(RestData)
        IndRest = RestData{RD};
        BaselineVars(RD) = var(IndRest(RestBuffer*RestFs:end));
    end
    BaselineVar = mean(BaselineVars);
    
    % Compare the CBV variance to the correlation for each trial
%     filenames = ls('*ProcData.mat');
%     Basefile = ls('*Baselines.mat');
%     load(Basefile);
    ProcFiles = dir('*ProcData.mat');
    filenames = {ProcFiles(:).name};
    Basefile = dir('*Baselines.mat');
    load(Basefile.name);
    rcell = cell(1,length(filenames));
    CBVVarCell = cell(1,length(filenames));
    for f = 1:length(filenames)
        clc
        display([num2str(a) ' of ' num2str(length(animals)) ' animals...' num2str(round(f/length(filenames)*100)) '%']);
        % Load individual filenames 
        load(filenames{f});
        [animal,hem,FileDate,FileID] = GetFileInfo(filenames{f});
        
        % Predict the CBV
        strdate = ConvertDate(FileDate);
        CurrentBaseline.(NeurType) = Baselines.(NeurType).(strdate);
        CurrentBaseline.(CBVType) = Baselines.(CBVType).(strdate);
        HRF = HRFStruct.(HRFBeh).GammaHRFs(a,:);
        [CBVPred,~] = CalculatePredictedCBV(ProcData,NeurType,CBVType,HRF,CurrentBaseline);
%         load([animal '_' hem '_' FileID '_CBVPred_' NeurType '.mat']); 
        
        % Filter out data associated with a puff
        AllSol = [0 ProcData.Sol.Contra ProcData.Sol.Ipsi ...
            ProcData.Sol.Tail ProcData.Sol.Control ProcData.TrialDur];
        dSol = diff(AllSol);
        SegAccept = gt(dSol,SegLength+PuffBuffer);
        
        filer = [];
        fileCBVVar = [];
        if any(SegAccept)
            strdate = ConvertDate(FileDate);
            Fs = ProcData.Fs.([CBVType '_fs']);
            [z,p,k] = butter(4,1/(Fs/2),'low');
            [sos,g] = zp2sos(z,p,k);
            eLen = ProcData.TrialDur*Fs;
            
            % Normalize and Filter the CBV
            nCBV = ProcData.(CBVType)(1:eLen)/mean(Baselines.(CBVType).(strdate).Means)-1;
            fCBV = filtfilt(sos,g,detrend(nCBV));
            
            
            segstarts = floor((AllSol([SegAccept false])+PuffBuffer)*Fs);
            filer = zeros(1,length(segstarts));
            fileCBVVar = zeros(1,length(segstarts));
            for seg = 1:length(segstarts) % Loop over each segment of non-sensory evoked data
                seginds = segstarts(seg)+1:segstarts(seg)+(SegLength*Fs);
%                 Pred = CBVPred.(HRFType)(seginds) - mean(CBVPred.(HRFType)(seginds));
                Pred = CBVPred(seginds) - mean(CBVPred(seginds));
                Act = fCBV(seginds)-mean(fCBV(seginds));
                r = corrcoef(Pred,Act);
                filer(seg) = r(2,1);
                fileCBVVar(seg) = var(Act)/BaselineVar;
%                 plot(Act); hold on; plot(Pred); 
%                 title(['Var = ' num2str(fileCBVVar(seg)) ...
%                     '; R = ' num2str(filer(seg))]); hold off;
%                 pause;
            end
        end
        rcell{f} = filer;
        CBVVarCell{f} = fileCBVVar;
    end
    cd(prevdir)
    r = [rcell{:}];
    CBVVar = [CBVVarCell{:}];
    IDArray = cell(size(r));
    IDArray(:) = {animal};
    
    animalRs{a} = r;
    animalVars{a} = CBVVar;
    animalIDArray{a} = IDArray;
    
    % Fit the scatter for each animal with a line
    [b,stats1] = robustfit(CBVVar,r);
    pvals(a) = stats1.p(2);
    slope(a) = b(2);
    intercept(a) = b(1);
    
%     h=figure; 
% %     [N,C] = hist3([CBVVar' r'],[20,20]);
%     scatter(CBVVar,r);
%     hold on; plot(CBVVar,polyval(flipud(b),CBVVar),'k');
%     scatter(MedVar,MedR(a),'k*');
% %     plot(min(CBVVar):0.1:max(CBVVar),tanhfit(tanhCoef,min(CBVVar):0.1:max(CBVVar)),'r');
%     hold off;
%     xlabel('CBV Variance/Baseline Variance')
%     ylabel('R')
%     ylim([-1 1])
%     title([animal ': slope = ' num2str(b(2)) '; p-val = ' num2str(stats1.p(2))])
%     SaveFig2Disk(h,['..\Results\Variance\' HRFType 'RvsCBVVariance']);
%     pause(0.001);
end
ActPredCorrs = [animalRs{:}];
CBVZscore = [animalVars{:}];
AnimalIDs = [animalIDArray{:}];

%% Fit the data with a mixed effects model
% Format data into a table
t = table(AnimalIDs',CBVZscore',ActPredCorrs','VariableNames',{'IDs','Zscore','R'});

% Process the data
t.IDs = categorical(t.IDs);

% Center the CBV Zscore to avoid correlations between slope and intercept
t.ZCentered = t.Zscore - mean(t.Zscore);

% Fit a random slope and intercept model using linear mixed effect models
clc
lme = fitlme(t,'R ~ 1 + ZCentered + (1+ZCentered|IDs)');
Coefficients = lme.Coefficients.Estimate;
correctedIntercept = -1*Coefficients(2)*mean(t.Zscore)+Coefficients(1);
CorrectedCoef = [Coefficients(2), correctedIntercept];
Regline = polyval(CorrectedCoef,0:10,'k');

%% Plot the result of the mixed effect fit
[N,C] = hist3([CBVZscore' ActPredCorrs'],[200,50]);
figure; imagesc(C{1},C{2},N'); axis xy;
axis square;
ylabel('Correlation Coefficient')
xlabel('CBV Var/Resting CBV Var')
caxis([0 10])
xlim([0 10])
hold on;
plot(0:10,Regline,'k')
pval = lme.Coefficients.pValue(2);
tst = lme.Coefficients.tStat(2);
DegF = lme.Coefficients.DF(2);
title(sprintf(['Based on ' NeurType ' neural activity and HRF calculated from '...
    HRFBeh '\nSlope=' num2str(Coefficients(2)) '; p=' num2str(pval) ...
    ', t=' num2str(tst) ', DF=' num2str(DegF)]))
