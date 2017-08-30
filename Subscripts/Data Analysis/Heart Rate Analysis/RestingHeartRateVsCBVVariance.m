function [Stats] = RestingHeartRateVsCBVVariance(animals,CBVType)
%   function [Stats] = RestingHeartRateVsCBVVariance(animals,CBVType)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Calculates the linear mixed effect regression model
%   between the heart rate and CBV variance.
%_______________________________________________________________
%   PARAMETERS:
%               animals - [cell array] animal IDs
%
%               CBVType - [string] designates the CBV ROI to be used
%_______________________________________________________________
%   RETURN:
%               Stats - [struct] results of the the statistical analysis
%_______________________________________________________________

RestCriteria.Fieldname = {'Duration','PuffDistance'};
RestCriteria.Comparison = {'gt','gt'};
RestCriteria.Value = {14,5};
RestBuffer = 4;
figure;
AllCBVVar = cell(1,length(animals));
AllHR = cell(1,length(animals));
AllIDs = cell(1,length(animals));
for a = 1:length(animals)
    animal = animals{a};
    %% Load the data
    prevdir = cd([animal filesep]);
    
    RestFile1 = dir(['*_RESTDATA_' CBVType '.mat']);
    load(RestFile1.name);
    AllRestData1 = RestData;
    DataField1 = fieldnames(RestData);
    DataField1 = DataField1{1};
    
    RestFile2 = dir('*_RESTDATA_HR.mat');
    load(RestFile2.name);
    DataField2 = fieldnames(RestData);
    DataField2 = DataField2{1};
    AllRestData2 = RestData;
    
    %% Throw out periods of rest according to the RestCriteria
    [RestFiltArray1] = FilterEvents(AllRestData1.(DataField1),RestCriteria);
    RestData1 = AllRestData1.(DataField1).NormData(RestFiltArray1);
    
    [RestFiltArray2] = FilterEvents(AllRestData2.(DataField2),RestCriteria);
    RestData2 = AllRestData2.(DataField2).Data(RestFiltArray2);
    
    %% Set up spectral filter for the data
    Fs = AllRestData1.(DataField1).Fs;
    
    CBVvar = zeros(1,length(RestData1));
    meanHR = zeros(1,length(RestData2));
    
    for RD = 1:length(RestData1)
        RestBuffer_Ind = RestBuffer*Fs;
        
        clippedRest1 = RestData1{RD}(RestBuffer_Ind:end);
        clippedRest2 = RestData2{RD}(RestBuffer_Ind:end);
        
        % Find any NaN in the HR data
        NaNind = not(isnan(clippedRest2));
        clippedRest1 = clippedRest1(NaNind);
        clippedRest2 = clippedRest2(NaNind);
        
        RestOffset1 = mean(clippedRest1);
        RestOffset2 = 0;
        
        Rest1 = detrend(clippedRest1-RestOffset1);
        Rest2 = clippedRest2-RestOffset2;
        
        CBVvar(RD) = var(Rest1);
        meanHR(RD) = mean(Rest2);
    end
    IDArray = cell(1,length(RestData1));
    IDArray(:) = {animal};
    AllIDs{a} = IDArray;
    AllCBVVar{a} = CBVvar;
    AllHR{a} = meanHR;
    cd(prevdir)
end
CBVVar = [AllCBVVar{:}];
HR = [AllHR{:}];
IDs = [AllIDs{:}];
scatter(HR,CBVVar,'.');
xlabel('Heart Rate (Hz)')
ylabel(sprintf('R.M.S of Resting\nCBV variance'))

%% Use a mixed effects model to fit the data
% Format the data into a table
t=table(IDs',CBVVar',HR','VariableNames',{'IDs','CBVVar','HR'});

t.IDs = categorical(t.IDs);
t.HRCentered = t.HR-mean(t.HR);
t.CBVCentered = t.CBVVar-mean(t.CBVVar);

clc
lme = fitlme(t,'CBVCentered ~ 1 + HRCentered + (1+HRCentered|IDs)');
Coefficients = lme.Coefficients.Estimate;
correctedIntercept = -1*Coefficients(2)*mean(t.HR)+Coefficients(1);
CorrectedCoef = [Coefficients(2), correctedIntercept];
Regline = polyval(CorrectedCoef,6:15)+mean(t.CBVVar);
hold on;
plot(6:15,Regline,'k')
xlim([6 15]);
ylim([0 7e-4])
Predicted = polyval(CorrectedCoef,t.HR)+mean(t.CBVVar);
Stats.R2 = CalculateRsquared(Predicted,t.CBVVar);
Stats.pval = lme.Coefficients.pValue(2);
Stats.tstat = lme.Coefficients.tStat(2);
Stats.df = lme.Coefficients.DF(2);
title(['Slope: ' num2str(lme.Coefficients.Estimate(2),2) ...
    '; R^2 = ' num2str(round(Stats.R2,3)) ' pvalue: ' ...
    num2str(round(lme.Coefficients.pValue(2),3)) '; t(' ...
    num2str(round(lme.Coefficients.DF(2),2)) ')='...
    num2str(round(lme.Coefficients.tStat(2),3)) '; n='...
    num2str(length(CBVVar))]);