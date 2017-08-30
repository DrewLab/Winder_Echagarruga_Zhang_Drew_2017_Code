function [XC, Reshuf95, lags] = RestingXCorr(animals,dataType1,dataType2)
%   [XC, Reshuf95, lags] = RestingXCorr(animals,dataType1,dataType2)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Calculates the cross correlation between dataType1 and
%   dataType2 during rest.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   animals - [cell array] animal IDs
%
%                   dataType1 - [string] fieldname of the RestData
%                   structure.
%
%                   dataType2 - [string] fieldname of the RestData
%                   structure.
%
%                       RestData structures were calculated from the
%                       function 'SingleAnimalProcessingMaster.m'
%_______________________________________________________________
%   RETURN:                     
%                   XC - [matrix] animals cross correlations, rows contain 
%                   the cross correlation for each animal. Columns
%                   correspond to correlation coefficients at times found
%                   in lags.
%
%                   Reshuf95 - [matrix] the [lower, upper] confidence
%                   intervals for each animal
%
%                   lags - [array] time vector containing the time lag of
%                   each correlation coefficient.
%_______________________________________________________________

%% Setup
RestCriteria.Fieldname = {'Duration','PuffDistance'};
RestCriteria.Comparison = {'gt','gt'};
RestCriteria.Value = {14,5};
RestBuffer = 4;
maxlags = 5;

for a = 1:length(animals)
    display([num2str(a) ' of ' num2str(length(animals)) ' animals...'])
    animal = animals{a};
    %% Load the data
    prevdir = cd([animal filesep]);
    
%     RestFile1 = ls(['*_RESTDATA_' dataType1 '.mat']);
%     load(RestFile1);
    RestFile1 = dir(['*_RESTDATA_' dataType1 '.mat']);
    load(RestFile1.name);
    AllRestData1 = RestData;
    DataField1 = fieldnames(RestData);
    DataField1 = DataField1{1};
    
%     RestFile2 = ls(['*_RESTDATA_' dataType2 '.mat']);
%     load(RestFile2);
    RestFile2 = dir(['*_RESTDATA_' dataType2 '.mat']);
    load(RestFile2.name);
    DataField2 = fieldnames(RestData);
    DataField2 = DataField2{1};
    AllRestData2 = RestData;
    
    %% Filter periods of rest according to the RestCriteria
    [RestFiltArray1] = FilterEvents(AllRestData1.(DataField1),RestCriteria);
    RestData1 = AllRestData1.(DataField1).Data(RestFiltArray1);
    
    [RestFiltArray2] = FilterEvents(AllRestData2.(DataField2),RestCriteria);
    RestData2 = AllRestData2.(DataField2).Data(RestFiltArray2);
    
    %% Set up spectral filter for the data
    [z,p,k] = butter(4,1/(AllRestData1.(DataField1).Fs/2),'low');
    [sos,g] = zp2sos(z,p,k);
    Fs = AllRestData1.(DataField1).Fs;
    
    %% Concatenate the resting events
    Rest1 = cell(1,length(RestData1));
    Rest2 = cell(1,length(RestData1));
    for RD = 1:length(RestData1)
        RestBuffer_Ind = RestBuffer*Fs;
        
        clippedRest1 = RestData1{RD}(RestBuffer_Ind:end);
        clippedRest2 = RestData2{RD}(RestBuffer_Ind:end);
        
        % Find any NaN in the data (for HR data)
        NaNind = not(isnan(clippedRest2));
        clippedRest1 = clippedRest1(NaNind);
        clippedRest2 = clippedRest2(NaNind);
        
        RestOffset1 = mean(clippedRest1);
        RestOffset2 = mean(clippedRest2);
        
        Rest1{RD} = filtfilt(sos,g,detrend(clippedRest1-RestOffset1));
        Rest2{RD} = filtfilt(sos,g,detrend(clippedRest2-RestOffset2));
    end
    numrest(a) = length(Rest1);
    JoinedRest1 = [Rest1{:}];
    JoinedRest2 = [Rest2{:}];
    
    %% Filter the concatenated events
%     FiltRest1 = filtfilt(sos,g,JoinedRest1);
%     FiltRest2 = filtfilt(sos,g,JoinedRest2);
    
    %% Calculate the xcorr for the filtered events
    [XC(a,:),lags] = xcorr(JoinedRest2,JoinedRest1,maxlags*Fs,'coeff');

    %% Calculate the 95% confidence interval on reshuffled data
    repetitions = 1000;
    CC = zeros(1,repetitions);
    for rep = 1:repetitions
        Random_inds = ceil(length(Rest1)*...
            rand(1,length(Rest1)));
        Reshuffled = randperm(length(Random_inds));
        Reshuf1 = [Rest1{Random_inds}];
        Reshuf2 = [Rest2{Random_inds(Reshuffled)}];
        
%         FiltReshuf1 = filtfilt(sos,g,Reshuf1);
%         FiltReshuf2 = filtfilt(sos,g,Reshuf2);
        
        r = corrcoef(Reshuf2,Reshuf1);
        CC(rep) = r(2,1);
    end
    
    Reshuf95(a,1) = quantile(CC,0.025);
    Reshuf95(a,2) = quantile(CC,0.975);
    cd(prevdir)
end
lags = lags/Fs;

