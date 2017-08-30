function [] = RestingHeartRate_LombSpectrum(animals)
%   function [] = RestingHeartRate_LombSpectrum(animals)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Calculates the power spectrum of the heart rate at rest
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   animals - [cell array] animal IDs                       
%_______________________________________________________________
%   RETURN:                     
%                               
%_______________________________________________________________

%% Setup
RestCriteria.Fieldname = {'Duration','PuffDistance'};
RestCriteria.Comparison = {'gt','gt'};
RestCriteria.Value = {14,5};
RestBuffer = 2;
AllRest = [];
for a = 1:length(animals)
    animal = animals{a};
    %% Load the data
    prevdir = cd([animal filesep]);
    RestFile = dir('*_RESTDATA_HR.mat');
    load(RestFile.name);
    DataField = fieldnames(RestData);
    DataField = DataField{1};
    AllRestData = RestData;
    
    %% Filter out periods of rest according to the RestCriteria    
    [RestFiltArray] = FilterEvents(AllRestData.(DataField),RestCriteria);
    RestData = AllRestData.(DataField).Data(RestFiltArray);
    
    %% Set up spectral filter for the data
    Fs = AllRestData.(DataField).Fs;
    
    %% Concatenate the resting events
    Rest = cell(1,length(RestData));
    for RD = 1:length(RestData)
        RestBuffer_Ind = RestBuffer*Fs;
        clippedRest = RestData{RD}(RestBuffer_Ind:end);
        
        % Find any NaN in the HR data
        NaNind = not(isnan(clippedRest));
        clippedRest = clippedRest(NaNind);
        RestOffset = mean(clippedRest);
        
        % Filter the resting heart rate
        Rest{RD} = [detrend(clippedRest-RestOffset) NaN];
    end
    JoinedRest = [Rest{:}];
    AllRest = [AllRest JoinedRest];
    cd(prevdir);
end

[pxx,f] = plomb(AllRest,Fs);
figure; plot(f,pxx);
xlim([0 0.5])
xlabel('Frequency (Hz)')
ylabel('Power_{Heart Rate}')