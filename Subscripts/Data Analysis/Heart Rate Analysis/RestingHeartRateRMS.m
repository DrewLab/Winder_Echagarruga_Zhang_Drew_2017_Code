function [] = RestingHeartRateRMS(animals)
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

RestCriteria.Fieldname = {'Duration','PuffDistance'};
RestCriteria.Comparison = {'gt','gt'};
RestCriteria.Value = {14,5};
RestBuffer = 4;

AllHR = NaN*ones(1,length(animals));
AllIDs = cell(1,length(animals));
for a = 1:length(animals)
    animal = animals{a};
    %% Load the data
    prevdir = cd([animal filesep]);
    
%     RestFile = ls('*_RESTDATA_HR.mat');
%     load(RestFile);
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
    
    HR_RMS = zeros(1,length(RestData));
    
    for RD = 1:length(RestData)
        RestBuffer_Ind = RestBuffer*Fs;
        clippedRest = RestData{RD}(RestBuffer_Ind:end);
        
        % Find any NaN in the HR data
        NaNind = not(isnan(clippedRest));
        clippedRest = clippedRest(NaNind);
        
        RestOffset = 0;
        
        Rest2 = clippedRest-RestOffset;
        
        HR_RMS(RD) = std(Rest2);
    end
    IDArray = cell(1,length(RestData));
    IDArray(:) = {animal};
    AllIDs{a} = IDArray;
    AllHR(a) = mean(HR_RMS);
    cd(prevdir)
end

figure; 
scatter(ones(size(AllHR)),AllHR,'MarkerEdgeColor','k');
ylim([0 1]);
ax = gca;
set(ax,'XTick',[]);
ylabel(sprintf('R.M.S. of Heart Rate\nat rest (Hz)'))