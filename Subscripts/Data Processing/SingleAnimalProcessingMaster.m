function [] = SingleAnimalProcessingMaster(animal,hem)
%   function [] = SingleAnimalAnalysisMaster(animal,hem)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Compiles all the standard analysis for a single animal
%   into a single script.
%   
%_______________________________________________________________
%   PARAMETERS:       
%                   animal - [string] animal ID
%
%                   hem - [string] hemisphere recorded
%                               
%_______________________________________________________________
%   RETURN:                 
%                   Nothing returned. Output of the script are plots and
%                   saved files.
%                               
%_______________________________________________________________

%% Setup
NeurTypes = {'Gam','MUpower','Beta','SubAlpha'};
CBVType = 'CrossCorrROI';

%% CATEGORIZE THE DATA ACCORDING TO BEHAVIOR
% Load each file and use the binarized (and linked) whisking to obtain
% information about each whisk. Identify solenoid firing times and use
% tracked whisking and body movement to classify air puffs. Identify
% periods of no detected whisker movement or puffing as periods of rest.
% Save the behavioral categorization as part of the ProcData.mat file

% Setup
ProcDataFileNames = dir('*ProcData.mat');

% Run
for f = 1:size(ProcDataFileNames,1)
    filename = ProcDataFileNames(f).name;
    CategorizeData(filename)
end


%% COMPILE ALL PERIODS OF REST
% Use categorization data added to each ProcData.mat file to extract the
% resting data from each file and save it in a separate structure. 

% Setup
Filenames = dir('*ProcData.mat');
ProcDataFileNames = {Filenames(:).name}';
dataTypes = [NeurTypes {CBVType}];

% Run
for dt = 1:length(dataTypes)
    dataType = dataTypes{dt};
    RestData = struct;  
    display(['Compiling rest data for ' dataTypes{dt}])
    [RestData] = GetAllRestingData(ProcDataFileNames,RestData,dataType);
    save([animal '_' hem '_RESTDATA_' dataType '.mat'],'RestData');
end

%% Calculate a daily resting baseline
% Use periods of rest (> 6 seconds duration, more than 5 seconds from
% puffs) to calculate a baseline value for a given day. 

% Setup
dataTypes = [NeurTypes {CBVType}];

% Run
for dt = 1:length(dataTypes)
    dataType = dataTypes{dt};
    RestDataFileName = dir(['*RESTDATA_' dataType '.mat']);
    load(RestDataFileName.name);
    CalculateRestingBaseline(RestData,dataType,animal,hem)
end

%% Normalize the resting periods to get a percent change from baseline

% Setup
dataTypes = [NeurTypes {CBVType}];

% Run
for dt = 1:length(dataTypes)
    dataType = dataTypes{dt};
    RestDataFileName = dir(['*RESTDATA_' dataType '.mat']);
    load(RestDataFileName.name);
    [RestData.(dataType).NormData, RestData.(dataType).ZNormData] = ...
        NormBehavioralDataStruct(RestData.(dataType),dataType);
    RestData.(dataType).Normalized = 1;
    save(RestDataFileName.name,'RestData');
end

%% COMPILE THE DATA SURROUNDING BEHAVIORAL EVENTS

% Setup
epoch.duration = 10; % seconds
epoch.offset = 4; % seconds
FileNames = dir('*ProcData.mat');
ProcDataFileNames = {FileNames(:).name}';
dataTypes = [NeurTypes {CBVType}];

% Run
for dt = 1:length(dataTypes)
    dataType = dataTypes{dt};
    EventData = struct;
    [EventData] = ExtractEventTriggeredData(ProcDataFileNames,EventData,dataType,epoch);
    save([animal '_' hem '_EVENTDATA_' dataType '.mat'],'EventData');
end

%% Normalize the behavioral data
% Setup
dataTypes = [NeurTypes {CBVType}];

% Run
for dt = 1:length(dataTypes)
    dataType = dataTypes{dt};
    EventDataFileName = dir(['*EVENTDATA_' dataType '.mat']);
    if not(isempty(EventDataFileName))
        display('EventData.mat structure already exists...loading.')
        load(EventDataFileName.name);
    else
        error('Cannot find Event Data in the current folder...')
    end
    
    EventTypes = fieldnames(EventData.(dataType));
    % Normalize the data
    for E = 1:length(EventTypes)
        [EventData.(dataType).(EventTypes{E}).NormData, EventData.(dataType).(EventTypes{E}).ZNormData] = ...
            NormBehavioralDataStruct(EventData.(dataType).(EventTypes{E}),dataType);
        EventData.(dataType).(EventTypes{E}).Normalized = 1;
    end
    
    % Save the structure
    save([animal '_' hem '_EVENTDATA_' dataType '.mat'],'EventData');
end