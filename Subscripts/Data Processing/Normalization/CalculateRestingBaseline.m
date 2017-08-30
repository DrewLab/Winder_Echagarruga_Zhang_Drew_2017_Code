function [] = CalculateRestingBaseline(RestData,dataType,animal,hem)
%   function [] = CalculateRestingBaseline(RestData,dataType,animal,hem)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Uses periods when animal is not being stimulated or moving
%   to establish a baseline for a given session of imaging.
%   
%_______________________________________________________________
%   PARAMETERS:    
%                   RestData - [Struct] contains all resting data to be
%                   used for the baseline.
%
%                   dataType = [string] data measurement desired for
%                   baseline
%
%                   animal = [string] animal ID
%
%                   hem = [string] hemisphere of animal recorded
%                           ('LH,'RH')
%                               
%_______________________________________________________________
%   RETURN:                     
%                   none, result of script is a saved file                   
%_______________________________________________________________

%% Load and setup
% Load previous Baselines
BaseFile = dir('*Baselines.mat');
if not(isempty(BaseFile.name))
    wraptext('CalculateRestingBaselines.m: Baselines file found...loading...');
    load(BaseFile.name);
else
    wraptext('CalculateRestingBaselines.m: No baselines file found...creating...');
    Baselines = struct();
end

% Check for "dataType" in RestData
if not(isfield(RestData,dataType))
    ErrorText = ['No Resting Data for ' dataType ...
        ' found in RestData.mat. Run GetAllRestingData.m for dataType'];
    error(ErrorText) 
end
% Determine if a baseline will be overwritten
if isfield(Baselines,dataType)
    wraptext(['Baselines for ' dataType ' already exist...overwriting'])
end

% Initialize dataType sub-structure
Baselines.(dataType) = [];

% Define Variables
target_duration = 6;
Min_PuffDistance = 5; % seconds between rest event and puff

[sessions,s_id] = GetUniqueDays(RestData.(dataType).FileID);
session_inds = [s_id' length(RestData.(dataType).FileID)];

% Get RestData Values for dataType
EventTimes = RestData.(dataType).EventTime;
PuffDistances = RestData.(dataType).PuffDistance;
Durations = RestData.(dataType).Duration;
    
for s = 1:length(sessions)
    
    % Convert session date to a string for use as a fieldname of structure
    session = sessions{s};
    fdname = ConvertDate(session);
    
    % Initialize fdname sub-structure
    Baselines.(dataType).(fdname) = [];
    
    % Identify periods corresponding to the session
    Session_filter = zeros(size(EventTimes));
    Session_filter(session_inds(s):session_inds(s+1)) = 1;
    
    % Identify periods isolated from puffs
    [PD_filter] = GetPuffDistanceFilter(PuffDistances,Min_PuffDistance);
    
    % Identify periods of sufficient length
    Len_filter = Durations > target_duration;
    
    % Combine filters to create final filter
    Baseline_filter = logical(Session_filter.*PD_filter.*Len_filter);
    
    % Keep values used to calculate baseline with the Baseline structure.
    Baselines.(dataType).(fdname).Vals = ...
        RestData.(dataType).Data(Baseline_filter);
    
    % Calculate the mean from each rest period, convert from cell array to
    % scalar array.
    % dimarray -> dimension of mean and variance functions for cellfun
    % syntax
    dimarray = mat2cell(2*ones(sum(Baseline_filter),1),...
        ones(sum(Baseline_filter),1));  
    BaseMeans = cellfun(@mean,RestData.(dataType).Data(Baseline_filter),...
        dimarray,'UniformOutput',0);
    Baselines.(dataType).(fdname).Means = [BaseMeans{:}]';
    
    % Calculate variance from each rest period, convert from cell array to
    % scalar array.
    BaseData = RestData.(dataType).Data(Baseline_filter);
%     Basemean = mean(Baselines.(dataType).(fdname).Means);
    for c = 1:length(BaseData)
%         NBase = BaseData{c}./(Basemean'*ones(1,size(BaseData{c},2)));
%         Baselines.(dataType).(fdname).StDev(c,:,:) = std(NBase,[],2);
        Baselines.(dataType).(fdname).StDev(c,:,:) = std(detrend(BaseData{c}),[],2);
    end
    
    % Save duration of rest period and distance of period from puff
    Baselines.(dataType).(fdname).Duration = ...
        RestData.(dataType).Duration(Baseline_filter);
    Baselines.(dataType).(fdname).PuffDistance = ...
        RestData.(dataType).PuffDistance(Baseline_filter);
    Baselines.(dataType).(fdname).FileID = ...
        RestData.(dataType).FileID(Baseline_filter);
end

save([animal '_' hem '_Baselines.mat'],'Baselines')
    
end
