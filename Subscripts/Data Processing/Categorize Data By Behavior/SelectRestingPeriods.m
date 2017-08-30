function [RestStruct] = SelectRestingPeriods(RestData,dataType,Criteria)
%   function [RestStruct] = SelectRestingPeriods(RestData,dataType,...
%       Min_Duration,Min_PuffDistance)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Selects periods of rest that satisfy duration and time
%   separation from puff conditions.
%   
%_______________________________________________________________
%   PARAMETERS:   
%                   RestData - [Struct] contains data for all periods of
%                   rest.
%
%                   dataType - [cell array] field names of Rest Data to be
%                   extracted
%
%                   Criteria - [struct] values used to identify the desired
%                   rest events:
%                       Min_Duration - [double] the minimum duration, in
%                       seconds, for accepted rest periods
%
%                       Min_PuffDist - [double] the minimum time between a 
%                       puff and accepted rest periods
%                               
%_______________________________________________________________
%   RETURN:                     
%                   RestStruct - [struct] contains data for resting periods
%                   that satisfy the criteria.
%_______________________________________________________________

%% Set thresholds for a period of rest
Min_Dur = Criteria.Min_Duration; % seconds
Min_PuffDist = Criteria.Min_PuffDist; % seconds

%% Load and setup
if not(isfield(RestData,dataType))
    error(sprintf(['No Resting Data for ' dataType ...
        ' found in RestData.mat. \nRun GetAllRestingData.m for dataType'])) 
end

Dur_Filter = RestData.(dataType).Duration > Min_Dur;
PD_Filter = GetPuffDistanceFilter(RestData.(dataType).PuffDistance,...
    Min_PuffDist);

RestFilter = logical(Dur_Filter.*PD_Filter);

fnames = fieldnames(RestData.(dataType));
for fn = 1:length(fnames)
    if size(RestData.(dataType).(fnames{fn}),1)==length(RestFilter)
        RestStruct.(fnames{fn}) = RestData.(dataType).(fnames{fn})(RestFilter);
    else
        RestStruct.(fnames{fn}) = RestData.(dataType).(fnames{fn});
    end
end
