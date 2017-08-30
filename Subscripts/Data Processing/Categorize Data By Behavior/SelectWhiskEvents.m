function [WhiskStruct] = SelectWhiskEvents(EventData,dataType,Criteria)
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

%% Set threshold values

Max_Dur = Criteria.MaxDur; % Seconds, this separates a volitional whisk from a struggle
Min_Dur = Criteria.MinDur; % Seconds, this is the minimum time for a movement to be considered a whisk
Min_RestTime = Criteria.Min_RestTime; % Seconds, minimum resting duration before start of event
Min_PuffDist = Criteria.Min_PuffDist; % Seconds, minimum time between rest and any puff.
Min_WhiskScore = Criteria.WhiskScore;

%% Load and setup

if not(isfield(EventData,dataType))
    error(sprintf(['No whisking Data for ' dataType ...
        ' found in RestData.mat. \nRun ExtractEventTriggeredData.m for dataType'])) 
end

Dur_Filter = and(EventData.(dataType).whisk.Duration<Max_Dur,...
    EventData.(dataType).whisk.Duration>Min_Dur);
RT_Filter = EventData.(dataType).whisk.RestTime>Min_RestTime;
PD_Filter = GetPuffDistanceFilter(EventData.(dataType).whisk.PuffDistance,...
    Min_PuffDist);
WS_Filter = EventData.(dataType).whisk.WhiskScore>Min_WhiskScore;


WhiskFilter = logical(Dur_Filter.*RT_Filter.*PD_Filter.*WS_Filter);

fnames = fieldnames(EventData.(dataType).whisk);
for fn = 1:length(fnames)
    if size(EventData.(dataType).whisk.(fnames{fn}),1)==length(WhiskFilter)
        WhiskStruct.(fnames{fn}) = EventData.(dataType).whisk.(fnames{fn})(WhiskFilter,:,:);
    else
        WhiskStruct.(fnames{fn}) = EventData.(dataType).whisk.(fnames{fn});
    end
end