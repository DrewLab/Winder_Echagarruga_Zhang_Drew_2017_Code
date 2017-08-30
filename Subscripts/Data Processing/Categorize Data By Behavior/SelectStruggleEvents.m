function [StrugStruct] = SelectStruggleEvents(EventData,dataType,Criteria)
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

Min_Duration = Criteria.Min_Dur; % Seconds, this separates a volitional whisk from a struggle
Min_RestTime = Criteria.Min_RestTime; % Seconds, minimum resting duration before start of event
Min_PuffDistance = Criteria.Min_PuffDistance; % Seconds, minimum distance to puff.

%% Load and setup

if not(isfield(EventData,dataType))
    error(sprintf(['No whisking Data for ' dataType ...
        ' found in RestData.mat. \nRun ExtractEventTriggeredData.m for dataType'])) 
end

Dur_Filter = EventData.(dataType).whisk.Duration>Min_Duration;
RT_Filter = EventData.(dataType).whisk.RestTime>Min_RestTime;
PD_Filter = GetPuffDistanceFilter(EventData.(dataType).whisk.PuffDistance,...
    Min_PuffDistance);

StrugFilter = logical(Dur_Filter.*RT_Filter.*PD_Filter);

fnames = fieldnames(EventData.(dataType).whisk);
for fn = 1:length(fnames)
    if size(EventData.(dataType).whisk.(fnames{fn}),1)==length(StrugFilter)
        StrugStruct.(fnames{fn}) = EventData.(dataType).whisk.(fnames{fn})(StrugFilter,:,:);
    else
        StrugStruct.(fnames{fn}) = EventData.(dataType).whisk.(fnames{fn});
    end
end