function [StimStruct] = SelectSensEvEvents(EventData,dataType,SolName,Criteria)
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

%% Set Threshold Values
Pre_WhiskScore_Max = Criteria.MaxWhiskScore_Pre; % 0.1
Pre_MovementScore_Max = Criteria.MaxMoveScore_Pre; %0.3
Min_PuffDistance = Criteria.MinPuffDistance; 

%% Load and setup
if not(isfield(EventData,dataType))
    error(wraptext(['No stimulation Data for ' dataType ...
        ' found in EventData.mat. Run ExtractEventTriggeredData.m']))
end

SolFilter = strcmp(EventData.(dataType).stim.Name,SolName);
PreWS_Filter = EventData.(dataType).stim.WhiskScore_Pre(:,1)<Pre_WhiskScore_Max;
PreMS_Filter = EventData.(dataType).stim.MoveScore_Pre(:,1)<...
    Pre_MovementScore_Max;
PD_Filter = GetPuffDistanceFilter(EventData.(dataType).stim.PuffDistance,...
    Min_PuffDistance);

StimFilter = logical(SolFilter.*PreWS_Filter.*PreMS_Filter.*PD_Filter);

fnames = fieldnames(EventData.(dataType).stim);
for fn = 1:length(fnames)
    if size(EventData.(dataType).stim.(fnames{fn}),1)==length(StimFilter)
        StimStruct.(fnames{fn}) = EventData.(dataType).stim.(fnames{fn})(StimFilter,:,:);
    else
        StimStruct.(fnames{fn}) = EventData.(dataType).stim.(fnames{fn});
    end
end