function [Elapsed] = ElapsedTime(FileTime,EventTime,StartTime)
%   function [Elapsed] = ElapsedTime(FileTime,EventTime,StartTime)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Calculates the time since a user-defined point until an
%   event time of a given trial.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   FileTime - [cell array] a cell array of fileIDs as
%                   output by the LabView acquisition program 
%                       [YYMMDD_HH_mm_SS]
%
%                   EventTime - [array, double] time in seconds since the
%                   beginning of the trial as a column vector
%
%                   StartTime - 
%_______________________________________________________________
%   RETURN:                     
%                   Elapsed - [array, double] time in minutes between the 
%                   event and the user defined file start.
%_______________________________________________________________

if isempty(FileTime)
    Elapsed = [];
    return
end

if ~iscell(FileTime)
    error('FileTime must be a cell array')
end

if nargin < 2 || isempty(EventTime)
    EventTime = zeros(size(FileTime,1),1);
end

if nargin < 3
    UniqueFiles = unique(FileTime);
    Selection = listdlg('ListString',UniqueFiles,'Name',...
        'Select Infusion time: ');
    StartTime = UniqueFiles{Selection};
end
StartTime_Tvec = ConvertTime(StartTime);

AllFiles_Tvec = ConvertTime(FileTime);
EventTime = round(EventTime);
EventTimes_Tvec = [zeros(length(FileTime),5) EventTime];
Combined_Tvec = AllFiles_Tvec+EventTimes_Tvec;
Elapsed = zeros(size(Combined_Tvec,1),1);
for CT = 1:size(Combined_Tvec,1)
    Elapsed(CT,:) = etime(Combined_Tvec(CT,:),StartTime_Tvec)/60; % Minutes
end