function [Norm] = Infusion_NormalizeData(Raw,FileDates,Baselines,dataType)
%   function [] =  ()
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Normalizes all data to pre-infusion amplitudes for
%   comparison across infusion conditions.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                               
%_______________________________________________________________
%   RETURN:                     
%                               
%_______________________________________________________________

Norm = Raw;
matDates = cell2mat(FileDates);
strdates = ConvertDate(matDates);
% Determine how many sessions occurred with resting data
UniqueRestingDates = fieldnames(Baselines.(dataType));

% Loop over each session to get the periods of rest which occur after
% infusion
for URD = 1:length(UniqueRestingDates)   
    Session_RestIDs = Baselines.(dataType).(UniqueRestingDates{URD}).FileID;
    MinutesSinceInfusion = Baselines.(dataType).(UniqueRestingDates{URD}).MinutesSinceInfusion;
    EventTimes = Baselines.(dataType).(UniqueRestingDates{URD}).RestingTime;
    
    % Filter all events that occur before infusion
    Baseline = Baselines.(dataType).(UniqueRestingDates{URD}).PreInfusionMeans;
    if Baseline==0
        error('No resting data before infusion')
    end
    if iscell(Raw)
        for R = 1:length(Raw)
            RawSnip = Raw{R};
            Norm{R} = RawSnip/Baseline-1;
        end
    else
        Basemat = Baseline*ones(size(Raw));
        Norm = Raw./Basemat-1;
    end
end