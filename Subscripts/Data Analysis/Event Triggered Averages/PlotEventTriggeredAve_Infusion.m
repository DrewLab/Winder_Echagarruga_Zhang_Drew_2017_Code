function PlotEventTriggeredAve_Infusion(EventData,dataType,Beh,...
    AverageThresh,colorcode)
%   function [] = PlotEventTriggeredAve_Infusion(EventData,dataType,Beh,...
%   AverageThresh,colorcode)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Plots the event triggered waveforms for each recording
%   session and for all sessions. Each plot will show an overlay of the
%   event triggered behavioral data for easy comparison in each imaging
%   session and for all data.
%_______________________________________________________________
%   PARAMETERS: 
%               animal - [string] animal ID
%
%               EventData - [struct] Info and data from around behavioral
%               events
%
%               dataTypes - [cell] data measures to plot
%
%               behaviors - [cell] behavior which corresponds to data
%_______________________________________________________________
%   RETURN: The output for this function is a set of figures showing the                  
% event triggered averages. Figures will be automatically saved.
%_______________________________________________________________

%% Setup
  
switch Beh
    case 'Contra'
        fdname = 'stim';
    case 'VW'
        fdname = 'whisk';
end

% Get data for the behavior
Criteria.Fieldname = {'Name','WhiskScore_Pre', 'MoveScore_Pre','PuffDistance'};
Criteria.Comparison = {'equal','lt','lt','gt'};
Criteria.Value = {'Contra',0.1, 0.3, 5};
DataStruct = EventData.(dataType).(fdname);
FiltArray = FilterEvents(DataStruct,Criteria);

% [~,FiltArray,Criteria] = SelectBehavioralEvents(EventData.(dataType),Beh);
FileDates = EventData.(dataType).(fdname).FileDate(FiltArray);
if strcmp(dataType,'RadiusROI')
    Data = EventData.(dataType).(fdname).NormData(FiltArray,:);
    PreInds = 1:EventData.(dataType).(fdname).epoch.offset*EventData.(dataType).(fdname).Fs;
    DC = mean(Data(:,PreInds),2)*ones(1,size(Data,2));
    Data = Data - DC;
else
    Data = EventData.(dataType).(fdname).Data(FiltArray,:);
end
MinutesSinceInfusion = EventData.(dataType).(fdname).MinutesSinceInfusion(FiltArray);

% Create time vector for triggered average
timevec = (0:1/EventData.(dataType).(fdname).Fs:EventData.(dataType).(fdname).epoch.duration)-...
    EventData.(dataType).(fdname).epoch.offset;

AllTimeFilter = false(size(MinutesSinceInfusion));
sessions = unique(EventData.(dataType).(fdname).FileDate(FiltArray));
for s = 1:length(sessions)
    session = sessions{s};
    
    % Identify all events that correspond to the current session
    EvFilter = strcmp(FileDates,session);
    
    % Convert the date to a string for use as a fieldname in
    % structures
    strdate = ConvertDate(session);
    
    % Identify all events that occur after the time threshold for
    % averaging
    if not(isfield(AverageThresh,strdate))
        display(['MUpower was never reduced below 0.2 for ' strdate])
        handles.(dataType)=[];
        continue;
    end
    TimeFilter = gt(MinutesSinceInfusion, AverageThresh.(strdate));
    % Track the events which are used for averaging
    AllTimeFilter = or(AllTimeFilter,TimeFilter);
end

plot(timevec,mean(Data(AllTimeFilter,:),1),'Color',colorcode);
xlabel('Time (s)');
ylabel('Normalized Amplitude (A.U)');
