function [EventData] = PlotEventTriggeredAve_Infusion_NormByPreInfusion(animal,EventData,...
    dataTypes,behaviors,AverageThresh,colorcode)
%   function [] = PlotEventTriggeredAve(animal,EventData,dataTypes,behaviors)
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

handles = struct(); %to aid in organizing and automating plots
legends.allsessions = []; %to correctly label the plot data
currdir = pwd;

if not(iscell(dataTypes))
    dataTypes = {dataTypes};
end

% Set the colormap for a variable number of behaviors and is robust in case
% a behavior does not occur.
h = plot(ones(10,length(behaviors)));
colors = get(h,'Color'); close(gcf);


for t = 1:length(dataTypes)
    dataType = dataTypes{t};
    for b = 1:length(behaviors)
        Beh = behaviors{b};
        
        switch Beh
            case 'Contra'
                fdname = 'stim';
            case 'VW'
                fdname = 'whisk';
        end
        % Get data for the behavior
        [~,FiltArray,Criteria] = SelectBehavioralEvents(EventData.(dataType),Beh);
        FileDates = EventData.(dataType).(fdname).FileDate(FiltArray);
        Data = EventData.(dataType).(fdname).NormData(FiltArray,:);
        MinutesSinceInfusion = EventData.(dataType).(fdname).MinutesSinceInfusion(FiltArray);
        
        % Create time vector for triggered average
        timevec = (0:1/EventData.(dataType).(fdname).Fs:EventData.(dataType).(fdname).epoch.duration)-...
            EventData.(dataType).(fdname).epoch.offset;
           
        % Loop over each session, get events after neural suppression
        AllPostTimeFilter = false(size(MinutesSinceInfusion));
        AllPreTimeFilter = false(size(MinutesSinceInfusion));
        sessions = unique(EventData.(dataType).(fdname).FileDate(FiltArray));
        for s = 1:length(sessions)
            session = sessions{s};
            
            % Identify all events that correspond to the current session
            SessionFilter = strcmp(FileDates,session);
            
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
            PostTimeFilter = gt(MinutesSinceInfusion, AverageThresh.(strdate));
            PreTimeFilter = lt(MinutesSinceInfusion,0);
            % Track the events which are used for averaging
            AllPostTimeFilter = or(AllPostTimeFilter,PostTimeFilter);
            AllPreTimeFilter = or(AllPreTimeFilter,PreTimeFilter);
            
            % Initiate field of handles structure if necessary
            if isfield(handles,dataType)==0
                handles.(dataType)=[];
            end
            
            % Create a new figure or change to figure corresponding to the
            % datatype and session using the handles structure
            if isfield(handles.(dataType),strdate)
                figure(handles.(dataType).(strdate))
            else
                handles.(dataType).(strdate) = figure;
            end            

            % Update the legends structure for the current data.
            legends.(strdate){b} = [Beh ' n= ' num2str(sum(SessionFilter&PostTimeFilter))];   
            
            % Control in case behavior is not seen in a session. Prevents
            % mislabelling the data
            emptyCells = cellfun(@isempty,legends.(strdate));
            legends.(strdate)(emptyCells) = [];
            l_len = length(legends.(strdate));
            
            % Plot the data for the individual session. Each plot will show
            % data from all behaviors.
            EventData.(dataType).(fdname).Averages.(strdate).Post = mean(Data(SessionFilter&PostTimeFilter,:),1);
            EventData.(dataType).(fdname).Averages.(strdate).Pre = mean(Data(SessionFilter&PreTimeFilter,:),1);
            plot(timevec,mean(Data(SessionFilter&PostTimeFilter,:),1),'Color',colorcode);
            hold on;
            xlabel('Time (s)');
            ylabel('Normalized Amplitude (A.U)');
            title([animal ': ' dataType ' _ ' session])
            legend(legends.(strdate)(1:l_len),'Location','southeast') %the index (1:b) avoids a warning when legend sizes are greater than plot sizes.
        end
        
        % Plot the maximum amplitudes over time
        if or(strcmp(dataType,'Gam'),strcmp(dataType,'MUpower'))
            Etime_Inds = timevec<0.5&timevec>-0.1;
            EvalFun = @(x)max(x,[],2);
        elseif strcmp(dataType,'RadiusROI')
            Etime_Inds = timevec<2.5&timevec>-0;
            EvalFun = @(x)min(x,[],2);
        end
        FiltData = Data(AllPostTimeFilter,Etime_Inds);
        IndAmps = EvalFun(FiltData-FiltData(:,1)*ones(1,size(FiltData,2)));
        figure; scatter(MinutesSinceInfusion(AllPostTimeFilter),IndAmps);
        close gcf;
        
        if isempty(handles)
            EventData.(dataType).(fdname).Averages.AllSessions = [];
            EventData.(dataType).(fdname).Averages.TimeVec = timevec;
            EventData.(dataType).(fdname).Averages.Criteria = Criteria;
            return
        end
        % Combine all the data into a cumulative plot for the data type
        
        % Create a new figure or change to figure corresponding to the
        % datatype using the handles structure
        if isfield(handles.(dataType),'alldays')
            figure(handles.(dataType).alldays)
        else
            handles.(dataType).alldays = figure;
        end
        
        % Build a legend
        legends.allsessions{b} = [Beh ' n= ' num2str(sum(AllPostTimeFilter))];
        emptyCells = cellfun(@isempty,legends.allsessions);
        legends.allsessions(emptyCells) = [];
        
        % Plot the combined data from a given behavior and data type
        EventData.(dataType).(fdname).Averages.AllSessions.Post = mean(Data(AllPostTimeFilter,:),1);
        EventData.(dataType).(fdname).Averages.AllSessions.Pre = mean(Data(AllPreTimeFilter,:),1);
        EventData.(dataType).(fdname).Averages.TimeVec = timevec;
        EventData.(dataType).(fdname).Averages.Criteria = Criteria;
        
        plot(timevec,mean(Data(AllPostTimeFilter,:),1),'Color',colorcode);
        hold on;
        xlabel('Time (s)');
        ylabel('Normalized Amplitude (A.U)');
        title([animal ': AllSessions: ' dataType]);
        legend(legends.allsessions(1:b),'Location','southeast') %the index (1:b) avoids a warning when legend sizes are greater than plot sizes.
    end
end