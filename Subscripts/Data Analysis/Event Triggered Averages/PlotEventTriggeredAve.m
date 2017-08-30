function [handles] = PlotEventTriggeredAve(animal,EventData,dataTypes,behaviors,handles)
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

legends.allsessions = []; %to correctly label the plot data
currdir = pwd;

if not(iscell(dataTypes))
    dataTypes = {dataTypes};
end

% Set the colormap for a variable number of behaviors and is robust in case
% a behavior does not occur.
figure(99);
h = plot(ones(10,length(behaviors)));
colors = get(h,'Color'); close(gcf);


for t = 1:length(dataTypes)
    dataType = dataTypes{t};
    for b = 1:length(behaviors)
        Beh = behaviors{b};
        
        % Get data for the behavior
        [DataStruct,FiltArray] = SelectBehavioralEvents(EventData.(dataType),Beh);
        BehData = DataStruct.NormData(FiltArray,:);
        
        % Create time vector for triggered average
        timevec = (0:1/DataStruct.Fs:DataStruct.epoch.duration)-...
            DataStruct.epoch.offset;
           
        sessions = unique(DataStruct.FileDate);
        for s = 1:length(sessions)
            session = sessions{s};
            
            % Identify all events that correspond to the current session
            EvFilter = strcmp(DataStruct.FileDate(FiltArray),session);
            
            % Convert the date to a string for use as a fieldname in
            % structures
            strdate = ConvertDate(session);
            
            % Initiate field of handles structure if necessary
            if isfield(handles,dataType)==0
                handles.(dataType)=[];
            end       

            % Update the legends structure for the current data.
            legends.(strdate){b} = [Beh ' n= ' num2str(sum(EvFilter))];   
            
            % Control in case behavior is not seen in a session. Prevents
            % mislabelling the data
            emptyCells = cellfun(@isempty,legends.(strdate));
            legends.(strdate)(emptyCells) = [];
            l_len = length(legends.(strdate));
            
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
        legends.allsessions{b} = [Beh ' n= ' num2str(size(BehData,1))];
        emptyCells = cellfun(@isempty,legends.allsessions);
        legends.allsessions(emptyCells) = [];
        
        % Plot the combined data from a given behavior and data type
        plot(timevec,mean(BehData,1)-1,'Color',colors{b});
        hold on;
        xlabel('Time (s)');
        ylabel('Normalized Amplitude (A.U)');
        title([animal ': AllSessions: ' dataType]);
        legend(legends.allsessions(1:b),'Location','southeast') %the index (1:b) avoids a warning when legend sizes are greater than plot sizes.
    end
end