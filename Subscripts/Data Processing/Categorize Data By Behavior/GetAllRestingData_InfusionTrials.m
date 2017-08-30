function [RestData] = GetAllRestingData_InfusionTrials(filenames,RestData,dataTypes,Thresholds)
%   [RestData] = GetAllRestingData_InfusionTrials(filenames,RestData,dataTypes)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Extracts data which correspond to periods of rest from
%   various data structures and compiles them into a single structure. The
%   duration of the resting period, time to nearest puff, and details of
%   the trial from which the data were taken are also added to the
%   structure. If resting data already exists for an element of dataTypes,
%   the resting data for that dataType is replaced. If no resting data
%   exists for an element of dataTypes, the resting data for that element
%   is added to the existing structure.
%
%_______________________________________________________________
%   PARAMETERS:
%                   filenames - [matrix] all filenames to be analyzed
%                   where filenames are organized as rows. Filenames can be
%                   of various types of data as long as the first set of
%                   fields contain fieldnames which correspond to elements
%                   of "dataTypes". Additionally, the structure must
%                   contain a field called "Flags" which identifies
%                   behaviors from the filename. The "Flags" field is
%                   obtained by running the script "CategorizeData.m".
%
%                   RestData - [Structure] Structure containing existing
%                   rest data
%
%                   dataTypes - [cell array] fields of the data structure
%                   to be extracted into the structure of resting data.
%_______________________________________________________________
%   RETURN:
%                   RestData - [structure] contains all resting data and
%                   information about each period of rest.
%_______________________________________________________________

%% Loop structure dataTypes>>filenames>>individual events
% (Speed up execution by putting dataTypes as the inner loop)

% Control for string input with length(dataTypes)=1
if not(iscell(dataTypes))
    dataTypes = {dataTypes};
end

for dt = 1:length(dataTypes)
    
    % Initialize cell arrays for resting data and other information.
    RestVals = cell(size(filenames,1),1);
    EventTimes = cell(size(filenames,1),1);
    Durations = cell(size(filenames,1),1);
    PuffDistances = cell(size(filenames,1),1);
    FileIDs = cell(size(filenames,1),1);
    FileDates = cell(size(filenames,1),1);
    
    for f = 1:length(filenames)
        % Load ProcData.mat file, call the loaded structure 
        filename = filenames{f};
        FileData = load(filename);
        fname1 = fieldnames(FileData);
        Data = FileData.(fname1{1});
        
        % Get the date and file identifier for the data to be saved with
        % each resting event
        [~,~,FileDate,FileID] = GetFileInfo(filename);
        
        % Sampling frequency for element of dataTypes
        Fs = Data.Fs.(dataTypes{dt});
        
        % Expected number of samples for element of dataType
        Expected_length = Data.Info.TrialDur*Fs;
        
        % Get information about periods of rest from the loaded file
        Trial_EventTimes = Data.Flags.rest.EventTime';
        Trial_PuffDistances = Data.Flags.rest.PuffDistance;
        Trial_Durations = Data.Flags.rest.Duration';
        
        % Initialize cell array for all periods of rest from the loaded
        % file
        Trial_RestVals = cell(size(Trial_EventTimes'));
        for ET = 1:length(Trial_EventTimes)
            
            % Extract the whole duration of the resting event. Coerce the 
            % start index to values above 1 to preclude rounding to 0.
            start_ind = max(floor(Trial_EventTimes(ET)*Fs),1);
            
            % Convert the duration from seconds to samples.
            dur = round(Trial_Durations(ET)*Fs); 
            
            % Get ending index for data chunk. If event occurs at the end of
            % the trial, assume animal whisks as soon as the trial ends and
            % give a 200ms buffer.
            stop_ind = min(start_ind+dur,Expected_length-round(0.2*Fs));
            
            % Extract data from the trial and add to the cell array for the
            % current loaded file
            fname = IdentifyStructureSubfield(rmfield(Data,'Fs'),dataTypes{dt});
            if isempty(fname)
                stop = 1;
            end
            Trial_RestVals{ET} = Data.(fname).(dataTypes{dt})(:,start_ind:stop_ind);
            
        end
        
        % Add all periods of rest to a cell array for all files
        RestVals{f} = Trial_RestVals';
        
        % Transfer information about resting periods to the new structure
        EventTimes{f} = Trial_EventTimes';
        Durations{f} = Trial_Durations';
        PuffDistances{f} = Trial_PuffDistances';
        FileIDs{f} = repmat({FileID},1,length(Trial_EventTimes));
        FileDates{f} = repmat({FileDate},1,length(Trial_EventTimes));
        
    end
    
    % Combine the cells from separate files into a single cell array of all
    % resting periods
    RestData.(dataTypes{dt}).Data = [RestVals{:}]';
    RestData.(dataTypes{dt}).EventTime = cell2mat(EventTimes);
    RestData.(dataTypes{dt}).Duration = cell2mat(Durations);
    RestData.(dataTypes{dt}).PuffDistance = [PuffDistances{:}]';
    RestData.(dataTypes{dt}).FileID = [FileIDs{:}]';
    RestData.(dataTypes{dt}).FileDate = [FileDates{:}]';
    RestData.(dataTypes{dt}).Fs = Fs;
    if strcmp(dataTypes{dt},'Specgram')
        RestData.(dataTypes{dt}).freqs = Data.(fname).freqs;
    end
    RestData.(dataTypes{dt}).Normalized = 0;
    
    RestData.(dataTypes{dt}).MinutesSinceInfusion = NaN*ones(size(RestData.(dataTypes{dt}).EventTime));
    UniqueDays = unique(RestData.(dataTypes{dt}).FileDate);
    for UD = 1:length(UniqueDays)
        CurrDay = UniqueDays{UD};
        StrDay = ConvertDate(CurrDay);
        DayFilt = strcmp(RestData.(dataTypes{dt}).FileDate,CurrDay);
        RestData.(dataTypes{dt}).MinutesSinceInfusion(DayFilt) = ...
            ElapsedTime(RestData.(dataTypes{dt}).FileID(DayFilt),...
            RestData.(dataTypes{dt}).EventTime(DayFilt),...
            Thresholds.Infusion_Start_Times.(StrDay));
    end
    
end
