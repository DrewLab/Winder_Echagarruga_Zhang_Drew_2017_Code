function [EventData] = ExtractEventTriggeredData_InfusionTrials(filenames,...
    EventData,dataTypes,epoch,InfusionTimes)
%   function [EventData] = ExtractEventTriggeredData(filenames,dataTypes,epoch)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Separates data corresponding to various behaviors into
%   structures
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   filenames - [cell array] names of files organized as rows.
%                   Files should already be processed using the script
%                   "ProcessRawDataFile.m" or "CalculatePredictedCBV.m".
%
%                   EventData - [struct] structure of existing data chunked
%                   around behavioral events
%
%                   dataTypes - [dataTypes] the measurements to be
%                   chunked into data epochs
%
%                   epoch - [struct] contains fields:
%                               duration - [double] the total time of the
%                               data epoch
%
%                               offset - [double] the amount of pre-event
%                               time (s) to include in the epoch.
%_______________________________________________________________
%   RETURN:                     
%                   EventData - [struct] data chunked around behavioral
%                   events
%_______________________________________________________________

% Control for dataTypes as string
if not(iscell(dataTypes))
    dataTypes = {dataTypes};
end

% Create temporary structure to work with before adding to the EventData
temp = struct();

for f = 1:length(filenames)   
    % Load ProcData File
    FileData = load(filenames{f});
    fname1 = fieldnames(FileData);
    StructData = FileData.(fname1{1});
    
    % Get the date and file ID to include in the EventData structure
    [~,~,FileDate,fileID] = GetFileInfo(filenames(f,:));
    
    % Get the types of behaviors present in the file (stim,whisk,rest)
    Beh_fields = fieldnames(StructData.Flags);
    
    for dT = 1:length(dataTypes)
        Data = [];
        dataType = dataTypes{dT};
        
        % Set the sampling frequency for the dataType
        Fs = StructData.Fs.(dataType);
        
        % Loop over the behaviors present in the file
        for Bf = 1:length(Beh_fields)
            
            %Preallocate space for unknown number of events using a
            %'temporary' structure of cells
            if not(isfield(temp,dataType))
                temp.(dataType) = [];
            end
            
            % Create behavioral subfields for the temp structure, if needed
            if not(isfield(temp.(dataType),Beh_fields{Bf}))
                Subfields = fieldnames(StructData.Flags.(Beh_fields{Bf}));
                blankcell = cell(1,size(filenames,1));
                StructVals = cell(size(Subfields));
                StructVals(:) = {blankcell};
                temp.(dataType).(Beh_fields{Bf}) = cell2struct(StructVals,...
                    Subfields, 1)';
                temp.(dataType).(Beh_fields{Bf}).FileID = blankcell;
                temp.(dataType).(Beh_fields{Bf}).FileDate = blankcell;
                temp.(dataType).(Beh_fields{Bf}).Data = blankcell;
            end
            
            % Assemble a structure to send to the sub-functions
            fname = IdentifyStructureSubfield(StructData,dataType);
            if isempty(fname)
                Data = StructData;
            else
                Data = StructData.(fname);
            end
            Data.Flags = StructData.Flags;
            Data.Info = StructData.Info;
            Data.Fs = StructData.Fs;
            
            % Extract the data from the epoch surrounding the event
            ScreenText = wraptext(['ExtractEventTriggeredData.m: Processing "'...
                Beh_fields{Bf} '" behavioral ' dataType ' data from file: '...
                fileID '...']);
            display(ScreenText)
            [chunk_data,EvFilter] = ExtractBehavioralData(Data,epoch,dataType,Beh_fields{Bf});
            
            % Add epoch details to temp struct
            [temp] = AddEpochInfo(Data,dataType,Beh_fields{Bf},temp,...
                fileID,FileDate,EvFilter,f,InfusionTimes);
            temp.(dataType).(Beh_fields{Bf}).Data{f} = chunk_data;
            
            % Add the sampling frequency, assume all Fs are the same for given
            % dataType
            temp.(dataType).(Beh_fields{Bf}).Fs = {Fs};
        end 
    end
end
% Convert the temporary stuct into a final structure
if isfield(StructData,'Neuro')
    if isfield(StructData.Neuro,'freqs')
        freqs = StructData.Neuro.freqs;
    else
        freqs = [];
    end
else
    freqs = [];
end
[EventData] = ProcessTempStruct(EventData,temp,epoch,freqs);


function [chunk_data,EvFilter] = ExtractBehavioralData(Data,epoch,dataType,Beh)
%   function [chunk_data,EvFilter] = ExtractBehavioralData(Data,epoch,dataType,Beh)
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

% Setup variables
EventTimes = Data.Flags.(Beh).EventTime;
TrialDur = Data.Info.TrialDur;
Fs = Data.Fs.(dataType);

% Get the content from Data.(dataType)
Data = getfield(Data, {}, dataType, {});

% Calculate start/stop times (seconds) for the events
all_epoch_starts = EventTimes - epoch.offset*ones(size(EventTimes));
all_epoch_ends = all_epoch_starts + epoch.duration*ones(size(EventTimes));

% Filter out events which are too close to the beginning or end of trials
start_filter = all_epoch_starts>0;
stop_filter = round(all_epoch_ends)<TrialDur; % Apply "round" to give an 
                                              % extra half second buffer 
                                              % and prevent indexing errors
EvFilter = logical(start_filter.*stop_filter);
ScreenText = wraptext(['ExtractEventTriggeredData>ExtractBehavioralData:' ...
    upper(Beh) ': Events at times: ' num2str(EventTimes(not(EvFilter))') ...
    ' seconds omitted. Too close to beginning/end of trial.']);
display(ScreenText);

% Convert the starts from seconds to samples, round down to the nearest
% sample, coerce the value above 1.
epoch_starts = max(floor(all_epoch_starts(EvFilter)*Fs),1);

% Calculate stops indices using the duration of the epoch, this avoids
% potential matrix dimension erros caused by differences in rounding when
% converting from seconds to samples.
sample_dur = round(epoch.duration*Fs);
epoch_stops = epoch_starts+sample_dur*ones(size(epoch_starts));

% Extract the chunk of data from the trial
chunk_data = zeros(sample_dur+1,length(epoch_starts),size(Data,1));
for st = 1:length(epoch_starts)
    chunk_inds = epoch_starts(st):epoch_stops(st);
    chunk_data(:,st,:) = Data(:,chunk_inds)';
end


function [temp] = AddEpochInfo(Data,dataType,Beh,temp,fileID,FileDate,EvFilter,f,InfusionTimes)
%   function [temp] = AddEpochInfo(Data,dataType,Beh,temp,fileID,FileDate,EvFilter,f)
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

% Get the field names for each behavior
fields = fieldnames(Data.Flags.(Beh));

% Filter out the events which are too close to the trial edge
for fd = 1:length(fields)
    field = fields{fd};
    temp.(dataType).(Beh).(field){f} = Data.Flags.(Beh).(field)(EvFilter,:)';
end

% Tag each event with the file ID, arrange cell array horizontally for
% later processing.
temp.(dataType).(Beh).FileID{f} = repmat({fileID},1,sum(EvFilter));
temp.(dataType).(Beh).FileDate{f} = repmat({FileDate},1,sum(EvFilter));

% Calculate the time since the infusion for each event
strdate = ConvertDate(fileID);
InfusionTime = InfusionTimes.(strdate);
temp.(dataType).(Beh).MinutesSinceInfusion{f} = ElapsedTime(temp.(dataType).(Beh).FileID{f},...
    temp.(dataType).(Beh).EventTime{f}',InfusionTime)';


function [EventData] = ProcessTempStruct(EventData,temp,epoch,freqs)
%   [EventData] = ProcessTempStruct(EventData,temp,epoch,freqs)
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
%______________________________________________________________

% Get the dataTypes from temp
dataTypes = fieldnames(temp);

for dT = 1:length(dataTypes)
    dataType = dataTypes{dT};
    
    % Get dataType names
    Beh_fields = fieldnames(temp.(dataType));
    
    % Intialize Behavior fields of the dataType sub-structure
    StructArray2 = cell(size(Beh_fields));
    EventData.(dataType) = cell2struct(StructArray2,Beh_fields, 1);
    
    for B = 1:length(Beh_fields)
        Beh = Beh_fields{B};
        
        % Get Behavior names
        Event_fields = fieldnames(temp.(dataType).(Beh));
        
        % Initialize Event fields for the Behavior sub-structure
        StructArray3 = cell(size(Event_fields));
        EventData.(dataType).(Beh) = cell2struct(StructArray3,Event_fields, 1);
        
        for Ef = 1:length(Event_fields);
            Evfield = Event_fields{Ef};
            TransferArray = [temp.(dataType).(Beh).(Evfield){:}];
            EventData.(dataType).(Beh).(Evfield) = permute(TransferArray,...
                unique([2,1,ndims(TransferArray)],'stable'));
        end
        EventData.(dataType).(Beh).epoch = epoch;
        if strcmp(dataType,'Specgram')
            EventData.Specgram.(Beh).freqs = freqs;
        end
        
        % Tag the data as not normalized
        EventData.(dataType).(Beh).Normalized = 0;
    end
end