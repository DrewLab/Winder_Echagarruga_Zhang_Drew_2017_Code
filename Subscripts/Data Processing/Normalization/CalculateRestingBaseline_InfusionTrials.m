function [] = CalculateRestingBaseline_InfusionTrials(RestData,dataType,animal,hem)
%   [] = CalculateRestingBaseline_InfusionTrials(RestData,dataType,animal,hem)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Uses periods when animal is not being stimulated or moving
%   to establish a baseline for a given session of imaging.
%   
%_______________________________________________________________
%   PARAMETERS:    
%                   RestData - [Struct] contains all resting data to be
%                   used for the baseline.
%
%                   dataType = [string] data measurement desired for
%                   baseline
%
%                   animal = [string] animal ID
%
%                   hem = [string] hemisphere of animal recorded
%                           ('LH,'RH')
%                               
%_______________________________________________________________
%   RETURN:                     
%                   none, result of script is a saved file                   
%_______________________________________________________________

%% Load and setup
% Set buffer size after beginning of rest
BufferTime = 4; % seconds
BufferSamples = BufferTime*RestData.(dataType).Fs;

% Load previous Baselines
BaseFile = dir('*Baselines.mat');
if not(isempty(BaseFile))
    wraptext('CalculateRestingBaselines_InfusionTrials.m: Baselines file found...loading...');
    load(BaseFile.name);
else
    wraptext('CalculateRestingBaselines_InfusionTrials.m: No baselines file found...creating...');
    Baselines = struct();
end

% Load the threshold file
ThreshFile = dir('*Thresholds.mat');
if isempty(ThreshFile)
    Thresholds = [];
    Thresholds.Infusion_Start_Times = [];
else
    load(ThreshFile.name);
    if not(isfield(Thresholds,'Infusion_Start_Times'))
        Thresholds.Infusion_Start_Times = [];
    end
end

% Check for "dataType" in RestData
if not(isfield(RestData,dataType))
    ErrorText = ['No Resting Data for ' dataType ...
        ' found in RestData.mat. Run GetAllRestingData.m for dataType'];
    error(ErrorText) 
end
% Determine if a baseline will be overwritten
if isfield(Baselines,dataType)
    wraptext(['Baselines for ' dataType ' already exist...overwriting'])
end

% Initialize dataType sub-structure
Baselines.(dataType) = [];

                % % Define Variables
                % target_duration = 14;
                % Min_PuffDistance = 5; % seconds between rest event and puff

% Get filter array for the resting aCSF events
Criteria.Fieldname = {'PuffDistance','Duration'};
Criteria.Comparison = {'gt','gt'};
Criteria.Value = {5,10+BufferTime};
FiltArray = FilterEvents(RestData.(dataType),Criteria);

% Apply the filter array to get info for each resting event
Fs = RestData.(dataType).Fs;
RestingDates = RestData.(dataType).FileDate(FiltArray);
RestingIDs = RestData.(dataType).FileID(FiltArray);
AllRestingData = RestData.(dataType).Data(FiltArray);
RestingTimes = RestData.(dataType).EventTime(FiltArray);
PuffDistances = RestData.(dataType).PuffDistance(FiltArray);
RestingDurations = RestData.(dataType).Duration(FiltArray);

% Calculate the variance and mean of each resting event
TrimFun = @(x,buffer)x(buffer:end);
BufferCell = num2cell(BufferSamples*ones(size(AllRestingData)));
TrimmedRest = cellfun(TrimFun,AllRestingData,BufferCell,'UniformOutput',0);
RestingMeans = cellfun(@mean,TrimmedRest);
RestingVars = cellfun(@var,TrimmedRest);

sessions = unique(RestingDates);

                % [sessions,s_id] = GetUniqueDays(RestData.(dataType).FileID);
                % session_inds = [s_id' length(RestData.(dataType).FileID)];

                % % Get RestData Values for dataType
                % EventTimes = RestData.(dataType).EventTime;
                % PuffDistances = RestData.(dataType).PuffDistance;
                % Durations = RestData.(dataType).Duration;
    
for s = 1:length(sessions)
    
    % Find the resting events from session{s}
    SessionFilt = strcmp(sessions{s},RestingDates);
    strdate = ConvertDate(sessions{s});
    
                %     % Convert session date to a string for use as a fieldname of structure
                %     session = sessions{s};
                %     strdate = ConvertDate(session);
    
    % Initialize fdname sub-structure
    Baselines.(dataType).(strdate) = [];
    
                %     % Identify periods corresponding to the session
                %     Session_filter = zeros(size(EventTimes));
                %     Session_filter(session_inds(s):session_inds(s+1)) = 1;
                %     
                %     % Identify periods isolated from puffs
                %     [PD_filter] = GetPuffDistanceFilter(PuffDistances,Min_PuffDistance);
                %     
                %     % Identify periods of sufficient length
                %     Len_filter = Durations > target_duration;
                %
                %     % Combine filters to create final filter
                %     Baseline_filter = logical(Session_filter.*PD_filter.*Len_filter);
    
                %     % Keep values used to calculate baseline with the Baseline structure.
                %     Baselines.(dataType).(strdate).Vals = ...
                %         RestData.(dataType).Data(Baseline_filter);
                %     % Calculate the mean from each rest period, convert from cell array to
                %     % scalar array.
                %     % dimarray -> dimension of mean and variance functions for cellfun
                %     % syntax
                %     dimarray = mat2cell(2*ones(sum(Baseline_filter),1),...
                %         ones(sum(Baseline_filter),1));  
                %     BaseMeans = cellfun(@mean,RestData.(dataType).Data(Baseline_filter),...
                %         dimarray,'UniformOutput',0);
                %     Baselines.(dataType).(strdate).Means = [BaseMeans{:}]';
                %     
                %     % Calculate variance from each rest period, convert from cell array to
                %     % scalar array.
                %     BaseData = RestData.(dataType).Data(Baseline_filter);
                % %     Basemean = mean(Baselines.(dataType).(fdname).Means);
                %     for c = 1:length(BaseData)
                % %         NBase = BaseData{c}./(Basemean'*ones(1,size(BaseData{c},2)));
                % %         Baselines.(dataType).(fdname).StDev(c,:,:) = std(NBase,[],2);
                %         Baselines.(dataType).(strdate).StDev(c,:,:) = std(detrend(BaseData{c}),[],2);
                %     end
    % Calculate the time since infusion for each event
    MinutesSinceInfusion = ElapsedTime(RestingIDs,RestingTimes,...
    Thresholds.Infusion_Start_Times.(strdate));
    PreInfusionFilt = MinutesSinceInfusion<=0;
    PreInfusionLength = round(RestingDurations(SessionFilt&PreInfusionFilt)*Fs);
    
    % Save duration of rest period and distance of period from puff
    Baselines.(dataType).(strdate).RestingDuration = RestingDurations(SessionFilt);
    Baselines.(dataType).(strdate).PuffDistance = PuffDistances(SessionFilt);
    Baselines.(dataType).(strdate).FileID = RestingIDs(SessionFilt);
    Baselines.(dataType).(strdate).RestingTime = RestingTimes(SessionFilt);
    Baselines.(dataType).(strdate).Means = RestingMeans(SessionFilt);
    Baselines.(dataType).(strdate).Vars = RestingVars(SessionFilt);
    Baselines.(dataType).(strdate).MinutesSinceInfusion = MinutesSinceInfusion(SessionFilt);
    Baselines.(dataType).(strdate).PreInfusionMeans = ...
        sum(RestingMeans(SessionFilt&PreInfusionFilt)...
            .*PreInfusionLength/sum(PreInfusionLength));
    Baselines.(dataType).(strdate).PreInfusionVars = ...
        sum(RestingVars(SessionFilt&PreInfusionFilt)...
            .*PreInfusionLength/sum(PreInfusionLength));
    Baselines.(dataType).(strdate).Number_PreInfusion = sum(SessionFilt&PreInfusionFilt);
    Baselines.Criteria = Criteria;
end

save([animal '_' hem '_Baselines.mat'],'Baselines')
    
end

% function [PD_Filter] = GetPuffDistanceFilter(PuffDistance, distance_thresh)
% %   function [PD_Filter] = GetPuffDistanceFilter(PuffDistance, distance_thresh)
% %
% %   Author: Aaron Winder
% %   Affiliation: Engineering Science and Mechanics, Penn State University
% %   https://github.com/awinde
% %
% %   DESCRIPTION: Creates a logical array which designates events at
% %   sufficient distance from sensory stimulation.
% %   
% %_______________________________________________________________
% %   PARAMETERS:             
% %                   PuffDistance - [nx1 Array of Cells] distances from
% %                   surrounding puffs.  
% %
% %                   distance_thresh - [double] minimum unacceptable
% %                   distance from surrounding puffs. Units should be the
% %                   same as that of PuffDistance.
% %_______________________________________________________________
% %   RETURN:                     
% %                 	PD_Filter - [Logical Array] designates events of
% %                   sufficient distance from surrounding puffs.
% %_______________________________________________________________
% 
% % Change to column vector, if needed
% if size(PuffDistance,2)>size(PuffDistance,1)
%     temp = PuffDistance';
%     PuffDistance = temp;
%     clear temp;
% end
% 
% % Identify cells containing NaN (signifies no puffs in trial)
% NaNTags = cellfun(@isnan,PuffDistance,'UniformOutput',0);
% NaNFilter = cellfun(@sum,NaNTags);
% 
% % Identify cells with sufficient distance from puffs in trial
% AbsDistance = cellfun(@abs,PuffDistance,'UniformOutput',0);
% temparray = mat2cell(distance_thresh*ones(size(PuffDistance)),...
%     ones(size(PuffDistance)));
% DistanceTags = cellfun(@gt,AbsDistance,temparray,'UniformOutput',0);
% % Use prod as logical "and" for all distances in a cell
% DistanceFilter = cellfun(@prod,DistanceTags);
% 
% % Combine the NaN Filter and Distance Filter, coerce to 0 or 1, convert to
% % a logical array for later indexing.
% PD_Filter = logical(min(NaNFilter+DistanceFilter,1));
% end
