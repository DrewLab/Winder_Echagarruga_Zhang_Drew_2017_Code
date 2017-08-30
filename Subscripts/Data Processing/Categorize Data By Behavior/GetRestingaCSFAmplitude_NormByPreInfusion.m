function [RestData] = GetRestingaCSFAmplitude_NormByPreInfusion(dataType,TimeThresh,Baselines,RestData)
%   [aCSFRMS] = GetRestingaCSFVAmplitude(animal,dataType)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Calculates the root mean squared (RMS) amplitude of a 
%   neural or CBV signal during periods of rest following aCSF infusion 
%   using the data in a Baselines structure.
%   
%_______________________________________________________________
%   PARAMETERS:             
%               dataType - [string] The data measurement of interest
%
%               TimeThresh - [double] number of minutes after infusion to
%               start calculating the amplitude.
%
%               Baselines - [structure] contains baseline data for various
%               data types.
%
%               RestData - [structure] repository for data and analysis of
%               resting data
%
%               Normalize - [logical] tells the program whether to
%               normalize the variance by the mean of the signal.
%_______________________________________________________________
%   RETURN:                     
%               RestData - [structure] repository for data and analysis of
%               resting data. Now with the resting amplitudes added.
%_______________________________________________________________


% Determine how many sessions occurred with resting data
UniqueRestingDates = unique(RestData.(dataType).FileDate);
BufferTime = 4; % Seconds
BufferSamples = BufferTime*RestData.(dataType).Fs;
Criteria.Fieldname = {'PuffDistance','Duration'};
Criteria.Comparison = {'gt','gt'};
Criteria.Value = {5,10+BufferTime};
FiltArray = FilterEvents(RestData.(dataType),Criteria);
TimeFilt = gt(RestData.(dataType).MinutesSinceInfusion,TimeThresh);
TrimFun = @(x,buffer)x(buffer:end);


% Loop over each session to get the periods of rest which occur after
% infusion
aCSFRestingAmps = [];
for URD = 1:length(UniqueRestingDates)
    Session = UniqueRestingDates{URD};
    SessionFilt = strcmp(RestData.(dataType).FileDate,Session);
    
    % Apply the filter to the Neuro data
    SessionData = RestData.(dataType).NormData(FiltArray&SessionFilt&TimeFilt);
    BufferCell = num2cell(BufferSamples*ones(size(SessionData)));
    TrimmedRest = cellfun(TrimFun,SessionData,BufferCell,'UniformOutput',0);
    sessionVars = cellfun(@var,TrimmedRest);
    RestingAmps = sqrt(sessionVars);
    aCSFRestingAmps = [aCSFRestingAmps; RestingAmps];
%     RestingData_IDs = [RestingData_IDs; Session_RestIDs(TimeFilt)];
%     TrialRestTimes = [TrialRestTimes; EventTimes(TimeFilt)];
%     AllMinutesSinceInfusion = [AllMinutesSinceInfusion; MinutesSinceInfusion(TimeFilt)];
end

%% Add the result to the Baselines.mat file

% Check for "dataType" in RestData
if not(isfield(RestData,dataType))
    RestData.(dataType) = [];
end

RestData.(dataType).TimeFilt = TimeFilt;
RestData.(dataType).SessionFilt = SessionFilt;
RestData.(dataType).CriteriaFilt = FiltArray;
RestData.(dataType).MeanRestingAmpFluctuation = aCSFRestingAmps;
Criteria.MinTimeSinceInfusion = TimeThresh;
RestData.(dataType).Criteria = Criteria;
