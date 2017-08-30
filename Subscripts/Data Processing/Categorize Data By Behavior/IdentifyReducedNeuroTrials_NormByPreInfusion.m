function [RestData] = IdentifyReducedNeuroTrials_NormByPreInfusion(NeuroThresh,aCSF_Amp,RestData,dataType)
%   function [] =  ()
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Identifies the trials where neural power (300-3000 Hz) is
%   less than that of aCSF controls. Then cycles through all trials
%   after the neural activity is reduced to allow the user to visually
%   inspect each one. The user can identify the end of the usable data. 
%   
%_______________________________________________________________
%   PARAMETERS:             
%               NeuroThresh - [double] the threshold of neural suppresion
%               as a ratio between 0-1
%
%               aCSF_Amp - [double] the resting aCSF amplitude
%
%               dataType - [string] the type of neural activity of interest
%_______________________________________________________________
%   RETURN:                     
%                               
%_______________________________________________________________

% Load and setup
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

% Filter the RestData
BufferTime = 4; % Seconds
BufferSamples = BufferTime*RestData.(dataType).Fs;
Criteria.Fieldname = {'PuffDistance','Duration'};
Criteria.Comparison = {'gt','gt'};
Criteria.Value = {5,10+BufferTime};
FiltArray = FilterEvents(RestData.(dataType),Criteria);


% Apply the filter array to get info for each resting event
RestingDates = RestData.(dataType).FileDate(FiltArray);
RestingIDs = RestData.(dataType).FileID(FiltArray);
AllRestingData = RestData.(dataType).NormData(FiltArray);
MinutesSinceInfusion = RestData.(dataType).MinutesSinceInfusion(FiltArray);
% RestTimes = RestData.(dataType).EventTime(FiltArray);

% Trim the beginnings from periods of rest
TrimFun = @(x,buffer)x(buffer:end);
BufferCell = num2cell(BufferSamples*ones(size(AllRestingData)));
TrimmedRest = cellfun(TrimFun,AllRestingData,BufferCell,'UniformOutput',0);

% Calculate the variance for each period of rest
RestingVars = cellfun(@var,TrimmedRest);

% Identify unique imaging sessions
sessions = unique(RestingDates);
for s = 1:length(sessions)
    % Find the resting events from session{s}
    SessionFilt = strcmp(sessions{s},RestingDates);
    strdate = ConvertDate(sessions{s});
    
    % Get the FileIDs for the current session
    sessionRestingIDs = RestingIDs(SessionFilt);
    
    % Get time since infusion for session
    SessionMinutesSinceInfusion = MinutesSinceInfusion(SessionFilt);
    
    % Calculate the standard deviation for periods of rest
    sessionVars = RestingVars(SessionFilt);
    RestingAmps = sqrt(sessionVars);
    
    % Calculate the ratio of resting amplitude deviations
    NeuroAmpRatio = RestingAmps/aCSF_Amp;
    
    % Plot and save the time since infusion against the amplitude ratios for
    % verification
    figure; scatter(SessionMinutesSinceInfusion,NeuroAmpRatio)
    ylabel('Resting Amp / aCSF Resting Amp')
    xlabel('Time since infusion (min)')
    title([strdate ': ' dataType])

    
    % Identify the first trial where the neural activity goes below
    % NeuroThresh
    AfterAve = zeros(1,(length(NeuroAmpRatio)-1));
    for NAR = 1:(length(NeuroAmpRatio)-1)
        AfterAve(NAR) = mean(NeuroAmpRatio(NAR:end));
    end
    hold on; plot(SessionMinutesSinceInfusion(1:end-1),AfterAve,'r');
    FileID_Ind = find(and(NeuroAmpRatio(1:end-1)'<NeuroThresh,AfterAve<=NeuroThresh),1,'first');
    Thresholds.(['Reduced' dataType '_Start']).(strdate) = ...
        sessionRestingIDs(FileID_Ind);
    
end
RestData.(dataType).RestingAmplitude = sqrt(RestingVars);
RestData.(dataType).RestingData_IDs = RestingIDs;
RestData.(dataType).TrialRestTimes = RestData.(dataType).EventTime(FiltArray);
RestData.(dataType).RestCriteria = Criteria;
RestData.(dataType).CriteriaFilt = FiltArray;
RestData.(dataType).SessionFilt = SessionFilt;
save(ThreshFile.name,'Thresholds')