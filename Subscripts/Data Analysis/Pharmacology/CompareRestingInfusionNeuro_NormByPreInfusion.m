function [RestData] = CompareRestingInfusionNeuro_NormByPreInfusion(aCSFNeuro,RestData,NeuroType)
%   function [] =  ()
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

% Load and setup
BufferTime = 0; % Seconds
BufferSamples = BufferTime*RestData.(NeuroType).Fs;

ThreshFile = dir('*Thresholds.mat');
if isempty(ThreshFile)
    Thresholds = [];
    Thresholds.Infusion_Start_Times = [];
else
    load(ThreshFile.name);
end

% Filter the RestData
Criteria.Fieldname = {'PuffDistance','Duration'};
Criteria.Comparison = {'gt','gt'};
Criteria.Value = {5,10+BufferTime};
FiltArray = FilterEvents(RestData.(NeuroType),Criteria);

% Apply the filter array to get info for each resting event
RestingDates = RestData.(NeuroType).FileDate(FiltArray);
RestingIDs = RestData.(NeuroType).FileID(FiltArray);
AllRestingData = RestData.(NeuroType).NormData(FiltArray);

% Trim the beginnings from periods of rest
TrimFun = @(x,buffer)x(buffer+1:end);
BufferCell = num2cell(BufferSamples*ones(size(AllRestingData)));
TrimmedRest = cellfun(TrimFun,AllRestingData,BufferCell,'UniformOutput',0);

% Calculate the variance for each period of rest
RestingVars = cellfun(@var,TrimmedRest);

% Identify unique imaging sessions
sessions = unique(RestingDates);
SessionVars = cell(1,length(sessions));
for s = 1:length(sessions)
    % Find the resting events from session{s}
    SessionFilt = strcmp(sessions{s},RestingDates);
    strdate = ConvertDate(sessions{s});
    
    % Get the FileIDs for the current session
    sessionRestingIDs = RestingIDs(SessionFilt);
    
    if isempty(Thresholds.(['ReducedMUpower_Start']).(strdate))
        SessionVars{s} = [];
        continue;
    end
    
    % Determine the time since infusion for verification
    MinutesSinceAttenuation = ElapsedTime(sessionRestingIDs,...
    zeros(size(sessionRestingIDs)),...
    Thresholds.(['ReducedMUpower_Start']).(strdate));
    
    % Filter out the resting events that occur before neural activity is
    % reduced
    TimeFilt = MinutesSinceAttenuation>=0;
    
    % Calculate the standard deviation for periods of rest
    AllSessionInfusionVars = RestingVars(SessionFilt);
    SessionVars{s} = AllSessionInfusionVars(TimeFilt)';
end    
% Calculate the confidence interval on the ratio of resting amplitude
% deviations
InfusionNeuroAmplitudeFluctuations = sqrt([SessionVars{:}]);
if length(InfusionNeuroAmplitudeFluctuations)<4
    RestData.(NeuroType).MeanRestingAmpFluctuations = InfusionNeuroAmplitudeFluctuations;
    RestData.(NeuroType).AmpRatio_SD = [];
    RestData.(NeuroType).AmpRatio_boot = [];
    RestData.(NeuroType).AmpRatio = [];
else
    [CBVAmpRatio_SD,CBVAmpRatio_boot] = Bootstrap_RatioAmplitudes(10000,...
        aCSFNeuro,InfusionNeuroAmplitudeFluctuations);
    
    % Add the confidence interval and bootstrap samples to RestData
    RestData.(NeuroType).MeanRestingAmpFluctuations = InfusionNeuroAmplitudeFluctuations;
    RestData.(NeuroType).AmpRatio_SD = CBVAmpRatio_SD;
    RestData.(NeuroType).AmpRatio_boot = CBVAmpRatio_boot;
    RestData.(NeuroType).AmpRatio = mean(CBVAmpRatio_SD);
end
end

function [AmpSD,AmpRatios] = Bootstrap_RatioAmplitudes(numreps,aCSFNeuro,InfusionNeuro)
%   function [] =  ()
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Generates a bootstrap sample of the CBV ratios following
%   aCSF and another infusion. Calculates the 95% confidence interval on
%   the ratio from the bootstrap sample.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   numreps - [integer] the number of bootstrap resamples
%
%                   aCSFCBV - [array] an array of the resting CBV
%                   amplitudes after aCSF infusion
%
%                   InfusionCBV - [array] an array of the resting CBV
%                   amplitudes after a pharmacological infusion
%_______________________________________________________________
%   RETURN:                     
%                   Ampci - [array] the [lower,upper] values of the 95% 
%                   confidence interval.
%
%                   AmpRatios - [array] the values of the bootstrap
%                   resampling.
%_______________________________________________________________

% Create a bootstrap sample of the CBV ratios
AmpRatios = NaN*ones(1,numreps);
for n = 1:numreps
    RandInds_aCSF = randsample(length(aCSFNeuro),length(aCSFNeuro),1);
    RandInds_Inf = randsample(length(InfusionNeuro),length(InfusionNeuro),1);
    aCSF_Mean(n) = mean(aCSFNeuro(RandInds_aCSF));
    Inf_Mean(n) = mean(InfusionNeuro(RandInds_Inf));
    AmpRatios(n) = Inf_Mean(n)/aCSF_Mean(n);
end

Mu = mean(AmpRatios);
Sigma = std(AmpRatios);
AmpSD(1) = Mu-Sigma;
AmpSD(2) = Mu+Sigma;
end