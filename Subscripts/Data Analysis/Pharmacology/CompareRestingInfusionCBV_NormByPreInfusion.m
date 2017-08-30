function [RestData] = CompareRestingInfusionCBV_NormByPreInfusion(aCSFCBV,RestData,CBVType)
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
BufferTime = 4; % Seconds
BufferSamples = BufferTime*RestData.(CBVType).Fs;
NeurType = 'MUpower'; % Sets MUA as the indicator of when infusion starts working

ThreshFile = dir('*Thresholds.mat');
if isempty(ThreshFile)
    Thresholds = [];
    Thresholds.Infusion_Start_Times = [];
else
    load(ThreshFile.name);
end

% if strcmp(InfusionType,'Muscimol')
%     NeurType = 'MUpower';
% elseif strcmp(InfusionType, 'Muscimol + CNQX + AP5')
%     NeurType = 'Gam';
% elseif strcmp(InfusionType, 'Muscimol + Adrenergic Blockers')
%     NeurType = 'Gam';
% else
%     error('Type of infusion not recognized, cannot identify the correct type of neural activity for analysis')
% end

% Filter the RestData
Criteria.Fieldname = {'PuffDistance','Duration'};
Criteria.Comparison = {'gt','gt'};
Criteria.Value = {5,10+BufferTime};
FiltArray = FilterEvents(RestData.(CBVType),Criteria);

% Apply the filter array to get info for each resting event
RestingDates = RestData.(CBVType).FileDate(FiltArray);
RestingIDs = RestData.(CBVType).FileID(FiltArray);
AllRestingData = RestData.(CBVType).NormData(FiltArray);

% Trim the beginnings from periods of rest
TrimFun = @(x,buffer)x(buffer:end);
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
    
    if isempty(Thresholds.(['Reduced' NeurType '_Start']).(strdate))
        SessionVars{s} = [];
        continue;
    end
    
    % Get the FileIDs for the current session
    sessionRestingIDs = RestingIDs(SessionFilt);
    
    % Determine the time since infusion for verification
    MinutesSinceAttenuation = ElapsedTime(sessionRestingIDs,...
    zeros(size(sessionRestingIDs)),...
    Thresholds.(['Reduced' NeurType '_Start']).(strdate));
    
    % Filter out the resting events that occur before neural activity is
    % reduced
    TimeFilt = MinutesSinceAttenuation>=0;
    
    % Calculate the standard deviation for periods of rest
    AllSessionInfusionVars = RestingVars(SessionFilt);
    SessionVars{s} = AllSessionInfusionVars(TimeFilt)';
end
% Calculate the confidence interval on the ratio of resting amplitude 
% deviations
InfusionCBVAmplitudeFluctuations = sqrt([SessionVars{:}]);
if length(InfusionCBVAmplitudeFluctuations)<4
    RestData.(CBVType).MeanRestingAmpFluctuation = InfusionCBVAmplitudeFluctuations;
    RestData.(CBVType).AmpRatio_SD = [];
    RestData.(CBVType).AmpRatio_boot = [];
    RestData.(CBVType).AmpRatio = [];
else
    [CBVAmpRatio_SD,CBVAmpRatio_boot] = Bootstrap_RatioAmplitudes(10000,...
        aCSFCBV,InfusionCBVAmplitudeFluctuations);
    
    % Add the confidence interval and bootstrap samples to RestData
    RestData.(CBVType).MeanRestingAmpFluctuation = InfusionCBVAmplitudeFluctuations;
    RestData.(CBVType).AmpRatio_SD = CBVAmpRatio_SD;
    RestData.(CBVType).AmpRatio_boot = CBVAmpRatio_boot;
    RestData.(CBVType).AmpRatio = mean(CBVAmpRatio_SD);
end

end

function [AmpSD,AmpRatios] = Bootstrap_RatioAmplitudes(numreps,aCSFCBV,InfusionCBV)
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
    RandInds_aCSF = randsample(length(aCSFCBV),length(aCSFCBV),1);
    RandInds_Inf = randsample(length(InfusionCBV),length(InfusionCBV),1);
    Inf_Mean(n) = mean(InfusionCBV(RandInds_Inf));
    aCSF_Mean(n) = mean(aCSFCBV(RandInds_aCSF));
    AmpRatios(n) = Inf_Mean(n)/aCSF_Mean(n);
end
Mu = mean(AmpRatios);
Sigma = std(AmpRatios);
AmpSD(1) = Mu-Sigma;
AmpSD(2) = Mu+Sigma;

end