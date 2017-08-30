function [Neuro, NeurFs] = ProcessNeuro_atw0(RawData,NeurType,animal,hem,numday)
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Nov 2013
%   Version 1
%
%   SUMMARY: Processes neural data according to variable NeurType
%_______________________________________________________________
%   INPUTS:
%                           RawData - [struct] a structure containing
%                           neural data (RawData.neuro) and sampling
%                           frequency (RawData.an_fs)
%
%                           NeurType - [string] the type of processing
%                           desired. This variable references parameters
%                           which are stored in Global variables and
%                           ensures that all neural analysis is done in a
%                           conserved manner.
%                               Acceptable inputs are:
%                                   -HiGam
%_______________________________________________________________
%   OUTPUTS:
%                           Neuro - [array] the processed neural data
%
%                           NeurFs - [integer] the sampling frequency of
%                           the new neural signal
%_______________________________________________________________
%   REQUIRED SCRIPTS:
%       [res_data1, res_data2] = MatchDataLengths(data1,data2)
%_______________________________________________________________
%   CALLED BY:
%                           Process RawDataFile.m
%_______________________________________________________________
%   FUTURE VERSIONS:
%_______________________________________________________________
%   CHANGES FROM PREV VERS:
%_______________________________________________________________

ExpectedLength = RawData.TrialDur*RawData.an_fs;
Trimmed_Neuro = RawData.neuro(1:min(ExpectedLength,length(RawData.neuro)));
ThreshFile = dir('*Thresholds.mat');
if isempty(ThreshFile)
    display('ProcessNeuro.m: Warning! Thresholds file not found... creating blank file');
    Thresholds = [];
else
    display('ProcessNeuro.m: Thresholds file found...Loading');
    load(ThreshFile.name);
end

switch NeurType
    case 'MUpower'
        fpass = [300 3000];
    case 'Gam'
        fpass = [40 100];
    case 'HiGam'
        fpass = [60 100];
    case 'LoGam'
        fpass = [40 60];
    case 'Beta'
        fpass = [10 30];
    case 'SubAlpha'
        fpass = [0.1 8];
end


% CALCULATE NEURAL POWER
if ismember(NeurType,[{'MUpower'},{'Gam'},{'HiGam'},{'LoGam'},{'Beta'},{'SubAlpha'}])
    display(['ProcessNeuro.m: Processing ' NeurType '; fpass = ' num2str(fpass)])
    tic;
    Fs = RawData.an_fs;
    [z,p,k] = butter(4,fpass/(Fs/2));
    [sos,g] = zp2sos(z,p,k);
    filtNeuro = filtfilt(sos,g,Trimmed_Neuro-mean(Trimmed_Neuro));
    [z1,p1,k1] = butter(4,10/(Fs/2),'low');
    [sos1,g1] = zp2sos(z1,p1,k1);
    Long_Neuro = filtfilt(sos1,g1,filtNeuro.^2);
    Neuro = max(resample(Long_Neuro,RawData.dal_fr,RawData.an_fs),0);
    NeurFs = RawData.dal_fr;
    elapsed = toc;
    
% FILTER WIDEBAND LFP    
elseif strcmp(NeurType,'Wideband_LFP')
    display(['ProcessNeuro.m: Processing ' NeurType '; fpass = < 300 HZ'])
    tic;
    Fs = RawData.an_fs;
    [z,p,k] = butter(4,300/(Fs/2),'low');
    [sos,g] = zp2sos(z,p,k);
    filtNeuro = filtfilt(sos,g,Trimmed_Neuro-mean(Trimmed_Neuro));
    NeurFs = 1500; % 300 Hz * 5 will give adequate oversampling
    Neuro = resample(filtNeuro,NeurFs,RawData.an_fs);
    elapsed = toc;
     


% CALCULATE THE MULTI-UNIT SPIKE RATE
elseif nonzeros(ismember(NeurType,'SR'))
    display(['ProcessNeuro.m: Processing ' NeurType ])
    tic;
    [b,a] = butter(2,300/(RawData.an_fs/2),'high');
    MU_data = filtfilt(b,a,Trimmed_Neuro);
    
    % CHECK FOR VOLTAGE STANDARD DEVIATION: identify spikes as
    % voltage amplitudes above some user-defined number of standard
    % deviations
    strday = ConvertDate(num2str(numday));
    ok1 = Check4Threshold(['MUA_StDev_' strday],animal,hem);
    if ok1 == 0
        Thresholds.(['MUA_StDev_' strday]) = GetMUAStDev(numday);
        save([animal '_' hem '_Thresholds.mat'],'Thresholds');
    end
    
    % Check the animal's shared variables structure for a user defined
    % spike threshold. ie: the number of standard deviations above the mean
    % voltage that should be considered a spike.
    if isfield(Thresholds, ['Spikes_' strday])
        display(['ProcessNeuro.m: Threshold Spikes_' strday ' found.']) 
        thresh = Thresholds.(['Spikes_' strday]);
    else
        display(['ProcessNeuro.m: Threshold Spikes_' strday ' not found...creating threshold.']) 
        thresh = CreateSpikeThresh(RawData,Thresholds,strday);
        Thresholds.(['Spikes_' strday]) = thresh;
        save([animal '_' hem '_Thresholds.mat'],'Thresholds');
    end
    
    % Calculate the spike rate
    Neuro = MUspikeRATE_gauss(MU_data,thresh,RawData.an_fs,Thresholds.(['MUA_StDev_' strday]));
    NeurFs = RawData.an_fs;
    elapsed = toc;
end

display(['ProcessNeuro.m: Done. Time elapsed: ' num2str(elapsed) ' seconds.']);
display('____________________________________________________________')