function [CC] = TrialCrossCorrelation_LFPvsCBV(Filenames,CBVType,SpecParams,...
    CBVFilterParams)
%   function [CC] = TrialCrossCorrelation_LFPvsCBV(Filenames,CBVType,...
%                   SpecParams, CBVFilterParams)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Calculates the cross correlation for data occurring
%   between sensory stimuli.
%   
%_______________________________________________________________
%   PARAMETERS:             
%              Inputs:
%               filenames - [cell array]
%
%               CBVType - [string] name of field in filenames containing
%               CBV data
%
%               SpecParams - [Struct] parameters for calculating
%               spectrogram
%                       .fpass - [array] [lower, upper] frequency.
%                       .tapers - [array] number of tapers.
%                       .movingwin - [array] [window size, window step]
%                           See chronux toolbox doc for more
%
%               CBVFilterParams - [Struct] parameters for the low-pass  
%               butterworth filter applied to the CBV data
%                       .cutoff - [int] frequency cutoff
%                       .order - [int] order of the butterworth filter
%_______________________________________________________________
%   RETURN:                     
%               CC - [struct] Contains the results of the cross correlation
%                       .LFP.vals - [matrix] (frequency x lag x period)
%                       .LFP.Lags - [array] time vector of the lags, sec.
%                       .LFP.Freqs - [array] frequencies of the specgram
%_______________________________________________________________


CBV_cutoff = CBVFilterParams.cutoff;
CBV_filt_ord = CBVFilterParams.order;
maxlags = 5;

% Parameters for filtering non-sensory data
TimeBuffer = 5; % Buffer between start of data and puff
TimeThresh = 10; % Minimum duration of the data between puffs

% BaseFile = ls('*Baselines.mat');
% load(BaseFile)
BaseFile = dir('*Baselines.mat');
load(BaseFile.name)
OutputIndex=1;
for f = 1:length(Filenames)
    load(Filenames{f})
    Fs = ProcData.Fs.([CBVType '_fs']);
    [animal,hem,FileDate,FileID] = GetFileInfo(Filenames{f});
    strdate = ConvertDate(FileDate);
    display(['Calculating cross-corr for: ' animal ', ' FileID '. '...
        num2str(f) ' of ' num2str(length(Filenames))]);

    % Calculate and process the neural spectrogram
    params.fpass = SpecParams.fpass;
    params.Fs = ProcData.Fs.Wideband_LFP_fs;
    movingwin = SpecParams.movingwin;
    params.tapers = SpecParams.tapers;
    [S,t,fr]=mtspecgramc(detrend(ProcData.Wideband_LFP),movingwin ,params);
    zpad1 = zeros(floor(t(1)/movingwin(2)),size(S,2));
    nNeur = S./(ones(size(S,1),1)*mean(S))-1;
    FiltNeur = [zpad1; nNeur; zpad1]';
    
    % Process the CBV
    nCBV = ProcData.(CBVType)/mean(Baselines.(CBVType).(strdate).Means);
    [z,p,k] = butter(CBV_filt_ord,CBV_cutoff/(Fs/2),'low');
    [sos,g] = zp2sos(z,p,k);
    FiltCBV = filtfilt(sos,g,nCBV-mean(nCBV));
    
    % Get all Puff times
    AllPuffs = [ProcData.Sol.Contra, ProcData.Sol.Ipsi, ProcData.Sol.Tail,...
        ProcData.Sol.Control];
    if not(isempty(AllPuffs))
        SortedPuffs = [0 sort(AllPuffs)];
        InterPuffTime = diff([SortedPuffs ProcData.TrialDur]);
        Spontinds = InterPuffTime>TimeThresh;
        Strtinds = max(round((SortedPuffs(Spontinds)+TimeBuffer*...
            ones(1,sum(Spontinds)))*Fs),1);
        Stpinds = min(Strtinds + round((InterPuffTime(Spontinds)-TimeBuffer)...
            *Fs),length(FiltNeur));
    else
        Strtinds = 1;
        Stpinds = ProcData.TrialDur*Fs;
    end
    for si = 1:length(Strtinds) % Looping over each inter-puff period
        for r = 1:size(FiltNeur,1) % Looping over each frequency of specgram
            [CC.LFP.vals(r,:,OutputIndex),lags] = xcorr(FiltCBV(Strtinds(si):Stpinds(si)),...
                FiltNeur(r,Strtinds(si):Stpinds(si)),maxlags*Fs,'coef');
        end
        OutputIndex=OutputIndex+1;
    end
end
CC.LFP.Lags = lags/Fs;
CC.LFP.Freqs = fr;
save([animal '_' hem '_CBVXCorr_TrialData.mat'],'CC');