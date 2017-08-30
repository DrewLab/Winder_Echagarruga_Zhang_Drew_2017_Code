function [] = BehavioralCrossCorrelation_LFPvsCBV()
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

% Set Variables
Beh = 'Rest';
CBV_cutoff = 1; % Low Pass Filter Threshold
CBV_filt_ord = 4; % Low Pass Filter Order
maxlags = 5; % Maximum Lags for cross correlation
CBVType = 'CrossCorrROI';
BufferTime = 4;

% Load any existing cross correlation file
CCfile = ls('*CBVXCorr.mat');
if isempty(CCfile)==0
    load(CCfile)
    CC.LFP.(Beh).vals = [];
end

% Load categorized behavioral data structure
if strcmp(Beh, 'Rest')
    RestFile = ls(['*RESTDATA_' CBVType '.mat']);
    load(RestFile)
    AllBehData = RestData.(CBVType);
    BehCriteria.Fieldname = {'Duration','PuffDistance'};
    BehCriteria.Comparison = {'gt','gt'};
    BehCriteria.Value = {14,5};
    [FiltArray] = FilterEvents(AllBehData,BehCriteria);
    DataStruct.FileID = AllBehData.FileID(FiltArray);
    DataStruct.EventTime = AllBehData.EventTime(FiltArray);
    DataStruct.Duration = AllBehData.Duration(FiltArray);
%     [RestStruct,FiltArray] = SelectBehavioralEvents(RestData.(CBVType),Beh);
%     DataStruct.FileID = RestStruct.FileID(FiltArray);
%     DataStruct.EventTime = RestStruct.EventTime(FiltArray);
%     DataStruct.Duration = RestStruct.Duration(FiltArray);
    clear RestData
end

filenames = unique(DataStruct.FileID);

EpochInd = 1;
for f = 1:size(filenames,1)
    ProcFile = ls(['*_' filenames{f} '_ProcData.mat']);
    load(ProcFile)
    
    Fs = ProcData.Fs.([CBVType '_fs']);
    [animal,hem,~,FileID] = GetFileInfo(ProcFile);
    FileMatch = strcmp(DataStruct.FileID,FileID);

    % Calculate the neural spectrogram
    params.fpass = [0.1 150];
    params.Fs = ProcData.Fs.Wideband_LFP_fs;
    movingwin = [1,1/30];
    params.tapers = [5 9];
    [S,t,fr]=mtspecgramc(detrend(ProcData.Wideband_LFP),movingwin ,params);
    zpad1 = zeros(floor(t(1)/movingwin(2)),size(S,2));
    nNeur = S./(ones(size(S,1),1)*mean(S))-1;
    FiltNeur = [zpad1; nNeur; zpad1]';
    
    % Process the CBV
    nCBV = ProcData.(CBVType)/mean(ProcData.(CBVType))-1;
    [z,p,k] = butter(CBV_filt_ord,CBV_cutoff/(Fs/2),'low');
    [sos,g] = zp2sos(z,p,k);
    FiltCBV = filtfilt(sos,g,nCBV-mean(nCBV));
    
    % Calculate cross correlation for each event
    EventTimes = DataStruct.EventTime(FileMatch);
    Durations = DataStruct.Duration(FileMatch);
    for ET = 1:length(EventTimes)
        Strtind = round(EventTimes(ET)+BufferTime*Fs);
        Stpind = min(Strtind+round(Durations(ET)*Fs),length(FiltCBV));
        for r = 1:size(FiltNeur,1)
            [CC.LFP.(Beh).vals(r,:,EpochInd),lags] = xcorr(FiltCBV(Strtind:Stpind),...
                FiltNeur(r,Strtind:Stpind),maxlags*Fs,'coef');
        end
        EpochInd=EpochInd+1;
    end
end
CC.LFP.(Beh).Lags = lags/Fs;
CC.LFP.(Beh).Freqs = fr;
save([animal '_' hem '_CBVXCorr.mat'],'CC')