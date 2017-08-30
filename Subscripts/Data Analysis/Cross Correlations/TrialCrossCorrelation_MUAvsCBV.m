function [CC] = TrialCrossCorrelation_MUAvsCBV(filenames, CBVType, CBVFilterParams)
%   function [CC] = TrialCrossCorrelation_MUAvsCBV(filenames, CBVType, CBVFilterParams)
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

CBV_cutoff = CBVFilterParams.cutoff;
CBV_filt_ord = CBVFilterParams.order;
maxlags = 5;

TimeBuffer = 5;
TimeThresh = 10;

BaseFile = dir('*Baselines.mat');
load(BaseFile.name)
i=1;
for f = 1:length(filenames)
    % Load ProcData and RawData files
    load(filenames{f})
    [~,~,FileDate,~] = GetFileInfo(filenames{f});
    strdate = ConvertDate(FileDate);
    
    % Process the Neural
    [z,p,k] = butter(4,1/(ProcData.Fs.MUpower_fs/2),'low');
    [sos,g] = zp2sos(z,p,k);
    MUpower = filtfilt(sos,g,ProcData.MUpower);
    nNeur = MUpower/mean(MUpower)-1;
    
    % Process the CBV
    Fs = ProcData.Fs.([CBVType '_fs']);
    nCBV = detrend(ProcData.(CBVType))/mean(Baselines.(CBVType).(strdate).Means);
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
            *Fs),length(nNeur));
    else
        Strtinds = 1;
        Stpinds = ProcData.TrialDur*Fs;
    end
    for si = 1:length(Strtinds)
        [CC.MUpower.vals(i,:),lags] = xcorr(FiltCBV(Strtinds(si):Stpinds(si)),...
            nNeur(Strtinds(si):Stpinds(si)),maxlags*Fs,'coef');
        i = i+1;
    end
end
CC.MUpower.Lags = lags/Fs;
CC.MUpower.vals = squeeze(CC.MUpower.vals);