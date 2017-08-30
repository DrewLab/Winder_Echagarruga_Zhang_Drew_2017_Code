function [Stats] = CompareBehavioralVariance(animals,Behaviors,CBVType)
%   function [Stats] = CompareBehavioralVariance(animals,Behaviors,CBVType)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Calculates the variance of the CBV data corresponding to
%   various behavioral stats. This code requires that the data have already
%   been categorized and separated according to behavior (see
%   SingleAnimalProcessingmaster.m)
%   
%_______________________________________________________________
%   PARAMETERS:             
%               animals - [cell array] containing the IDs for all animals 
%               be analyzed.
%
%               Behaviors - [cell array] containing designations for the
%               behavioral categories to be tested.
%
%               CBVType - [string] name of the CBV ROI to be used for the
%               analysis.
%_______________________________________________________________
%   RETURN:                     
%               Stats - [structure] contains the pvalues and statistics for
%               the tests performed here.
%_______________________________________________________________

Vars = [];

% Display
clc
display('Calculating CBV variance for each behavior...')
for a = 1:length(animals)
    display([num2str(a) ' of ' num2str(length(animals)) ' animals...'])
    animal = animals{a};
%     EventFile = ls([animal filesep '*EventData_' CBVType '.mat']);
%     RestFile = ls([animal filesep '*RestData_' CBVType '.mat']);
    EventFile = dir([animal filesep '*EVENTDATA_' CBVType '.mat']);
    RestFile = dir([animal filesep '*RESTDATA_' CBVType '.mat']);
    load([animal filesep EventFile.name]);
    load([animal filesep RestFile.name]);
    display('Getting rest data...')
    Criteria.Min_Duration = 8;
    Criteria.Min_PuffDist = 5;
    [DataStruct] = SelectRestingPeriods(RestData,CBVType,Criteria);
    Rest = DataStruct.NormData;
    StrtInd = 2*RestData.(CBVType).Fs+1; % Buffer for CBV to return to baseline
    StpInd = 1*RestData.(CBVType).Fs; % Buffer for movement preparation
    RestVars = NaN*ones(1,length(Rest));
    for c = 1:length(Rest)
        RestVars(c) = var(Rest{c}(StrtInd:(end-StpInd)));
    end
    Vars.Rest(a) = mean(RestVars);
    FileIDs.Rest = DataStruct.FileID;
    EventStart.Rest = round(DataStruct.EventTime*...
        DataStruct.Fs);
    EventStop.Rest = EventStart.Rest+...
        round(DataStruct.Duration*DataStruct.Fs);
    
    % Volitional Whisking
    display('Getting voluntary movement data...')
    [DataStruct,FiltArray,~] = SelectBehavioralEvents(EventData.(CBVType),'VW');
    VWData = DataStruct.NormData(FiltArray,:);
    strt = (DataStruct.epoch.offset-1)*EventData.(CBVType).stim.Fs;
    stp = strt + 6*EventData.(CBVType).stim.Fs;
    inds = strt:stp;
    Vars.VW(a) = mean(var(VWData(:,inds),[],2));
    FileIDs.VW = DataStruct.FileID(FiltArray);
    EventStart.VW = round(DataStruct.EventTime(FiltArray)*...
        DataStruct.Fs);
    EventStop.VW = EventStart.VW+...
        round((DataStruct.epoch.duration-DataStruct.epoch.offset)*DataStruct.Fs);

    % Extended Movement
    [DataStruct,FiltArray,~] = SelectBehavioralEvents(EventData.(CBVType),'Str');
    StrData = DataStruct.NormData(FiltArray,:);
    strt = (DataStruct.epoch.offset-1)*EventData.(CBVType).stim.Fs;
    stp = size(StrData,2);
    inds = strt:stp;
    Vars.Str(a) = mean(var(StrData(:,inds),[],2));
    FileIDs.Str = DataStruct.FileID(FiltArray);
    EventStart.Str = round(DataStruct.EventTime(FiltArray)*...
        DataStruct.Fs);
    EventStop.Str = EventStart.Str+...
        round((DataStruct.epoch.duration-DataStruct.epoch.offset)*DataStruct.Fs);

    % Contralateral Puff
    display('Getting stimulus evoked data...')
    [DataStruct,FiltArray,~] = SelectBehavioralEvents(EventData.(CBVType),'Contra');
    ContraData = DataStruct.NormData(FiltArray,:);
    strt = (DataStruct.epoch.offset)*EventData.(CBVType).stim.Fs;
    stp = size(ContraData,2);
    inds = strt:stp;
    Vars.Contra(a) = mean(var(ContraData(:,inds),[],2));
    FileIDs.Contra = DataStruct.FileID(FiltArray);
    EventStart.Contra = round(DataStruct.EventTime(FiltArray)*...
        DataStruct.Fs);
    EventStop.Contra = EventStart.Contra+...
        round((DataStruct.epoch.duration-DataStruct.epoch.offset)*DataStruct.Fs);
    
    % Ipsilateral Puff
    [DataStruct,FiltArray,~] = SelectBehavioralEvents(EventData.(CBVType),'Ipsi');
    IpsiData = DataStruct.NormData(FiltArray,:);
    strt = (DataStruct.epoch.offset)*EventData.(CBVType).stim.Fs;
    stp = size(IpsiData,2);
    inds = strt:stp;
    Vars.Ipsi(a) = mean(var(IpsiData(:,inds),[],2));
    FileIDs.Ipsi = DataStruct.FileID(FiltArray);
    EventStart.Ipsi = round(DataStruct.EventTime(FiltArray)*...
        DataStruct.Fs);
    EventStop.Ipsi = EventStart.Ipsi+...
        round((DataStruct.epoch.duration-DataStruct.epoch.offset)*DataStruct.Fs);
    
    % Auditory Puff
    [DataStruct,FiltArray,~] = SelectBehavioralEvents(EventData.(CBVType),'Control');
    CtlData = DataStruct.NormData(FiltArray,:);
    strt = (DataStruct.epoch.offset)*EventData.(CBVType).stim.Fs;
    stp = size(CtlData,2);
    inds = strt:stp;
    Vars.Control(a) = mean(var(CtlData(:,inds),[],2));
    FileIDs.Control = DataStruct.FileID(FiltArray);
    EventStart.Control = round(DataStruct.EventTime(FiltArray)*...
        DataStruct.Fs);
    EventStop.Control = EventStart.Control+...
        round((DataStruct.epoch.duration-DataStruct.epoch.offset)*DataStruct.Fs);
    
    
%     % Get the variance from the rest of the data
    display('Gathering all other data...')
%     Basefile = ls([animal filesep '*Baselines.mat']);
%     ProcFiles = ls([animal filesep '*ProcData.mat']);
%   load([animal filesep Basefile]);
    Basefile = dir([animal filesep '*Baselines.mat']);
    ProcFiles = dir([animal filesep '*ProcData.mat']);
    load([animal filesep Basefile.name]);
    for PF = 1:size(ProcFiles,1)
%         load([animal filesep ProcFiles(PF,:)])
        load([animal filesep ProcFiles(PF).name])
%         [~,~,FileDate,FileID] = GetFileInfo(ProcFiles(PF,:));
        [~,~,FileDate,FileID] = GetFileInfo(ProcFiles(PF).name);
        strdate = ConvertDate(FileDate);
        TrialData = ProcData.(CBVType)/mean(Baselines.(CBVType).(strdate).Means)-1;
        % Setup index for dataType
        inds = ones(size(ProcData.(CBVType)));
        % Behavior at the beginning/end of the trial cannot be definitively
        % classified, omit from analysis
        inds(1:round(DataStruct.epoch.offset*DataStruct.Fs)) = 0;
        inds(end-((DataStruct.epoch.duration-DataStruct.epoch.offset)*DataStruct.Fs):end) = 0;
        for B = 1:length(Behaviors)
            filefilt = strcmp(FileIDs.(Behaviors{B}), FileID);
            strts = EventStart.(Behaviors{B})(filefilt);
            stps = EventStop.(Behaviors{B})(filefilt);
            for st = 1:length(strts)
                inds(strts(st):stps(st)) = 0;
            end
        end
        ups = nonzeros(find(diff(inds)==1));
        downs = nonzeros(find(diff(inds)==-1));
        peakdurs = downs-ups;
        durationfilter = find(peakdurs < (6*DataStruct.Fs));
        for df = 1:length(durationfilter)
            inds(ups(durationfilter(df)):downs(durationfilter(df))) = 0;
        end
        if sum(inds) == 0
            continue
        else
            OtherStruct.Data{PF} = TrialData(logical(inds));
        end
    end
    
    OtherVars = NaN*ones(size(ProcFiles,1),1);
    loopind = 1;
    for c = 1:length(OtherStruct.Data)
        if not(isempty(OtherStruct.Data{c}))
            OtherVars(loopind) = var(OtherStruct.Data{c});
            loopind=loopind+1;
        end
    end
    NaNFilt = not(isnan(OtherVars));
    Vars.Other(a) = mean(OtherVars(NaNFilt));
end

%% Statistics
Bonferroni = 6; % Number of comparisons
nContra = Vars.Contra./Vars.Rest;
[Stats.Contra.pval,~,stat] = signrank(Vars.Contra,Vars.Rest,...
    'method','approximate');
Stats.Contra.Pcorrected = Stats.Contra.pval*Bonferroni;
Stats.Contra.Zval = stat.zval;

nIpsi = Vars.Ipsi./Vars.Rest;
[Stats.Ipsi.pval,~,stat] = signrank(Vars.Ipsi,Vars.Rest,...
    'method','approximate');
Stats.Ipsi.Pcorrected = Stats.Ipsi.pval*Bonferroni;
Stats.Ipsi.Zval = stat.zval;

nVW = Vars.VW./Vars.Rest;
[Stats.Whisk.pval,~,stat] = signrank(Vars.VW,Vars.Rest,...
    'method','approximate');
Stats.Whisk.Pcorrected = Stats.Whisk.pval*Bonferroni;
Stats.Whisk.Zval = stat.zval;

nStr = Vars.Str./Vars.Rest;
[Stats.ExtMov.pval,~,stat] = signrank(Vars.Str,Vars.Rest,...
    'method','approximate');
Stats.ExtMov.Pcorrected = Stats.ExtMov.pval*Bonferroni;
Stats.ExtMov.Zval = stat.zval;

nOther = Vars.Other./Vars.Rest;
[Stats.BehTrans.pval,~,stat] = signrank(Vars.Other,Vars.Rest,...
    'method','approximate');
Stats.BehTrans.Pcorrected = Stats.BehTrans.pval*Bonferroni;
Stats.BehTrans.Zval = stat.zval;

NaNFilt = not(isnan(Vars.Control));
nCtl = Vars.Control(NaNFilt)./Vars.Rest(NaNFilt);
[Stats.Aud.pval,~,stat] = signrank(Vars.Control(NaNFilt),Vars.Rest(NaNFilt),...
    'method','approximate');
Stats.Aud.Pcorrected = Stats.Aud.pval*Bonferroni;
Stats.Aud.Zval = stat.zval;

figure;
scatter([ones(size(nStr)), 2*ones(size(nContra)), 3*ones(size(nIpsi)), ...
    4*ones(size(nOther)), 5*ones(size(nCtl)), 6*ones(size(nVW))],...
    [nStr, nContra, nIpsi, nOther, nCtl, nVW],'o','linewidth',1.5,...
    'MarkerEdgeColor',[0.3 0.3 0.3]);

hold on; scatter([1,2,3,4,5,6],...
    [mean(nStr),mean(nContra),mean(nIpsi),mean(nOther),mean(nCtl),mean(nVW)],...
    'ko','linewidth',2,'MarkerFaceColor','k');
hold off;
ylabel('Variance/Resting Variance');
xlim([0.8 6.2]);
set(gca,'XTick',[1;2;3;4;5;6],'XTickLabel',...
    {'Ext. movement';'Stim-Contra';'Stim-Ipsi';'Beh.Transition';'Stim-Auditory';'Whisk'})
set(gcf,'name','Figure 1e','numbertitle','off')
ylim([0 5])

