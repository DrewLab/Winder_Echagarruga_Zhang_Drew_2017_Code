function [] = TotalBehavioralTimes(animals,CBVType)
%   [] = TotalBehavioralTimes(animals,CBVType)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Calculates the amount of data which falls under various
%   behavioral categories
%   
%_______________________________________________________________
%   PARAMETERS:             
%               animals - [cell array] contains the animal IDs
%
%               CBVType - [string] fieldname which corresponds to the CBV
%               ROI.
%_______________________________________________________________
%   RETURN:                     
%                               
%_______________________________________________________________

% Preallocate
TrialLen = NaN*ones(1,length(animals));
epochLen = NaN*ones(1,length(animals));
RestLen = NaN*ones(1,length(animals));
RestNum = NaN*ones(1,length(animals));
WhiskLen = NaN*ones(1,length(animals));
WhiskNum = NaN*ones(1,length(animals));
WhiskChunkLen = NaN*ones(1,length(animals));
StrugLen = NaN*ones(1,length(animals));
StrugNum = NaN*ones(1,length(animals));
StrugChunkLen = NaN*ones(1,length(animals));
numContra = NaN*ones(1,length(animals));
ContraChunkLen = NaN*ones(1,length(animals));
numIpsi = NaN*ones(1,length(animals));
IpsiChunkLen = NaN*ones(1,length(animals));
numControl = NaN*ones(1,length(animals));
ControlChunkLen = NaN*ones(1,length(animals));

% Output display
clc
display ('Figure 1b: Calculating amount of behavioral data...')

for a = 1:length(animals)
    display([num2str(a) ' of ' num2str(length(animals)) ' animals done...']);
    prevdir = cd([animals{a} filesep]);
%     EventFile = ls(['*EVENTDATA_' CBVType '.mat']);
%     RestFile = ls(['*RESTDATA_' CBVType '.mat']);
%     ProcFiles = ls('*ProcData.mat');
%     load(strtrim(EventFile));
%     load(strtrim(RestFile));
%     load(strtrim(ProcFiles(1,:)))
    EventFile = dir(['*EVENTDATA_' CBVType '.mat']);
    RestFile = dir(['*RESTDATA_' CBVType '.mat']);
    ProcFiles = dir('*ProcData.mat');
    load(EventFile.name);
    load(RestFile.name);
    load(ProcFiles(1).name);
    TrialLen(a) = ProcData.TrialDur*size(ProcFiles,1);
    epochLen(a) = EventData.(CBVType).stim.epoch.duration - ...
        EventData.(CBVType).stim.epoch.offset;
    [RestStruct,RestFilt] = SelectBehavioralEvents(RestData.(CBVType),'Rest');
    RestLen(a) = sum(RestStruct.Duration(RestFilt));
    RestNum(a) = length(RestStruct.Duration(RestFilt));
    
    [DataStruct,VWFilt] = SelectBehavioralEvents(EventData.(CBVType),'VW');
    WhiskLen(a) = sum(DataStruct.Duration(VWFilt));
    WhiskNum(a) = length(DataStruct.Duration(VWFilt));
    WhiskChunkLen(a) = WhiskNum(a)*epochLen(a);
    
    [DataStruct,StrFilt] = SelectBehavioralEvents(EventData.(CBVType),'Str');
    StrugLen(a) = sum(DataStruct.Duration(StrFilt));
    StrugNum(a) = length(DataStruct.Duration(StrFilt));
    StrugChunkLen(a) = StrugNum(a)*epochLen(a);
    
    [DataStruct,ContraFilt] = SelectBehavioralEvents(EventData.(CBVType),'Contra');
    numContra(a) = size(DataStruct.Data(ContraFilt,:),1);
    ContraChunkLen(a) = numContra(a)*epochLen(a);
    
    [DataStruct,IpsiFilt] = SelectBehavioralEvents(EventData.(CBVType),'Ipsi');
    numIpsi(a) = size(DataStruct.Data(IpsiFilt,:),1);
    IpsiChunkLen(a) = numIpsi(a)*epochLen(a);
    
    [DataStruct,CtlFilt] = SelectBehavioralEvents(EventData.(CBVType),'Control');
    numControl(a) = size(DataStruct.Data(CtlFilt,:),1);
    ControlChunkLen(a) = numControl(a)*epochLen(a);
    cd(prevdir)
end

OtherChunkLen = TrialLen-(RestLen+WhiskChunkLen+StrugChunkLen+ContraChunkLen+...
    IpsiChunkLen+ControlChunkLen);

figure; plot([1,2,3,4,5,6,7], ...
    [OtherChunkLen', RestLen', WhiskChunkLen', ContraChunkLen', StrugChunkLen', ...
    IpsiChunkLen', ControlChunkLen', ]/60, 'Marker','o','LineStyle','none',...
    'MarkerEdgeColor','k','MarkerSize',8)
ylabel(sprintf('Total Time of\nBehavioral Events (min)'))
xlim([0.9 7.1]);
set(gca,'XTickLabel',{'Behavioral Transition','Rest','Volitional Whisking',...
    'Stim. - Contra','Extended Movement','Stim. - Ipsi.','Stim. - Auditory'})
set(gcf,'name','Figure 1b','numbertitle','off')

TRest = sum(RestLen);
TWhisk = sum(WhiskChunkLen);
TStrug = sum(StrugChunkLen);
TContra = sum(ContraChunkLen);
TIpsi = sum(IpsiChunkLen);
TControl = sum(ControlChunkLen);
TOther = sum(OtherChunkLen);
TTotal = sum(TrialLen);
labels = {'Beh. Transition','Rest','Stim. Contra','Vol. Whisk','Ext. Movement','Stim. Ipsi','Stim. Aud'};
axes('Position',[0.5, 0.5, 0.35, 0.35]); box on;
pie([TOther TRest TContra TWhisk TStrug TIpsi TControl]/TTotal,labels);