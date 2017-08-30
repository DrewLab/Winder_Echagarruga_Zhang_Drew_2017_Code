function [AveR2,Med_IndR2,CBVPred] = EvaluateCBVPredictionAccuracy(dataType1,dataType2,Beh,HRFs,CBVPred)
%   function [AveR2,Med_IndR2,CBVPred] = EvaluateCBVPredictionAccuracy(dataType1,dataType2,Beh,HRFs,CBVPred)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Evaluates the goodness of fit of the HRF predictions of
%   averaged and individual CBV.
%   
%_______________________________________________________________
%   PARAMETERS:             
%               dataType1 - [string] fieldname of neural measure
%
%               dataType2 - [string] fieldname of the CBV ROI
%
%               Beh - [string] behavioral category to be used
%
%               HRFs - [struct] contains the HRF array as a field
%
%               CBVPred - [struct] contains the CBV predictions
%_______________________________________________________________
%   RETURN:                     
%               AveR2 - [struct] contains the coefficient of determination
%               for the averaged data
%
%               IndR2 - [struct] contains the coefficient of determination
%               for the individual data
%
%               CBVPred - [struct] contains the CBV predictions with the
%               current predictions appended
%_______________________________________________________________

%% Setup
Event_Inds.CalcStart = 1;
Event_Inds.TestStart = 2;
if strcmp(Beh,'Str')
    Event_Inds.Increment = 1;
else
    Event_Inds.Increment = 2;
end

%% Load in the data structure
if strcmp(Beh,'Rest')
%     RestDataFileName1 = ls(['*RESTDATA_' dataType1 '.mat']);
%     RestDataFileName2 = ls(['*RESTDATA_' dataType2 '.mat']);
    RestDataFileName1 = dir(['*RESTDATA_' dataType1 '.mat']);
    RestDataFileName2 = dir(['*RESTDATA_' dataType2 '.mat']);
    if or(isempty(RestDataFileName1),isempty(RestDataFileName2))
        error(wraptext(['Missing RestData structure in current directory. '...
            'Change directory to location of file, move file to current '...
            'directory, or run "ExtractEventTriggeredData.m" to create RestData.']));
    else
%         LoadData = load(RestDataFileName1);
%         RestData.(dataType1) = LoadData.RestData.(dataType1);
%         LoadData = load(RestDataFileName2);
%         RestData.(dataType2) = LoadData.RestData.(dataType2);
%         BehData = RestData;
        LoadData = load(RestDataFileName1.name);
        RestData.(dataType1) = LoadData.RestData.(dataType1);
        LoadData = load(RestDataFileName2.name);
        RestData.(dataType2) = LoadData.RestData.(dataType2);
        BehData = RestData;
        clear RestData;
    end
else
%     EventDataFileName1 = ls(['*EVENTDATA_' dataType1 '.mat']);
%     EventDataFileName2 = ls(['*EVENTDATA_' dataType2 '.mat']);
    EventDataFileName1 = dir(['*EVENTDATA_' dataType1 '.mat']);
    EventDataFileName2 = dir(['*EVENTDATA_' dataType2 '.mat']);    
    if or(isempty(EventDataFileName1),isempty(EventDataFileName2))
        error(wraptext(['Missing EventData structure in current directory. '...
            'Change directory to location of file, move file to current '...
            'directory, or run "ExtractEventTriggeredData.m" to create RestData.']));
    else
%         LoadData = load(EventDataFileName1);
%         EventData.(dataType1) = LoadData.EventData.(dataType1);
%         LoadData = load(EventDataFileName2);
%         EventData.(dataType2) = LoadData.EventData.(dataType2);
%         BehData = EventData;
        LoadData = load(EventDataFileName1.name);
        EventData.(dataType1) = LoadData.EventData.(dataType1);
        LoadData = load(EventDataFileName2.name);
        EventData.(dataType2) = LoadData.EventData.(dataType2);
        BehData = EventData;
        clear EventData;
    end
end

%% Get the arrays for the calculation
[DataStruct1,FiltArray1] = SelectBehavioralEvents(BehData.(dataType1),Beh);
NormData1 = DataStruct1.NormData(FiltArray1,:);
[DataStruct2,FiltArray2] = SelectBehavioralEvents(BehData.(dataType2),Beh);
NormData2 = DataStruct2.NormData(FiltArray2,:);


%% Setup the data
test_inds = Event_Inds.TestStart:Event_Inds.Increment:size(NormData1,1);
if strcmp(Beh,'Rest')
    NormData1 = NormData1(test_inds);
    % Mean subtract the data
    Processed1 = cell(size(NormData1));
    for c = 1:length(NormData1)
        template = zeros(size(NormData1{c}));
        strt = 2*DataStruct1.Fs;
        stp = size(template,2);
        template(:,strt:stp) = NormData1{c}(:,strt:stp)-mean(NormData1{c}(:,strt:stp));
        Processed1{c} = template;
    end
    Data1 = Processed1;
    clear Processed1;
elseif strcmp(Beh,'VW')
    Data1_end = 5; 
    strt = round((DataStruct1.epoch.offset-1)*DataStruct1.Fs);
    stp = strt + round(Data1_end*DataStruct1.Fs);
    Data1 = zeros(size(NormData1(test_inds,:)));
    offset1 = mean(NormData1(test_inds,1:strt),2)*ones(1,stp-strt+1);
    Data1(:,strt:stp) = NormData1(test_inds,strt:stp)-offset1;
elseif strcmp(Beh,'Str')
    Data1_end = 7;
    strt = round((DataStruct1.epoch.offset-1)*DataStruct1.Fs);
    stp = strt + round(Data1_end*DataStruct1.Fs);
    Data1 = zeros(size(NormData1(test_inds,:)));
    offset1 = mean(NormData1(test_inds,1:strt),2)*ones(1,stp-strt+1);
    Data1(:,strt:stp) = NormData1(test_inds,strt:stp)-offset1;
else
    Data1_end = 1.5; 
    strt = round((DataStruct1.epoch.offset)*DataStruct1.Fs); 
    stp = strt + (Data1_end*DataStruct1.Fs);
    Data1 = zeros(size(NormData1(test_inds,:)));
    offset1 = mean(NormData1(test_inds,1:strt),2)*ones(1,stp-strt+1);
    Data1(:,strt:stp) = NormData1(test_inds,strt:stp)-offset1;
end

if strcmp(Beh,'Rest')
    NormData2 = NormData2(test_inds);
    % Mean subtract the data
    Processed2 = cell(size(NormData2));
    for c = 1:length(NormData2)
        template = zeros(size(NormData2{c}));
        strt = 2*DataStruct2.Fs;
        stp = size(template,2);
        offset = mean(NormData2{c})*ones(1,stp-strt+1);
        template(:,strt:stp) = detrend(NormData2{c}(:,strt:stp)-offset);
        Processed2{c} = template;
    end
    Data2 = Processed2;
    clear Processed2
elseif strcmp(Beh,'VW')
    Data2_end = 7;
    strt = round((DataStruct2.epoch.offset-1)*DataStruct2.Fs);
    stp = strt + round(Data2_end*DataStruct2.Fs);
    Data2 = zeros(size(NormData2(test_inds,:)));
    offset2 = mean(NormData2(test_inds,1:strt),2)*ones(1,stp-strt+1);
    Data2(:,strt:stp) = NormData2(test_inds,strt:stp)-offset2;
elseif strcmp(Beh,'Str')
    Data2_end = 7;
    strt = round((DataStruct2.epoch.offset-1)*DataStruct2.Fs);
    stp = strt + round(Data2_end*DataStruct2.Fs);
    Data2 = zeros(size(NormData2(test_inds,:)));
    offset2 = mean(NormData2(test_inds,1:strt),2)*ones(1,stp-strt+1);
    Data2(:,strt:stp) = NormData2(test_inds,strt:stp)-offset2;
else
    Data2_end = 3;
    strt = round(DataStruct2.epoch.offset*DataStruct2.Fs);
    stp = strt + (Data2_end*DataStruct2.Fs);
    Data2 = zeros(size(NormData2(test_inds,:)));
    offset2 = mean(NormData2(test_inds,1:strt),2)*ones(1,stp-strt+1);
    Data2(:,strt:stp) = NormData2(test_inds,strt:stp)-offset2;
end

%% Calculate R-squared on average data
if strcmp(Beh,'Rest')
    AveR2 = NaN;
else
    [Act,Pred] = ConvolveHRF(HRFs.HRF,mean(Data1),mean(Data2),0);
    mPred = Pred(strt:stp)-mean(Pred(strt:stp));
    mAct = Act(strt:stp)-mean(Act(strt:stp));
    AveR2 = CalculateRsquared(mPred,mAct);
end

%% Calculate R-squared on individual data
IndR2 = NaN*ones(1,size(Data2,1));
if strcmp(Beh,'Rest')
    if not(isfield(CBVPred,Beh))
        CBVPred.(Beh) = cell(size(Data2));
    end
    for tc = 1:length(Data2)
        strt = 2*DataStruct2.Fs;
        stp = length(Data2{tc});
        [Act,Pred] = ConvolveHRF(HRFs.HRF,detrend(Data1{tc}),...
            detrend(Data2{tc}),0);
        CBVPred.(Beh){tc} = Pred;
        mPred = Pred(strt:stp)-mean(Pred(strt:stp));
        mAct = Act(strt:stp)-mean(Act(strt:stp));
        IndR2(tc) = CalculateRsquared(mPred,mAct);
    end
    Med_IndR2 = median(IndR2);
else
    if not(isfield(CBVPred,Beh))
        CBVPred.(Beh) = NaN*ones(size(Data2));
    end
    for tc = 1:size(Data2,1)
        [Act,Pred] = ConvolveHRF(HRFs.HRF,Data1(tc,:),Data2(tc,:),0);
        CBVPred.(Beh)(tc,:) = Pred;
        mPred = Pred(strt:stp)-mean(Pred(strt:stp));
        mAct = Act(strt:stp)-mean(Act(strt:stp));
        IndR2(tc) = CalculateRsquared(mPred,mAct);
    end
    Med_IndR2 = median(IndR2);
end