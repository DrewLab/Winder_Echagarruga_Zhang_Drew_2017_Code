function [HRFs] = CalculateHRF_Deconvolution(dataType1,dataType2, Beh)
%   function [HRFs] = CalculateHRF_Deconvolution(dataType1,dataType2, Beh)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Finds the hemodynamic response function between two data
%   measurements during a given behavioral condition. HRF is found using
%   the deconvolution method. It also fits the HRF using a gamma-variate
%   distribution (based on fminsearch).
%   
%_______________________________________________________________
%   PARAMETERS:   
%               dataType1 - [string] the fieldname of the first data
%               measurement
%
%               dataType2 - [string] the fieldname of the second data
%               measurement
%
%               Beh - [string] the behavioral condition                               
%_______________________________________________________________
%   RETURN:                     
%               HRFs - [struct] a structure of the HRF and various 
%               pertinent parameters
%_______________________________________________________________

%% Load and Setup
HRFParams.offset = 3; 
HRFParams.dur = 10; 
Event_Inds.CalcStart = 1;
Event_Inds.TestStart = 2;
Event_Inds.Increment = 2;


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

%% Separate events for HRF calculation from events used for later testing.
% Insert padding of zeros with size equal to the HRF between individual
% events.
zpad_inds = Event_Inds.TestStart:Event_Inds.Increment:size(NormData1,1);
calc_inds = Event_Inds.CalcStart:Event_Inds.Increment:size(NormData1,1);
zpad1 = zeros(1,HRFParams.dur*DataStruct1.Fs);
zpad2 = zeros(1,HRFParams.dur*DataStruct2.Fs);


if strcmp(Beh,'Rest')
    NormData1(zpad_inds) = {zpad1};
    % Mean subtract the data
    Processed1 = cell(size(NormData1));
    for c = 1:length(NormData1)
        template = zeros(size(NormData1{c}));
        strt = 2*DataStruct1.Fs;
        stp = size(template,2);
        template(:,strt:stp) = detrend(NormData1{c}(:,strt:stp)-mean(NormData1{c}(:,strt:stp)));
        Processed1{c} = template;
    end
    Data1 = [Processed1{:}];
    clear Processed1;
elseif strcmp(Beh,'VW')
    Data1_end = 4;
    strt = round((DataStruct1.epoch.offset-1)*DataStruct1.Fs); 
    stp = strt + round(Data1_end*DataStruct1.Fs);
    template = zeros(size(NormData1(calc_inds,:)));
    offset1 = mean(NormData1(calc_inds,1:strt),2)*ones(1,stp-strt+1);
    template(:,strt:stp) = NormData1(calc_inds,strt:stp)-offset1;
    Data1Pad = [template ones(length(calc_inds),1)*zpad1];
    Data1 = reshape(Data1Pad',1,numel(Data1Pad));
else
    Data1_end = 1.5;
    strt = round((DataStruct1.epoch.offset)*DataStruct1.Fs); 
    stp = strt + (Data1_end*DataStruct1.Fs);
    template = zeros(size(NormData1(calc_inds,:)));
    offset1 = mean(NormData1(calc_inds,1:strt),2)*ones(1,stp-strt+1);
    template(:,strt:stp) = NormData1(calc_inds,strt:stp)-offset1;
    Data1Pad = [template ones(length(calc_inds),1)*zpad1];
    Data1 = reshape(Data1Pad',1,numel(Data1Pad));
end


if strcmp(Beh,'Rest')
    NormData2(zpad_inds) = {zpad2};
    Processed2 = cell(size(NormData2));
    for c = 1:length(NormData2)
        template = zeros(size(NormData2{c}));
        strt = 2*DataStruct1.Fs;
        stp = size(template,2);
        offset = mean(NormData2{c})*ones(1,stp-strt+1);
        template(:,strt:stp) = detrend(NormData2{c}(:,strt:stp)-offset);
        Processed2{c} = template;
    end
    Data2 = [Processed2{:}];
    clear Processed2
elseif strcmp(Beh,'VW')
    Data2_end = 6;
    strt = round((DataStruct2.epoch.offset-1)*DataStruct2.Fs);
    stp = strt + round(Data2_end*DataStruct2.Fs);
    template = zeros(size(NormData2(calc_inds,:)));
    offset2 = mean(NormData2(calc_inds,1:strt),2)*ones(1,stp-strt+1);
    template(:,strt:stp) = NormData2(calc_inds,strt:stp)-offset2;
    Data2Pad = [template ones(length(calc_inds),1)*zpad2];
    Data2 = reshape(Data2Pad',1,numel(Data2Pad));
else
    Data2_end = 3;
    strt = round(DataStruct2.epoch.offset*DataStruct2.Fs);
    stp = strt + (Data2_end*DataStruct2.Fs);
    template = zeros(size(NormData2(calc_inds,:)));
    offset2 = mean(NormData2(calc_inds,1:strt),2)*ones(1,stp-strt+1);
    template(:,strt:stp) = NormData2(calc_inds,strt:stp)-offset2;
    Data2Pad = [template ones(length(calc_inds),1)*zpad2];
    Data2 = reshape(Data2Pad',1,numel(Data2Pad));
end

%% Calculate HRF based on deconvolution
Fs = DataStruct2.Fs;
IR_est=IR_analytic(Data1',Data2',HRFParams.offset*Fs,HRFParams.dur*Fs);

HRF = sgolayfilt(IR_est.IR',3,Fs+1);
timevec = (1:length(HRF))/Fs-HRFParams.offset;
num_events = length(calc_inds);
currdate = date;

HRFs.HRF = HRF;
HRFs.timevec = timevec;
HRFs.Fs = Fs;
HRFs.HRFParams = HRFParams;
HRFs.num_calc_events = num_events;
HRFs.Event_Inds = Event_Inds;
HRFs.calc_date = currdate;

%% Calculate the gamma HRF
options = optimset('MaxFunEvals',2e3,'MaxIter',2e3,'TolFun',1e-7,'TolX',1e-7);
initvals = [-1e-3, 0.75, 1];
HRFDur = 5; % seconds
[gam_params,~,~] = fminsearch(@(x)gammaconvolve(x,Data1,Data2,...
    DataStruct2.Fs,HRFDur),initvals,options);
t = 0:1/DataStruct2.Fs:HRFDur;
a = ((gam_params(2)/gam_params(3))^2*8*log10(2));
beta = ((gam_params(3)^2)/gam_params(2)/8/log10(2));
gamma = gam_params(1)*(t/gam_params(2)).^a.*exp((t-gam_params(2))/(-1*beta));
HRFs.GammaHRF = gamma;
HRFs.GammaTime = t;
