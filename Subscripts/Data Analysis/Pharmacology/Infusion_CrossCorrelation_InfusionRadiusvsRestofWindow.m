function [] = Infusion_CrossCorrelation_InfusionRadiusvsRestofWindow(animal,DataType,Beh)
%   [] = Infusion_CrossCorrelation_MuscDiffROIvsRestofWindow(animal,DataType,Beh)
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

%% Setup variables
TimeThresh = 45; %min
orig_file_dir = '..\..\Original Files\';
maxlags = 10;
Threshfile = dir('*Thresholds.mat');
load(Threshfile.name);
MultiplicationFactor = 1e5; % Used to allow a large double array to be saved as integer, to conserve memory


%% Identify behavior, load the behavioral data, create data filter criteria
if strcmp(Beh,'Rest')
    RestFile = [animal '_LH_RestData_' DataType '.mat'];
    if ~isempty(RestFile)
        load(RestFile)
    else
        error(['SpatialCrossCorrelation_CBVROIVsCBV: ' animal ...
            '_LH_RestData_' DataType '.mat - No such file found'])
    end
    AllBehData = RestData.(DataType);
    clear RestData;
    
    % Create criteria to filter the behavioral data
    BehCriteria.Fieldname = {'Duration','PuffDistance'};
    BehCriteria.Comparison = {'gt','gt'};
    BehCriteria.Value = {14,5};
    Buffer = 4;
else    
    % Load the eventfile
    AllBehFileID = [animal '_LH_EventData_ ' DataType '.mat'];
    if ~isempty(AllBehFileID)
        load(AllBehFileID)
    else
        error(['SpatialCrossCorrelation_CBVROIVsCBV: ' animal ...
            '_LH_EventData_' DataType '.mat - No such file found'])
    end
    AllBehData = EventData.(DataType).DataStructCategory;
    clear EventData;
    
    % Create criteria to filter the behavioral data
    BehCriteria.Fieldname = {'WhiskScore_Pre','WhiskScore_Post',...
        'MoveScore_Pre','MoveScore_Post'};
    BehCriteria.Comparison = {'lt','lt','lt','lt'};
    BehCriteria.Value = {0.1,0.33,0.1,0.33};
    Buffer = 0;
end

%% Get criteria filter for the behavioral data
[BehFiltArray] = FilterEvents(AllBehData,BehCriteria);

%% Calculate the time elapsed since infusion for each behavior and get filter
Fs = AllBehData.Fs;
FileDates = AllBehData.FileDate(BehFiltArray);
BehOnset = AllBehData.EventTime(BehFiltArray);
FileIDs = AllBehData.FileID(BehFiltArray);
Data = AllBehData.Data(BehFiltArray,:);
UniqueDates = unique(FileDates);
TimeFilt = false(size(FileDates));
ROInames = cell(size(FileDates));
for UD = 1:length(UniqueDates)
    DateMatch = strcmp(FileDates,UniqueDates{UD});
    strdate = ConvertDate(UniqueDates{UD});
    InfusionTime = Thresholds.Infusion_Start_Times.(strdate);
    TimeSinceInfusion = ElapsedTime(FileIDs(DateMatch),BehOnset(DateMatch),InfusionTime);
    TimeFilt(DateMatch) = gt(TimeSinceInfusion,TimeThresh);
    ROInames(DateMatch) = {[DataType '_' strdate]};
end

% Fs = AllBehData.Fs;
% AllBehFileID = AllBehData.FileID;
% AllBehTimes = round(AllBehData.EventTime);
% TimeSinceInfusion = ElapsedTime(AllBehFileID,AllBehTimes,InfusionTime);
% BehTimeFilt = gt(TimeSinceInfusion,TimeThresh);

%% Apply the filters
% BehData = AllBehData.Data(BehFiltArray&BehTimeFilt,:);
% BehTimes = AllBehData.EventTime(BehFiltArray&BehTimeFilt)';
% BehFileID = AllBehFileID(BehFiltArray&BehTimeFilt)';
BehData = Data(TimeFilt,:);
BehTimes = BehOnset(TimeFilt)';
BehFileID = FileIDs(TimeFilt)';
BehDates = FileDates(TimeFilt)';

%% Process the structure data
if iscell(BehData) 
    for BD = 1:length(BehData)
        ShortData = BehData{BD}(Buffer*Fs:end);
        BehData{BD} = detrend(ShortData);
    end
    BufferedDurs = cellfun('length',BehData);
    CombinedBehData = [BehData{:}];
else
    BehData = BehData(:,Buffer*Fs+1:end)';
    BufferedDurs = size(BehData,2)*ones(size(BehData,1),1);
    CombinedBehData = BehData(:);
end

[z,p,k] = butter(4,3/(Fs/2),'low');
[sos,g] = zp2sos(z,p,k);
InsideCBV = filtfilt(sos,g,CombinedBehData);
DataLen = length(InsideCBV);
clear CombinedBehData;

ROIfile = dir('*ROIs.mat');
load(ROIfile.name);

UniqueDates = unique(BehDates);
% PreAllocate an array of frames
OutsideCBV = zeros(1,DataLen);
AllCBVStart = 1;
for UD = 1:length(UniqueDates);
    DayFilt = strcmp(BehDates,UniqueDates{UD});
    DayIDs = BehFileID(DayFilt);
    DayROI = ROInames(DayFilt);
    % Create a library of raw files and indexes
    [uniqueFileIDs,~,Uinds] = unique(DayIDs);

    % Create ROI for the rest of the window
    prevdir = cd(orig_file_dir);
    IndFrame = GetSingleCBVFrame([uniqueFileIDs{1} '_dalsa.bin'], 256, 256);
    figure; imagesc(IndFrame);
    colormap('gray');
    impoly(gca, [ROIs.(DayROI{1}).xi ROIs.(DayROI{1}).yi]);
    display('Draw ROI around rest of window: ')
    mask = roipoly;
    close gcf;

    %% Compile and process the camera frames
    for uFI = 1:length(uniqueFileIDs)
        FileID = [uniqueFileIDs{uFI} '_dalsa.bin'];
        
        % Match the behavior times and duration to the FileID
        EventMatch = Uinds==uFI;
        FileBehTimes = floor(BehTimes(EventMatch)*Fs);
        FileBehDurs = BufferedDurs(EventMatch);
        
        % Create a list of frame numbers to import from file
        FrameInds_cell = cell(length(FileBehTimes),1);
        for FBT = 1:length(FileBehTimes)
            StartFrame1 = round(FileBehTimes(FBT)+Buffer*Fs);
            StopFrame1 = round(StartFrame1+FileBehDurs(FBT)-1);
            FrameInds_cell{FBT} = StartFrame1:StopFrame1;
        end
        FrameInds = [FrameInds_cell{:}];
        
        % Import the frames from the dalsa file
        t1 = tic;
        Frames = GetCBVFrameSubset(FileID, 256, 256, FrameInds);
        elapsed = toc(t1);
        display([num2str(elapsed) ' seconds to fetch ' ...
            num2str(length(FrameInds)) ' frames.'])
        
%         CBVother = zeros(1,size(Frames,3));
%         FileBehDurs = [0; FileBehDurs];
        NormFile = cell(1,length(FileBehDurs));
        for D = 1:length(FileBehDurs)
            StartFrame2 = sum(FileBehDurs(1:D-1))+1;
            StopFrame2 = sum(FileBehDurs(1:D));
            PeriodFrames = Frames(:,:,StartFrame2:StopFrame2);
            PeriodPixelMean = mean(PeriodFrames,3);
            PeriodROIMean = mean(PeriodPixelMean(mask));
            NormPeriod = zeros(1,size(PeriodFrames,3));
            for PF = 1:size(PeriodFrames,3)
                indFrame = PeriodFrames(:,:,PF);
                NormPeriod(PF) = mean(indFrame(mask))./PeriodROIMean-1;
            end
            %Detrend the Resting period
            NormFile{D} = detrend(NormPeriod);
        end
        clear PeriodFrames;
        clear Frames;
        
        % Low pass filter the detrended/normalized frames
        FiltCBV = filtfilt(sos,g,[NormFile{:}]);
        
        % Add the processed frames to cell array
        AllCBVStop = AllCBVStart+sum(FileBehDurs)-1;
        OutsideCBV(AllCBVStart:AllCBVStop) = FiltCBV;
        AllCBVStart = AllCBVStop+1;
    end
end

display('Acquiring Camera Frames: Done')


%% Calculate cross correlation
[CrossROI,lags] = xcorr(InsideCBV,OutsideCBV,maxlags*Fs,'coeff');
figure; plot(lags/30,CrossROI);
%% Save the data
cd(prevdir)
save([animal '_CBVROIXCorr_' Beh '.mat'],'CrossROI','lags','mask',...
    'IndFrame','DataLen');
