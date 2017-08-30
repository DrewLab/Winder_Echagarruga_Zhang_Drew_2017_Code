function [] = SpatialCrossCorrelation(animal,DataType,Beh)
%   function [] = SpatialCrossCorrelation_CBVROIVsCBV(animal,ROIname,Beh)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: This function calculates the cross correlation between
%   "DataType" and each pixel of the CBV image. The CBV signal
%   across the entire window (the global signal) is subtracted from each
%   pixel. The cross correlation is calculated between "DataType" and this
%   subtracted signal.
%   
%_______________________________________________________________
%   PARAMETERS:             
%               animal - [string] the animal ID
%
%               DataType - [string] the fieldname of the *_RestData.mat
%               structure. These data are used for the cross correlation
%               with the pixel-wise CBV signal
%
%               Beh - [string] the behavioral category corresponding to the
%               data.
%_______________________________________________________________
%   RETURN:                     
%               The output of this function is a structure which is saved 
%               to the animal directory.
%_______________________________________________________________

%% Setup variables
orig_file_dir = uigetdir(pwd,'Select the location of the raw .bin files...');
maxlags = 10;
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
Fs = AllBehData.Fs;

%% Get criteria filter for the behavioral data
[BehFiltArray] = FilterEvents(AllBehData,BehCriteria);

%% Apply the filters
AllBehFileID = AllBehData.FileID;
BehData = AllBehData.Data(BehFiltArray,:);
BehTimes = AllBehData.EventTime(BehFiltArray)';
BehFileID = AllBehFileID(BehFiltArray)';
BehDates = AllBehData.FileDate(BehFiltArray)';
BehDurs = AllBehData.Duration(BehFiltArray);

%% Process the structure data
if iscell(BehData) 
    for BD = 1:length(BehData)
        ShortData = BehData{BD}(Buffer*Fs:end);
%         ShortData = BehData{BD}(Buffer*Fs:(Buffer*Fs+min(floor(BehDurs(BD)*Fs),300)));
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
FiltBehData = filtfilt(sos,g,CombinedBehData);
clear CombinedBehData;

%% Create a library of raw files and indexes
[uniqueDates,~,DateInds] = unique(BehDates);

%% PreAllocate and setup
AllCameraFramesStart = 1;
GlobalSignal = zeros(1,sum(BufferedDurs));


IndFrame = zeros(256,256,length(uniqueDates));
mask = zeros(size(IndFrame));
prevdir = cd(orig_file_dir);
for uD = 1:length(uniqueDates)
    UniqueDateInds = DateInds==uD;
    DateIDs = BehFileID(UniqueDateInds);
    DayBehTimes = BehTimes(UniqueDateInds);
    DayDurs = BufferedDurs(UniqueDateInds);
    [uniqueFileIDs,~,Uinds] = unique(DateIDs);
    
    % Create mask and register images
    IndFrame(:,:,uD) = GetSingleCBVFrame([uniqueFileIDs{1} '_dalsa.bin'], 256, 256);
    
    if uD == 1
        % Create mask of the window for the first day
        figure; imagesc(IndFrame(:,:,uD));
        colormap('gray');
        display('Draw ROI around the thinned skull window:')
        mask(:,:,uD) = roipoly;
        close gcf;
        [row,col] = find(mask(:,:,uD));
        numrows = max(row)-min(row)+1;
        numcols = max(col)-min(col)+1;
        AllCameraFrames = zeros(numrows,numcols,length(FiltBehData),'int16');
    else
        % Register the image to uD=1
        [~,RA,RB,~] = RegisterCBVImage(IndFrame(:,:,1),IndFrame(:,:,uD));
        
        % Shift the mask based on the registered images
        % xdimension: RB-RA>0, shift mask left; RB-RA<0, shift mask right.
        % ydimension: RB-RA>0, shift mask down; RB-RA<0, shift mask up.
        xshift = round(RB.XWorldLimits(1)-RA.XWorldLimits(1));
        yshift = round(RB.YWorldLimits(1)-RA.YWorldLimits(1));
        shiftrow = row-yshift;
        shiftcol = col-xshift;
        for sr = 1:length(shiftrow)
            mask(shiftrow(sr),shiftcol(sr),uD) = 1;
        end
    end
    
    [maskrow,maskcol] = find(mask(:,:,uD));
    minrow = min(maskrow);
    maxrow = max(maskrow);
    mincol = min(maskcol);
    maxcol = max(maskcol);
    
    % Match the behavior times and durations to events from each day
    
    
    %% Compile and process the camera frames
    for uFI = 1:length(uniqueFileIDs)
        FileID = [uniqueFileIDs{uFI} '_dalsa.bin'];
        
        % Match the behavior times and duration to the FileID
        EventMatch = Uinds==uFI;
        FileBehTimes = floor(DayBehTimes(EventMatch)*Fs);
%         FileBehDurs = min(DayDurs(EventMatch),300);
        FileBehDurs = DayDurs(EventMatch);
        
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
        
        % Clip the Frames around the ROI
        clippedFrames = Frames(minrow:maxrow,mincol:maxcol,:);
        
        % Normalize, mean subtract, and filter the frames
        NormFrames = zeros(size(clippedFrames));
        FileBehDurs = [0; FileBehDurs];
        for D = 2:length(FileBehDurs)
            StartFrame2 = sum(FileBehDurs(1:D-1))+1;
            StopFrame2 = sum(FileBehDurs(1:D));
            PeriodFrames = clippedFrames(:,:,StartFrame2:StopFrame2);
            PeriodMean = mean(PeriodFrames,3);
            NormPeriodFrames = zeros(size(PeriodFrames));
            for PF = 1:size(PeriodFrames,3)
                NormPeriodFrames(:,:,PF) = PeriodFrames(:,:,PF)./PeriodMean-1;
            end
            %Detrend the Resting period
            for Col = 1:size(NormPeriodFrames,2)
                NormPeriodFramesColumn = squeeze(NormPeriodFrames(:,Col,:));
                DTColumn = detrend(NormPeriodFramesColumn');
                NormPeriodFrames(:,Col,:) = DTColumn';
            end
            NormFrames(:,:,StartFrame2:StopFrame2) = NormPeriodFrames;
        end
        clear NormPeriodFrames;
        clear PeriodFrames;
        clear PeriodMean;
        clear Frames;
        
        % Low pass filter the detrended/normalized frames
        FiltNormFrames = zeros(size(NormFrames));
        for Col = 1:size(NormFrames,2)
            NormFramesColumn = squeeze(NormFrames(:,Col,:));
            FiltColumn = filtfilt(sos,g,NormFramesColumn');
            FiltNormFrames(:,Col,:) = FiltColumn';
        end
        clear NormFrames;
        
        % Subtract the Global Signal
        mask2 = double(mask(minrow:maxrow,mincol:maxcol,uD));
        GS_Subtracted = zeros(size(FiltNormFrames));
        PeriodGlobalSignal = zeros(size(FiltNormFrames,3),1);
        for FNF = 1:size(FiltNormFrames,3)
            IndFiltFrame = nonzeros(FiltNormFrames(:,:,FNF).*mask2);
            PeriodGlobalSignal(FNF) = mean(IndFiltFrame);
            GS_Subtracted(:,:,FNF) = FiltNormFrames(:,:,FNF)-...
                PeriodGlobalSignal(FNF).*mask2;
        end
        clear FiltNormFrames
        
        % Add the processed frames to cell array
        AllCameraFramesStop = AllCameraFramesStart+sum(FileBehDurs)-1;
        AllCameraFrames(:,:,AllCameraFramesStart:AllCameraFramesStop)...
            = int16(round(GS_Subtracted*MultiplicationFactor));
        GlobalSignal(AllCameraFramesStart:AllCameraFramesStop) = ...
            PeriodGlobalSignal;
        AllCameraFramesStart = AllCameraFramesStop+1;
        clear GS_Subtracted
        clear PeriodGlobalSignal
    end
end
display('Acquiring Camera Frames: Done')


%% Calculate cross correlation
[row,col] = find(mask(:,:,1));
clippedmask = mask(min(row):max(row),min(col):max(col),1);
[cliprow,clipcol] = find(clippedmask);
XC_Map = zeros(size(clippedmask,1),size(clippedmask,2),2*maxlags*Fs+1);
for cr = 1:length(cliprow)
    [XC_Map(cliprow(cr),clipcol(cr),:),lags] = ...
        xcorr(double(squeeze(AllCameraFrames(cliprow(cr),clipcol(cr),:))),FiltBehData,...
        maxlags*Fs,'coeff');
end
display('Calculating Cross Correlation: Done')

%% Save the data
cd(prevdir)
save([animal '_GlobalSubtractSpatialXCorr_' DataType '_' Beh '.mat'],'XC_Map','lags','mask',...
    'IndFrame');



