function [SpatStruct] = EventTriggeredCBVSpatialMap(animal,CBVType,Beh,Day,Poly,FileDir)
%   function [] = EventData2SpatialMap_atw1(CBVType,Beh)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Averages frames around a behavioral event and plots the
%   mean pixel intensity over a user-defined time interval as an image. 
%_______________________________________________________________
%   PARAMETERS:             
%               animal - [string] animal ID, used to find t
%
%               CBVType - [string] name of CBV region of interest
%
%               Beh - [string] behavioral category
%_______________________________________________________________
%   RETURN:                     
%               No output, function saves a structure and figure                  
%_______________________________________________________________

%% Load and Setup
prevdir = cd([animal filesep]);
% EFile = ls(['*_EVENTDATA_' CBVType '.mat']);
% load(EFile);
EFile = dir(['*_EVENTDATA_' CBVType '.mat']);
load(EFile.name);
[DataStruct,FiltArray] = SelectBehavioralEvents(EventData.(CBVType),Beh);
FileDates = DataStruct.FileDate(FiltArray);
FileIDs = DataStruct.FileID(FiltArray);
EventTimes = DataStruct.EventTime(FiltArray);
Fs = DataStruct.Fs;
% [animal,hem,~,~] = GetFileInfo(EFile);

% Find individual days
days = unique(FileDates);
if isempty(Day)
    [selection,~] = listdlg('ListString',days,'SelectionMode','single','Name',...
        'Choose Day to analyze:');
    ind_day = days{selection};
else
    ind_day = Day;
end
% strday = ConvertDate(ind_day);
matches = strcmp(FileDates,{ind_day});
inds = find(matches);
prevfile = [];

% Set Color map bounds
if strcmp(Beh,'Contra')
    mean_cmap_bounds = [0.97 1.03];
elseif strcmp(Beh,'VW')
    mean_cmap_bounds = [0.99 1.01];
end

if isempty(FileDir)
    orig_file_dir = uigetdir(pwd,'Choose the location of the original CBV camera files...');
else
    orig_file_dir = FileDir;
end

time_intervals = [0.5,1.5; 3,4];

for i = 1:sum(matches)
    currfile = FileIDs{inds(i)};
    if not(strcmp(prevfile,currfile))
        % convert the .bin file to frames
        binname = [orig_file_dir '\' currfile '_dalsa.bin'];
        [frames]=ReadDalsaBinary_Matrix(binname, 256, 256);
    end
    
    if i==1 && isempty(Poly)
        display('Outline the edges of the window...');
        ReferenceFrame = frames(:,:,1);
        figure; imagesc(ReferenceFrame); axis square; colormap('gray')
        [mask, xi, yi] = roipoly;
        display('ROI Entered')
        NormIntervalMeans = zeros(size(frames,1),size(frames,2),...
            size(time_intervals,1),sum(matches));
        close(gcf);
    elseif i==1 && not(isempty(Poly))
        mask = poly2mask(Poly.x,Poly.y,256,256);
        xi = Poly.x;
        yi = Poly.y;
    end
        
    startframe = round(EventTimes(inds(i))*Fs);
    Baseline_Start_Frame = startframe - round((DataStruct.epoch.offset*Fs));
    Baseline_Stop_Frame = startframe;
    Baseline_Image = mean(frames(:,:,Baseline_Start_Frame:Baseline_Stop_Frame),3);
    % Chunk Data and frames
    for ti = 1:size(time_intervals,1)
        sample_intervals = time_intervals(ti,:)*Fs;
        Interval_Start = startframe+sample_intervals(1);
        Interval_Stop = startframe+sample_intervals(2);
        frame_inds = Interval_Start:Interval_Stop;
        IntervalMean = mean(frames(:,:,frame_inds),3);
        NormIntervalMeans(:,:,ti,i) = IntervalMean./Baseline_Image;
    end
    prevfile = currfile;
end
SpatStruct.Labels = {'PixelRows','PixelColumns','TimeInterval','DailyEventNumber'};
SpatStruct.Data = NormIntervalMeans;
SpatStruct.Mask = mask;
SpatStruct.xi = xi;
SpatStruct.yi = yi;
SpatStruct.FileDir = orig_file_dir;
SpatStruct.Date = ind_day;
SpatStruct.TimeIntervals = time_intervals;
% save(['..\Processed Data\' animal '_' hem '_' strday '_' Beh '_SpatStruct.mat'],...
%     'SpatStruct','baseline','SpatStruct_Labels','ReferenceFrame','time_intervals','mold');



cd(prevdir);