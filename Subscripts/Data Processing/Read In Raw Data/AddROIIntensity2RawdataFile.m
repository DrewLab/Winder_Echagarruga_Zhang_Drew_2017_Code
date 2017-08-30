function [] = AddROIIntensity2RawdataFile(ROIname)
%   [] = AddROIIntensity2RawdataFile(ROIname)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Create an ROI from camera frames and average the
%   reflectance from the ROI to get a timeseries. Add the result into the
%   *rawdata.mat file.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   ROIname - [string] Description of ROI                        
%_______________________________________________________________
%   RETURN:                     
%                               
%_______________________________________________________________
close all
% Setup
filenames = uigetfile('*rawdata.mat','multiselect','on');
orig_dir = uigetdir(pwd,'Select the directory for the original files');
for f = 1:length(filenames)
    filename = filenames{f};
    [animal,hem,filedate,fileID] = GetFileInfo(filename);
    ROIFile = ls('*ROIs.mat');
    if not(isempty(ROIFile))
        load(ROIFile)
    else
        ROIs = [];
    end
    
    strday = ConvertDate(filedate);
    load(filename)
    
    [isok] = CheckROI([ROIname '_' strday] ,animal,hem);
    prevdir = cd(orig_dir);
    [Frames]=ReadDalsaBinary([fileID '_dalsa.bin'],256,256);
    cd(prevdir);
    if not(isok)
        mask = CreateROI(Frames{2},[ROIname '_' strday],animal,hem);
    elseif isok
        mask = GetROI(Frames{1},[ROIname '_' strday],animal,hem);
    end
    
    mean_intensity = Bin2Intensity([fileID '_dalsa.bin'],mask,Frames);
    RawData.Data.CBV.(ROIname) = mean_intensity;
    save(filename,'RawData')
end

