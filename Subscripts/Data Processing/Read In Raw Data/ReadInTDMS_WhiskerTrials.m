function [TDMSFile]=ReadInTDMS_WhiskerTrials(filename)
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University
%   Adapted from "ReadInAcquireTDMS.m", written by Patrick Drew. Jan 2016.
%
%   SUMMARY: Reads in .TDMS files from a LABVIEW acquisition program using
%   the convertTDMS.m script (acquired through Mathworks file exchange)
%   and organizes the acquired data into a structure.
%
%http://www.mathworks.com/matlabcentral/fileexchange/44206-converttdms--v10-
%_______________________________________________________________
%   PARAMETER TYPE:             
%                               filename - [string] list of filames
%                               with the extension _dalsa.bin
%_______________________________________________________________
%   RETURN:                     
%                               TDMSFile - [struct] contains measured
%                               analog data and trial notes from the
%                               LabVIEW acquisition program
%_______________________________________________________________

%load in the struct
[TempStruct,~]=convertTDMS(0,filename);

% Extract Whisker Camera info and delete from TempStruct
TDMSFile.Number_Dropped_WhiskerCam_Frames = 0;
TDMSFile.Dropped_WhiskerCam_Frame_Numbers = [];
for stct = length(TempStruct.Data.MeasuredData):-1:1 %Count backward to erase from the end and keep proper indexes
    if strcmp(TempStruct.Data.MeasuredData(stct).Name,'WhiskerCamBufferNumber')
        TDMSFile.Dropped_WhiskerCam_Frame_Numbers = TempStruct.Data.MeasuredData(stct).Data;
        TempStruct.Data.MeasuredData(stct) = [];
    elseif strcmp(TempStruct.Data.MeasuredData(stct).Name,'WhiskerCamBuffersDropped')
        TDMSFile.Number_Dropped_WhiskerCam_Frames = TempStruct.Data.MeasuredData(stct).Data;
        TempStruct.Data.MeasuredData(stct) = [];
    end
end
        
        

TDMSFile.Data=zeros(length(TempStruct.Data.MeasuredData),length(TempStruct.Data.MeasuredData(1).Data));
for k=1:length(TempStruct.Data.MeasuredData)
    if not(strcmp(TempStruct.Data.MeasuredData(k).Name(1:6),'Analog'))%!!!!!
        continue;
    end
    hold_data=TempStruct.Data.MeasuredData(1,k).Data;
    TDMSFile.Data(k,:)=hold_data;
end

% Add trial properties to the structure <- FUTURE, DISPLAY TABLE OF TRIAL
% PROPERTIES
TDMSFile.Experimenter=TempStruct.Data.Root.Experimenter;
TDMSFile.Animal_Name=TempStruct.Data.Root.Animal_Name;
TDMSFile.Hemisphere_Recorded=TempStruct.Data.Root.Hemisphere_Recorded;
TDMSFile.Solenoid_Pressure_psi = TempStruct.Data.Root.Solenoid_Pressure_psi;
TDMSFile.Isoflurane_Time=str2double(TempStruct.Data.Root.Isoflurane_Time);
if isfield(TempStruct.Data.Root,'Session_ID')
    TDMSFile.Session_ID=TempStruct.Data.Root.Session_ID;
elseif isfield(TempStruct.Data.Root,'Session_Number')
    TDMSFile.Session_ID=TempStruct.Data.Root.Session_Number;
end
TDMSFile.Flash_Trial = TempStruct.Data.Root.Flash_Trial;
if isfield(TempStruct.Data.Root,'Amplifier_Gain');
    TDMSFile.Amplifier_Gain = TempStruct.Data.Root.Amplifier_Gain;
else
    TDMSFile.Amplifier_Gain = 10000;
end
if isfield(TempStruct.Data.Root,'Ai0');
    TDMSFile.Pinout = {TempStruct.Data.Root.Ai0, TempStruct.Data.Root.Ai1,...
        TempStruct.Data.Root.Ai2, TempStruct.Data.Root.Ai3, ...
        TempStruct.Data.Root.Ai4, TempStruct.Data.Root.Ai5, ...
        TempStruct.Data.Root.Ai6, TempStruct.Data.Root.Ai7};
end
if isfield(TempStruct.Data.Root,'CBVCam_Frame_Rate')
    TDMSFile.CBVCam_Frame_Rate=str2double(TempStruct.Data.Root.CBVCam_Frame_Rate);
elseif isfield(TempStruct.Data.Root,'Dalsa_Frame_Rate')
    TDMSFile.CBVCam_Frame_Rate=30; %!!!!str2double(TempStruct.Data.Root.Dalsa_Frame_Rate);
end
if isfield(TempStruct.Data.Root,'WhiskerCam_Frame_Rate')
    TDMSFile.WhiskerCam_Frame_Rate=str2double(TempStruct.Data.Root.WhiskerCam_Frame_Rate);
elseif isfield(TempStruct.Data.Root,'Basler_Frame_Rate');
    TDMSFile.WhiskerCam_Frame_Rate=150; %!!!!str2double(TempStruct.Data.Root.Basler_Frame_Rate);
end
TDMSFile.Analog_Sampling_Freq=20000; %!!!!str2double(TempStruct.Data.Root.Analog_Sampling_Freq);
TDMSFile.WebCam_Frame_Rate=15; %!!!str2double(TempStruct.Data.Root.WebCam_Frame_Rate);
if isfield(TempStruct.Data.Root,'Trial_Duration')
    TDMSFile.Trial_Duration=str2double(TempStruct.Data.Root.Trial_Duration);
else
    TDMSFile.Trial_Duration = 300; %!!!size(TDMSFile.Data,2)/TDMSFile.Analog_Sampling_Freq;
end

