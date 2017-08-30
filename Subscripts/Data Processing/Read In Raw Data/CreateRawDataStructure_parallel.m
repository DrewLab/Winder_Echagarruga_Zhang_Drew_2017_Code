function CreateRawDataStructure_parallel(filenames,TrackWhiskers)
%
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University
%
%   SUMMARY: Data acquired during trials must be in a form that Matlab can
%   work with easily. This code converts the various forms of data listed 
%   below into a matlab structure that can be easily manipulated. This code
%   requires a GPU to process the whisker images quickly.
%
%       .bin - Cameras
%       .tdms - Digital and Analog Data
%_______________________________________________________________
%   PARAMETER TYPE:             filenames - [cell array] list of filames
%                               with the extension _dalsa.bin
%
%                               TrackWhiskers - [binary] tells code whether
%                               to track the whiskers or not.
%_______________________________________________________________
%   RETURN:                     Output is a file saved to the current
%                               directory
%_______________________________________________________________

%% Load and setup

% Identify Files
if nargin == 0;
    filenames = uigetfile('*_fwire.bin','MultiSelect', 'on');
end

% Default
if nargin < 2;
    TrackWhiskers = 1;
end

% Control for single file instead of a list
if iscell(filenames) == 0;
    fileid = filenames;
    filenames = 1;
end

for fil = 1:length(filenames)
    close all;
    % Adapt to list or single file
    if iscell(filenames) == 1;
        indfile = filenames{fil};
    else
        indfile = fileid;
    end
    file_exist = ls(['*' indfile(1:17) '_rawdata.mat']);
    if not(isempty(file_exist));
        display('File already exists, continuing...'); continue;
    end
    %% Import .tdms data (All channels);
    trialdata = ReadInTDMS_WhiskerTrials([indfile(1:17) '.tdms']);
    
    % Transcribe Trial Notes
    Experimenter = trialdata.Experimenter;
    animal = trialdata.Animal_Name;
    hem = trialdata.Hemisphere_Recorded;
    Sol_pressure_psi = trialdata.Solenoid_Pressure_psi;
    Isoflurane_Time = trialdata.Isoflurane_Time;
    SessionID = trialdata.Session_ID;
    FlashTrial = trialdata.Flash_Trial;
    AmpGain = str2double(trialdata.Amplifier_Gain);
    Pinout = trialdata.Pinout;
    if strcmp(trialdata.Flash_Trial,'y')
        CBV_fr = trialdata.CBVCam_Frame_Rate/2;
        display(['Dalsa Frame Rate has been divided by 2. Current Frame Rate is: ' num2str(CBV_fr)]); %<- Temporary
    else
        CBV_fr = trialdata.CBVCam_Frame_Rate;
    end
    WhiskerCam_fr = trialdata.WhiskerCam_Frame_Rate;
    an_fs = trialdata.Analog_Sampling_Freq;
    WebCam_fr = trialdata.WebCam_Frame_Rate;
    Trial_Duration = trialdata.Trial_Duration;
   
    fswf = trialdata.Data(2,:); % force sensor, see pinout
    neuro = trialdata.Data(1,:);
    Sol = trialdata.Data(3,:);
    
    %% Start Whisker tracker
    if TrackWhiskers
        [WhiskerAngle] = whiskertracker_parallel(indfile(1:17));
        inds = isnan(WhiskerAngle)==1;
        WhiskerAngle(inds) = [];
    else
        WhiskerAngle = [];
    end

    %% Evaluate Data and Save
    RawData.Notes.Experimenter = Experimenter;
    RawData.Notes.animal = animal;
    RawData.Notes.hem = hem;
    RawData.Notes.Sol_pressure = Sol_pressure_psi;
    RawData.Notes.Isoflurane_Time = Isoflurane_Time;
    RawData.Notes.SessionID = SessionID;
    RawData.Notes.FlashTrial = FlashTrial;
    RawData.Notes.AmpGain = AmpGain;
    RawData.Notes.TrialDur = Trial_Duration;
    RawData.Notes.Pinout = {'Ai0','Ai1','Ai2','Ai3','Ai4','Ai5','Ai6','Ai7';...
        Pinout{:}};
    
    RawData.Data.Neuro = neuro/AmpGain; % Convert to volts
    RawData.Data.Sol = Sol;
    RawData.Data.fswf = fswf;
    RawData.Data.WhiskerAngle = WhiskerAngle;
    
    RawData.Fs.CBVCam = CBV_fr;
    RawData.Fs.WhiskCam = WhiskerCam_fr;
    RawData.Fs.Analog = an_fs;
    RawData.Fs.WebCam = WebCam_fr;

    save([animal '_' hem '_' indfile(1:17) '_rawdata'],'RawData')
end