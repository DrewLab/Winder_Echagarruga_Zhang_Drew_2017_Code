function ProcessRawDataFile_atw2(filename)
%%  [] = ProcessRawDataFile_atw2(filename)
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   SUMMARY: Saves processed forms of measured data to the _rawdata.mat
%   file in order to reduce run times in later scripts and increase
%   consistency. No data in the appended data will be normalized.
%_______________________________________________________________
%   PARAMETERS:
%                           filename: [str] name of _rawdata.mat file to be
%                           processed
%
%   RETURN:
%                           Function output is a saved structure called
%                           'ProcData.mat'
%_______________________________________________________________

%% LOAD AND SETUP
display(['Processing RawData File: ' filename]);
load(filename);
underscore_indexes = strfind(filename,'_');
animal = filename(1:(underscore_indexes(1)-1));
hem = filename(underscore_indexes(1)+1:underscore_indexes(2)-1);
date_index = underscore_indexes(2)+1;
numday = filename(date_index:date_index+5);
strday = ConvertDate(numday);

% Desired frequencies of processed data
SR_fs = RawData.Fs.CBVCam;                                                                % <- Ln 34: Final sampling frequency of the spike rate                   
fswf_fs = RawData.Fs.CBVCam;                                                               % <- Ln 35: Final sampling frequency of the force sensor
wwf_fs = RawData.Fs.CBVCam;


ProcData.Info.TrialDur = RawData.Notes.TrialDur;

% Filtering parameters for whisking and force sensor
whisk_filt_thresh = 20;
whisk_filt_ord = 2;
fs_filt_thresh = 20;
fs_filt_ord = 2;
ExpectedLength = RawData.Notes.TrialDur*RawData.Fs.Analog;

%% SAVE SOLENOID TIMES (IN SECONDS)

% Identify the solenoids by amplitude
[pp_id,sp_id,tp_id,cp_id] = IdentifySolenoids(RawData.Notes,hem);

ProcData.Sol.Contra = find(diff(RawData.Data.Sol) == pp_id)/RawData.Fs.Analog;
ProcData.Sol.Ipsi = find(diff(RawData.Data.Sol) == sp_id)/RawData.Fs.Analog;
ProcData.Sol.Tail = find(diff(RawData.Data.Sol) == tp_id)/RawData.Fs.Analog;
ProcData.Sol.Control = find(diff(RawData.Data.Sol) == cp_id)/RawData.Fs.Analog;
ProcData.Fs.Sol_fs = RawData.Fs.Analog;

%% CBV -> 

% CBV from ROIs
CBVfields = fieldnames(RawData.Data.CBV);
for field = 1:length(CBVfields)
    ProcData.CBV.([CBVfields{field}]) = RawData.Data.CBV.(CBVfields{field});
    ProcData.Fs.([CBVfields{field}]) = RawData.Fs.CBVCam;
end

%% NEURAL DATA -> PROCESS INTO VARIOUS FORMS:

if isfield(RawData.Data,'Neuro')
    
    % Wide-Band LFP
    [ProcData.Neuro.Wideband_LFP, ProcData.Fs.Wideband_LFP] = ...
        ProcessNeuro(RawData,'Wideband_LFP',animal,hem,numday);
    
    % MU Band [300-3000]
    [ProcData.Neuro.MUpower, ProcData.Fs.MUpower] = ProcessNeuro(RawData,...
        'MUpower',animal,hem,numday);
    
    % Gamma Band [40-100]
    [ProcData.Neuro.Gam, ProcData.Fs.Gam] = ProcessNeuro(RawData,'Gam',...
        animal,hem,numday);
    
    % Beta [10-30 Hz]
    [ProcData.Neuro.Beta, ProcData.Fs.Beta] = ProcessNeuro(RawData,'Beta',...
        animal,hem,numday);
    
    % Sub-Alpha Frequencies
    [ProcData.Neuro.SubAlpha, ProcData.Fs.SubAlpha] = ProcessNeuro(RawData,...
        'SubAlpha',animal,hem,numday);
end
%% WHISKER ANGLE -> BINARIZE, SET RESTING ANGLE TO ZERO DEGREES

% Track dropped basler frames
ExpectedFrames = RawData.Notes.TrialDur*RawData.Fs.WhiskCam;
ProcData.Info.DroppedBaslerFrames = ExpectedFrames-length(RawData.Data.WhiskerAngle);
if ProcData.Info.DroppedBaslerFrames>15
    display([filename ' Dropped ' num2str(ProcData.Info.DroppedBaslerFrames) ' Frames.'])
end
% Trim any additional frames, for resample
WhiskerAngle = RawData.Data.WhiskerAngle(1:min(ExpectedFrames,...
    length(RawData.Data.WhiskerAngle)));

% Create filter for whisking/movement
[z,p,k] = butter(whisk_filt_ord,whisk_filt_thresh/(RawData.Fs.WhiskCam/2),'low');
[sos,g] = zp2sos(z,p,k);
filt_wwf = filtfilt(sos,g,WhiskerAngle-mean(WhiskerAngle));
res_wwf = resample(filt_wwf,wwf_fs,RawData.Fs.WhiskCam);

% Binarize the whisker waveform (wwf)
Threshfile = dir('*Thresholds.mat');
if not(isempty(Threshfile));
    load(Threshfile.name)
end
[ok] = Check4Threshold(['BinWWF_Lower_' strday], animal, hem);
if ok == 0;
    [wwf_thresh1,wwf_thresh2] = CreateWhiskThreshold(res_wwf, wwf_fs, strday);
    Thresholds.(['BinWWF_Lower_' strday]) = wwf_thresh1;
    Thresholds.(['BinWWF_Upper_' strday]) = wwf_thresh2;
    save([animal '_' hem '_Thresholds.mat'],'Thresholds');
end
load([animal '_' hem '_Thresholds.mat']);
bin_wwf = binarize_wwf(res_wwf,wwf_fs,Thresholds.(['BinWWF_Lower_' strday]),...
    Thresholds.(['BinWWF_Upper_' strday]));
[linked_bin_wwf] = link_binary_events(bin_wwf,[round(wwf_fs/3), ...
    0]);

inds = linked_bin_wwf == 0;
rest_angle = mean(res_wwf(inds));

ProcData.Beh.wwf = res_wwf - rest_angle;
ProcData.Fs.wwf = wwf_fs;
ProcData.Beh.Bin_wwf = bin_wwf;
ProcData.Fs.Bin_wwf = wwf_fs;

%% FORCE SENSOR -> BINARIZE, DOWNSAMPLE TO DESIRED FREQUENCY

% Trim any additional data points for resample
trimmedpswf = RawData.Data.fswf(1:min(ExpectedLength,length(RawData.Data.fswf)));

% Filter then downsample the force sensory waveform to desired frequency 
[z,p,k] = butter(fs_filt_ord,fs_filt_thresh/(RawData.Fs.Analog/2),'low');
[sos,g] = zp2sos(z,p,k);
filtpswf = filtfilt(sos,g,trimmedpswf);
ProcData.Beh.pswf = resample(filtpswf,fswf_fs,RawData.Fs.Analog);
ProcData.Fs.pswf = fswf_fs;

% Binarize the force sensor waveform
[ok] = Check4Threshold(['BinPSWF_' strday], animal, hem);
if ok == 0;
    [fswf_thresh] = CreatePswfThreshold(ProcData.Beh.pswf);
    Thresholds.(['BinPSWF_' strday]) = fswf_thresh;
    save([animal '_' hem '_Thresholds.mat'],'Thresholds');
end
ProcData.Beh.Bin_pswf = binarize_pswf(ProcData.Beh.pswf,Thresholds.(['BinPSWF_' strday]));
ProcData.Fs.Bin_pswf = fswf_fs;

%% SAVE THE PROCESSED DATA STRUCTURE

save([filename(1:(underscore_indexes(end)-1)) '_ProcData.mat'],'ProcData');