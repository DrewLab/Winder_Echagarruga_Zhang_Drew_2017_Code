function ProcessRawDataFile_V1(filename)
%%  [] = ProcessRawData()
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Jul 2014
%   Version 1
%
%   SUMMARY: Saves processed forms of measured data to the _rawdata.mat
%   file in order to reduce run times in later scripts and increase
%   consistency. No data in the appended data will be normalized.
%_______________________________________________________________
%   INPUTS:
%                           filename: [str] name of _rawdata.mat file to be
%                           processed
%
%
%_______________________________________________________________
%   OUTPUTS:
%_______________________________________________________________
%   REQUIRED SCRIPTS:
%_______________________________________________________________
%   CALLED BY:
%_______________________________________________________________
%   FUTURE VERSIONS:
%_______________________________________________________________
%   CHANGES FROM PREV VERS:
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
SR_fs = RawData.dal_fr;                                                                % <- Ln 34: Final sampling frequency of the spike rate                   
pswf_fs = RawData.dal_fr;                                                               % <- Ln 35: Final sampling frequency of the force sensor
wwf_fs = RawData.dal_fr;
ProcData.TrialDur = RawData.TrialDur;
whisk_filt_thresh = 20;
whisk_filt_ord = 2;
ps_filt_thresh = 20;
ps_filt_ord = 2;
ExpectedLength = RawData.TrialDur*RawData.an_fs;

%% SAVE SOLENOID TIMES (IN SECONDS)

% Identify the solenoids by amplitude
if RawData.TrialDur> 60
    % Extended acquisition program (includes upgrades that changed solenoid
    % numbers)
    [pp_id,sp_id,tp_id,cp_id] = IdentifySolenoids2(hem);
else
    % Previous acquisition program (short trials, ~ 40 seconds)
    [pp_id,sp_id,tp_id,cp_id] = IdentifySolenoids_atw0(hem); 
end

ProcData.Sol.Contra = find(diff(RawData.Sol) == pp_id)/RawData.an_fs;
ProcData.Sol.Ipsi = find(diff(RawData.Sol) == sp_id)/RawData.an_fs;
ProcData.Sol.Tail = find(diff(RawData.Sol) == tp_id)/RawData.an_fs;
ProcData.Sol.Control = find(diff(RawData.Sol) == cp_id)/RawData.an_fs;
ProcData.Fs.Sol_fs = RawData.an_fs;

%% CBV -> 

% CBV from whisker barrel ROI
ProcData.CrossCorrROI = RawData.CrossCorrROI;
ProcData.Fs.CrossCorrROI_fs = RawData.dal_fr;

% Heart Rate
% ProcData.HR = GetHeartRate(filename);
% ProcData.Fs.HR_fs = RawData.dal_fr;

%% NEURAL DATA -> PROCESS INTO VARIOUS FORMS:

if isfield(RawData,'neuro')
    
    % Wide-Band LFP
    [ProcData.Wideband_LFP, ProcData.Fs.Wideband_LFP_fs] = ...
        ProcessNeuro_atw0(RawData,'Wideband_LFP',animal,hem,numday);
    
    % MU Band [300-3000]
    [ProcData.MUpower, ProcData.Fs.MUpower_fs] = ProcessNeuro_atw0(RawData,...
        'MUpower',animal,hem,numday);
    
    % Gamma Band [40-100]
    [ProcData.Gam, ProcData.Fs.Gam_fs] = ProcessNeuro_atw0(RawData,'Gam',...
        animal,hem,numday);
    
    % Beta [10-30 Hz]
    [ProcData.Beta, ProcData.Fs.Beta_fs] = ProcessNeuro_atw0(RawData,'Beta',...
        animal,hem,numday);
    
    % Sub-Alpha Frequencies
    [ProcData.SubAlpha, ProcData.Fs.SubAlpha_fs] = ProcessNeuro_atw0(RawData,...
        'SubAlpha',animal,hem,numday);
    
end
%% WHISKER ANGLE -> BINARIZE, SET RESTING ANGLE TO ZERO DEGREES

% Track dropped basler frames
ExpectedFrames = RawData.TrialDur*RawData.bas_fr;
ProcData.DroppedBaslerFrames = ExpectedFrames-length(RawData.angle);
if ProcData.DroppedBaslerFrames>15
    display([filename ' Dropped ' num2str(ProcData.DroppedBaslerFrames) ' Frames.'])
end
% Trim any additional frames, for resample
angl = RawData.angle(1:min(ExpectedFrames,length(RawData.angle)));

% Create filter for whisking/movement
[z,p,k] = butter(whisk_filt_ord,whisk_filt_thresh/(RawData.bas_fr/2),'low');
[sos,g] = zp2sos(z,p,k);
filt_wwf = filtfilt(sos,g,angl-mean(angl));
res_wwf = resample(filt_wwf,wwf_fs,RawData.bas_fr);

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
[linked_bin_wwf] = link_binary_events(gt(bin_wwf,0),[round(wwf_fs/3), ...
    0]);

inds = linked_bin_wwf == 0;
rest_angle = mean(res_wwf(inds));

ProcData.wwf = res_wwf - rest_angle;
ProcData.Fs.wwf_fs = wwf_fs;
ProcData.Bin_wwf = bin_wwf;
ProcData.Fs.Bin_wwf_fs = wwf_fs;

%% FORCE SENSOR -> BINARIZE, DOWNSAMPLE TO DESIRED FREQUENCY

% Trim any additional data points for resample
trimmedpswf = RawData.pswf(1:min(ExpectedLength,length(RawData.pswf)));

% Filter then downsample the force sensory waveform to desired frequency 
[z,p,k] = butter(ps_filt_ord,ps_filt_thresh/(RawData.an_fs/2),'low');
[sos,g] = zp2sos(z,p,k);
filtpswf = filtfilt(sos,g,trimmedpswf);
ProcData.pswf = resample(filtpswf,pswf_fs,RawData.an_fs);
ProcData.Fs.pswf_fs = pswf_fs;

% Binarize the force sensor waveform
[ok] = Check4Threshold(['BinPSWF_' strday], animal, hem);
if ok == 0;
    [pswf_thresh] = CreatePswfThreshold(ProcData.pswf);
    Thresholds.(['BinPSWF_' strday]) = pswf_thresh;
    save([animal '_' hem '_Thresholds.mat'],'Thresholds');
end
ProcData.Bin_pswf = binarize_pswf(ProcData.pswf,Thresholds.(['BinPSWF_' strday]));
ProcData.Fs.Bin_pswf_fs = pswf_fs;

%% SAVE THE PROCESSED DATA STRUCTURE

save([filename(1:(underscore_indexes(end)-1)) '_ProcData.mat'],'ProcData');