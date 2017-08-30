function [] = CategorizeData(filename)
%   [] = CategorizeData(filename)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Identifies periods of sensory stimulation, volitional
%   movement, and rest. Calculates relevant details for each behavioral
%   period:
%           Stimulation:    Whisk Score - measure of the intensity of 
%                           whisking before and after onset of a puff. 
%                           A 0 indicates no whisking, a 1 indicates 
%                           maximum whisking over a 1 second period.
%                           Movement Score - Same as whisk score except
%                           uses the force sensor beneath the animal to
%                           detect body movment.
%
%           Whisking:       Duration - the time, in seconds, from onset to
%                           cessation of a whisking event.
%                           Rest Time - the duration of resting behavior
%                           prior to onset of the volitional whisk
%                           Whisk Score - a measure of the intensity of
%                           whisking for the duration of the whisk a
%                           maximum whisk for the whole duration will give
%                           a score of 1. No whisking will give a score of
%                           0.
%                           Movement Score - Same as whisk score exept uses
%                           the force sensor beneath the animal to detect
%                           body movement over the duration of the
%                           volitional whisk.
%                           Puff Distance - The time, in seconds, between
%                           the onset of each whisk an every puff
%                           administered during the trial.
%
%
%          Rest:            Duration - the time, in seconds without any 
%                           detected whisking or body movement.
%                           Start Time - the trial time corresponding to
%                           the cessation of all volitional movement.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                           filename - [string] file identifier                      
%_______________________________________________________________
%   RETURN:                     
%                           None, output of the script is additions to the
%                           ProcData structure.
%_______________________________________________________________

%% Load and Setup
display(['Loading filename: ' filename])
load(filename)
wwf_fs = ProcData.Fs.wwf_fs;
link_thresh = 0.5; % seconds
break_thresh = 0;% seconds

%% Further data processing

% Assume that whisks at the beginning/end of trial continue outside of the
% trial time. This will link any event occurring within "link_thresh"
% seconds to the beginning/end of the trial rather than assuming that it is
% a new/isolated event.
mod_Bin_wwf = ProcData.Bin_wwf;
mod_Bin_wwf([1,end]) = 1;
Bin_wwf = link_binary_events(gt(mod_Bin_wwf,0),...
    [link_thresh break_thresh]*wwf_fs);
SuperWhisk_wwf = ProcData.Bin_wwf==1;

%% Categorize data by behavior

% Retrieve details on whisking events
[ProcData.Flags.whisk] = GetWhiskingData(ProcData,Bin_wwf,SuperWhisk_wwf);

% Retrieve details on puffing events
[ProcData.Flags.stim] = GetStimData(ProcData);

% Identify and separate resting data
[ProcData.Flags.rest] = GetRestData(ProcData);

% Save ProcData structure
save(filename,'ProcData');

function [Puff_Times] = GetPuffTimes(ProcData)
%   function [Puff_Times] = GetPuffTimes(ProcData)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Gets the time in seconds of all puffs administered during
%   a trial.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                       ProcData - [struct] structure obtained using the 
%                       function ProcessRawDataFile.
%_______________________________________________________________
%   RETURN:                     
%                       Puff_Times - [array] time in seconds of all puffs              
%_______________________________________________________________
SolNames = fieldnames(ProcData.Sol);
Puff_List = cell(1,length(SolNames));
for SN = 1:length(SolNames)
    Puff_List{SN} = ProcData.Sol.(SolNames{SN});
end
Puff_Times = cell2mat(Puff_List);

function [Stim] = GetStimData(ProcData)
%   function [Stim] = GetStimData(ProcData)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Returns details on puffs administered during a trial.
%   Including: 
%                           Whisk Score - measure of the intensity of 
%                           whisking before and after onset of a puff. 
%                           A 0 indicates no whisking, a 1 indicates 
%                           maximum whisking over a 1 second period.
%                           Movement Score - Same as whisk score except
%                           uses the force sensor beneath the animal to
%                           detect body movment.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                       ProcData - [struct] structure obtained using the 
%                       function ProcessRawDataFile.    
%_______________________________________________________________
%   RETURN:                     
%                       Stim - [struct] structure containing a nested 
%                       structure for each puff administered. Each nested 
%                       structure contains details about puffs from a
%                       single solenoid.
%_______________________________________________________________

% Setup
wwf_fs = ProcData.Fs.wwf_fs;
pswf_fs = ProcData.Fs.pswf_fs;
Puff_Times = GetPuffTimes(ProcData);
Trial_Dur = ProcData.TrialDur;

% Set time intervals for calculation of the whisk scores
pre_time = 1;
post_time = 1;

% Get puffer IDs
SolNames = fieldnames(ProcData.Sol);
Stim.Name = cell(length(Puff_Times),1);
Stim.EventTime = zeros(length(Puff_Times),1);
Stim.WhiskScore_Pre = zeros(length(Puff_Times),1);
Stim.WhiskScore_Post = zeros(length(Puff_Times),1);
Stim.MoveScore_Pre = zeros(length(Puff_Times),1);
Stim.MoveScore_Post = zeros(length(Puff_Times),1);
i=1;
for SN = 1:length(SolNames)
    sol_puff_times = ProcData.Sol.(SolNames{SN});
    for p = 1:length(sol_puff_times) 
        if Trial_Dur-sol_puff_times(p)<=post_time;
            display(['Puff at time: ' num2str(sol_puff_times(p)) ' is too close to trial end'])
            continue;
        end
        % Set indexes for pre and post periods
        w_puff_ind = round(sol_puff_times(p)*wwf_fs);
        m_puff_ind = round(sol_puff_times(p)*pswf_fs);
        w_pre_start = max(round((sol_puff_times(p)-pre_time)*wwf_fs),1);
        m_pre_start = max(round((sol_puff_times(p)-pre_time)*pswf_fs),1);
        w_post_end = round((sol_puff_times(p)+post_time)*wwf_fs);
        m_post_end = round((sol_puff_times(p)+post_time)*pswf_fs);
        
        % Calculate the percent of the pre-stim time that the animal moved
        % or whisked
        whisk_score_pre = sum(ProcData.Bin_wwf(w_pre_start:w_puff_ind))...
            /(pre_time*wwf_fs);
        whisk_score_post = sum(ProcData.Bin_wwf(w_puff_ind:w_post_end))...
            /(post_time*wwf_fs);
        move_score_pre = sum(ProcData.Bin_pswf(m_pre_start:m_puff_ind))...
            /(pre_time*pswf_fs);
        move_score_post = sum(ProcData.Bin_pswf(m_puff_ind:m_post_end))...
            /(post_time*pswf_fs);
        
        % Add to Stim structure
        Stim.Name{i} = SolNames{SN};
        Stim.EventTime(i) = sol_puff_times(p)';
        Stim.WhiskScore_Pre(i) = whisk_score_pre';
        Stim.WhiskScore_Post(i) = whisk_score_post';
        Stim.MoveScore_Pre(i) = move_score_pre'; 
        Stim.MoveScore_Post(i) = move_score_post';
        i=i+1;
    end
end

% Calculate the time to the closest puff, omit comparison of puff to itself
% (see nonzeros)
Puff_mat = ones(length(Puff_Times),1)*Puff_Times;
TimeElapsed = abs(nonzeros(Puff_mat - Puff_mat'));

% If no other puff occurred during the trial, store 0 as a place holder.
if isempty(TimeElapsed)
    PuffTimeElapsed = 0;
else
% if not empty, Reshape the array to compensate for nonzeros command
    PuffTimeElapsed = reshape(TimeElapsed,numel(Puff_Times)-1,...
        numel(Puff_Times));
end

% Convert to cell and add to struct, if length of Puff_Times = 0, coerce to
% 1 to accommodate the NaN entry.
PuffTimeCell = mat2cell(PuffTimeElapsed',ones(max(length(Puff_Times),1),1));
Stim.PuffDistance = PuffTimeCell;

function [Whisk] = GetWhiskingData(ProcData,Bin_wwf,SuperWhisk_Thresh)
%   function [Whisk] = GetWhiskingData(ProcData,Bin_wwf)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Returns details on whisks which occurred during a trial.
%   Including:
%                           Duration - the time, in seconds, from onset to
%                           cessation of a whisking event.
%                           Rest Time - the duration of resting behavior
%                           prior to onset of the volitional whisk
%                           Whisk Score - a measure of the intensity of
%                           whisking for the duration of the whisk a
%                           maximum whisk for the whole duration will give
%                           a score of 1. No whisking will give a score of
%                           0.
%                           Movement Score - Same as whisk score exept uses
%                           the force sensor beneath the animal to detect
%                           body movement over the duration of the
%                           volitional whisk.
%                           Puff Distance - The time, in seconds, between
%                           the onset of each whisk an every puff
%                           administered during the trial.
%_______________________________________________________________
%   PARAMETERS:             
%                       ProcData - [struct] structure obtained using the 
%                       function ProcessRawDataFile.    
%_______________________________________________________________
%   RETURN:                     
%                       Whisk - [struct] structure containing a nested 
%                       structure for each whisk performed.
%_______________________________________________________________

%% Setup
wwf_fs = ProcData.Fs.wwf_fs;
pswf_fs = ProcData.Fs.pswf_fs;

%% Get Puff Times
[Puff_Times] = GetPuffTimes(ProcData);

%% Find the starts of whisking
Whisk_edge = diff(Bin_wwf);
whisk_samples = find(Whisk_edge>0);
whisk_starts = whisk_samples/wwf_fs;

%% Classify each whisking event by duration, whisking intensity, rest durations
samplevec = 1:length(Bin_wwf); 

% Identify periods of whisking/resting, include beginning and end of trial
% if needed (hence unique command) for correct interval calculation
high_samples = unique([0, samplevec(Bin_wwf), samplevec(end)+1]); 
low_samples = unique([0, samplevec(not(Bin_wwf)), samplevec(end)+1]);

% Calculate the number of samples between consecutive high/low samples.
d_high = diff(high_samples);
d_low = diff(low_samples);

% Identify skips in sample numbers which correspond to rests/whisks,
% convert from samples to seconds.
rest_len = d_high(d_high>1);
whisk_len = d_low(d_low>1);
rest_durs = rest_len/wwf_fs;
whisk_durs = whisk_len/wwf_fs;

% Control for the beginning/end of the trial to correctly map rests/whisks
% onto the whisk_starts.
if Bin_wwf(1)
    whisk_durs(1) = [];
    whisk_len(1) = [];
elseif and(not(Bin_wwf(1)),Bin_wwf(2))
    whisk_samples(1) = [];
    whisk_starts(1) = [];
    whisk_durs(1) = [];
    whisk_len(1) = [];
end
    
if not(Bin_wwf(end))
    rest_durs(end) = [];
end

if Bin_wwf(end)
    whisk_samples(end) = [];
    whisk_starts(end) = [];
    whisk_durs(end) = [];
    whisk_len(end) = [];
    rest_durs(end) = [];
end

% Tag the whisking periods with accelerations above threshold for tracking
ST_inds = 2:length(SuperWhisk_Thresh);
ST_Whisk_starts = ST_inds(diff(SuperWhisk_Thresh)==1);

% Remove all periods of detected whisking without a threshold crossing
if not(isempty(ST_Whisk_starts))
    WhiskLogic = zeros(length(ST_Whisk_starts),length(whisk_durs));
    for ST = 1:length(ST_Whisk_starts)
        WhiskLogic(ST,:) = ST_Whisk_starts(ST)>=whisk_samples & ...
            ST_Whisk_starts(ST)<=(whisk_samples+whisk_len);
    end
    Thresh_Filter = sum(WhiskLogic)==0;
    whisk_samples(Thresh_Filter) = [];
    whisk_starts(Thresh_Filter) = [];
    whisk_durs(Thresh_Filter) = [];
    whisk_len(Thresh_Filter) = [];
    rest_durs(Thresh_Filter) = [];
end
    
% Calculate the whisking intensity -> sum(ProcData.Bin_wwf)/sum(Bin_wwf)
% over the duration of the whisk. Calculate the movement intensity over the
% same interval.
whisk_int = zeros(size(whisk_starts));
movement_int = zeros(size(whisk_starts));
for ws = 1:length(whisk_samples)
    % Whisking intensity
    wh_inds = whisk_samples(ws):whisk_samples(ws)+whisk_len(ws);
    whisk_int(ws) = sum(ProcData.Bin_wwf(wh_inds))/numel(wh_inds);
    
    % Movement intensity
    mo_start = round(whisk_starts(ws)*pswf_fs);
    mo_dur = round(whisk_durs(ws)*pswf_fs);
    mo_inds = max(mo_start,1):min(mo_start+mo_dur,length(ProcData.Bin_pswf));
    movement_int(ws) = sum(ProcData.Bin_pswf(mo_inds))/numel(mo_inds);
end

% Code Testing: Plot to verify correct whisk onset times
% plot((1:length(SuperWhisk_Thresh))/30,SuperWhisk_Thresh,'c');
% hold on; plot((1:length(Bin_wwf))/30,Bin_wwf,'r');
% hold on; scatter(whisk_starts,ones(size(whisk_starts)),'g^');
% hold on; scatter(whisk_starts+whisk_durs,ones(size(whisk_starts)),'rv');
% hold off;
% ylim([-0.1 1.1])
% pause;

% Calculate the time to the closest puff
% If no puff occurred during the trial, store 0 as a place holder.
if isempty(Puff_Times)
    Puff_Times = 0;
end
Puff_mat = ones(length(whisk_samples),1)*Puff_Times;
whisk_mat = whisk_samples'*ones(1,length(Puff_Times))/wwf_fs;
PuffTimeElapsed = abs(whisk_mat - Puff_mat);

% Convert to cell
PuffTimeCell = mat2cell(PuffTimeElapsed,ones(length(whisk_starts),1));

%% Compile into final structure

Whisk.EventTime = whisk_starts';
Whisk.Duration = whisk_durs';
Whisk.RestTime = rest_durs';
Whisk.WhiskScore = whisk_int';
Whisk.MovementScore = movement_int';
Whisk.PuffDistance = PuffTimeCell;

function [Rest] = GetRestData(ProcData)
%   function [Rest] = GetRestData(ProcData)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Returns details on periods of rest during a trial.
%   Including:
%          Rest:            Duration - the time, in seconds without any 
%                           detected whisking or body movement.
%                           Start Time - the trial time corresponding to
%                           the cessation of all volitional movement.   
%_______________________________________________________________
%   PARAMETERS:             
%                       ProcData - [struct] structure obtained using the 
%                       function ProcessRawDataFile.    
%_______________________________________________________________
%   RETURN:                     
%                       Rest - [struct] structure containing a nested 
%                       structure for each period of rest.
%_______________________________________________________________

% Setup
link_thresh = 1; % seconds
break_thresh = 0.03;
wwf_fs = ProcData.Fs.wwf_fs;
pswf_fs = ProcData.Fs.pswf_fs;
[Puff_Times] = GetPuffTimes(ProcData);

% Assume that whisks at the beginning/end of trial continue outside of the
% trial time. This will link any event occurring within "link_thresh"
% seconds to the beginning/end of the trial rather than assuming that it is
% a new/isolated event.
mod_Bin_wwf = ProcData.Bin_wwf;
mod_Bin_wwf([1,end]) = 1;

mod_Bin_pswf = ProcData.Bin_pswf;
mod_Bin_pswf([1,end]) = 1;

% Recalculate linked binarized wwf without omitting any possible whisks
Bin_wwf = link_binary_events(gt(mod_Bin_wwf,0),...
    [link_thresh break_thresh]*wwf_fs);
Bin_pswf = link_binary_events(mod_Bin_pswf,...
    [link_thresh break_thresh]*pswf_fs);

% Downsample bin_wwf to match length of bin_pswf
samplevec = 1:length(Bin_wwf); 
whisk_high = samplevec(Bin_wwf)/wwf_fs;
ds_Bin_wwf = zeros(size(Bin_pswf));
% Map 'high' index into pswf sampling frequency, coerce between 1 and
% length(Bin_pswf) to avoid dimension errors. Take only unique values.
ds_inds = min(max(round(whisk_high*pswf_fs),1),length(Bin_pswf));
ds_Bin_wwf(unique(ds_inds)) = 1;

% Combine Bin_wwf and Bin_pswf, new sampling rate is pswf_fs
Bin_wf = logical(min(ds_Bin_wwf+Bin_pswf,1));
Fs = pswf_fs;

% figure; plot(ProcData.wwf);
% hold on; scatter(find(Bin_wf),45*ones(size(find(Bin_wf))),'k.')
% ylim([-50 50]);
% pause;

% Add puff times into the Bin_wf
puff_inds = round(Puff_Times*Fs);
Bin_wf(puff_inds) = 1;

% Find index for end of whisking event
edge = diff(Bin_wf);
samples = find([not(Bin_wf(1)) edge<0]);
stops = samples/Fs;

% Identify periods of whisking/resting, include beginning and end of trial
% if needed (hence unique command) for correct interval calculation
samplevec = 1:length(logical(Bin_wf));
high_samples = unique([1, samplevec(Bin_wf), samplevec(end)]); 
low_samples = unique([1, samplevec(not(Bin_wf)), samplevec(end)]); 

% Calculate the number of samples between consecutive high/low samples.
d_high = diff(high_samples);
d_low = diff(low_samples);

% Identify skips in sample numbers which correspond to rests/whisks,
% convert from samples to seconds.
rest_len = d_high(d_high>1);
rest_durs = rest_len/Fs;
whisk_len = d_low(d_low>1);
whisk_durs = whisk_len/Fs;

% Control for the beginning/end of the trial to correctly map rests/whisks
% onto the whisk_starts. Use index 2 and end-1 since it is assumed that the
% first and last indexes of a trial are the end/beginning of a volitional
% movement.
if not(Bin_wf(2)) 
    whisk_durs = [NaN whisk_durs];
end
if Bin_wf(end-1)
    whisk_durs(end) = [];
end

% Calculate the time to the closest puff
% If no puff occurred during the trial, store 0 as a place holder.
if isempty(Puff_Times)
    Puff_Times = 0;
end
Puff_mat = ones(length(samples),1)*Puff_Times;
rest_mat = samples'*ones(1,length(Puff_Times))/Fs;
PuffTimeElapsed = abs(rest_mat - Puff_mat);

% Convert to cell
PuffTimeCell = mat2cell(PuffTimeElapsed,ones(length(samples),1));


% Compile into a structure
Rest.EventTime = stops';
Rest.Duration = rest_durs';
Rest.PuffDistance = PuffTimeCell;
Rest.WhiskDurs = whisk_durs';