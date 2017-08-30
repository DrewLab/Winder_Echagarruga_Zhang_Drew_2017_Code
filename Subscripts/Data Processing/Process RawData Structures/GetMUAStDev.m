function [StDev] = GetMUAStDev(fileday)
%   [StDev] = GetMUAStDev(fileday)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Calculates the standard deviation of multi-unit activity 
%       for a single day.
%   
%_______________________________________________________________
%   PARAMETERS:     
%                   fileday - [string] six digit designation of
%                    the day (yymmdd)                       
%_______________________________________________________________
%   RETURN:                     
%                   StDev - [double] MUA standard deviation for 'fileday'       
%_______________________________________________________________

strday = ConvertDate(fileday);
display(['GetMUAStDev.m: Calculating MUA standard deviation for: ' strday '...'])
tic;

% Load animal thresholds if exist
Threshfile = ls('*Thresholds.mat');
if not(isempty(Threshfile))
    load(Threshfile)
else
    Thresholds = [];
end

% Check for MUA standard deviation for current day
if isfield(Thresholds,['MUA_StDev_' strday]);
    StDev = Thresholds.(['MUA_StDev_' (strday)]);
else
    % Get all filenames that match current day
    all_files = ls('*_rawdata.mat');
    underscore_inds = strfind(all_files(1,:),'_');
    ind = underscore_inds(2)+1:underscore_inds(3)-1;
    [filenames] = GetFilenames('*_rawdata.mat',fileday, ind);
    
    % Calculate the variance for the current day
    MUA_var = zeros(length(filenames),1);
    for f = 1:length(filenames)
        name = filenames{f};
        load(name)
        [b,a] = butter(2,300/(RawData.Fs.Analog/2),'high');
        MU_data = filtfilt(b,a,RawData.Data.Neuro);
        MUA_var(f) = var(MU_data);
    end
    
    %% Add Standard Deviation to Thresholds structure
    StDev = sqrt(mean(MUA_var));
end
elapsed = toc;
display(['GetMUAstDev.m: Done... Time Elapsed: ' num2str(elapsed) ' seconds...'])