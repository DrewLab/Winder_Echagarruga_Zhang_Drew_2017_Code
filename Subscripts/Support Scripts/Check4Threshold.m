function [ok] = Check4Threshold(sfield, animal, hem)
%
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Nov 2013
%   Version 1
%
%   SUMMARY: This function checks shared variable for baselines that can 
%               be used to normalize data.
%_______________________________________________________________
%   INPUTS:
%                           sfield - a string giving the subfield name
%                           (second level) of the shared variables
%                           structure
%                                       BinWWF - for binarizing whisker
%                                       angle measurements
%                                       BinPSWF - for binarizing pressure
%                                       sensor measurements
%                                       Spike - for detecting a spike
%
%                           animal - animal name as a string
%
%                           hem - recorded hemisphere as a string
%                                   [RH/LH]
%_______________________________________________________________
%   OUTPUTS:
%                           ok - output of 1 or 0 to indicate whether
%                           hemo baseline exist for the given day
%
%_______________________________________________________________
%   REQUIRED SCRIPTS:
%                           [Q] = DetectMachine(prevpath)
%_______________________________________________________________
%   CALLED BY:
%                           HemoBaseline.m
%_______________________________________________________________
%   FUTURE VERSIONS: Make the code generalizable to other shared variables
%_______________________________________________________________
%   CHANGES FROM PREV VERS:
%_______________________________________________________________

% Navigate to Shared Variables folder

% Begin Check
ok = 0;
if exist([animal '_' hem '_Thresholds.mat'],'file') == 2;
    load([animal '_' hem '_Thresholds.mat']);
    if isfield(Thresholds, sfield)
        ok = 1;
    end
end