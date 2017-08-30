function [ok] = Check4Baseline(sfield, day, animal, hem)

%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Nov 2013
%   Version 1

%   SUMMARY: This function checks shared variable for baselines that can 
%               be used to normalize data.
%_______________________________________________________________
%   INPUTS:
%                           sfield - a string giving the subfield name
%                           (second level) of the shared variables
%                           structure
%                                       Hemo - CBV baseline
%                                       High_Gam - 70-110 Hz
%                                       Low_Gam - 30-50 Hz
%                                       SR - Spike Rate

%                           day - date of data as a string [format: mmmdd]

%                           animal - animal name as a string

%                           hem - recorded hemisphere as a string
%                                   [RH/LH]
%_______________________________________________________________
%   OUTPUTS:
%                           ok - output of 1 or 0 to indicate whether
%                           hemo baseline exist for the given day

%_______________________________________________________________
%   REQUIRED SCRIPTS:
%                           [Q] = DetectMachine(prevpath)
%_______________________________________________________________
%   CALLED BY:
%_______________________________________________________________
%   FUTURE VERSIONS: Make the code generalizable to other shared variables
%_______________________________________________________________
%   CHANGES FROM PREV VERS:
%_______________________________________________________________


% Begin Check
ok = 0;
if exist([animal '_' hem '_.mat'],'file') == 2;
    load([animal '_' hem '_SharVars.mat']);
    if isfield(SharVars,'Baselines');
        if isfield(SharVars.Baselines, sfield);
            if isfield(SharVars.Baselines.(sfield),day);
                ok = 1;
            end
        else
            SharVars.Baselines.(sfield) = [];
        end
    else
        SharVars.Baselines = [];
    end
end
cd(prevdir);
