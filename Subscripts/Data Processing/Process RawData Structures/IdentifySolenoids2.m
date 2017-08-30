function [pp_id,sp_id,t_id,c_id] = IdentifySolenoids2(hem)

%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, April
%   2015
%   Version 2

%   SUMMARY: Solenoids from whisker trials are identified by their
%   amplitude after being processed into "_RawData.mat". This script gives
%   a standard conversion from hemisphere recorded into amplitude of
%   primary, secondary, tail, and control solenoids

%   Solenoid Map: 
%           Direction                           Amplitude
%           ---------                           ---------
%        Left Whisker Pad                           2
%       Right Whisker Pad                           3
%             Tail                                  4
%            Control                                5
%_______________________________________________________________
%   INPUTS:
%               hem - the recorded brain hemisphere during whisker trial
%_______________________________________________________________
%   OUTPUTS:
%               pp_id - amplitude of the primary solenoid in "_RawData.m"
%               sp_id - amplitude of the secondary solenoid in "_RawData.m"
%               t_id - amplitude of the tail solenoid in "_RawData.m"
%               c_id - amplitude of the control solenoid in "_RawData.m"
%_______________________________________________________________
%   REQUIRED SCRIPTS: None
%_______________________________________________________________
%   CALLED BY: 
%               ChunkData.m
%_______________________________________________________________
%   FUTURE VERSIONS:
%_______________________________________________________________
%   CHANGES FROM PREV VERS:
%       Compatible with LabVIEW vi: WhiskerTrials_Stream2Disk
%_______________________________________________________________

if strcmp(hem,'RH') == 1;
    pp_id = 2;
    sp_id = 3;
elseif strcmp(hem,'LH') == 1;
    pp_id = 3;
    sp_id = 2;
end
t_id = 4;
c_id = 5;