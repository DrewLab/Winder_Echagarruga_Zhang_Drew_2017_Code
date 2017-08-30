function [pp_id,sp_id,t_id,c_id] = IdentifySolenoids_atw0(hem)
%
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Oct 2013
%   Version 1
%
%   SUMMARY: Solenoids from whisker trials are identified by their
%   amplitude after being processed into "_RawData.mat". This script gives
%   a standard conversion from hemisphere recorded into amplitude of
%   primary, secondary, tail, and control solenoids

%   Solenoid Map: 
%           Direction                           Amplitude
%           ---------                           ---------
%        Left Whisker Pad                           1
%       Right Whisker Pad                           2
%             Tail                                  3
%            Control                                4
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
%_______________________________________________________________

if strcmp(hem,'RH') == 1;
    pp_id = 1;
    sp_id = 2;
elseif strcmp(hem,'LH') == 1;
    pp_id = 2;
    sp_id = 1;
end
t_id = 3;
c_id = 4;