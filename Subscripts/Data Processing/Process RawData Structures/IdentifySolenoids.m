function [contra_id,ipsi_id,tail_id,auditory_id] = IdentifySolenoids(RDNotes,hem)
%%  [contra_id,ipsi_id,tail_id,auditory_id] = IdentifySolenoids(RDNotes,hem)
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   SUMMARY: Solenoids from whisker trials are identified by their
%   amplitude after being processed into "_RawData.mat". This script gives
%   a standard conversion from hemisphere recorded into amplitude of
%   primary, secondary, tail, and auditory solenoids
%_______________________________________________________________
%   PARAMETERS:
%                           filename: [str] name of _rawdata.mat file to be
%                           processed
%
%   RETURN:
%                           Function output is a saved structure called
%                           'ProcData.mat'
%_______________________________________________________________
%

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
%        Left Whisker Pad                           3
%       Right Whisker Pad                           4
%             Tail                                  5
%            Control                                6
%_______________________________________________________________
%   INPUTS:
%               RDNotes - [structure] contains the notes from a imaging
%               session. Corresponds to the notes subfield of a RawData
%               structure.
%
%               hem - [string] the recorded brain hemisphere during whisker 
%               trial
%_______________________________________________________________
%   OUTPUTS:
%               contra_id - amplitude of the solenoid in contralateral to the
%               imaged hemisphere
%               ispi_id - amplitude of the secondary solenoid in "_RawData.m"
%               tail_id - amplitude of the tail solenoid in "_RawData.m"
%               auditory_id - amplitude of the control solenoid in "_RawData.m"
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
    contraSolmatch = strcmp(RDNotes.Notes.SolenoidName,'Left Pad');
    ipsiSolmatch = strcmp(RDNotes.Notes.SolenoidName,'Right Pad');
    contra_id = RDNotes.Notes.SolenoidWaveformAmplitude(contraSolmatch);
    ipsi_id = RDNotes.Notes.SolenoidWaveformAmplitude(ipsiSolmatch);
    tailSolmatch = strcmp(RDNotes.SolenoidName,'Tail');
    tail_id = RDNotes.SolenoidWaveformAmplitude(tailSolmatch);
    AuditorySolmatch = strcmp(RDNotes.SolenoidName,'Auditory');
    auditory_id = RDNotes.SolenoidWaveformAmplitude(AuditorySolmatch);
elseif strcmp(hem,'LH') == 1;
    contraSolmatch = strcmp(RDNotes.SolenoidName,'Right Pad');
    ipsiSolmatch = strcmp(RDNotes.SolenoidName,'Left Pad');
    contra_id = RDNotes.SolenoidWaveformAmplitude(contraSolmatch);
    ipsi_id = RDNotes.SolenoidWaveformAmplitude(ipsiSolmatch);
    tailSolmatch = strcmp(RDNotes.SolenoidName,'Tail');
    tail_id = RDNotes.SolenoidWaveformAmplitude(tailSolmatch);
    AuditorySolmatch = strcmp(RDNotes.SolenoidName,'Auditory');
    auditory_id = RDNotes.SolenoidWaveformAmplitude(AuditorySolmatch);
else
    contra_id = [];
    ipsi_id = [];
    tail_id = [];
    auditory_id = [];
end
