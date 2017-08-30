function [bin_pswf] = binarize_pswf(pswf,thresh)
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Sep 2013
%   Version 1
%
%   SUMMARY: Binarizes the pressure sensor based of the force applied
%________________________________________________________________________
%   INPUTS:                     pswf - the array of voltage values from the
%                                       pressure sensor.
%________________________________________________________________________
%   OUTPUTS:                    bin_pswf - the binarized pressure sensor as
%                                       an array
%________________________________________________________________________
%   REQUIRED SCRIPTS:           None
%________________________________________________________________________
%   CALLED BY:                  ChunkData.m
%________________________________________________________________________
%   FUTURE VERSIONS:
%________________________________________________________________________
%   CHANGES FROM PREV VERS:     v2 - Now uses the force applied rather than
%                                   its derivative.


% Low Pass filter
% [b,a] = butter(2,20/fs/2,'low');
% filt_pswf = filtfilt(b,a,pswf-mean(pswf));
% Half-wave rectify, threshold above 0.25 (approximate resting is 0.1)
% identify values greater than 0.
% bin_pswf = gt(max(pswf,thresh)-thresh,0);

y = hilbert(diff(pswf));
env = abs(y);
bin_pswf = gt(env,thresh);