function  [bin_wwf] = binarize_wwf(angl,fs,thresh1,thresh2)
%% function [bin_wwf] = binarize_wwf_v1(angl,fs,thresh1,thresh2)

% Written by Aaron Winder, Drew Lab, ESM, Penn State University, Sept 2013
% Version 1

% SUMMARY: Converts a timeseries of whisker angles into a binarized
% waveform. Binarization is based on the second derivative of the whisker
% angle so that forces on the whiskers may be considered when defining a
% threshold.
%________________________________________________________________________
% INPUTS:          
%                   angl - a whisker angle timeseries given as a row vector

%                   fs - the sampling rate of the whisker angle timeseries

%                   thresh - a user-defined threshold. 
%________________________________________________________________________
% OUTPUTS:         
%                   bin_wwf - the binarized waveform as a row vector
%________________________________________________________________________

% Differentiate, rectify, subtract off noise, rectify
dd_wwf = abs((diff(angl,2)))*fs^2;
bin_wwf1 = gt(dd_wwf,thresh1);                                              % Acceleration exceeds lower threshold
bin_wwf2 = gt(dd_wwf,thresh2);                                              % Acceleration exceeds upper threshold

% Combine the two waveforms
bin_wwf = (bin_wwf1+bin_wwf2)/2;
