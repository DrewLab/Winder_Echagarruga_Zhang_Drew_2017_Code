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

[z,p,k] = butter(4,2/(fs/2),'high');
[sos,g] = zp2sos(z,p,k);
Filtwwf = filtfilt(sos,g,angl);
Pos = abs(Filtwwf);

bin_wwf1 = gt(Pos,thresh1);
bin_wwf2 = gt(Pos,thresh2);

% Combine the two waveforms
bin_wwf = (bin_wwf1+bin_wwf2)/2;
