function [routp,convout] = ConvolveHRF(HRF, inp, outp, offset)
%   function [routp,convout] = ConvolveHRF(HRF, inp, outp, offset)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Convolves neural activity with the HRF to predict CBV
%   data.
%
%_______________________________________________________________
%   PARAMETERS:
%                   HRF - [Array] the hemodynamic response function
%
%                   inp - [Array] the neural data
%
%                   outp - [Array] the CBV data
%
%                   offset - [int] indicates whether the kernel is shifted
%                   from time = 0.
%_______________________________________________________________
%   RETURN:
%                   routp - [array] the resampled outp vector
%
%                   convout - [array] the result of the convolution of
%                   neural activity and the HRF
%_______________________________________________________________

% Match Data Lengths
outdc = outp(1);
indc = inp(1);
[rinp,routp] = MatchDataLengths(inp-indc, outp-outdc);
%routp = routp+outdc;
rinp = rinp+indc;

% Perform Convolution
convl = conv(rinp,HRF);

% Add DC offset to result
convout = convl(1-offset:length(routp)-offset);


