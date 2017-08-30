function [sR2] = SlidingRSquared(data1,data2,winsize)
%   function [sR2] = SlidingRSquared(data1,data2,winsize)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Calculates the coefficient of determination within a
%   sliding window of size=winsize
%_______________________________________________________________
%   PARAMETERS:             
%                   data1 - [array]
%
%                   data2 - [array]
%
%                   winsize - [int] the length of the window in samples
%_______________________________________________________________
%   RETURN:                     
%                   sR2 - [array] coefficients of determination for each
%                   window as an array
%_______________________________________________________________
R2 = zeros(1,length(data1)-winsize);
for si = 1:(length(data1)-winsize)
    snip1 = data1(si:si+winsize);
    snip2 = data2(si:si+winsize);
    R2(si) = CalculateRsquared(snip1,snip2);
end
sR2 = [NaN*ones(1,winsize/2), R2, NaN*ones(1,winsize/2)];