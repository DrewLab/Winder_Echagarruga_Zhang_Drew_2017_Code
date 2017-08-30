function [res_data1, res_data2] = MatchDataLengths(data1,data2)
%   [res_data1, res_data2] = MatchDataLengths(data1,data2)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Given 2 arrays of the same duration but different sampling
%   rates, resample the larger of the arrays to match the sampling rate of
%   the lesser
%
%_______________________________________________________________
%   PARAMETERS:
%                   data1/data2 - [array] Any data array as a row vector
%_______________________________________________________________
%   RETURN:
%                   res_data1 - [array] If data1 > data2, this will be the 
%                   resampled data1 vector. If data1 <= data2, this will be
%                   the same as data1.
%
%                   res_data2 - [array] If data2 > data1, this will be the 
%                   resampled data1 vector. If data2 <= data1, this will be
%                   the same as data1.
%_______________________________________________________________

% Convert to column vectors
data1 = data1';
data2 = data2';
if size(data1,1) ~= size(data2,1)
    [~,larger] = max([size(data1,1) size(data2,1)]);
    [~,lesser] = min([size(data1,1) size(data2,1)]);
    P = eval(sprintf('size(data%d,1)',lesser));
    Q = eval(sprintf('size(data%d,1)',larger));
    eval(sprintf('res_data%d  = resample(data%d,P,Q);',larger,larger));
    eval(sprintf('res_data%d = data%d;',lesser,lesser));
    res_data1 = res_data1';
    res_data2 = res_data2';
else
    res_data1 = data1';
    res_data2 = data2';
end