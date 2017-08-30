function [tvec] = ConvertTime(filename)
%   [tvec] = ConvertTime(filename)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Converts the date format output from the LabVIEW
%   acquisition program to a time vector
%   
%_______________________________________________________________
%   PARAMETERS:            
%                   filename - [            
%_______________________________________________________________
%   RETURN:                     
%                               
%_______________________________________________________________

if iscell(filename)
    % Change to cell columns
    if size(filename,1)<size(filename,2)
        filename = filename';
    end
    filename = cell2mat(filename);
end
tvec = zeros(size(filename,1),6);
for f = 1:size(filename,1)
    tvec(f,:) = [2000+str2double(filename(f,1:2)) str2double(filename(f,3:4)) ...
        str2double(filename(f,5:6)) str2double(filename(f,8:9))...
        str2double(filename(f,11:12)) str2double(filename(f,14:15))];
end