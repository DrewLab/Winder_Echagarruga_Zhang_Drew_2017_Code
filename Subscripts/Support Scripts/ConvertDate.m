function [days] = ConvertDate(DateTag)
%   [days] = ConvertDate(DateTag)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Converts the date format output by the LabVIEW acquisition
%   program into a date string.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   DateTag - [string] date vector output by the LabVIEW
%                       acquisition
%_______________________________________________________________
%   RETURN:                     
%                   days - [String] date as a string MMMDD                    
%_______________________________________________________________
days = cell(size(DateTag,1),1);
for f = 1:size(DateTag,1)
    days{f} = datestr([2000+str2double(DateTag(f,1:2)) str2double(DateTag(f,3:4)) ...
        str2double(DateTag(f,5:6)) 00 00 00],'mmmdd');
end

if length(days) == 1
    days = days{1};
end