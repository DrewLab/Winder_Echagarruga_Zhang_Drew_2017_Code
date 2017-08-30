function [UniqueDays,Day_ID] = GetUniqueDays(DateList)
%   [UniqueDays,Day_ID] = GetUniqueDays(DateList)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Identifies the unique entries of a list containing
%   timestamp identifiers.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   DateList - [matrix] list of dates or filenames            
%_______________________________________________________________
%   RETURN:                     
%                   UniqueDays - [cell array] list of session dates.                    
%_______________________________________________________________

if iscellstr(DateList)
    temp = cell2mat(DateList);
    DateList = temp;
end

filebreaks = strfind(DateList(1,:),'_');
if isempty(filebreaks)
    AllDates = DateList;
elseif or(length(filebreaks)==3,length(filebreaks)==4)
    AllDates = DateList(:,1:filebreaks(1)-1);
elseif length(filebreaks)==6
    date_ind = filebreaks(2)+1:filebreaks(3)-1;
    AllDates = DateList(:,date_ind);
else
    error('Format of the list of dates not recognized...')
end
All_days = mat2cell(AllDates,ones(1,size(AllDates,1)));
[UniqueDays,Day_ID] = unique(All_days);