function [ID,hem,filedate,fileID] = GetFileInfo(filename)
%   function [ID,hem,filedate,fileID] = GetFileInfo(filename)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Uses file name in a standard format to get the animal ID
%   and hemisphere recorded.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                  filename - [string] filename in a standard format.
%                       'SubjectID_Hemisphere_Date_HH_mm_ssdd'
%_______________________________________________________________
%   RETURN:                     
%                   Animal_ID - [string] the subject identifier
%
%                   hem - [string] the hemisphere recorded
%
%                   filedate - [string] the date the file was recorded
%
%                   fileID - [string] the timestamp of the file
%_______________________________________________________________

% Handle cell arrays so that the output of uigetfile can be used.
if iscell(filename)
    arraylen = length(filename);
    celllen = length(filename{1});
    CellContents = cell2mat(filename);
    filename = reshape(CellContents,celllen,arraylen)';
end

% Identify the extension
Ext_ind = strfind(filename(1,:),'.');
extension = filename(1,Ext_ind+1:end);

% Identify the underscores
filebreaks = strfind(filename(1,:),'_');

switch extension
    case 'bin'
        ID = [];
        hem = [];
        filedate = filename(:,1:filebreaks(1)-1);
        fileID = filename(:,1:filebreaks(4)-1);
    case 'tdms'
        ID = [];
        hem = [];
        filedate = filename(:,1:filebreaks(1)-1);
        fileID = filename(:,1:Ext_ind-1);
    case 'mat'       
        % Use the known format to parse
        ID = filename(:,1:filebreaks(1)-1);
        hem = filename(:,filebreaks(1)+1:filebreaks(2)-1);
        if numel(filebreaks)>3
            filedate = filename(:,filebreaks(2)+1:filebreaks(3)-1);
            fileID = filename(:,filebreaks(2)+1:filebreaks(6)-1);
        else
            filedate = [];
            fileID = [];
        end
end