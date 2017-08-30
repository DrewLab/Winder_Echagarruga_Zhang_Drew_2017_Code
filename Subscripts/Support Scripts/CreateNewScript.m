function [] = CreateNewScript(Fun_Name)
%   function [] = CreateNewScript(Fun_Name)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Creates and saves a new script and inserts a script 
%   template into it.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                     Fun_Name - [string] the name of the function as it 
%                           will be saved.
%_______________________________________________________________
%   RETURN:                     
%                     None, the result of this script is that a new MATLAB
%                           script will appear and be saved to the
%                           user-defined location.
%_______________________________________________________________
if(~isdeployed)
	prevdir = cd(fileparts(which('ScriptTemplate.m')));
end

targetdir = uigetdir(userpath,'Choose Save Location for Script: ');

if nargin == 0
    Fun_Name = input('Input Function Name: ','s');
end

copyfile('ScriptTemplate.m',[targetdir '\' Fun_Name])
edit(Fun_Name)
cd(prevdir)