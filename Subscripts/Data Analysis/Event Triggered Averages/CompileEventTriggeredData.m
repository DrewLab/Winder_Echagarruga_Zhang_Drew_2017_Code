function [BehData,ZData,timevec] = CompileEventTriggeredData(animals,dataType,Behavior)
%   function [BehData,ZData,timevec] = CompileEventTriggeredData(animals,dataType,Behavior)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: This function gathers and combines data surrounding a
%   behavioral event.
%   
%_______________________________________________________________
%   PARAMETERS:             
%               animals - [cell array] animal IDs
%
%               dataType - [string] fieldname of the data structures to be
%               compiled
%
%               Behavior - [string] behavioral category
%_______________________________________________________________
%   RETURN:                     
%               BehData - [matrix] the baseline normalized event triggered
%               data of size (length(animals) x length(data))
%
%               Zdata - [matrix] the zscore normalized event triggered data
%               of size (length(animals) x length(data))
%
%               timevec - [array] vector containing the peri-stimulus time
%_______________________________________________________________

for a = 1:length(animals)
    animal = animals{a};
%     EventFile = ls([animal filesep '*_EVENTDATA_' dataType '.mat']);
%     load([animal filesep EventFile])
    EventFile = dir([animal filesep '*_EVENTDATA_' dataType '.mat']);
    load([animal filesep EventFile.name])
    [DataStruct,FiltArray] = SelectBehavioralEvents(EventData.(dataType),Behavior);
    if a == 1
        BehData = NaN*ones(length(animals),size(DataStruct.NormData,2));
        ZData = NaN*ones(length(animals),size(DataStruct.NormData,2));
    end
    BehData(a,:) = mean(DataStruct.NormData(FiltArray,:)-1);
    ZData(a,:) = mean(DataStruct.ZNormData(FiltArray,:));
    timevec = ((1:size(BehData,2))/DataStruct.Fs)-DataStruct.epoch.offset;
end