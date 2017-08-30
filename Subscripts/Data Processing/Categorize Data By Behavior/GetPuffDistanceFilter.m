function [PD_Filter] = GetPuffDistanceFilter(PuffDistance, distance_thresh)
%   function [PD_Filter] = GetPuffDistanceFilter(PuffDistance, distance_thresh)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Creates a logical array which designates events at
%   sufficient distance from sensory stimulation.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   PuffDistance - [nx1 Array of Cells] distances from
%                   surrounding puffs.  
%
%                   distance_thresh - [double] minimum unacceptable
%                   distance from surrounding puffs. Units should be the
%                   same as that of PuffDistance.
%_______________________________________________________________
%   RETURN:                     
%                 	PD_Filter - [Logical Array] designates events of
%                   sufficient distance from surrounding puffs.
%_______________________________________________________________

% Change to column vector, if needed
if size(PuffDistance,2)>size(PuffDistance,1)
    temp = PuffDistance';
    PuffDistance = temp;
    clear temp;
end

% Identify cells containing NaN (signifies no puffs in trial)
NaNTags = cellfun(@isnan,PuffDistance,'UniformOutput',0);
NaNFilter = cellfun(@sum,NaNTags);

% Identify cells with sufficient distance from puffs in trial
AbsDistance = cellfun(@abs,PuffDistance,'UniformOutput',0);
temparray = mat2cell(distance_thresh*ones(size(PuffDistance)),...
    ones(size(PuffDistance)));
DistanceTags = cellfun(@gt,AbsDistance,temparray,'UniformOutput',0);
% Use prod as logical "and" for all distances in a cell
DistanceFilter = cellfun(@prod,DistanceTags);

% Combine the NaN Filter and Distance Filter, coerce to 0 or 1, convert to
% a logical array for later indexing.
PD_Filter = logical(min(NaNFilter+DistanceFilter,1));