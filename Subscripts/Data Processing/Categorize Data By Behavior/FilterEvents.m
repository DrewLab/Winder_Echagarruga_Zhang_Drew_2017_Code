function [FiltArray] = FilterEvents(DataStruct,Criteria)
%   function [] =  ()
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Filters a data structure according to a set of
%   user-defined criteria.
%   
%_______________________________________________________________
%   PARAMETERS:             
%               DataStruct - [structure] contains the data to be filtered.
%
%               Criteria - [structure] contains fieldnames with 
%                   instructions on how to filter DataStruct:
%                       Required fields:
%                           Fieldname - [cells of strings] the fieldnames of
%                           DataStruct to be used for filtering
%
%                           Comparison - [cells of strings] instruction of how to
%                           filter the fieldnames. This input in restricted
%                           to the three commands: 'gt','lt','equal'.
%
%                           Value - [cell of doubles] the value that the data in 
%                           Criteria.Fieldnames should be compared to using
%                           the instruction in Criteria.Comparison.
%_______________________________________________________________
%   RETURN:                     
%               FiltArray - [logical array] an array for filtering the 
%                   data in "DataStruct" according the instructions in
%                   "Fieldnames".
%_______________________________________________________________
FName = Criteria.Fieldname;
Comp = Criteria.Comparison;
Val = Criteria.Value;

if length(FName)~=length(Comp)
    error(' ')
elseif length(FName)~=length(Val)
    error(' ')
end

FiltArray = true(size(DataStruct.Data,1),1);
for FN = 1:length(FName)
    if ~isfield(DataStruct,FName{FN})
        error('Criteria field not found')
    end
    switch Comp{FN}
        case 'gt'
            if iscell(DataStruct.(FName{FN}))
                if ischar(DataStruct.(FName{FN}){1})
                    error(' ')
                else
                    IndFilt = false(size(FiltArray));
                    for c = 1:length(DataStruct.(FName{FN}))
                        IndFilt(c) = all(gt(abs(DataStruct.(FName{FN}){c}), Val{FN}));
                    end
                end
            else
                IndFilt = gt(DataStruct.(FName{FN}), Val{FN});
            end
        case 'lt'
             if iscell(DataStruct.(FName{FN}))
                if ischar(DataStruct.(FName{FN}){1})
                    error(' ')
                else
                    IndFilt = false(size(FiltArray));
                    for c = 1:length(DataStruct.(FName{FN}))
                        IndFilt(c) = all(lt(DataStruct.(FName{FN}){c}, Val{FN}));
                    end
                end
            else
                IndFilt = lt(DataStruct.(FName{FN}), Val{FN});
            end
        case 'equal'
            if iscell(DataStruct.(FName{FN}))
                IndFilt = strcmp(DataStruct.(FName{FN}),Val{FN});
            else
                IndFilt = DataStruct.(FName{FN}) == Val{FN};
            end
        otherwise
            error(' ')
    end
    
    % This is a patch for old data, only consider the pre-event movement
    if or(strcmp(FName{FN},'WhiskScore'),strcmp(FName{FN},'MoveScore'))
        IndFilt = IndFilt(:,1);
    end
    
    FiltArray = and(FiltArray,IndFilt);
end