function [fname] = IdentifyStructureSubfield(Structure,field)
%   function [fname] = IdentifyStructureSubfield(Structure,field)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Finds a the subfield of a structure that contains the
%   requested field.
%   
%_______________________________________________________________
%   PARAMETERS:   
%                   Structure - [struct] the structure to search
%
%                   field - [string] the name of the desired field
%                               
%_______________________________________________________________
%   RETURN:   
%                   fname = [string] the subfield of "Structure" that
%                   contains "field"
%                               
%_______________________________________________________________

fnames = fieldnames(Structure);
fname = [];
for fn = 1:length(fnames)
    if isfield(Structure.(fnames{fn}),field)
        fname = fnames{fn};
    end
end