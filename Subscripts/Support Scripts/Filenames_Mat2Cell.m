function [Filenames_Cell] = Filenames_Mat2Cell(Filenames_Matrix)
%   [Filenames_Cell] = Filenames_Mat2Cell(Filenames_Matrix)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Converts a matrix list of filenames into a list within a
%   cell array. This improves compatibility between scripts that use 'ls'
%   and 'uigetfile' 
%   
%_______________________________________________________________
%   PARAMETERS:             
%               Filenames_Matrix - [matrix] individual filenames contained
%               in each row.
%_______________________________________________________________
%   RETURN:                     
%               Filenames_Cell - [cell array] individual filenames
%               contained in each cell
%_______________________________________________________________

Filenames_Cell = mat2cell(Filenames_Matrix,...
    ones(1,size(Filenames_Matrix,1)),...
    size(Filenames_Matrix,2));