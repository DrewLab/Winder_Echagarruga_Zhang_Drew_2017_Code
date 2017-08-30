function [tform,RA,RB] = RegisterCBVImage(Fixed,Moving)
%   function [tform,RA,RB] = RegisterCBVImage(Fixed,Moving)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Uses translation to align a moving image to a fixed image
%   and gets the world coordinates for the moving image.
%   
%_______________________________________________________________
%   PARAMETERS:             
%               Fixed - [array, integer] an image to be used as a reference
%               for the moving image
%
%               Moving - [array, integer] an image to be aligned to the
%               fixed image
%_______________________________________________________________
%   RETURN:                     
%               tform - [array, double] the translation matrix that aligns
%               the moving image to the fixed image.
%
%               RA - [object] imref2d object containing the world
%               coordinates for the fixed image.
%
%               RB - [object] imref2d object containing the world
%               coordinates for the moving image.
%
%               h - [handle] figure handle for the alignment image.
%_______________________________________________________________

% Register the Moving image
[optimizer,metric] = imregconfig('multimodal');
tform = imregtform(Moving,Fixed,'Translation',optimizer,metric);

% Get the world coordinates for the image after the transformation
[~,RB] = imwarp(Moving,tform);

% Plot the image regeistration in false color for verification
RA = imref2d(size(Fixed));
% h = figure; 
% imshowpair(Moving,RB,Fixed,RA,'Scaling','joint');