function [Frame] = GetSingleCBVFrame(filename, image_width, image_height)
%   function [Frame] = GetSingleCBVFrame(filename, image_width, image_height)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Reads in a single frame of a binary image
%   
%_______________________________________________________________
%   PARAMETERS:      
%               filename - [string] binary file name with extension
%
%               image_width - [double] number of pixels in width of image
%
%               image_height - [double] number of pixels in height of image
%                               
%_______________________________________________________________
%   RETURN:                     
%               Frame - [array] single image of the binary file           
%_______________________________________________________________

% Calculate the number of pixels in a single frame
pixels_per_frame=image_width*image_height;

% Open the Binary File
fid=fopen(filename);

% Read the image from the binary file
FramePix=fread(fid, pixels_per_frame,'*int16','b');

% Reshape the image into rows and columns
img=reshape(FramePix,image_height,image_width);

% Orient the frame so that rostral is up
Frame = rot90(img',2);