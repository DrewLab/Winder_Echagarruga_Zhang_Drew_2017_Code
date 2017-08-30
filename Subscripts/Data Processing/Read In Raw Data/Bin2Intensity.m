function [Refl] = Bin2Intensity(filename,ROImask,Frames)
%   function [Refl] = Bin2Intensity(filename,ROIx,ROIy)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Converts .bin file into a timeseries of average pixel
%   intensity values of a user defined ROI.
%   
%_______________________________________________________________
%   PARAMETERS:    
%                   filename - [string] name of the '*dalsa.bin' camera 
%                   containing the reflectance data from cranial
%                   window.
%
%                   ROImask - [matrix] same size as the camera images where
%                   '1' occupy the region within the roi and '0' values
%                   occur everywhere else. Create 'ROImask' using roipoly
%                   function.
%                               
%_______________________________________________________________
%   RETURN:                     
%                   Refl - [array] a timeseries of the reflectance data            
%_______________________________________________________________

% Import camera frames from dalsa file
if nargin<3
    [Frames]=ReadDalsaBinary(filename,256,256);
end

nFrames = length(Frames);


Refl = zeros(1,nFrames);
for n = 1:nFrames
    mask = ROImask.*double(Frames{n});
    Refl(n) = mean(nonzeros(mask));
end