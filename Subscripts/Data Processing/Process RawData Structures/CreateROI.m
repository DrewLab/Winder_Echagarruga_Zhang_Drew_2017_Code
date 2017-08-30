function [mask] = CreateROI(img,ROIname,animal,hem)
%   [mold] = CreateROI(img,ROIname,animal,hem)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Prompts user to identify a region of interest on an image
%   and saves the selection to the shared variables.
%_______________________________________________________________
%   PARAMETERS:      
%                   img - [matrix] the as output by ReadDalsaBinary.m
%
%                   ROIname - [string] a designation for the region of
%                   interest, usually should have some description and a
%                   date (i.g. 'Barrels_May20')
%
%                   animal - [string] ID for the animal
%
%                   hem - [string] hemisphere recorded
%_______________________________________________________________
%   RETURN:       
%                   mask - [matrix] same size as image where pixels within
%                   the ROI are designated with a '1'. All other elements
%                   are '0'.
%_______________________________________________________________

display('Please select your region of interest.')
figure(99); imagesc(img); colormap(gray); axis image;
xlabel('Caudal');
if strcmp(hem,'LH')
    ylabel('Lateral');
elseif strcmp(hem,'RH')
    ylabel('Medial')
end
[mask, xi, yi] = roipoly;
ROIs.(ROIname).xi = xi;
ROIs.(ROIname).yi = yi;
asksave = input('Save this ROI for future? (y/n) ','s');
if asksave == 'y';
    save([animal '_' hem '_ROIs.mat'], 'ROIs');
end
impoly(gca,[xi,yi]);