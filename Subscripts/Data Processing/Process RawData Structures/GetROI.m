function [mask] = GetROI(img, ROIname, animal, hem)
%   [mold] = GetROI(img,ROIname,animal,hem)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Retrieves coordinates for a region of interest from shared
%   variables.
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

isok = 'n';
ROIfile = ls('*ROIs.mat');
if not(isempty(ROIfile))
    load(ROIfile)
else
    ROIs = [];
end

if isfield(ROIs,ROIname);
    xi = ROIs.(ROIname).xi;
    yi = ROIs.(ROIname).yi;
    display('Previous ROI found for this animal and hemisphere.')
    figure(99); imagesc(img); colormap(gray); axis image; 
    xlabel('Caudal');
    if strcmp(hem,'LH')
        ylabel('Lateral');
    elseif strcmp(hem,'RH');
        ylabel('Medial')
    end
    impoly(gca,[xi,yi]);
    if strcmp(isok,'y');
        mask = roipoly(img,xi,yi);
    end
    while strcmp(isok,'n') == 1;
        display('Previous ROI found... abort function and delete ROI file to create new.');
        isok = 'y';
        if strcmp(isok,'y')
            mask = roipoly(img,xi,yi);
            continue;
        elseif strcmp(isok,'n')
            [mask] = CreateROI(img,animal,hem);
        end
    end
else
    [mask] = CreateROI(img,ROIname,animal,hem);
end
end
