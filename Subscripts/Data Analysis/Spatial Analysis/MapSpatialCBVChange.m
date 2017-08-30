function [ChangeMap,World_Fixed,World_Moving] = MapSpatialCBVChange(Fixed,RA,Moving,RB)
%   function [] =  ()
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Aligns images in world coordinates and returns the
%   intensity difference between them.
%   
%_______________________________________________________________
%   PARAMETERS:    
%                   Fixed - [Integer Array] Fixed image
%
%                   RA - [Object] imref2d object containing the world
%                   coordinates for the fixed image.
%
%                   Moving - [Integer Array] Moving image
%
%                   RB - [Object] imref2d object containing the world
%                   coordinates for the fixed image.                       
%_______________________________________________________________
%   RETURN:   
%                   ChangeMap - [Integer Array] difference image between
%                   the moving and fixed images.
%                               
%_______________________________________________________________

% [Fused,RC] = imfuse(Moving,RB,Fixed,RA,'Diff');
% World_Moving = NaN*ones(size(Fused));
% World_Fixed = NaN*ones(size(Fused));
% Moving_Inds = zeros(2);
% [Moving_Inds(1,1),Moving_Inds(1,2)] = worldToSubscript(RC,RB.XWorldLimits(1),...
%     RB.YWorldLimits(1));
% [Moving_Inds(2,1),Moving_Inds(2,2)] = worldToSubscript(RC,RB.XWorldLimits(2),...
%     RB.YWorldLimits(2));
% Fixed_Inds = zeros(2);
% [Fixed_Inds(1,1),Fixed_Inds(1,2)] = worldToSubscript(RC,RA.XWorldLimits(1),...
%     RA.YWorldLimits(1));
% [Fixed_Inds(2,1),Fixed_Inds(2,2)] = worldToSubscript(RC,RA.XWorldLimits(2),...
%     RA.YWorldLimits(2));
% Moving_Xinds = Moving_Inds(1,1):Moving_Inds(2,1)-1;
% Moving_Yinds = Moving_Inds(1,2):Moving_Inds(2,2)-1;
% World_Moving(Moving_Xinds, Moving_Yinds) = Moving;
% Fixed_Xinds = Fixed_Inds(1,1):Fixed_Inds(2,1)-1;
% Fixed_Yinds = Fixed_Inds(1,2):Fixed_Inds(2,2)-1;
% World_Fixed(Fixed_Xinds, Fixed_Yinds) = Fixed;

[World_Fixed,World_Moving] = AlignRegisteredImages(Fixed,RA,Moving,RB);

ChangeMap = World_Moving-World_Fixed;