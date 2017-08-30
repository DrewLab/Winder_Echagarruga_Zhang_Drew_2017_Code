function [angle] = whiskertracker_parallel(filename)
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University
%
%   SUMMARY: Tracks the approximate movement of whiskers using the radon
%           transform of an image of the whiskers. Identifies the whisker
%           angles by analyzing the column variance of the transformed
%           image. Columns with the highest variance will correspond to the
%           correct angle where the highest value of the radon transform
%           will correspond to the correct angle and y-position in the
%           image, while the lowest value will correspond to the correct
%           angle but an incorrect image position.
%
%           This version of the code uses an onboard GPU to speed up the
%           calculation of the whisker angles.
%_______________________________________________________________
%   PARAMETER TYPE:             
%                               filename - [string] list of filames
%                               without the extension.
%_______________________________________________________________
%   RETURN:                     
%                               TDMSFile - [struct] contains measured
%                               analog data and trial notes from the
%                               LabVIEW acquisition program
%_______________________________________________________________

% Variable Setup
theta = -40:80; % Angles used for radon

% Import whisker movie
import_start = tic;
basler_frames = ReadBinFile_U8MatrixGradient([filename '_fwire.bin'],350,30);
import_time = toc(import_start);
display(['whiskertracker_parallel: Binary file import time was ' ...
    num2str(import_time) ' seconds.']);

% Transfer the images to the GPU
gpu_trans1 = tic;
gpu_frame = gpuArray(basler_frames);
gpu_transfer = toc(gpu_trans1);
display(['whiskertracker_parallel: GPU transfer time was ' ...
    num2str(gpu_transfer) ' seconds.']);

% PreAllocate array of whisker angles, use NaN as a place holder
angle = NaN*ones(1,length(basler_frames));
radon_time1 = tic;
for f = 1:(length(basler_frames)-1);
    % Radon on individual frame
    [R,~] = radon(gpu_frame(:,:,f),theta);
    % Get transformed image from GPU and calculate the variance
    col_var = var(gather(R));
    % Sort the columns according to variance
    ord_var = sort(col_var);
    % Choose the top 0.1*number columns which show the highest variance
    thresh = round(numel(ord_var)*0.9);
    sieve = gt(col_var,ord_var(thresh));
    % Associate the columns with the corresponding whisker angle
    angles = nonzeros(theta.*sieve);
    % Calculate the average of the whisker angles
    angle(f) = mean(angles);
end
radon_time = toc(radon_time1);
display(['whiskertracker_parallel: Whisker Tracking time was ' ...
    num2str(radon_time) ' seconds.']);

inds = isnan(angle)==1;
angle(inds) = [];


