function [CBVPred,Fs] = CalculatePredictedCBV(ProcData,root,CBVField,HRF,Baselines)
%   [CBVPred,Fs] = CalculatePredictedCBV(ProcData,root,CBVField,HRF,Baselines)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Uses hemodynamic response functions for a given datatype
%   and behavior to predict the measured CBV.
%
%_______________________________________________________________
%   PARAMETERS:
%                   ProcData - [Struct] contains the processed trial data
%                   for a single file
%
%                   root - [string] datatype used to predict CBV
%
%                   CBVField - [string] The CBV ROI used for the HRF
%                   calculation.
%
%                   HRF - [array] 
%_______________________________________________________________
%   RETURN:
%                  CBVPred - [array] the predicted CBV, resulting from the
%                  convoluation of "root" with the HRF.
%
%                  Fs - [double] the sampling frequency of CBVPred
%_______________________________________________________________
 
    root_fs = ProcData.Fs.([root '_fs']);
    CBV_fs = ProcData.Fs.([CBVField '_fs']);
    
    % Normalize the data measurement
    root_baseline = mean(Baselines.(root).Means);
    norm_root = detrend(ProcData.(root)/root_baseline-1);
    CBV_baseline = mean(Baselines.(CBVField).Means);
    normCBV = detrend(ProcData.(CBVField)/CBV_baseline-1);
    
    % Match the lengths of data
    [Root,CBV] = MatchDataLengths(norm_root,normCBV);
    
    % Apply the HRF to the measured data
    [~,CBVPred] = ConvolveHRF(HRF,Root,CBV,0);
    Fs = min([root_fs, CBV_fs]);
end
