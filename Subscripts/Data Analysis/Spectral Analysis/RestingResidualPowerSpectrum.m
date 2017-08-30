function [S,f] = RestingResidualPowerSpectrum(Actual,Predicted,Fs)
%   function [] =  ()
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: 
%   
%_______________________________________________________________
%   PARAMETERS:             
%                               
%_______________________________________________________________
%   RETURN:                     
%                               
%_______________________________________________________________

% RestBuffer = 2;
params.tapers = [1 1];
params.fpass = [0.1 2];
params.Fs = Fs;
params.err = [1 0.05];
        
CellRest = cell(size(Actual));
CellSeg = cell(size(Actual));
CellPred = cell(size(Predicted));
for a = 1:length(Actual)
    % Combine measured data
%     RestBuffer_Ind = RestBuffer*Fs;
%     clippedRest = Actual{a}(RestBuffer_Ind:end);
%     RestOffset = mean(clippedRest);
%     CellRest{a} = detrend(clippedRest-RestOffset);
    CellRest{a} = detrend(Actual{a}-mean(Actual{a}));

    % Combine predicted data
%     RestBuffer_Ind = RestBuffer*Fs;
%     clippedPred = Predicted{a}(RestBuffer_Ind:end);
%     PredOffset = mean(clippedPred);
%     CellPred{a} = detrend(clippedPred-PredOffset);
    CellPred{a} = detrend(Predicted{a}-mean(Predicted{a}));
    
    % Track the segment starts/end
    CellSeg{a} = zeros(size(CellRest{a}));
    CellSeg{a}(1) = 1;
    CellSeg{a}(end) = -1;
end

AllRest = [CellRest{:}];
AllPred = [CellPred{:}];
SegMarkers = [CellSeg{:}];
SegStarts = find(SegMarkers==1);
SegStops = find(SegMarkers==-1);

% Calculate residuals
Residual = AllRest - AllPred;

[S,f]= mtspectrumc_unequal_length_trials(Residual',[10 10],params,...
    [SegStarts', SegStops']);
% Calculate Power Spectrum of residuals

% [CBVPow.Act,Freqs.Act] = mtspectrumc(CellRest,params);
% [CBVPow.Pred,Freqs.Pred] = mtspectrumc(CellPred,params);
% [CBVPow.resid,Freqs.resid] = mtspectrumc(Residual-mean(Residual),params);