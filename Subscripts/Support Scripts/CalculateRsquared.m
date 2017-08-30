function r2 = CalculateRsquared(pred,act)
%   function r2 = CalculateRsquared(pred,act)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Calculates the coefficient of determination between two
%   signals. Formula was obtained from 
%       http://en.wikipedia.org/wiki/Coefficient_of_determination
%   and is 1 - [(sum of squares of residuals)/(total sum of squares)]
%
%_______________________________________________________________
%   PARAMETERS:
%                   pred - [array] column vector of the model which 
%                   predicts the actual signal
%
%                   act - [array] the actual signal as a column
%                   vector
%_______________________________________________________________
%   RETURN:
%                   r2 - [double] the coefficient of determination
%_______________________________________________________________

% Error Variance
SSE = sum((act-pred).^2);

% Total Variance
SST = sum((act-(ones(size(act,1),1)*mean(act))).^2);

% Check that the sum of the residuals is small compared to SSE + SSR
r2 = ones(size(SSE))-SSE./SST;