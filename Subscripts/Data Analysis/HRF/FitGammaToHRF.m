function [gammafit,coef,t,HRF] = FitGammaToHRF(RawHRF,Fs)
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

options = optimset('MaxFunEvals',2e3,'MaxIter',2e3,'TolFun',1e-7,'TolX',1e-7);
x = [3e-3,1.25,0.5];
m2fun = @(raw,fit) sum((raw-fit).^2);
for r = 1:size(RawHRF.timevec,1)
    HRFTime = RawHRF.timevec(r,:);
    timeinds = and(HRFTime>=0, HRFTime<=5);
    fit_inds = and(HRFTime>=0, HRFTime<=2.5);
    zeroind = HRFTime(timeinds)==0;
    t = HRFTime(timeinds);
    HRF = RawHRF.HRF(r,timeinds);
    HRFDC = HRF(zeroind);
    HRF = HRF-HRFDC;
    FitHRF = RawHRF.HRF(r,fit_inds)-HRFDC;
    fittime = HRFTime(fit_inds);
    coef = lsqcurvefit(gammafun,x,fittime,FitHRF,[],[],options);
    gammafit = gammafun(coef,t);
end

function [m2error]=CompareFit(gam_params,RawHRF,Fs,HRFDur)
t = 0:1/Fs:HRFDur;
a = ((gam_params(2)/gam_params(3))^2*8*log10(2));
beta = ((gam_params(3)^2)/gam_params(2)/8/log10(2));
gamma = gam_params(1)*(t/gam_params(2)).^a.*exp((t-gam_params(2))/(-1*beta));
m2error = mean((gamma-HRF).^2);