function [S,fr] = TrialPowerSpectrum(Data,Fs,dataType,fpass)
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

params.Fs = Fs;
params.tapers = [5 9];
params.fpass = fpass;
SpecData = detrend(Data);
[S,fr] = mtspectrumc(SpecData,params);