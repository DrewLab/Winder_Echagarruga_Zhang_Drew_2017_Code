function [XC,lags,Reshuf95] = RestingXCorr_CBVvsHR(CBVType)
%   function [XC,lags,Reshuf95] = RestingXCorr_CBVvsHR(CBVType)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Calculates the cross correlation between HR and CBV
%   fluctuations during periods of rest.
%
%_______________________________________________________________
%   PARAMETERS:
%               CBVType - [string] designates the CBV ROI to be used
%_______________________________________________________________
%   RETURN:
%               XC - [array] the correlation coefficients at various
%               temporal lags
%
%               lags - [array] time vector of the temporal lags.
%
%               Reshuf95 - [array] 95% confidence interval calculated from
%               reshuffled data
%_______________________________________________________________

filt_ord = 3;
filt_cutoff = 2; % Hz
maxlags = 5;
% RestFile = ls('*RestData_HR.mat');
% load(RestFile)
RestFile = dir('*RESTDATA_HR.mat');
load(RestFile.name)
RestCriteria.Fieldname = {'Duration','PuffDistance'};
RestCriteria.Comparison = {'gt','gt'};
RestCriteria.Value = {14,5};
[FiltArray] = FilterEvents(RestData.HR,RestCriteria);
RestDataCriteria = RestData.HR.Data(FiltArray);

RestFile = dir(['*RESTDATA_' CBVType '.mat']);
load(RestFile.name)
CBVDataCriteria = RestData.(CBVType).Data(FiltArray);

Fs = RestData.(CBVType).Fs;
[z,p,k] = butter(filt_ord,filt_cutoff/(Fs/2),'low');
[sos,g] = zp2sos(z,p,k);

% Remove NaN from the Rest Struct
HRcell = cell(1,length(RestDataCriteria));
CBVcell = cell(1,length(RestDataCriteria));
for c = 1:length(RestDataCriteria)
    NaNFilt = not(isnan(RestDataCriteria{c}));
    HRcell{c} = RestDataCriteria{c}(NaNFilt)-mean(RestDataCriteria{c}(NaNFilt));
%     DataStruct.HR.Data{c} = HR;
    CBVcell{c} = filtfilt(sos,g,CBVDataCriteria{c}(NaNFilt) - mean(CBVDataCriteria{c}));
%     DataStruct.CBV.Data{c} = filtfilt(sos,g,CBV);
end
AllHR = [HRcell{:}];
AllCBV = [CBVcell{:}];
[XC,lags] = xcorr(AllCBV,AllHR,maxlags*Fs,'coeff');
lags = lags/Fs;

%% Calculate the 95% confidence interval for the cross correlation
repetitions = 1000;
CC = zeros(1,repetitions);
for rep = 1:repetitions
    Random_inds = ceil(length(HRcell)*...
        rand(1,length(HRcell)));
    Reshuffled = randperm(length(Random_inds));
    Reshuf1 = [HRcell{Random_inds}];
    Reshuf2 = [CBVcell{Random_inds(Reshuffled)}];
    r = corrcoef(Reshuf2,Reshuf1);
    CC(rep) = r(2,1);
end

Reshuf95(1) = quantile(CC,0.025);
Reshuf95(2) = quantile(CC,0.975);