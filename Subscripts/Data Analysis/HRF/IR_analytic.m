% IR_analytic.m
% Bingxing Huo, May 2013
% Modified from:
% IOS_analytic.m, August 2012
% The script aims to perform direct calculation to solve for impulse
% response in a linear time invariant system.
% Input:
%   datain: input vector. In form of samples-by-1. 
%   dataout: output matrix. In form of samples-by-channels.
% Output:
%   IR_est: a structure containing information of estimated impulse response 
%       - IR: impulse response vector.
%       - tshift: copy the time shift as defined in the input.
function IR_est=IR_analytic(datain,dataout,tshift,IR_length)
if size(datain,1)~=size(dataout,1)
    error('Input and output should be of same lengths.')
end
datain=datain(tshift+1:end); % shift input time forward, assuming non-causal system
dataout=dataout(1:end-tshift); % truncate output time to keep same lengths
N=length(datain);
TL=[1:N];
if nargin<4
    N_IR=round(N/5); % Let the length of IR be 1/5 of the length of the input data
else
    if isempty(IR_length==1)
        N_IR=round(N/5); 
    else
    N_IR=IR_length;
    end
end
channels=size(dataout,2);
% 1. solving for impulse response (IR)
datain_col=zeros(2*N-1,1); % allocate the first column of the Toeplitz matrix
datain_col(TL)=datain; % first column of the Toeplitz matrix
datain_row=zeros(N_IR,1); % allocate the first row of the Toeplitz matrix
datain_row(1)=datain_col(1); % first row of the Toeplitz matrix, all zeros except first entry
datainM=toeplitz(datain_col,datain_row); % generate the Toeplitz matrix
% datainM=[ones(2*N-1,1),datainM]; % a constant term
impulseM_sym=datainM'*datainM;
% impulseM_inv=inv(impulseM_sym);
response=zeros(2*N-1,channels);
response(TL,:)=dataout(TL,:);
IR_est.IR = impulseM_sym\datainM'*response;
IR_est.tshift=tshift;
end