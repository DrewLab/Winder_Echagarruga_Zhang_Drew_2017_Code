function [thresh] = CreatePswfThreshold(pswf)
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Dec 2013
%   Version 1

%   SUMMARY:
%_______________________________________________________________
%   INPUTS:
%_______________________________________________________________
%   OUTPUTS:
%_______________________________________________________________
y = hilbert(diff(pswf));
force = abs(y);
figure;
isok='n';
while strcmp(isok,'y') == 0
    plot(force);
    thresh = input('No Threshold to binarize pressure sensor found. Please enter a threshold: ');
    bin_pswf = binarize_pswf(pswf,thresh);
    bin_inds = find(bin_pswf);
    subplot(211); plot(pswf); axis tight;
    hold on; scatter(bin_inds,max(pswf)*ones(size(bin_inds)),'r.');
    subplot(212); plot(force); axis tight;
    hold on; scatter(bin_inds,max(force)*ones(size(bin_inds)),'r.');
    isok = input('Is this threshold okay? (y/n) ','s');
    hold off;
end