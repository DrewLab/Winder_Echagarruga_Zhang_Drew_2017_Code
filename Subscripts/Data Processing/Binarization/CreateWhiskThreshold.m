function [thresh1,thresh2] = CreateWhiskThreshold(angl, fs, strday)
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Dec 2013
%   Version 1

%   SUMMARY:
%_______________________________________________________________
%   INPUTS:
%_______________________________________________________________
%   OUTPUTS:
%_______________________________________________________________


isok='n';

% High pass filter the signal to pull out individual whisks
[z,p,k] = butter(4,2/(fs/2),'high');
[sos,g] = zp2sos(z,p,k);
Filtwwf = filtfilt(sos,g,angl);
Pos = abs(Filtwwf);
while strcmp(isok,'y') == 0
    close all; plot(Pos);
    title(strday);
    thresh2 = input('No Threshold for volitional whisks found. Please enter a threshold: ');
    thresh1 = input('No Threshold for resting behavior found. Please enter a threshold: ');
    bin_wwf = binarize_wwf(angl,fs,thresh1,thresh2);
    
    Th1 = find(bin_wwf==1);
    Th2 = find(bin_wwf==0.5);
    ax1 = subplot(211);
    plot(angl); axis tight;
    hold on;
    scatter(Th2,15*ones(size(Th2)),'c.');
    scatter(Th1,15*ones(size(Th1)),'k.');
    ylim([-30 30])
    ylabel('Position')
    title(strday);
    hold off;
    ax2 = subplot(212);
    plot(Pos);
    hold on;
    scatter(Th1,15*ones(size(Th1)),'b.');
    scatter(Th2,15*ones(size(Th2)),'c.');
    ylim([0 10]);
    ylabel('Filtered Position')
    hold off;
    linkaxes([ax1,ax2],'x');
    isok = input('Is this threshold okay? (y/n) ','s');
end
