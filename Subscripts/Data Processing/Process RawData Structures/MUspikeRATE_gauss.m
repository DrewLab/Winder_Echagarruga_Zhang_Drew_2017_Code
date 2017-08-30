function [Spike_Rate,Spiketimes] = MUspikeRATE_gauss(MU_data,thresh,Fs,stdev)
%%      [Spike_Rate,Spiketimes] = MUspikeRATE_gauss(MU_data,thresh,Fs,stdev)
%
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Aug 2013
%
%   SUMMARY:    Calculates the multi-unit spike rate by convolving a train
%   of spikes with a gaussian function.
%_________________________________________________________________________
%   INPUTS:                 MU_data - Vector of raw multi-unit recording
%
%                           thresh - # of standard deviations above which
%                           the signal can be considered a spike. Thresh
%                           will be applied as both a negative and positive
%                           value
%
%                           Fs - [double] Sampling frequency during neural
%                           acquisition
%
%                           stdev - [double] standard deviation of MU data
%_________________________________________________________________________
%   OUTPUTS:                Spike_Rate - [array] spike rates as a function
%                           of time
%
%                           Spiketimes - [array] spike times in samples
%_________________________________________________________________________

%% Get binary spikes

    % Subtract mean, half-wave rectify signal, compare to thresh
    % Per Dayan/Abbot, spikes represented as delta functions -> integral
    % equals 1. Multiply height by Fs.
    
spikes1 = lt(MU_data,stdev*(-1)*thresh)*Fs;
spikes2 = gt(MU_data,stdev*thresh)*Fs;
all_spikes = spikes1+spikes2;

%% Set Minimum time between spikes

    % find spike times, take derivative, eliminate all "spikes" within 
    % the minimum time of each other. Avoids counting spikes more than
    % once.
    
min_time = 0.001; % 1 ms minimum time between spikes 

All_spiketimes = find(all_spikes);
Spiketimes = All_spiketimes;
for s = 1:length(All_spiketimes)
    if ismember(All_spiketimes(s),Spiketimes)
        eraser = and(Spiketimes>All_spiketimes(s),...
            Spiketimes<(All_spiketimes(s)+round(min_time*Fs)));
        Spiketimes(eraser) = [];
    end
end
spikes = zeros(1,length(all_spikes));
spikes(Spiketimes) = 1;

%% Create Gaussian

    % Formula for gaussian taken from Dayan/Abbott p.13. The variable "wid"
    % gives the amount of time to be contained in the middle ~99% of the
    % distribution. "tau" gives the time
    
wid = 0.1; % Lowest frequency resolvable is wid*Fs Hz.
sig = round(wid/6*Fs); % ~99 percent of weight given to values within 100 ms of center
tau = -round(wid*Fs):round(wid*Fs);
gau = 1/(sig*sqrt(2*pi()))*exp(-1*(tau.^2)/(2*sig^2));

%% Get Multi_unit Firing Rate

    % Convolve spike train with gaussian, extract only valid portion,
    % divide by width of the window.
    
Spike_Rate = conv(double(spikes*Fs),gau,'same');


