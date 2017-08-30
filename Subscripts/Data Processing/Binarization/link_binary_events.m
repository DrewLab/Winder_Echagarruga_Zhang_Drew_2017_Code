function [linked_wf] = link_binary_events(bin_wf,dCrit)
%   [linked_wf] = link_binary_events(bin_wf,dCrit)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Takes a binary waveform and links the peaks which
%   occur within a close period of time of each other creating a single 
%   peak.   
%_______________________________________________________________
%   PARAMETERS:             
%                   bin_wf - [Array] The binary waveform
%
%                   dCrit - [1x2 Array] The distances, in samples, between 
%                   the falling edge of the previous event and the rising 
%                   edge of the current waveform. 
%                       Input should be given as a 2D array: 
%                           [dCrit for breaks ,dCrit for peaks]. 
%                       Any events with breaks between them less than dCrit
%                           for breaks will be merged. 
%                       Any events that have duration less than dCrit for 
%                           peaks will be removed. 
%                       Either dCrit value can be zero to avoid merging or 
%                           erasing events.                    
%_______________________________________________________________
%   RETURN:                     
%                   linked_wf - [Array] The new, linked binary waveform as 
%                   an array.
%_______________________________________________________________

%% Identify Edges, control for trial start/stop
d_bin_wf = diff(gt(bin_wf,0));
up_ind = find(d_bin_wf == 1)+1;
down_ind = find(d_bin_wf == -1);
if bin_wf(end)>0
    down_ind = [down_ind length(bin_wf)];
end
if bin_wf(1)>0
    up_ind = [1 up_ind];
end

%% Link periods of bin_wf==0 together if less than dCrit(1)
% Calculate time between events
brktimes = up_ind(2:length(up_ind))-down_ind(1:(length(down_ind)-1));
% Identify times less than user-defined period
sub_dCrit_Downs = find(lt(brktimes,dCrit(1)));

% Link any identified breaks together
if isempty(sub_dCrit_Downs) == 0;
    for d = 1:length(sub_dCrit_Downs);
        start = down_ind(sub_dCrit_Downs(d));
        stop = up_ind(sub_dCrit_Downs(d)+1);
        bin_wf(start:stop) = 1;
    end
end

% Recalculate the edges
d_bin_wf = diff(gt(bin_wf,0));
up_ind = find(d_bin_wf == 1);
down_ind = find(d_bin_wf == -1);
if bin_wf(end)>0
    down_ind = [down_ind length(bin_wf)];
end
if bin_wf(1)>0
    up_ind = [1 up_ind];
end


%% Link periods of bin_wf==1 together if less than dCrit(2)
hitimes = down_ind - up_ind;
blips = find(lt(hitimes,dCrit(2))==1);
if isempty(blips) == 0;
    for b = 1:length(blips);
        start = up_ind(blips(b));
        stop = down_ind(blips(b));
        bin_wf(start:stop) = 0;
    end
end
linked_wf = bin_wf;