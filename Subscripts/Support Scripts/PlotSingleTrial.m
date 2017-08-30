function PlotSingleTrial(filenames,CBVType)
%   function [] = PlotSingleTrial(filenames)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Plots the data summary for a single trial. Top panel shows
%   the whisker angle overlaid with indicators of puffing times, locations
%   of detected whisks using a threshold, and locations of detected body
%   movement using a threshold of the force sensor. Second panel shows the
%   normalized light reflection intensity (CBV). If data has already been
%   separated by behavior (stored in *_Baselines.mat structures) the resting
%   data will be used for normalization. Otherwise, it is normalized by the
%   mean intensity for that trial. The third panel shows a spectrogram of
%   the local field potentials. The bottom panel displays the MUA spike
%   rate.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   filenames - [cell array] list of all filenames   
%
%                   CBVType - [string] fieldname corresponding to the
%                   desired CBV ROI.
%_______________________________________________________________
%   RETURN:                     
%                   None. Function outputs a plot
%_______________________________________________________________
%% Load and Setup
% Control for non-cell input
if iscell(filenames)==0
    rows = size(filenames,1);
    if rows == 1;
        temp = {filenames};
        filenames = temp;
    else
        temp = mat2cell(filenames,ones(1,size(filenames,1)));
        filenames = temp;
    end
end

Sol_cmap_numeric = [0.05,0.45,0.72;0.5,0.5,0.5;0,0,0;0.6,0.79,0.24];
CBV_filt_thresh = 2;
CBV_filt_order = 5;
% BasFile = ls('*Baselines.mat');
BasFile = dir('*Baselines.mat');
if isempty(BasFile)
    Baselines = [];
else
%     load(BasFile);
    load(BasFile.name);
end
params.tapers = [5 9];
params.fpass = [5 150];
movingwin = [1 1/5];
for fil = 1:length(filenames);
    filename = filenames{fil};
    [animal,hem,FileDate,FileID] = GetFileInfo(filename);
    strdate = ConvertDate(FileDate);
    load(filename)
    Whisk_Fs = ProcData.Fs.wwf_fs;
    CBV_Fs = ProcData.Fs.BarrelCBV_fs;
    MU_Fs = ProcData.Fs.MUpower_fs;
    PS_Fs = ProcData.Fs.pswf_fs;
    Spec_Fs = ProcData.Fs.Wideband_LFP_fs;

    % Create time vectors
    Whisk_time = 1/Whisk_Fs:1/Whisk_Fs:ProcData.TrialDur;
    CBV_time = 1/CBV_Fs:1/CBV_Fs:ProcData.TrialDur;
    Sol_times = [ProcData.Sol.Contra, ProcData.Sol.Ipsi, ...
        ProcData.Sol.Tail, ProcData.Sol.Control];
    MU_time = 1/MU_Fs:1/MU_Fs:ProcData.TrialDur;
    % Set plot color scheme
    cmap = [ones(length(ProcData.Sol.Contra),1)*Sol_cmap_numeric(1,:);...
        ones(length(ProcData.Sol.Ipsi),1)*Sol_cmap_numeric(2,:);...
        ones(length(ProcData.Sol.Tail),1)*Sol_cmap_numeric(3,:);...
        ones(length(ProcData.Sol.Control),1)*Sol_cmap_numeric(4,:)];
    % Set y-values of scatter points
    sol_y = max(ProcData.wwf);
    whisk_y = max(ProcData.wwf)+5;
    move_y = max(ProcData.wwf)+10;
    % Identify locations of detected whisks and movements
    whisktimes = find(ProcData.Bin_wwf)/Whisk_Fs;
    movetimes = find(ProcData.Bin_pswf)/PS_Fs;
    % Whisker Position, overlay solenoid puff times
    WhiskAngle = ProcData.wwf(1:ProcData.TrialDur*Whisk_Fs);
    
    figure;
    set(gcf,'name','Figure 1f','numbertitle','off')
    subplot(411);
    scatter(Sol_times,sol_y*ones(1,length(Sol_times)),50,cmap,'filled','v');
    hold on;
    scatter(whisktimes,whisk_y*ones(1,length(whisktimes)),'.',...
        'MarkerEdgeColor',[235, 176, 32]/255);
    scatter(movetimes,move_y*ones(1,length(movetimes)),'.',...
        'MarkerEdgeColor',[0.5 0.5 0.5]);
    plot(Whisk_time,WhiskAngle,'k');
    ylim([min(WhiskAngle)-10 max(WhiskAngle)+10]);
    ylabel(sprintf('Whisker\nAngle'));
    xlim([1 ProcData.TrialDur-1]);
    ax = gca;
    set(ax,'XTickLabel',[]);
    hold off;
    
    % Plot Normalized CBV
    subplot(412);
    CBVcolors = {'r','k','b','m'};
    CBV = ProcData.(CBVType)(1:ProcData.TrialDur*CBV_Fs);
    if isempty(Baselines)
        CBVBaseline = mean(CBV);
    else
        CBVBaseline = mean(Baselines.(CBVType).(strdate).Means);
    end
    nCBV = detrend(CBV)/CBVBaseline-1;
    [z,p,k] = butter(CBV_filt_order,CBV_filt_thresh/(CBV_Fs/2),'low');
    [sos,g] = zp2sos(z,p,k);
    filt_nCBV = filtfilt(sos,g,nCBV-mean(nCBV));
    plot(CBV_time,filt_nCBV,'k','LineWidth',2); hold on;
    hold off;
    labl = sprintf('Normalized\nReflectance');
    ylim([-0.06 0.05])
    ylabel(labl);
    xlim([1 ProcData.TrialDur-1]);
    ax = gca;
    set(ax,'XTickLabel',[]);
    
    % Plot the power in the Multi-Unit Band.
    subplot(414);
    MUpower = ProcData.MUpower;
    plot(MU_time,MUpower/mean(MUpower)-1,'k');
    xlim([1 ProcData.TrialDur-1]);
    xlabel('Time (s)')
    ylabel(sprintf('\\DeltaP/P\n(0.3-3kHz'))
    ylim([-1 5])
    
    % Plot Neural Spectrogram
    Specax = subplot(413);
    params.Fs = Spec_Fs;
    [S,t,f] = mtspecgramc(ProcData.Wideband_LFP,movingwin,params);
%     zpad1 = zeros(floor(t(1)/movingwin(2)),size(S,2));
    SNorm = S./(ones(size(S,1),1)*mean(S,1));
    imagesc(t,f,SNorm'); axis xy;
    xlim([1 ProcData.TrialDur-1]);
    ylabel(sprintf('Frequency\n(Hz)'));
    caxis([0 2])
    ax = gca;
    set(ax,'XTickLabel',[]);
    colormap('parula')
    
    % Add the colorbar
    SpecPosition = get(Specax,'Position');
    ypos = SpecPosition(2);
    yheight = SpecPosition(4);
    SpecXBoundary = SpecPosition(1)+SpecPosition(3)+0.01;
    cbar = colorbar('Position',[SpecXBoundary, ypos, 0.02,...
        yheight]);
    ylabel(cbar,'\DeltaP/P')
end