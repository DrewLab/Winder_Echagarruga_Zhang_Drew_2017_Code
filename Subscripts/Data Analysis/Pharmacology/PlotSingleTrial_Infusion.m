function PlotSingleTrial_Infusion(filenames,CBVType)
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
%   separated by behavior (stored in *_Chunked.mat structures) the resting
%   data will be used for normalization. Otherwise, it is normalized by the
%   mean intensity for that trial. The third panel shows a spectrogram of
%   the local field potentials. The bottom panel displays the MUA spike
%   rate.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   filenames - [cell array] list of all filenames
%
%                   CBVType = [string] name of the CBV ROI
%_______________________________________________________________
%   RETURN:                     
%                   None. Function outputs a plot and saves the figure if 
%                   desired.
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
BasFile = dir('*Baselines.mat');
if isempty(BasFile)
    Baselines = [];
else
    load(BasFile.name);
end

for fil = 1:length(filenames);
    
    filename = filenames{fil};
    [animal,hem,FileDate,FileID] = GetFileInfo(filename);
    strdate = ConvertDate(FileDate);
    load(filename)

    Whisk_Fs = ProcData.Fs.wwf;
    CBV_Fs = ProcData.Fs.(CBVType);
    MU_Fs = ProcData.Fs.MUpower;
    PS_Fs = ProcData.Fs.pswf;
    
    % Create time vectors 
    Whisk_time = 1/Whisk_Fs:1/Whisk_Fs:ProcData.Info.TrialDur;
    CBV_time = 1/CBV_Fs:1/CBV_Fs:ProcData.Info.TrialDur;
    Sol_times = [ProcData.Sol.Contra, ProcData.Sol.Ipsi, ...
        ProcData.Sol.Tail, ProcData.Sol.Control];
    MU_time = 1/MU_Fs:1/MU_Fs:ProcData.Info.TrialDur;
    
    % Set plot color scheme
    cmap = [ones(length(ProcData.Sol.Contra),1)*Sol_cmap_numeric(1,:);...
        ones(length(ProcData.Sol.Ipsi),1)*Sol_cmap_numeric(2,:);...
        ones(length(ProcData.Sol.Tail),1)*Sol_cmap_numeric(3,:);...
        ones(length(ProcData.Sol.Control),1)*Sol_cmap_numeric(4,:)];
    
    % Set y-values of scatter points
    sol_y = max(ProcData.Beh.wwf);
    whisk_y = max(ProcData.Beh.wwf)+5;
    move_y = max(ProcData.Beh.wwf)+10;
    
    % Identify locations of detected whisks and movements
    whisktimes = find(ProcData.Beh.Bin_wwf)/Whisk_Fs;
    movetimes = find(ProcData.Beh.Bin_pswf)/PS_Fs;
    
    % Whisker Position, overlay solenoid puff times
    WhiskAngle = ProcData.Beh.wwf;
    subplot(411);
    scatter(Sol_times,sol_y*ones(1,length(Sol_times)),50,cmap,'filled','v');
    hold on;
    scatter(whisktimes,whisk_y*ones(1,length(whisktimes)),'.',...
        'MarkerEdgeColor',[235, 176, 32]/255);
    scatter(movetimes,move_y*ones(1,length(movetimes)),'.',...
        'MarkerEdgeColor',[0.5 0.5 0.5]);
    plot(Whisk_time(1:length(WhiskAngle)),WhiskAngle,'k');
    ylim([min(WhiskAngle)-10 max(WhiskAngle)+10]);
    ylabel(sprintf('Whisker\nAngle'));
    xlim([1 ProcData.Info.TrialDur-1]);
    ax = gca;
    set(ax,'XtickLabel',[])
    hold off;
    
    % Plot Normalized CBV
    subplot(412);
    CBV = ProcData.CBV.(CBVType)(1:ProcData.Info.TrialDur*CBV_Fs);
    if ~isfield(Baselines,CBVType)
        CBVBaseline = mean(CBV);
    else
        if ~isfield(Baselines.(CBVType),strdate)
            CBVBaseline = mean(CBV);
        else
            CBVBaseline = mean(Baselines.(CBVType).(strdate).Means);
        end
    end
    nCBV = CBV/CBVBaseline-1;
    [z,p,k] = butter(CBV_filt_order,CBV_filt_thresh/(CBV_Fs/2),'low');
    [sos,g] = zp2sos(z,p,k);
    filt_nCBV = filtfilt(sos,g,nCBV-mean(nCBV));
    plot(CBV_time,filt_nCBV,'LineWidth',1.5); hold on;
    hold off;
    labl = sprintf('\\DeltaR/R');
    ylim([-0.1 0.1])
    ylabel(labl); 
    xlim([1 ProcData.Info.TrialDur-1]);
    ax = gca;
    set(ax,'XTickLabel',[]);

    % Plot the power in the Gamma Band.
    subplot(413);
    Gam = ProcData.Neuro.Gam;
    plot(MU_time,Gam,'color',[0.5 0.5 0.5]);
    ylim([0 5e-9])
    xlim([0 ProcData.Info.TrialDur]);
    ylabel(sprintf('\\DeltaP/P\n(40-100Hz)'))
    ax = gca;
    set(ax,'XTickLabel',[]);

    
    % Plot the power in the Multi-Unit Band.
    subplot(414);
    MUpower = ProcData.Neuro.MUpower;
    plot(MU_time,MUpower,'color',[0.5 0.5 0.5]);
    xlim([1 ProcData.Info.TrialDur-1]);
    ylim([0 1e-9])
    xlabel('Time (s)')
    ylabel(sprintf('\\DeltaP/P\n(300-3000Hz)'))
end