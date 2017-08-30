function [] = PlotSpatialXCorrExample(animal)
%   function [] = PlotSpatialXCorrExample(animal)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Plots a spatial map of the correlation coefficients from a
%   cross correlation between MUA and lagged CBV for a single animal. The
%   CCFile variable was created using the SpatialCrossCorrelation.m script
%   and a description of the methods are detailed in the Methods section.
%   
%_______________________________________________________________
%   PARAMETERS:             
%               animal - [string] the animal ID              
%_______________________________________________________________
%   RETURN:                     
%                               
%_______________________________________________________________

%
clc
display('Plotting Figure 1d...')

% Load and setup
prevdir = cd([animal filesep]);
% CCFile = ls([animal '_GlobalSubtractSpatialXCorr_MUpower_Rest.mat']);
CCFile = dir([animal '_GlobalSubtractSpatialXCorr_MUpower_Rest.mat']);
load(CCFile.name)
Fs = 30; % Sampling Frequency = 30 Hz;
CCLagTimes = [0.5 1.5];

% Plot the median spatial map
timelags = lags/Fs; % Get lags in seconds
LagInds = timelags>=CCLagTimes(1) & timelags<=CCLagTimes(2);
Frame = median(XC_Map(:,:,LagInds),3);

% Plot example image
figure; 
set(gcf,'name','Figure 1d','numbertitle','off')
subplot(221)
imagesc(IndFrame(:,:,1).*mask(:,:,1)); 
axis image; colormap(gca,'gray');
[rows,cols] = find(mask(:,:,1));
xlim([min(cols),max(cols)])
ylim([min(rows),max(rows)])
ax = gca;
set(ax,'XTick',[]);
set(ax,'YTick',[]);

ax = subplot(222);
imagesc(Frame);
axis image; caxis([-0.1 0.1]); colormap(gca,'parula');
SpecPosition = get(ax,'Position');
xpos = SpecPosition(1);
xwidth = SpecPosition(3);
cbar = colorbar('location','southoutside','position',[xpos, 0.53, xwidth,...
    0.04]);
ylabel(cbar,sprintf('Median Corr. Coeff.\n(\\DeltaR/R vs MUA)'));
ax = gca;
set(ax,'XTick',[],'YTick',[])
impoly(gca,[ROI.Barrel.xi, ROI.Barrel.yi]);
impoly(gca,[ROI.Forepaw.xi, ROI.Forepaw.yi]);


BarrelCC = zeros(1,size(XC_Map,3),'double');
for fr = 1:size(XC_Map,3)
    indframe = XC_Map(:,:,fr);
    BarrelCC(fr) = mean(indframe(ROI.Barrel.Mask));
end

ForePawCC = zeros(1,size(XC_Map,3),'double');
for fr = 1:size(XC_Map,3)
    indframe = XC_Map(:,:,fr);
    ForePawCC(fr) = mean(indframe(ROI.Forepaw.Mask));
end

subplot('position',[xpos, 0.1, xwidth, 0.3]);
plot(timelags,BarrelCC,'Color',[0,0.447,0.741],'linewidth',1.5);
hold on; plot(timelags,ForePawCC,'Color',[0.929,0.694,0.125],'linewidth',1.5);
plot(CCLagTimes(1)*ones(size(-1:0.1:1)),-1:0.1:1,'k:')
plot(CCLagTimes(2)*ones(size(-1:0.1:1)),-1:0.1:1,'k:')
ylim([-0.2 0.1])
xlim([0 max(timelags)])
xlabel('Lags (s)')
ylabel('Corr. Coeff.')
legend({'Barrels','ForeLimb'},'Location','northoutside',...
    'Orientation','horizontal')

cd(prevdir);