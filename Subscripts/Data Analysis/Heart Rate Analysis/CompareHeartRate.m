function [] = CompareHeartRate(animals,CBVType)
%   function [] = CompareHeartRate(animals,CBVType)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Generates figures S2c-e. 
%
%_______________________________________________________________
%   PARAMETERS:
%               animals - [cell array] animal IDs
%
%               CBVType - [string] name of the CBV ROI
%_______________________________________________________________
%   RETURN:
%
%_______________________________________________________________
RestDataStD = zeros(1,length(animals));
for a = 1:length(animals)
    animal = animals{a};
    prevdir = cd([animal filesep]);

    RestFile = dir('*_RESTDATA_HR.mat');
    load(RestFile.name)
    RestCriteria.Fieldname = {'Duration','PuffDistance'};
    RestCriteria.Comparison = {'gt','gt'};
    RestCriteria.Value = {14,5};
    [FiltArray] = FilterEvents(RestData.HR,RestCriteria);
    RestDataFilt = RestData.HR.Data(FiltArray);
    RestLength = zeros(size(RestDataFilt));
    RestDataSqDev = zeros(size(RestDataFilt));
    for RDF = 1:length(RestDataFilt)
        RestSnip = RestDataFilt{RDF};
        NaNFilt = not(isnan(RestSnip));
        RestLength(RDF) = sum(NaNFilt);
        RestDataSqDev(RDF) = sum((RestSnip(NaNFilt)-mean(RestSnip(NaNFilt))).^2);
    end
    RestDataStD(a) = sum(RestDataSqDev)/sum(RestLength);
    
    EventFile = dir('*_EVENTDATA_HR.mat');
    load(EventFile.name);
    if a==1
        AllHR1 = zeros(length(animals),size(EventData.HR.stim.Data,2));
        AllHR2 = zeros(length(animals),size(EventData.HR.stim.Data,2));
    end
    NameFilter1 = strcmp(EventData.HR.stim.Name,'Contra');
    HRData1 = EventData.HR.stim.Data(NameFilter1,:);
    NaNFilter1 = not(isnan(sum(HRData1,2)));
    AllHR1(a,:) = mean(HRData1(NaNFilter1,:));
    NameFilter2 = strcmp(EventData.HR.stim.Name,'Control');
    HRData2 = EventData.HR.stim.Data(NameFilter2,:);
    NaNFilter2 = not(isnan(sum(HRData2,2)));
    AllHR2(a,:) = mean(HRData2(NaNFilter2,:));
    cd(prevdir)
end
% NaNinds1 = not(isnan(sum(AllHR1,2)));
NaNinds2 = not(isnan(sum(AllHR2,2)));
timevec = (1:size(EventData.HR.stim.Data,2))/EventData.HR.stim.Fs-...
    EventData.HR.stim.epoch.offset;
pretime = timevec<0;
posttime = timevec>0 & timevec<3;
onset = EventData.HR.stim.epoch.offset*EventData.HR.stim.Fs;

AllHR1_Sub = AllHR1(NaNinds2,:)-(mean(AllHR1(NaNinds2,pretime),2)*ones(1,size(AllHR1,2)));
AllHR2_Sub = AllHR2(NaNinds2,:)-(mean(AllHR2(NaNinds2,pretime),2)*ones(1,size(AllHR2,2)));
AllHRZ1 = AllHR1_Sub./(RestDataStD(NaNinds2)'*ones(1,size(AllHR1,2)));
AllHRZ2 = AllHR2_Sub./(RestDataStD(NaNinds2)'*ones(1,size(AllHR2,2)));

PreHRAvg1 = mean(AllHR1(NaNinds2,pretime),2);
PostHRAvg1 = mean(AllHR1(NaNinds2,posttime),2);

PreHRAvg2 = mean(AllHR2(NaNinds2,pretime),2);
PostHRAvg2 = mean(AllHR2(NaNinds2,posttime),2);

NormHRDifference1 = PostHRAvg1./PreHRAvg1-1;
NormHRDifference2 = PostHRAvg2./PreHRAvg2-1;

PostHRMax_Z1 = max(AllHRZ1(:,onset:end),[],2);
PostHRMax_Z2 = max(AllHRZ2(:,onset:end),[],2);

%% Figure S2C
figure; 
subplot(121);
plot(timevec,mean(AllHR1(NaNinds2,:)));
hold on; plot(timevec,mean(AllHR2(NaNinds2,:)),'Color',[0.7 0.7 0.7]);
plot(timevec,mean(AllHR1(NaNinds2,:))+std(AllHR1(NaNinds2,:)),'k:');
plot(timevec,mean(AllHR1(NaNinds2,:))-std(AllHR1(NaNinds2,:)),'k:');
legend({'Contalateral Stim.','Auditory Stim.'})
ylim([7 15]);
xlim([min(timevec) max(timevec)])
ylabel('Heart Rate (Hz)')
xlabel('Peristimulus Time (s)')
set(gcf,'name','Figure S2c','numbertitle','off')
title(['HR inc.: Contra ' num2str(round(mean(NormHRDifference1)*1000)/10) ...
    '\pm' num2str(round(std(NormHRDifference1)*1000)/10) '% Aud. '...
    num2str(round(mean(NormHRDifference2)*1000)/10) ...
    '\pm' num2str(round(std(NormHRDifference2)*1000)/10)]);

subplot(144);
scatter(ones(1,sum(NaNinds2)),PostHRMax_Z1);
hold on; 
scatter(2*ones(1,sum(NaNinds2)),PostHRMax_Z2,'MarkerEdgeColor',[0.7 0.7 0.7]);
scatter([1,2],[mean(PostHRMax_Z1), mean(PostHRMax_Z2)],'ko','MarkerFaceColor','k');
ylabel('Post-stimulus max HR (zscore)')
xlim([0.5 2.5])
ylim([0 2.5])
ax=gca;
set(ax,'XTick',[1,2],'XTickLabel',{'Contra.','Auditory'});
[~,p,~,stat] = ttest(PostHRMax_Z1,PostHRMax_Z2);
title(['p=' num2str(p) ', t(' num2str(stat.df) ')=' num2str(stat.tstat)])


%% Plot CBV
for a = 1:length(animals)
    animal = animals{a};
    prevdir = cd([animal filesep]);
    EventFile = dir(['*_EVENTDATA_' CBVType '.mat']);
    load(EventFile.name);
    NameFilter1 = strcmp(EventData.(CBVType).stim.Name,'Contra');
    CBVData1 = EventData.(CBVType).stim.NormData(NameFilter1,:)-1;
    if a==1
        AllCBV1 = zeros(length(animals),size(EventData.(CBVType).stim.Data,2));
        AllCBV2 = zeros(length(animals),size(EventData.(CBVType).stim.Data,2));
    end
    AllCBV1(a,:) = mean(CBVData1);
    NameFilter2 = strcmp(EventData.(CBVType).stim.Name,'Control');
    CBVData2 = EventData.(CBVType).stim.NormData(NameFilter2,:)-1;
    AllCBV2(a,:) = mean(CBVData2);
    cd(prevdir)
end
timevec = (1:size(EventData.(CBVType).stim.Data,2))/EventData.(CBVType).stim.Fs-...
    EventData.(CBVType).stim.epoch.offset;
pretime = timevec<0;
posttime = timevec>0 & timevec<3;
NaNinds2 = not(isnan(sum(AllCBV2,2)));

PreCBVAvg1 = mean(AllCBV1(NaNinds2,pretime),2);
PostCBVMin1 = min(AllCBV1(NaNinds2,posttime)-PreCBVAvg1*ones(1,sum(posttime)),[],2);

PreCBVAvg2 = mean(AllCBV2(NaNinds2,pretime),2);
PostCBVMin2 = min(AllCBV2(NaNinds2,posttime)-PreCBVAvg2*ones(1,sum(posttime)),[],2);

figure; 
set(gcf,'name','Figure S2d','numbertitle','off')
subplot(121); 
plot(timevec,mean(AllCBV1(NaNinds2,:)));
hold on; plot(timevec,mean(AllCBV2(NaNinds2,:)),'Color',[0.7 0.7 0.7]);
plot(timevec,mean(AllCBV1(NaNinds2,:))+std(AllCBV1(NaNinds2,:)),'k:');
plot(timevec,mean(AllCBV1(NaNinds2,:))-std(AllCBV1(NaNinds2,:)),'k:');
legend({'Contra','Auditory'},'location','southwest')
xlim([min(timevec) max(timevec)])
ylabel('\DeltaR/R')
xlabel('Peristimulus Time (s)')

subplot(144);
scatter(ones(1,sum(NaNinds2)),PostCBVMin1);
hold on; 
scatter(2*ones(1,sum(NaNinds2)),PostCBVMin2,'MarkerEdgeColor',[0.7 0.7 0.7]);
scatter([1,2],[mean(PostCBVMin1), mean(PostCBVMin2)],'ko','MarkerFaceColor','k');
xlim([0 3])
ylim([-0.04 0])
ylabel('\DeltaR/R')
xlim([0.5 2.5])
ax=gca;
set(ax,'XTick',[1,2],'XTickLabel',{'Contra.','Auditory'});
[~,p2,~,stat] = ttest(PostCBVMin1,PostCBVMin2);
title(['p=' num2str(p2) ', t(' num2str(stat.df) ')=' num2str(stat.tstat)])

%% Figure S2e

% NormHRDifference = PostHRAvg2./PreHRAvg2-PostHRAvg1./PreHRAvg2;
% NormCBVMaxDifference = PostCBVMin2-PostCBVMin1;
% figure; scatter(NormHRDifference,NormCBVMaxDifference);
% [b,stats] = robustfit(NormHRDifference,NormCBVMaxDifference);
% linefit = polyval(flipud(b),NormHRDifference);
% hold on; plot(NormHRDifference,linefit);
% SDline1 = polyval(flipud(b)+flipud(stats.se),NormHRDifference);
% SDline2 = polyval(flipud(b)-flipud(stats.se),NormHRDifference);
% plot(NormHRDifference,SDline1,'r:')
% plot(NormHRDifference,SDline2,'r:')
% hold off
% xlim([-0.04 0.03])
% ylim([-0.04 0.03])

NormHRDifference1 = PostHRMax_Z1-PostHRMax_Z2;
NormCBVMaxDifference = abs(PostCBVMin1)-abs(PostCBVMin2);

figure; scatter(NormHRDifference1,NormCBVMaxDifference,'ko');
[b,stats] = robustfit(NormHRDifference1,NormCBVMaxDifference);
pval = stats.p(2);
df = stats.dfe;
tstat = stats.t(2);
linefit = polyval(flipud(b),NormHRDifference1);
hold on; plot(NormHRDifference1,linefit,'k');
SDline1 = polyval(flipud(b)+flipud(stats.se),NormHRDifference1);
SDline2 = polyval(flipud(b)-flipud(stats.se),NormHRDifference1);
plot(NormHRDifference1,SDline1,'k:')
plot(NormHRDifference1,SDline2,'k:')
hold off
xlim([-2 2])
ylim([0 0.02])
xlabel('Max Contra HR - Max Auditory HR')
ylabel('Max Contra CBV - Max Auditory CBV')
title(['Slope: ' num2str(round(b(2)*1000)/1000) '\pm' ...
    num2str(round(stats.se(2)*1000)/1000) '; p=' ...
    num2str(round(pval*1000)/1000) '; t(' num2str(df) ')=' num2str(tstat)]);
set(gcf,'name','Figure S2e','numbertitle','off')