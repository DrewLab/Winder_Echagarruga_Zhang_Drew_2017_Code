function [Stats] = CompareResidualRMSE(animals,NeurType,CBVType,CBVPred)
%   function [Stats] = CompareResidualRMSE(animals,NeurType,CBVType,CBVPred)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Compares the variance of the prediction residuals among
%   behavioral categories
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   animals - [cell array] animal IDs
%
%                   NeurType - [string] Type of neural activity used for
%                   the prediction
%
%                   CBVType - [string] fieldname corresponding to the CBV
%                   ROI for the actual CBV
%
%                   CBVPred - [struct] contains the CBV predictions
%_______________________________________________________________
%   RETURN:                     
%                   Stats - [struct] results of the statistical analysis                
%_______________________________________________________________

teststart = 2;
increment = 2;

for a = 1:length(animals)
    
    % CD to location of files
    prevdir = cd([animals{a} filesep]);
    
    % Load data chunked around behavioral events
%     EvFile = ls(['*EVENTDATA_' (CBVType) '.mat']);
%     load(EvFile);
    EvFile = dir(['*EVENTDATA_' (CBVType) '.mat']);
    load(EvFile.name);
    
    % Process volitional whisking data
    Fs = EventData.(CBVType).whisk.Fs;
    StrtInd = (EventData.(CBVType).whisk.epoch.offset-1)*Fs;
    StpInd = StrtInd + 5*Fs;
    EventInds = StrtInd:StpInd;
    VWPred = CBVPred.(animals{a}).(NeurType).VW(:,EventInds);
    [~,FiltArray,~] = SelectBehavioralEvents(EventData.(CBVType),'VW');
    VWAct = EventData.(CBVType).whisk.NormData(FiltArray,:);
    VWAct = VWAct(teststart:increment:size(VWAct,1),EventInds);
    mVWAct = mean(VWAct,2)*ones(1,size(VWAct,2));
    mVWPred = mean(VWPred,2)*ones(1,size(VWPred,2));
    VWResid = (VWAct-mVWAct)-(VWPred-mVWPred);
    ReshapedVW = reshape(VWResid',1,numel(VWResid));
    Vars.VW(a) = var(ReshapedVW);
    
    Fs = EventData.(CBVType).stim.Fs;
    StrtInd = EventData.(CBVType).stim.epoch.offset*Fs;
    StpInd = StrtInd + 3*Fs;
    EventInds = StrtInd:StpInd;
    ContraPred = CBVPred.(animals{a}).(NeurType).Contra(:,EventInds);
    [~,FiltArray,~] = SelectBehavioralEvents(EventData.(CBVType),'Contra');
    ContraAct = EventData.(CBVType).stim.NormData(FiltArray,:);
    ContraAct = ContraAct(teststart:increment:size(ContraAct,1),EventInds);
    mContraAct = mean(ContraAct,2)*ones(1,size(ContraAct,2));
    mContraPred = mean(ContraPred,2)*ones(1,size(ContraPred,2));
    ContraResid = (ContraAct-mContraAct)-(ContraPred-mContraPred);
    ReshapedContra = reshape(ContraResid',1,numel(ContraResid));
    Vars.Contra(a) = var(ReshapedContra);   
    
%     RestFile = ls(['*RestData_' CBVType '.mat']);
%     load(RestFile);
    RestFile = dir(['*RESTDATA_' CBVType '.mat']);
    load(RestFile.name);
    RestPred = CBVPred.(animals{a}).(NeurType).Rest;
    [~,FiltArray,~] = SelectBehavioralEvents(RestData.(CBVType),'Rest');
    RestAct = RestData.(CBVType).NormData(FiltArray);
    RestAct = RestAct(teststart:increment:size(RestAct,1),:);
    resid = cell(1,length(RestAct));
    for c = 1:length(RestAct)
        resid{c} = (RestAct{c}-mean(RestAct{c}))-(RestPred{c}-mean(RestPred{c}));
    end
    JoinedRest = [resid{:}];
    Vars.Rest(a) = var(JoinedRest);
    clear RestData PredRestCBV
    cd(prevdir)
end

nContra = sqrt(Vars.Contra)./sqrt(Vars.Rest);
nVW = sqrt(Vars.VW)./sqrt(Vars.Rest);

figure;
bar([1,2],[mean(nContra),mean(nVW)],0.1,'EdgeColor','k',...
    'FaceColor',[1 1 1]);
hold on;
scatter([ones(size(nContra)), 2*ones(size(nVW))],...
    [nContra, nVW],'ko');
% ContraVar = 2.1; % Result of figure 1e
% VWVar = 1.6; % Result of figure 1e
% plot(0:3, ContraVar*ones(size(0:3)),'--','Color',[13 115 184]/255);
% plot(0:3, VWVar*ones(size(0:3)),'--','Color',[235, 176, 32]/255);
hold off;
xlim([0.9 2.1]);
ylim([0 2.5])
set(gca,'XTick',[1;2],'XTickLabel',{'Sensory Evoked';'Whisk'})

[h,p] = adtest(nContra);
[h,p] = adtest(nVW);

Bonferroni = 3;
if strcmp(NeurType,'Gam')
    [~,Stats.StimVRest.pval,~,stat] = ttest(nContra,1);
    Stats.StimVRest.pCorrected = Stats.StimVRest.pval*Bonferroni;
    Stats.StimVRest.tstat = stat.tstat;
    Stats.StimVRest.df = stat.df;
    [~,Stats.WhiskVRest.pval,~,stat] = ttest(nVW,1);
    Stats.WhiskVRest.Pcorrected = Stats.WhiskVRest.pval*Bonferroni;
    Stats.WhiskVRest.tstat = stat.tstat;
    Stats.WhiskVRest.df = stat.df;
    [~,Stats.StimVWhisk.pval,~,stat] = ttest(nContra,nVW);
    Stats.StimVWhisk.Pcorrected = Stats.StimVWhisk.pval*Bonferroni;
    Stats.StimVWhisk.tstat = stat.tstat;
    Stats.StimVWhisk.df = stat.df;
elseif strcmp(NeurType,'MUpower')
    [Stats.StimVRest.pval,~,stat] = signrank(nContra,ones(size(nContra)),...
        'method','approximate');
    Stats.StimVRest.pCorrected = Stats.StimVRest.pval*Bonferroni;
    Stats.StimVRest.Zval = stat.zval;
    [Stats.WhiskVRest.pval,~,stat] = signrank(nVW,ones(size(nVW)),...
        'method','approximate');
    Stats.WhiskVRest.pCorrected = Stats.WhiskVRest.pval*Bonferroni;
    Stats.WhiskVRest.Zval = stat.zval;
    [Stats.StimVWhisk.pval,~,stat] = signrank(nContra,nVW,...
        'method','approximate');
    Stats.StimVWhisk.pCorrected = Stats.StimVWhisk.pval*Bonferroni;
    Stats.StimVWhisk.Zval = stat.zval;
end