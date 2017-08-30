function [Stats] = CompareInterTrialVariance(animals,CBVType)
%   function [Stats] = CompareInterTrialVariance(animals,CBVType)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Compares the variance of the behavior-evoked CBV to the
%   CBV variance during rest.
%_______________________________________________________________
%   PARAMETERS:             
%                   animals - [cell array] contains animal IDs
%
%                   CBVType - [string] fieldname corresponding to the CBV
%                   ROI
%_______________________________________________________________
%   RETURN:                     
%                   Stats - [struct] results of the statistical analysis                  
%_______________________________________________________________
EventDur = 4;
for a = 1:length(animals)
    ProcDir = [animals{a} filesep];
    prevdir = cd(ProcDir);
%     EventFile = ls(['*EventData_' CBVType '.mat']);
%     RestFile = ls(['*RestData_' CBVType '.mat']);
%     load(EventFile)
%     load(RestFile)
    EventFile = dir(['*EVENTDATA_' CBVType '.mat']);
    RestFile = dir(['*RESTDATA_' CBVType '.mat']);
    load(EventFile.name)
    load(RestFile.name)
    [~,VWfilt] = SelectBehavioralEvents(EventData.(CBVType),'VW');
    [~,Contrafilt] = SelectBehavioralEvents(EventData.(CBVType),'Contra');
    [~,Restfilt] = SelectBehavioralEvents(RestData.(CBVType),'Rest');
    
    VW.Data = EventData.(CBVType).whisk.NormData(VWfilt,:);
    Contra.Data = EventData.(CBVType).stim.NormData(Contrafilt,:);
    Rest.Data = RestData.(CBVType).NormData(Restfilt);
    strt = (EventData.(CBVType).whisk.epoch.offset)*...
        EventData.(CBVType).whisk.Fs;
    stp = strt + EventDur*EventData.(CBVType).whisk.Fs;
    VWData = VW.Data(:,strt:stp);
    mVW = mean(VWData,2)*ones(1,size(VWData,2));
    Vars.VW(a) = mean(var(VWData-mVW));
    strt = (EventData.(CBVType).stim.epoch.offset)*...
        EventData.(CBVType).stim.Fs;
    stp = strt + EventDur*EventData.(CBVType).stim.Fs;
    ContraData = Contra.Data(:,strt:stp);
    mContra = mean(ContraData,2)*ones(1,size(ContraData,2));
    Vars.Contra(a) = mean(var(ContraData-mContra));
    RestVars = zeros(size(Rest.Data));
    meanSubRest = cell(size(Rest.Data));
    for c = 1:length(Rest.Data)
        meanSubRest{c} = Rest.Data{c}-mean(Rest.Data{c});
        RestVars(c) = var(meanSubRest{c});
    end
    Vars.Rest(a) = var([meanSubRest{:}]);
    cd(prevdir)
end
nContra = sqrt(Vars.Contra)./sqrt(Vars.Rest);
nVW = sqrt(Vars.VW)./sqrt(Vars.Rest);

[~,Stats.StimVRest.pval,~,stat] = ttest(nContra,1);
Stats.StimVRest.tstat = stat.tstat;
Stats.StimVRest.df = stat.df;
[~,Stats.WhiskVRest.pval,~,stat] = ttest(nVW,1);
Stats.WhiskVRest.tstat = stat.tstat;
Stats.WhiskVRest.df = stat.df;
[~,Stats.StimVWhisk.pval,~,stat] = ttest(nContra,nVW);
Stats.StimVWhisk.tstat = stat.tstat;
Stats.StimVWhisk.df = stat.df;

figure;
bar([1,2],[mean(nContra),mean(nVW)],0.1,'EdgeColor','k',...
    'FaceColor',[1 1 1]);
hold on; 
scatter([ones(size(nContra)), 2*ones(size(nVW))],...
    [nContra, nVW],'ko');

hold off;
ylabel('Behavioral CBV Variance/Resting Variance');
xlim([0.9 2.1]);
ylim([0 2.5])
set(gca,'XTick',[1;2],'XTickLabel',{'Sensory Evoked';'Whisk'})