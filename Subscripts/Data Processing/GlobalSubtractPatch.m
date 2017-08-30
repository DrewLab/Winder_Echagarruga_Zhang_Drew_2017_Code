clear
animal = 'w327';
hem = 'LH';
GlobFile = ls('*RESTDATA_GlobalROI.mat');
load(GlobFile)
Global = RestData.GlobalROI.NormData;
CCFile = ls('*RESTDATA_CrossCorrROI.mat');
load(CCFile)
GlobalSubtract = RestData.CrossCorrROI;
for D = 1:length(RestData.CrossCorrROI.Data)
%     plot(RestData.CrossCorrROI.NormData{D}-1);
%     hold on; 
%     plot(Global{D}-1);
    GlobalSubtract.NormData{D} = RestData.CrossCorrROI.NormData{D}-...
        Global{D};
%     plot(GlobalSubtract.Data{D});
%     hold off;
%     pause;
end
clear RestData
RestData.GlobalSubtract = GlobalSubtract;
clear GlobalSubtract
save([animal '_' hem '_RESTDATA_GlobalSubtract.mat'],'RestData')

GlobFile = ls('*EVENTDATA_GlobalROI.mat');
load(GlobFile)
fnames = fieldnames(EventData.GlobalROI);
for f = 1:length(fnames)
    fname = fnames{f};
    Global.(fname) = EventData.GlobalROI.(fname).NormData;
end
CCFile = ls('*EVENTDATA_CrossCorrROI.mat');
load(CCFile)
GlobalSubtract = EventData.CrossCorrROI;
fnames = fieldnames(EventData.CrossCorrROI);
for f = 1:length(fnames)
    fname = fnames{f};
    GlobalSubtract.(fname).NormData = EventData.CrossCorrROI.(fname).NormData-...
        Global.(fname);
end
clear EventData
EventData.GlobalSubtract = GlobalSubtract;
clear GlobalSubtract
save([animal '_' hem '_EVENTDATA_GlobalSubtract.mat'],'EventData')