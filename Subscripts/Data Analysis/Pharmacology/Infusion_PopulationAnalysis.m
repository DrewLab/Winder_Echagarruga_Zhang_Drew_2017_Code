function [] = Infusion_PopulationAnalysis(animals, InfusionTypes)
%   function [] =  ()
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: 
%   
%_______________________________________________________________
%   PARAMETERS:             
%                       
% Muscimol animals = {'CAN08','CAN10','CAN11','CAN12','CAN13','CAN15','CAN16','CAN18','CAN24'};
% 
%_______________________________________________________________
%   RETURN:                     
%                               
%_______________________________________________________________

%% Load and setup
cd('F:\Infusion Data\')

%% Gather all the CBV and neural fluctuation amplitudes for population comparison
for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    if strcmp(InfusionType,'aCSF')
        continue;
    end
    if exist(['G:\Infusion Data\Results\' InfusionType '\PharmStruct.mat'],'file')==2
        load(['G:\Infusion Data\Results\' InfusionType '\PharmStruct.mat'])
    elseif exist(['G:\Infusion Data\Results\' InfusionType '\PharmStruct.mat'],'file')==0
        PharmData = struct;
    end
    for a = 1:length(animals)
            
        animal = animals{a};
        if not(isdir(['G:\Infusion Data\' animal '\' InfusionType '\Processed Data\']))
            continue;
        end
        prevdir = cd(['G:\Infusion Data\' animal '\' InfusionType '\Processed Data\']);
        % Load the CBV data
        RestFile = ls(['*_RESTDATA_' CBVType '.mat']);
        if isempty(RestFile)
            error(['No rest data found for ' CBVType])
        else
            load(RestFile)
        end
        if isempty(RestData.(CBVType).AmpRatio)
            continue;
        end
        PharmData.(CBVType).(animal) =...
            RestData.(CBVType).AmpRatio_SD;
        PharmData.(CBVType).Means(a) = RestData.(CBVType).AmpRatio;
        
        for NT = 1:length(NeurTypes)
            NeurType = NeurTypes{NT};
            % Not all animals received all infusions, check if animal
            % received infusion
            if not(isdir(['G:\Infusion Data\' animal '\' InfusionType '\Processed Data\']))
                continue;
            end
           
            % Load the neural rest data
            RestFile = ls(['*_RESTDATA_' NeurType '.mat']);
            if isempty(RestFile)
                error(['No rest data found for ' NeurType])
            else
                load(RestFile)
            end
            PharmData.(NeurType).(animal) = RestData.(NeurType).AmpRatio_SD;
            PharmData.(NeurType).Means(a) = RestData.(NeurType).AmpRatio;
            
        end
        cd(prevdir)
    end
    CreateSaveLocation(['G:\Infusion Data\Results\' InfusionType '\'])
    save(['G:\Infusion Data\Results\' InfusionType '\PharmStruct.mat'],'PharmData')
end

%% Plot the CBV vs Neuro amplitude comparison
cd('G:\Infusion Data\')
close all;
for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    if strcmp(InfusionType,'aCSF')
        continue;
    end
    if exist(['Results\' InfusionType '\PharmStruct.mat'],'file')==2
        load(['Results\' InfusionType '\PharmStruct.mat'])
    elseif exist(['Results\' InfusionType '\PharmStruct.mat'],'file')==0
        error(['No PharmStruct exists for ' InfusionType])
    end
    for NT = 1:length(NeurTypes)
        NeurType = NeurTypes{NT};
        handle = figure(10*IT+NT+1); % Create figure number
        LegendFilt = true(1,length(animals));
        AnimalNeur = [];
        AnimalCBV = [];
        ii=1;
        for a = 1:length(animals)
            animal = animals{a};
            if not(isfield(PharmData.(NeurType),animal))
                LegendFilt(a) = 0;
                continue;
            end
            if isempty(PharmData.(NeurType).(animal))
                LegendFilt(a) = 0;
                continue;
            elseif isempty(PharmData.(CBVType).(animal))
                LegendFilt(a) = 0;
                continue;
            end
            handle = xy_error_ellipses(handle,mean(PharmData.(NeurType).(animal)),...
                mean(PharmData.(CBVType).(animal)),PharmData.(NeurType).(animal),...
                PharmData.(CBVType).(animal));
            AnimalNeur(ii) = mean(PharmData.(NeurType).(animal));
            AnimalCBV(ii) = mean(PharmData.(CBVType).(animal));
            ii = ii+1;
        end
        hold on; 
        scatter(1.55,mean(AnimalCBV),'ko');
        errorbar(1.55,mean(AnimalCBV),std(AnimalCBV));
        herrorbar(mean(AnimalNeur),1.55,std(AnimalNeur),std(AnimalNeur),'ko');
        ylim([0 1.6])
        xlim([0 1.6])
        title(['Pharmacology effects compared to aCSF: ' InfusionType])
        xlabel([InfusionType ' ' NeurType '/' InfusionType ' aCSF'])
        ylabel([InfusionType ' CBV/' InfusionType ' aCSF'])
        legend(animals(LegendFilt));
        axis square;
        pause;
        SaveFig2Disk(handle,['Results\' InfusionType ...
            '\Pharmacology effects on REST compared to aCSF_' NeurType '.fig'])
    end
end


%% Compile the population average CBV and neural activity
dataTypes = {'RadiusROI','MUpower','Gam'};
BehTypes = {'Contra','VW'};
AveStruct = [];
for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    InfusionField = strrep(InfusionType, ' + ','_');
    AveStruct.(InfusionField) = [];
    for dT = 1:length(dataTypes)
        dataType = dataTypes{dT};
        AveStruct.(InfusionField).(dataType) = [];
        AnimalInd = 1;
        for a = 1:length(animals)
            animal = animals{a};
            FileLocation = ['G:\Infusion Data\' animal '\' InfusionType '\Processed Data\' ...
                    animal '_LH_EVENTDATA_' dataType '.mat'];
            if exist(FileLocation,'file')==2
                load(FileLocation)
            elseif exist(FileLocation,'file')==0
                display(['No EVENTDATA Structure exists for ' animal '; ' ...
                    dataType '; ' InfusionType]);
                continue
            end
            
            for BT = 1:length(BehTypes)
                BehType = BehTypes{BT};
                if strcmp(BehType,'Contra')
                    fdname = 'stim';
                elseif strcmp(BehType,'VW')
                    fdname = 'whisk';
                end
                if AnimalInd==1
                    AveStruct.(InfusionField).(dataType).(BehType).PostAve = ...
                        NaN*ones(length(animals),...
                        length(EventData.(dataType).(fdname).Averages.TimeVec));
                    AveStruct.(InfusionField).(dataType).(BehType).PreAve = ...
                        NaN*ones(length(animals),...
                        length(EventData.(dataType).(fdname).Averages.TimeVec));
                    AveStruct.(InfusionField).(dataType).(BehType).AnimalID = ...
                        cell(1,length(animals));
                end
                if isfield(EventData.(dataType).(fdname),'Averages')
                    AveStruct.(InfusionField).(dataType).(BehType).PostAve(AnimalInd,:)...
                        = EventData.(dataType).(fdname).Averages.AllSessions.Post;
                    AveStruct.(InfusionField).(dataType).(BehType).PreAve(AnimalInd,:)...
                        = EventData.(dataType).(fdname).Averages.AllSessions.Pre;
                    AveStruct.(InfusionField).(dataType).(BehType).AnimalID{AnimalInd}...
                        = animal;
                    AveStruct.(InfusionField).(dataType).(BehType).timevec = ...
                        EventData.(dataType).(fdname).Averages.TimeVec;
                else
                    continue;
                end
            end
            AnimalInd = min(AnimalInd+1,length(animals));
        end
    end
end

for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    InfusionField = strrep(InfusionType, ' + ','_');
    for dT = 1:length(dataTypes)
        dataType = dataTypes{dT};
        for BT = 1:length(BehTypes)
            BehType = BehTypes{BT};
            NaNFilt = isnan(sum(AveStruct.(InfusionField).(dataType).(BehType).PreAve,2));
            AveStruct.(InfusionField).(dataType).(BehType).PostAve(NaNFilt,:) = [];
            AveStruct.(InfusionField).(dataType).(BehType).PreAve(NaNFilt,:) = [];
            AveStruct.(InfusionField).(dataType).(BehType).AnimalID = ...
                AveStruct.(InfusionField).(dataType).(BehType).AnimalID(not(NaNFilt));
        end
    end
end
save('Results\Averages\AveStruct.mat','AveStruct')

%% Plot the averages

cd('G:\Infusion Data\')
load('Results\Averages\AveStruct.mat')
close all;
BehaviorTypes = {'Contra'};
for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    InfusionType = strrep(InfusionType, ' + ','_');
    if strcmp(InfusionType,'aCSF')
        continue;
    end
    for NT = 1:length(NeurTypes)
        NeurType = NeurTypes{NT};
        for BT = 1:length(BehaviorTypes)
            handle = figure(10*IT+NT+BT+1); % Create figure number
            BehType = BehaviorTypes{BT};
            
            % Set the time periods of interest
            NeurInds = AveStruct.(InfusionType).(NeurType).(BehType).timevec>-0.1 & ...
                AveStruct.(InfusionType).(NeurType).(BehType).timevec<0.5;
            CBVInds = AveStruct.(InfusionType).RadiusROI.(BehType).timevec>0 & ...
                AveStruct.(InfusionType).RadiusROI.(BehType).timevec<3;
            PreInds = AveStruct.(InfusionType).RadiusROI.(BehType).timevec<0;
            
            % Normalize everything to the PreStim average of PreInfusion
            % Puffs
            SumPharmNeur = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
            SumaCSFNeur = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
            SumaCSFCBV = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
            SumPharmCBV = zeros(1,size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1));
            for a = 1:size(AveStruct.(InfusionType).(NeurType).(BehType).PostAve,1)
                animal = AveStruct.(InfusionType).(NeurType).(BehType).AnimalID{a};
                
                % Normalize the post infusion Pharm neural activity to
                % preinfusion amplitudes
                PrePharmNeur = mean(AveStruct.(InfusionType).(NeurType).(BehType).PreAve(a,PreInds));
                NormPharmPreNeur = AveStruct.(InfusionType).(NeurType).(BehType).PreAve(a,:)/PrePharmNeur;
                NormPharmPostNeur = AveStruct.(InfusionType).(NeurType).(BehType).PostAve(a,:)/PrePharmNeur;
                SumPharmPreNeur = sum(NormPharmPreNeur(NeurInds));
                SumPharmPostNeur = sum(NormPharmPostNeur(NeurInds));
                SumPharmNeur(a) = SumPharmPostNeur/SumPharmPreNeur;
                
                % Normalize the post infusion aCSF neural activity to
                % preinfusion amplitudes
                AnimalInd = strcmp(AveStruct.aCSF.(NeurType).(BehType).AnimalID,animal);
                PreaCSFNeur = mean(AveStruct.aCSF.(NeurType).(BehType).PreAve(AnimalInd,PreInds));
                NormaCSFPreNeur = AveStruct.aCSF.(NeurType).(BehType).PreAve(AnimalInd,:)/PreaCSFNeur;
                NormaCSFPostNeur = AveStruct.aCSF.(NeurType).(BehType).PostAve(AnimalInd,:)/PreaCSFNeur;
                SumaCSFPreNeur = sum(NormaCSFPreNeur(NeurInds));
                SumaCSFPostNeur = sum(NormaCSFPostNeur(NeurInds));
                SumaCSFNeur(a) = SumaCSFPostNeur/SumaCSFPreNeur;
                
                % Normalize the post infusion Pharm CBV to preinfusion
                % amplitudes
                PrePharmCBV = mean(AveStruct.(InfusionType).(CBVType).(BehType).PreAve(a,PreInds));
                NormPharmPreCBV = AveStruct.(InfusionType).(CBVType).(BehType).PreAve(a,:)/PrePharmCBV-1;
                NormPharmPostCBV = AveStruct.(InfusionType).(CBVType).(BehType).PostAve(a,:)/PrePharmCBV-1;
                SumPharmPreCBV = sum(NormPharmPreCBV(CBVInds));
                SumPharmPostCBV = sum(NormPharmPostCBV(CBVInds)...
                    -mean(NormPharmPostCBV(PreInds)));
                SumPharmCBV(a) = SumPharmPostCBV/SumPharmPreCBV;
                
                PreaCSFCBV = mean(AveStruct.aCSF.(CBVType).(BehType).PreAve(a,PreInds));
                NormaCSFPreCBV = AveStruct.aCSF.(CBVType).(BehType).PreAve(a,:)/PreaCSFCBV-1;
                NormaCSFPostCBV = AveStruct.aCSF.(CBVType).(BehType).PostAve(a,:)/PreaCSFCBV-1;
                SumaCSFPreCBV = sum(NormaCSFPreCBV(CBVInds));
                SumaCSFPostCBV = sum(NormaCSFPostCBV(CBVInds)...
                    -mean(NormaCSFPostCBV(PreInds)));
                SumaCSFCBV(a) = SumaCSFPostCBV/SumaCSFPreCBV;
                
            end
            NormNeur = SumPharmNeur./SumaCSFNeur;
            NormCBV = SumPharmCBV./SumaCSFCBV;
            scatter(NormNeur,NormCBV);
            ylim([0 2])
            xlim([0 2])
            title(['Pharmacology effects on ' BehType ' compared to aCSF: ' strrep(InfusionType,'_',' ')])
            xlabel([NeurType ' ' strrep(InfusionType,'_',' ') '/' NeurType ' aCSF'])
            ylabel(['CBV ' strrep(InfusionType,'_',' ') '/CBV aCSF']);
            axis square;
            hold on;
            scatter(1.95,mean(NormCBV),'ko');
            errorbar(1.95,mean(NormCBV),std(NormCBV));
            herrorbar(mean(NormNeur),1.95,std(NormNeur),std(NormNeur),'ko');
            hold off;
            pause;
            SaveFig2Disk(handle,['G:\Infusion Data\Results\' InfusionType ...
                '\Pharmacology effects on ' BehType ' compared to aCSF_' NeurType '.fig'])
            MaxStruct.(InfusionType).(NeurType).(BehType).NormNeur = NormNeur;
            MaxStruct.(InfusionType).(CBVType).(BehType).NormCBV = NormCBV;
        end
    end
end
save('G:\Infusion Data\Results\EventTriggeredMaxComparison','MaxStruct')
