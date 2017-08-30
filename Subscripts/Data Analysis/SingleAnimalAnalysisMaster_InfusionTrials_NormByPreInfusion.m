function [] = SingleAnimalAnalysisMaster_InfusionTrials_NormByPreInfusion(animal,hem,InfusionTypes)
%   [] = SingleAnimalAnalysisMaster_InfusionTrials(animal,hem)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Compiles all the standard analysis for infusion data of a
%   single animal into a single script.
%   
%_______________________________________________________________
%   PARAMETERS:       
%                   animal - [string] animal ID
%
%                   hem - [string] hemisphere recorded
%
%                   InfusionTypes - [cell array of strings] designates type
%                   of infusions to analyze. The contents
%                   of the cell array must match the folder where the data
%                   is found.
%                               
%_______________________________________________________________
%   RETURN:                 
%                   Nothing returned. Output of the script are plots and
%                   saved files.
%                               
%_______________________________________________________________

%% CATEGORIZE THE DATA ACCORDING TO BEHAVIOR
% Load each file and use the binarized (and linked) whisking to obtain
% information about each whisk. Identify solenoid firing times and use
% tracked whisking and body movement to classify air puffs. Identify
% periods of no detected whisker movement or puffing as periods of rest.
% Save the behavioral categorization as part of the ProcData.mat file

for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    prevdir = cd([InfusionType filesep]);
    ProcDataFileNames = dir('*ProcData.mat');
    for f = 1:size(ProcDataFileNames,1)
        filename = ProcDataFileNames(f).name;
        CategorizeData_InfusionTrials(filename)
    end
    cd(prevdir)
end








%% IDENTIFY THE INFUSION START TIME
% Commented out since these times have already been entered.
% for IT = 1:length(InfusionTypes)
%     InfusionType = InfusionTypes{IT};
%     prevdir = cd([InfusionType '\Processed Data\']);
%     ThreshFile = ls('*Thresholds.mat');
%     if ~isempty(ThreshFile)
%         load(ThreshFile)
%         if ~isfield(Thresholds,'Infusion_Start_Times')
%             Thresholds.Infusion_Start_Times = [];
%         end
%     else
%         Thresholds = [];
%         Thresholds.Infusion_Start_Times = [];
%     end
%     
%     ProcDataFileNames = ls('*ProcData.mat');
%     [~,~,FileDate,~] = GetFileInfo(ProcDataFileNames);
%     DateCell = mat2cell(FileDate,ones(1,size(FileDate,1)),size(FileDate,2));
%     UniqueDates = unique(DateCell);
%     
%     for UD = 1:length(UniqueDates)
%         FileInds = strcmp(DateCell,UniqueDates{UD});
%         DayFiles = ProcDataFileNames(FileInds,:);
%         strdate = ConvertDate(UniqueDates{UD});
%         
%         % Log the infusion time if hasn't been recorded
%         if not(isfield(Thresholds.Infusion_Start_Times,strdate))
%             InfusionStartFile = uigetfile(['*' UniqueDates{UD} '*_rawdata.mat'],...
%                 ['Select the first file after infusion start on ' strdate ':']);
%             [~,~,~,InfusID] = GetFileInfo(InfusionStartFile);
%             Thresholds.Infusion_Start_Times.(strdate) = InfusID;
%             save(ThreshFile,'Thresholds');
%         end
%     end
%     save([animal '_' hem '_Thresholds.mat'],'Thresholds')
%     cd(prevdir)
% end








%% COMPILE ALL PERIODS OF REST
% Use categorization data added to each ProcData.mat file to extract the
% resting data from each file and save it in a separate structure. 

for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    prevdir = cd([InfusionType filesep]);
    ProcDataFileNames = dir('*ProcData.mat');
    ProcDataFileNames = {ProcDataFileNames(:).name}';
    
    ThreshFile = dir('*Thresholds.mat');
    if ~isempty(ThreshFile)
        load(ThreshFile.name)
    end
    
    % User select type of data to be analyzed
    if IT==1
        load(ProcDataFileNames{1})
        AllDataTypes = [fieldnames(ProcData.CBV); fieldnames(ProcData.Neuro); ...
            fieldnames(ProcData.Beh)];
        selections = listdlg('ListString',AllDataTypes,'SelectionMode','multiple',...
            'PromptString','Get resting data, select the data types: ');
        dataTypes = AllDataTypes(selections);
    end
    for dt = 1:length(dataTypes)
        dataType = dataTypes{dt};
        wraptext(['Overwriting RestData.mat structure  for ' dataType])
        RestData = struct;
        
        % Identify periods of rest from the ProcData files
        [RestData] = GetAllRestingData_InfusionTrials(ProcDataFileNames,RestData,dataType,Thresholds);
        save([animal '_' hem '_RESTDATA_' dataType '.mat'],'RestData');
    end
    cd(prevdir)
end




%% Calculate the daily baselines, as periods of rest longer than 14 seconds,
% for neural and hemodynamics

dataTypes = {};
for IT = 1:length(InfusionTypes)
    prevdir = cd([InfusionTypes{IT} filesep]);
    if IT==1
        ProcDataFileNames = dir('*ProcData.mat');
        load(ProcDataFileNames(1).name)
        AllDataTypes = [fieldnames(ProcData.CBV); fieldnames(ProcData.Neuro)];
        selections = listdlg('ListString',AllDataTypes,'SelectionMode','multiple',...
            'PromptString','Get resting data, select the data types: ');
        dataTypes = AllDataTypes(selections);
    end
    for dt = 1:length(dataTypes)
        dataType = dataTypes{dt};
        RestFile = dir(['*_RESTDATA_' dataType '.mat']);
        if isempty(RestFile)
            error(['No resting ' dataType ' found for ' animal ])
        else
            load(RestFile.name)
        end
        CalculateRestingBaseline_InfusionTrials(RestData,dataType,animal,hem);
    end
    cd(prevdir)
end

%% Normalize the data by pre-infusion rest from the same session

dataTypes = {};
for IT = 1:length(InfusionTypes)
    prevdir = cd([InfusionTypes{IT} filesep]);
    if IT==1
        ProcDataFileNames = dir('*ProcData.mat');
        load(ProcDataFileNames(1).name)
        AllDataTypes = [fieldnames(ProcData.CBV); fieldnames(ProcData.Neuro)];
        selections = listdlg('ListString',AllDataTypes,'SelectionMode','multiple',...
            'PromptString','Get resting data, select the data types: ');
        dataTypes = AllDataTypes(selections);
    end
    BaseFile = dir('*_Baselines.mat');
    if isempty(BaseFile)
        Baselines = [];
    else
        load(BaseFile.name)
    end
    for dt = 1:length(dataTypes)
        dataType = dataTypes{dt};
        RestFile = dir(['*_RESTDATA_' dataType '.mat']);
        if isempty(RestFile)
            error(['No resting ' dataType ' found for ' animal ])
        else
            load(RestFile.name)
        end
        RestData.(dataType).NormData = Infusion_NormalizeData(RestData.(dataType).Data,...
            RestData.(dataType).FileDate,Baselines,dataType);
        save(RestFile.name,'RestData')
    end
    cd(prevdir)
end
        
%% Calculate the resting aCSF NEURAL fluctuation amplitude 20 minutes after infusion start
prevdir = cd(['aCSF' filesep]);
ProcDataFileNames = dir('*ProcData.mat');
load(ProcDataFileNames(1).name)
AllDataTypes = fieldnames(ProcData.Neuro);
selections = listdlg('ListString',AllDataTypes,'SelectionMode','multiple',...
    'PromptString','Calculate Baselines, select the data types: ');
dataTypes = AllDataTypes(selections);

BaseFile = dir('*_Baselines.mat');
if isempty(BaseFile)
    Baselines = [];
else
    load(BaseFile.name)
end

for dt = 1:length(dataTypes)
    RestFile = dir(['*_RESTDATA_' dataTypes{dt} '.mat']);
    if isempty(RestFile)
        error(['No resting ' dataTypes{dt} ' found for ' animal ])
    else
        load(RestFile.name)
    end
    RestData = GetRestingaCSFAmplitude_NormByPreInfusion(dataTypes{dt},20,Baselines,RestData);
    display([num2str(length(RestData.(dataTypes{dt}).MeanRestingAmpFluctuation)) ' resting events used for aCSF amplitude...'])
    save([animal '_' hem '_RESTDATA_' dataTypes{dt} '.mat'],'RestData')
end
cd(prevdir)


%% Find the point when the resting neural activity drops below 0.2*aCSF resting MUA
dataTypes = [];
for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    if strcmp(InfusionType,'aCSF')
        continue;
    end
    display(['----' InfusionType '-----'])
    prevdir1 = cd([InfusionTypes{IT} filesep]);
    
    if isempty(dataTypes)
        ProcDataFileNames = dir('*ProcData.mat');
        load(ProcDataFileNames(1).name)
        AllDataTypes = fieldnames(ProcData.Neuro);
        selections = listdlg('ListString',AllDataTypes,'SelectionMode','multiple',...
            'PromptString','Select the neural data types for monitoring infusion: ');
        dataTypes = AllDataTypes(selections);
    end
    
    for dt = 1:length(dataTypes)
        dataType = dataTypes{dt};
        
        % Get the aCSF resting data
        prevdir = cd([prevdir1 filesep 'aCSF' filesep]);
        RestFile = dir(['*_RESTDATA_' dataType '.mat']);
        if isempty(RestFile)
            error('No Resting aCSF Amplitudes found...')
        else
            load(RestFile.name)
        end
        aCSFAmp = mean(RestData.(dataType).MeanRestingAmpFluctuation);
        cd(prevdir)

        RestFile = dir(['*_RESTDATA_' dataType '.mat']);
        if isempty(RestFile)
            error(['No Resting Amplitudes found for ' InfusionTypes{IT} '...'])
        else
            load(RestFile.name)
        end
        % Make sure that the aCSF resting amplitude can be found
        RestData = IdentifyReducedNeuroTrials_NormByPreInfusion(0.2,aCSFAmp,RestData,dataTypes{dt});
        save(RestFile.name,'RestData')
    end
    cd(prevdir1)
end


%% Calculate the resting aCSF CBV amplitude beginning 20 minutes after infusion
prevdir = cd(['aCSF' filesep]);
ProcDataFileNames = dir('*ProcData.mat');
load(ProcDataFileNames(1).name)
AllDataTypes = fieldnames(ProcData.CBV);
selections = listdlg('ListString',AllDataTypes,'SelectionMode','multiple',...
    'PromptString','Calculate Baselines, select the data types: ');
dataTypes = AllDataTypes(selections);

BaseFile = dir('*_Baselines.mat');
if isempty(BaseFile)
    Baselines = [];
else
    load(BaseFile.name)
end

for dt = 1:length(dataTypes)
    RestFile = dir(['*_RESTDATA_' dataTypes{dt} '.mat']);
    if isempty(RestFile)
        error(['No resting ' dataType{dt} ' found for ' animal ])
    else
        load(RestFile.name)
    end
    RestData = GetRestingaCSFAmplitude_NormByPreInfusion(dataTypes{dt},20,Baselines,RestData);
    save([animal '_' hem '_RESTDATA_' dataTypes{dt} '.mat'],'RestData')
end
cd(prevdir)




%% Compare the pharmacological infusion CBV to aCSF
dataTypes = [];
for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    if strcmp(InfusionType,'aCSF')
        continue;
    end
    display(['----' InfusionType '-----'])
    prevdir1 = cd([InfusionTypes{IT} filesep]);
    
    if isempty(dataTypes)
        ProcDataFileNames = dir('*ProcData.mat');
        load(ProcDataFileNames(1).name)
        AllDataTypes = [fieldnames(ProcData.CBV)];
        selections = listdlg('ListString',AllDataTypes,'SelectionMode','multiple',...
            'PromptString','Calculate Baselines, select the data types: ');
        dataTypes = AllDataTypes(selections);
    end
    
    for dt = 1:length(dataTypes)
        CBVdataType = dataTypes{dt};
        
        % Get the aCSF resting data
        prevdir = cd([prevdir1 filesep 'aCSF' filesep]);
        RestFile = dir(['*_RESTDATA_' CBVdataType '.mat']);
        if isempty(RestFile)
            error('No Resting aCSF Amplitudes found...')
        else
            load(RestFile.name)
        end
        aCSFCBV = RestData.(CBVdataType).MeanRestingAmpFluctuation;
        cd(prevdir)
        
        RestFile = dir(['*_RESTDATA_' CBVdataType '.mat']);
        if isempty(RestFile)
            error(['No Resting Amplitudes found for ' InfusionTypes{IT} '...'])
        else
            load(RestFile.name)
        end
        % Make sure that the aCSF resting amplitude can be found
        RestData = CompareRestingInfusionCBV_NormByPreInfusion(aCSFCBV,RestData,CBVdataType);
        save(RestFile.name,'RestData');
    end
    cd(prevdir1)
end



%% Compare the pharmacological infusion Neuro to aCSF
dataTypes = [];
for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    if strcmp(InfusionType,'aCSF')
        continue;
    end
    display(['----' InfusionType '-----'])
    prevdir1 = cd([InfusionTypes{IT} filesep]);
    
    if isempty(dataTypes)
        ProcDataFileNames = dir('*ProcData.mat');
        load(ProcDataFileNames(1).name)
        AllDataTypes = fieldnames(ProcData.Neuro);
        selections = listdlg('ListString',AllDataTypes,'SelectionMode','multiple',...
            'PromptString','Calculate Baselines, select the data types: ');
        dataTypes = AllDataTypes(selections);
    end
    
    for dt = 1:length(dataTypes)
        NeurodataType = dataTypes{dt};
        
        % Get the aCSF resting data
        prevdir = cd([prevdir1 filesep 'aCSF' filesep]);
        RestFile = dir(['*_RESTDATA_' NeurodataType '.mat']);
        if isempty(RestFile)
            error('No Resting aCSF Amplitudes found...')
        else
            load(RestFile.name)
        end
        aCSFNeuro = RestData.(NeurodataType).MeanRestingAmpFluctuation;
        cd(prevdir)
        
        RestFile = dir(['*_RESTDATA_' NeurodataType '.mat']);
        if isempty(RestFile)
            error(['No Resting Amplitudes found for ' InfusionTypes{IT} '...'])
        else
            load(RestFile.name)
        end
        % Make sure that the aCSF resting amplitude can be found
        RestData = CompareRestingInfusionNeuro_NormByPreInfusion(aCSFNeuro,RestData,NeurodataType);
        save(RestFile.name,'RestData');
    end
    cd(prevdir1)
end



    


%% COMPILE THE DATA SURROUNDING BEHAVIORAL EVENTS
% Setup
epoch.duration = 8; % seconds
epoch.offset = 2; % seconds

for IT = 1:length(InfusionTypes)
    InfusionType = InfusionTypes{IT};
    prevdir = cd([InfusionType filesep]);
    ProcDataFileNames = dir('*ProcData.mat');
    ThreshFile = dir('*Thresholds.mat');
    load(ThreshFile.name);
    
    if IT==1
    load(ProcDataFileNames(1).name)
    AllDataTypes = [fieldnames(ProcData.CBV); fieldnames(ProcData.Neuro); ...
        fieldnames(ProcData.Beh)];
    selections = listdlg('ListString',AllDataTypes,'SelectionMode','multiple',...
        'PromptString','Get event related data, select the data types: ');
    dataTypes = AllDataTypes(selections);
    end
    
    ProcDataFileNames = {ProcDataFileNames(:).name}';
    for dt = 1:length(dataTypes)
        dataType = dataTypes{dt};
        % Load previous resting structure if it exists
        EventData = struct;
        [EventData] = ExtractEventTriggeredData_InfusionTrials(ProcDataFileNames,...
            EventData,dataType,epoch,Thresholds.Infusion_Start_Times);
        
        % Save the structure
        save([animal '_' hem '_EVENTDATA_' dataType '.mat'],'EventData');
    end
    cd(prevdir)
end



%% Normalize the Behavioral Structure

dataTypes = {};
for IT = 1:length(InfusionTypes)
    prevdir = cd([InfusionTypes{IT} filesep]);
    if IT==1
        ProcDataFileNames = dir('*ProcData.mat');
        load(ProcDataFileNames(1).name)
        AllDataTypes = [fieldnames(ProcData.CBV); fieldnames(ProcData.Neuro)];
        selections = listdlg('ListString',AllDataTypes,'SelectionMode','multiple',...
            'PromptString','Get resting data, select the data types: ');
        dataTypes = AllDataTypes(selections);
    end
    BaseFile = dir('*_Baselines.mat');
    if isempty(BaseFile)
        Baselines = [];
    else
        load(BaseFile.name)
    end
    for dt = 1:length(dataTypes)
        dataType = dataTypes{dt};
        EventFile = dir(['*_EVENTDATA_' dataType '.mat']);
        if isempty(EventFile)
            error(['No event ' dataType ' found for ' animal ])
        else
            load(EventFile.name)
        end
        EventData.(dataType).stim.NormData = Infusion_NormalizeData(EventData.(dataType).stim.Data,...
            EventData.(dataType).stim.FileDate,Baselines,dataType);
        EventData.(dataType).whisk.NormData = Infusion_NormalizeData(EventData.(dataType).whisk.Data,...
            EventData.(dataType).whisk.FileDate,Baselines,dataType);
        save(EventFile.name,'EventData')
    end
    cd(prevdir)
end





% %% CALCULATE AND PLOT THE TRIGGERED AVERAGE FOR EACH TYPE OF BEHAVIOR
% % User select the data types to measure
% for IT = 1:length(InfusionTypes)
%     InfusionType = InfusionTypes{IT};
%     prevdir = cd([InfusionType filesep]);
%     ProcDataFileNames = dir('*ProcData.mat');
%     ThreshFile = dir('*Thresholds.mat');
%     if not(isempty(ThreshFile))
%         load(ThreshFile.name)
%     else
%         error('No Threshold file found...')
%     end
%     
%     % Get the minutes after infusion to start averaging events
%     if strcmp(InfusionType,'aCSF')
%         Threshfields = fieldnames(Thresholds.Infusion_Start_Times);
%         for Tf = 1:length(Threshfields)
%             AverageThresh.(Threshfields{Tf}) = 20; 
%         end
%     else
%         AverageThresh = [];
%         Threshfields = fieldnames(Thresholds.ReducedMUpower_Start);
%         for Tf = 1:length(Threshfields)
%             if isempty(Thresholds.ReducedMUpower_Start.(Threshfields{Tf}))
%                 continue;
%             end
%             AverageThresh.(Threshfields{Tf}) = ElapsedTime(Thresholds.ReducedMUpower_Start.(Threshfields{Tf}),...
%                 0,Thresholds.Infusion_Start_Times.(Threshfields{Tf}));
%         end
%     end
%     
%     if IT==1
%         load(ProcDataFileNames(1).name)
%         AllDataTypes = [fieldnames(ProcData.CBV); fieldnames(ProcData.Neuro); ...
%             fieldnames(ProcData.Beh)];
%         selections = listdlg('ListString',AllDataTypes,'SelectionMode','multiple',...
%             'PromptString','Plot event data, select the data types: ');
%         dataTypes = AllDataTypes(selections);
%     end
%     
%     % Plot the event triggered waveforms
%     behaviors = {'Contra','VW'};
%     for dt = 1:length(dataTypes)
%         dataType = dataTypes{dt};
%         EventDataFileName = dir(['*EVENTDATA_' dataType '.mat']);
%         load(EventDataFileName.name);
%         EventData = PlotEventTriggeredAve_Infusion_NormByPreInfusion(animal,EventData,dataType,...
%             behaviors,InfusionType,AverageThresh);
%         
%         % Save the structure
%         save([animal '_' hem '_EVENTDATA_' dataType '.mat'],'EventData');
%     end
%     close all
%     cd(prevdir)
% end

%%