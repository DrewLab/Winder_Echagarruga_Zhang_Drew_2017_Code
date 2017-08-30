function [] = SingleAnimalAnalysisMaster(animal,hem)
%   function [] = SingleAnimalAnalysisMaster_atw2(animal,hem)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Compiles all the standard analysis for a single animal
%   into a single script.
%   
%_______________________________________________________________
%   PARAMETERS:       
%                   animal - [string] animal ID
%
%                   hem - [string] hemisphere recorded
%                               
%_______________________________________________________________
%   RETURN:                 
%                   Nothing returned. Output of the script are plots and
%                   saved files.
%                               
%_______________________________________________________________

%% DECLARE VARIABLES

NeurTypes = {'Gam','MUpower'};
CBVType = 'CrossCorrROI';

%% CALCULATE AND PLOT THE TRIGGERED AVERAGE FOR EACH TYPE OF BEHAVIOR

% Setup Analysis
dataTypes = [NeurTypes CBVType];
behaviors = {'Contra','Str','VW','Control','Ipsi'};
handles = struct();


for dt = 1:length(dataTypes)
    dataType = dataTypes{dt};
    EventDataFileName = ls(['*EVENTDATA_' dataType '.mat']);
    load(EventDataFileName);
    handles = PlotEventTriggeredAve(animal,EventData,dataType,behaviors,handles);
end
% 
% % close all



%% CALCULATE ANALYTIC HRF

% Setup Analysis
behaviors = {'Contra','Str','VW','Rest'};
% behaviors = {'Contra'};
% behaviors = {'VW'};
HRFLims = [0 5]; % Truncate the HRF between t=0,5 seconds

for nt = 1:length(NeurTypes)
    NeurType = NeurTypes{nt};
    figure;
    for b = 1:length(behaviors)
        Beh = behaviors{b};
        [HRFs.(NeurType).(Beh)] = ...
            CalculateHRF_Deconvolution(NeurType,CBVType, Beh);
        figure(nt);
        plot(HRFs.(NeurType).(Beh).timevec,...
            HRFs.(NeurType).(Beh).HRF); 
        hold on; 
        figure(nt+10);
        plot(HRFs.(NeurType).(Beh).GammaTime,...
            HRFs.(NeurType).(Beh).GammaHRF);
        hold on;
    end
    figure(nt);
    hold off;
    xlabel('HRF time (s)');
    ylabel('HRF Amplitude (A.U.)')
    title([animal ' HRF: ' NeurType])
    LegendBeh = behaviors;
    for b = 1:length(behaviors)
        Beh = behaviors{b};
        HRFs.(NeurType).(Beh).Lims = HRFLims;
        LegendBeh{b} = [behaviors{b} ' n=' ...
            num2str(HRFs.(NeurType).(Beh).num_calc_events)];
    end
    legend(LegendBeh)
end
save([animal '_' hem '_HRFs.mat'],'HRFs');

% close all

%% FIT ANALYTIC HRF WITH A GAMMA FUNCTION - Obsolete

% % Setup
% % behaviors = {'Contra','Str','VW','Rest'};
% behaviors = {'Contra'};
% HRFfile = ls('*HRFs.mat');
% if isempty(HRFfile)==0;
%     load(HRFfile);
% end
% 
% for nt = 1:length(NeurTypes)
%     NeurType = NeurTypes{nt};
%     for b = 1:length(behaviors)
%         Beh = behaviors{b};
%         IndHRF = HRFs.(NeurType).(Beh).Analytic;
%         [GammaHRF,gammacoef,t,HRF] = FitGammaToHRF(IndHRF);
%         plot(t,HRF,'k','Linewidth',3);
%         hold on;
%         plot(t,GammaHRF,'c--');
%         hold off;
%         ylim([-30e-4 3e-4])
%         R2 = CalculateRsquared(GammaHRF,HRF);
%         title([': R2 = ' num2str(R2)])
%         pause(0.1);
% %         SaveFig2Disk(gcf,[DriveName Group '\' animal ...
% %             '\Results\HRFs\' animal '_' dataType '_' Beh '_GammaHRFFit.fig']);
%         HRFs.(NeurType).(Beh).Gamma.HRF = GammaHRF;
%         HRFs.(NeurType).(Beh).Gamma.coef = gammacoef;
%         HRFs.(NeurType).(Beh).Gamma.timevec = t;
%         HRFs.(NeurType).(Beh).Gamma.Fs = HRFs.(NeurType).(Beh).Analytic.Fs;
%         HRFs.(NeurType).(Beh).Gamma.num_calc_events = ...
%             HRFs.(NeurType).(Beh).Analytic.num_calc_events;
%         HRFs.(NeurType).(Beh).Gamma.Event_Inds = ...
%             HRFs.(NeurType).(Beh).Analytic.Event_Inds;
%         HRFs.(NeurType).(Beh).Gamma.calc_date = date;
%         HRFs.(NeurType).(Beh).Gamma.HRFParams.offset = -1*t(1);
%         HRFs.(NeurType).(Beh).Gamma.HRFParams.dur = t(end)-t(1);
%         HRFs.(NeurType).(Beh).Gamma.Lims = HRFs.(NeurType).(Beh).Analytic.Lims;
%         HRFs.(NeurType).(Beh).Gamma.R2 = R2;
%         save(HRFfile,'HRFs')
%     end
% end


%% Evaluate the goodness of fit of the HRF
load([animal '_' hem '_HRFs.mat'])
behaviors = {'Contra'};
GOF = struct;
for nt = 1:length(NeurTypes)
    for b = 1:length(behaviors)
        Beh = behaviors{b};
        HRF = HRFs.(NeurTypes{nt}).(Beh).Gamma.HRF;
        GOF = CalculateHRFGoodnessOfFit_NEW(NeurTypes{nt},CBVType,Beh,HRF,GOF);
    end
end
save([animal '_' hem '_GOF.mat'],'GOF')









%% USE HRF TO PREDICT CBV

% Setup
behaviors = {'Contra','Str','VW','Rest'};

HRF_file = ls('*HRFs.mat');
if isempty(HRF_file)
    error('Cannot find HRF data...move HRFs.mat file to current directory, change directory, or run CompileHRFs.m')
end
load(HRF_file);

Base_file = ls('*Baselines.mat');
if isempty(Base_file)
    error(['Cannot find baseline data for the selections...move Baselines.mat'...
        'file to current directory, change directory, or run CalculateDailyBaselines.m']);
end
load(Base_file);

ProcDataFileNames = ls('*ProcData.mat');

for FileIter = 1:length(ProcDataFileNames)
    filename = ProcDataFileNames(FileIter,:);
    display(['CalculatePredictedCBV.m: Loading ' filename]);
    load(filename)
    [animal,hem,datename,FileID] = GetFileInfo(filename);
    strdate = ConvertDate(datename);
    
    for nt = 1:length(NeurTypes)
        CBVPred = [];
        NeurType = NeurTypes{nt};
        if exist([animal '_' hem '_' FileID '_CBVPRED_' NeurType '.mat'],'file')==2
            display(['CalculatePredictedCBV.m: Loading ' animal '_' hem '_'...
                FileID '_CBVPRED_' NeurType '.mat.']);
            load([animal '_' hem '_' FileID '_CBVPRED_' NeurType '.mat']);
        end
        CurrentBaseline.(NeurType) = Baselines.(NeurType).(strdate);
        CurrentBaseline.CrossCorrROI = Baselines.CrossCorrROI.(strdate);
        for b = 1:length(behaviors)
            Beh = behaviors{b};
            HRFType = [NeurType Beh];
            % Use Gamma-based HRF to predict
            CurrentHRF = HRFs.(NeurType).(Beh).GammaHRF;
            [CBVPred.(HRFType),CBVPred.Fs.([HRFType '_fs'])] = ...
                CalculatePredictedCBV(ProcData,NeurType,CBVType,CurrentHRF,CurrentBaseline);
        end
        CBVPred.Flags = ProcData.Flags;
        CBVPred.TrialDur = ProcData.TrialDur;
        CurrentBaseline = rmfield(CurrentBaseline,NeurType);
        % Save file in structure format for easy loading
        display(['Saving CBVPRED_' NeurType ...
            ProcDataFileNames(FileIter,1:25) 'HRF: ' NeurType ', ' Beh]);
        save([animal '_' hem '_' FileID '_CBVPRED_' NeurType '.mat'], 'CBVPred');
    end
end








%% CATEGORIZE PREDICTED CBV ACCORDING TO BEHAVIOR
% 
% % Setup
% behaviors = {'Contra','Str','VW','Rest'};
% epoch.duration = 8; % seconds
% epoch.offset = 2; % seconds
% 
% for nt = 1:length(NeurTypes)
%     NeurType = NeurTypes{nt};
%     CBVPredFileNames = ls(['*CBVPRED_' NeurType '.mat']);
%     PredEventCBV = struct;
%     PredRestCBV = struct;
%     for b = 1:length(behaviors)
%         %------------------------------------------------------------------
%         % EVENT RELATED BEHAVIORS
%         %------------------------------------------------------------------
%         % Load and setup
%         Beh = behaviors{b};
%         HRFType = [NeurType Beh];
%         
%         % Separate predictions according to event types
%         [PredEventCBV] = ExtractEventTriggeredData(CBVPredFileNames,...
%             PredEventCBV,HRFType,epoch);
%                
%         %------------------------------------------------------------------
%         % RESTING BEHAVIORS
%         %------------------------------------------------------------------
%         [PredRestCBV] = GetAllRestingData(CBVPredFileNames,PredRestCBV,HRFType);
%     end
%     
%     % Save the Event Structure
%     save([animal '_' hem '_EVENTPRED_' NeurType '.mat'],'PredEventCBV')
%     save([animal '_' hem '_RESTPRED_' NeurType '.mat'],'PredRestCBV')
% end
% 
% 
% 
% 
% 
% 
% 
% 
% %% CALCULATE R2 BETWEEN MEASURED AND PREDICTED CBV
% 
% % Setup
% HRFfile = ls('*_HRFs.mat');
% load(HRFfile)
% 
% EventFile = ls(['*EventData_' CBVType '.mat']);
% load(EventFile)
% RestFile = ls(['*RestData_' CBVType '.mat']);
% load(RestFile);
% 
% behaviors = {'Contra','Str','VW','Rest'};
% 
% 
% % Run
% for nt = 1:length(NeurTypes)
%     NeurType = NeurTypes{nt};
%     PredEventFileName = ls(['*EVENTPRED_' NeurType '.mat']);
%     load(PredEventFileName);
%     PredRestFileName = ls(['*RESTPRED_' NeurType '.mat']);
%     load(PredRestFileName);
%     
%     for b = 1:length(behaviors)
%         Beh = behaviors{b};
%         HRFType = [NeurType Beh];
%         
%         %------------------------------------------------------------------
%         % EVENT DRIVEN BEHAVIOR
%         %------------------------------------------------------------------
%         
%         % Calculate the R-squared between HRF-predicted CBV and actual CBV
%         [PredEventCBV] = CalculateEventTriggeredRsquared(PredEventCBV,...
%             EventData,HRFType,CBVType); 
%     
% %         % Calculate the R-squared using the corrected R-squared
% %         [PredEventCBV] = CalculateEventTriggeredRsquared(PredEventCBV,...
% %             EventData,HRFType,NonLin.(dataType).(Beh)); !!!
%         
%         %------------------------------------------------------------------
%         % RESTING BEHAVIOR
%         %------------------------------------------------------------------
%         [PredRestCBV] = CalculateRestingRsquared(PredRestCBV,RestData.(CBVType),...
%             HRFType);
%     end
%     save([animal '_' hem '_EVENTPRED_' NeurType '.mat'],'PredEventCBV')
%     save([animal '_' hem '_RESTPRED_' NeurType '.mat'],'PredRestCBV')
% end
% display('CALCULATE R2 BETWEEN MEASURED AND PREDICTED CBV: Done')
% 
% 
% 
% 
% 
% 
% 
% 
% %% CALCULATE CORRELATION COEFFICIENT BETWEEN MEASURED AND PREDICTED CBV
% 
% % Load and Setup
% EventDataFileName = ls(['*EVENTDATA_' CBVType '.mat']);
% load(EventDataFileName);
% 
% RestDataFileName = ls(['*RESTDATA_' CBVType '.mat']);
% load(RestDataFileName);
% 
% behaviors = {'Contra','Str','VW','Rest'};
% 
% % Loop over each HRF
% for nt = 1:length(NeurTypes)
%     NeurType = NeurTypes{nt};
%     PredEventFileName = ls(['*EVENTPRED_' NeurType '.mat']);
%     load(PredEventFileName);
%     PredRestFileName = ls(['*RESTPRED_' NeurType '.mat']);
%     load(PredRestFileName);
%     for b = 1:length(behaviors)
%         Beh = behaviors{b};
%         HRFType = [NeurType Beh];
%         
%         %------------------------------------------------------------------
%         % EVENT RELATED BEHAVIOR
%         %------------------------------------------------------------------
%         [PredEventCBV] = CalculateEventTriggeredCC(PredEventCBV,EventData,...
%             HRFType,CBVType);
% 
%         
%         %------------------------------------------------------------------
%         % RESTING BEHAVIOR
%         %------------------------------------------------------------------
%         [PredRestCBV] = CalculateRestingCC(PredRestCBV,RestData,HRFType,CBVType);
%     end
%     save([animal '_' hem '_EVENTPRED_' NeurType '.mat'],'PredEventCBV')
%     save([animal '_' hem '_RESTPRED_' NeurType '.mat'],'PredRestCBV')
% end
% display('CALCULATE R BETWEEN MEASURED AND PREDICTED CBV: Done')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% COMPILE THE GOODNESS OF FIT (GOF)
% 
% behaviors = {'Contra','Str','VW','Rest'};
% 
% EventDataFileName = ls(['*EVENTDATA_' CBVType '.mat']);
% load(EventDataFileName);
% RestDataFileName = ls(['*EVENTDATA_' CBVType '.mat']);
% load(RestDataFileName);
% 
% GOF = [];
% 
% for nt = 1:length(NeurTypes)
%     NeurType = NeurTypes{nt};
%     PredEventFileName = ls(['*EVENTPRED_' NeurType '.mat']);
%     load(PredEventFileName);
%     PredRestFileName = ls(['*RESTPRED_' NeurType '.mat']);
%     load(PredRestFileName);
%     
%     % Create HRF names
%     HRFTypes = cell(1,length(behaviors));
%     i=1;
%     for b = 1:length(behaviors)
%         Beh = behaviors{b};
%         HRFTypes{i} = [NeurType Beh];
%         i=i+1;
%     end
%     
%     for b = 1:length(behaviors)
%         Beh = behaviors{b};
%         for HT = 1:length(HRFTypes)
%             HRFType = HRFTypes{HT};
%             if strcmp(Beh,'Rest')
%                 [GOF] = CompileIndividualBehavioralRsquared(PredRestCBV,...
%                     GOF,HRFType,Beh);
%                 [GOF] = CompileIndividualBehavioralCC(PredRestCBV,...
%                     GOF,HRFType,Beh);
%             else
%                 [GOF] = CompileIndividualBehavioralRsquared(PredEventCBV,...
%                     GOF,HRFType,Beh);
%                 [GOF,AveCBV,PredCBV] = BehaviorAverageRsquared(GOF,PredEventCBV,...
%                     EventData,HRFType,Beh,CBVType);
%                 [GOF] = CompileIndividualBehavioralCC(PredEventCBV,...
%                     GOF,HRFType,Beh);
%                 [GOF,AveCBV,PredCBV] = BehaviorAverageCC(GOF,PredEventCBV,...
%                     EventData,HRFType,Beh,CBVType);
%             end
%             Fs = EventData.(CBVType).stim.Fs; 
%             timevec = (1:length(AveCBV))/Fs-EventData.(CBVType).stim.epoch.offset;
% %             
%             if not(strcmp(Beh,'Rest'))
% %                 figure(1); plot(timevec,[(AveCBV-mean(AveCBV))' (PredCBV-mean(PredCBV))']); axis tight;
% %                 xlabel('Time (s)'); ylabel('delR/R');
% % %                 title([HRFType ' Prediction of Average ' Beh ' CBV R^2=' ...
% % %                     num2str(GOF.(Beh).(HRFType).AveR2)]);
% %                 title([HRFType ' Prediction of Average ' Beh ' CBV R^2=' ...
% %                     num2str(GOF.(Beh).(HRFType).AveR2) '; R=' ...
% %                     num2str(GOF.(Beh).(HRFType).AveR)]);
% %                 legend({'Actual','Predicted'})
% % %                 SaveFig2Disk(gcf,['..\Results\GOF\' ...
% % %                     HRFType ' Prediction of Average ' Beh ' CBV.fig']);
%             end
% % %             
% %             figure(2);
% %             hist(GOF.(Beh).(HRFType).IndR2,...
% %                 max(10,round(numel(GOF.(Beh).(HRFType).IndR2)/4)));
% %             xlabel('R^2'); ylabel('Count');
% %             title([HRFType ' Prediction of ' Beh ' CBV Trials: Median R^2=' ...
% %                 num2str(median(GOF.(Beh).(HRFType).IndR2))]);
% % %             SaveFig2Disk(gcf,['..\Results\GOF\' ...
% % %                 HRFType ' Prediction of individual ' Beh ' CBV.fig']);
%             
% %             figure;
% %             hist(GOF.(Beh).(HRFType).IndR,...
% %                 max(10,round(numel(GOF.(Beh).(HRFType).IndR)/4)));
% %             xlabel('R'); ylabel('Count');
% %             title([HRFType ' Prediction of ' Beh ' CBV Trials: Median R=' ...
% %                 num2str(median(GOF.(Beh).(HRFType).IndR))]);
% %             SaveFig2Disk(gcf,['..\Results\GOF\' ...
% %                 HRFType ' CORRCOEF_Prediction of individual ' Beh ' CBV.fig'])
% % pause;
%         end
%         close all;
%     end
% end
% save([animal '_' hem '_GOF.mat'],'GOF')
% display('COMPILE THE GOODNESS OF FIT (GOF): Done')
