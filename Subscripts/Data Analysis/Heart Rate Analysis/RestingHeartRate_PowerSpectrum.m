function [] = RestingHeartRate_PowerSpectrum(animals)
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
%_______________________________________________________________
%   RETURN:                     
%                               
%_______________________________________________________________

%% Setup
RestCriteria.Fieldname = {'Duration','PuffDistance'};
RestCriteria.Comparison = {'gt','gt'};
RestCriteria.Value = {14,5};
RestBuffer = 0;
AllRest = [];
for a = 1:length(animals)
    animal = animals{a};
    %% Load the data
    prevdir = cd([animal filesep]);
    
    RestFile = ls('*_RESTDATA_HR.mat');
    load(RestFile);
    DataField = fieldnames(RestData);
    DataField = DataField{1};
    AllRestData = RestData;
    
    %% Throw out periods of rest according to the RestCriteria    
    [RestFiltArray] = FilterEvents(AllRestData.(DataField),RestCriteria);
    RestData = AllRestData.(DataField).Data(RestFiltArray);
    
    %% Set up spectral filter for the data
%     [z,p,k] = butter(4,2/(AllRestData.(DataField).Fs/2),'low');
%     [sos,g] = zp2sos(z,p,k);
    Fs = AllRestData.(DataField).Fs;
    
    %% Concatenate the resting events
    Rest = cell(1,length(RestData));
    CellSeg = cell(size(RestData));
    for RD = 1:length(RestData)
        RestBuffer_Ind = RestBuffer*Fs+1;
        
        clippedRest = RestData{RD}(RestBuffer_Ind:end);
        
        % Find any NaN in the HR data
        NaNind = not(isnan(clippedRest));
        clippedRest = clippedRest(NaNind);
        RestOffset = mean(clippedRest);
        
        % Filter the resting heart rate
        Rest{RD} = detrend(clippedRest-RestOffset);
        
        % Track the segment starts/end
        CellSeg{RD} = zeros(size(clippedRest));
        CellSeg{RD}(1) = 1;
        CellSeg{RD}(end) = -1;
    end
    JoinedRest = [Rest{:}];
    SegMarkers = [CellSeg{:}];
    SegStarts = find(SegMarkers==1);
    SegStops = find(SegMarkers==-1);
    params.tapers = [1 1];
    params.fpass = [0.1 2];
    params.Fs = Fs;
    params.err = [1 0.05];
    [S,f]= mtspectrumc_unequal_length_trials(JoinedRest',[14 14],params,...
    [SegStarts', SegStops']);
    if a == 1
        Spectra = NaN*ones(length(animals),length(S));
    end
    Spectra(a,:) = S;
    
    cd(prevdir);
end

figure; plot(f,[mean(Spectra)', (mean(Spectra)+std(Spectra))', ...
    (mean(Spectra)-std(Spectra))'],'k');
[pxx,f] = plomb(AllRest,Fs);