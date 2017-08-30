function [MeanTrigSpecgram,timevec,Freqs] = TriggeredSpecgram(EventData,Behavior)
%   function [MeanTrigSpecgram,timevec,Freqs] = TriggeredSpecgram(EventData,Behavior)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Calculates and compiles the spectrogram surrounding a
%   behavioral event.
%   
%_______________________________________________________________
%   PARAMETERS:             
%               EventData - [Struct] Chunked data structure containing the
%               data surrounding behavioral events. This structure is used
%               to compile the event related spectrograms.
%
%               Behavior - [string] the behavioral category
%_______________________________________________________________
%   RETURN:                     
%               MeanTrigSpecgram - [matrix] contains the spectrograms
%               surrounding behavioral events
%
%               timevec - [array] vector containing the peri-event time
%
%               Freqs - [array] the frequency bins within the spectrogram
%_______________________________________________________________

[DataStruct,FiltArray] = SelectBehavioralEvents(EventData,Behavior);

FileIDs = DataStruct.FileID(FiltArray);
EventTimes = DataStruct.EventTime(FiltArray);
UniqueIDs = unique(FileIDs);

EventInd = 1;
for UID = 1:length(UniqueIDs)
%     SpFile = ls(['*' UniqueIDs{UID} '_Specgram.mat']);
%     load(SpFile);
    SpFile = dir(['*' UniqueIDs{UID} '_Specgram.mat']);
    load(SpFile.name);
    % Preallocate array for all events
    if UID == 1
        ChunkLength = DataStruct.epoch.duration*Specgram.Fs+1;
        TrigSpecgrams = NaN*ones(length(Specgram.Freqs),ChunkLength,length(FileIDs));
    end
    
    FileInds = strcmp(FileIDs,UniqueIDs{UID});
    FilePuffTimes = EventTimes(FileInds);
    
    for FPT = 1:length(FilePuffTimes)
        [~,SpecInd] = min(abs(Specgram.Time-FilePuffTimes(FPT)));
        strt = SpecInd - DataStruct.epoch.offset*Specgram.Fs;
        if strt <= 0
            continue;
        end
        stp = strt + DataStruct.epoch.duration*Specgram.Fs;
        
        % Normalize to preEvent time
        NormData = ones(ChunkLength,1)*mean(Specgram.Power(strt:SpecInd,:));
        TrigSpecgrams(:,:,EventInd) = (Specgram.Power(strt:stp,:)./NormData-1)';
        EventInd = EventInd+1;
    end
end
TrigSpecgrams(:,:,EventInd:end) = [];
timevec = (0:1/Specgram.Fs:DataStruct.epoch.duration)-DataStruct.epoch.offset;
Freqs = Specgram.Freqs;
MeanTrigSpecgram = mean(TrigSpecgrams,3);