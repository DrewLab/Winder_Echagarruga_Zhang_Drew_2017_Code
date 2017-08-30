function [] = GetHeartRate()
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

Fs = 30;
params.Fs = Fs;
params.fpass = [5 15];
params.tapers = [3 5];
SpGramWindow = [4 1/Fs];
lpfreq = (1/SpGramWindow(1));
[z,p,k] = butter(2,lpfreq/(Fs/2),'low');
[sos,g] = zp2sos(z,p,k);

filenames = ls('*_ProcData.mat');
[~,~,FileDates,FileIDs] = GetFileInfo(filenames);
CellDates = mat2cell(FileDates,ones(1,size(FileDates,1)),...
    size(FileDates,2));
[UniqueDates,~] = unique(CellDates);
Orig_Files_Dir = uigetdir(pwd,'Select the folder containing the raw .bin files...');

IndFrame = zeros(256,256,length(UniqueDates));
mask = zeros(256,256,length(UniqueDates));
for UD = 1:length(UniqueDates)
    ROIFileMatch = strcmp(CellDates,UniqueDates{UD});
    ProcFileNames = filenames(ROIFileMatch,:);
    ProcFileIDs = FileIDs(ROIFileMatch,:);
    prevdir = cd(Orig_Files_Dir);
    IndFrame(:,:,UD) = GetSingleCBVFrame([FileIDs(1,:) '_dalsa.bin'],256,256);
    if UD == 1
        % Create mask of the window for the first day
        figure; imagesc(IndFrame(:,:,UD));
        colormap('gray');
        display('Draw ROI around the thinned skull window:')
        mask(:,:,UD) = roipoly;
        [FixedRows,FixedCols] = find(mask(:,:,UD));
        close gcf;
    else
        % Register the image to uD=1
        [~,RA,RB,~] = RegisterCBVImage(IndFrame(:,:,1),IndFrame(:,:,UD));
        
        % Shift the mask based on the registered images
        % xdimension: RB-RA>0, shift mask left; RB-RA<0, shift mask right.
        % ydimension: RB-RA>0, shift mask down; RB-RA<0, shift mask up.
        xshift = round(RB.XWorldLimits(1)-RA.XWorldLimits(1));
        yshift = round(RB.YWorldLimits(1)-RA.YWorldLimits(1));
        shiftrow = FixedRows-yshift;
        shiftcol = FixedCols-xshift;
        for sr = 1:length(shiftrow)
            mask(shiftrow(sr),shiftcol(sr),UD) = true;
        end
    end
    
    for PFI = 1:size(ProcFileIDs,1)
        load([prevdir filesep ProcFileNames(PFI,:)])
        filename = [ProcFileIDs(PFI,:) '_dalsa.bin'];
        
        % Calculate the MEDIAN intensity across the window
        [Frames]=ReadDalsaBinary_Matrix(filename,256,256);
        mFrames = mean(Frames,3);
        BrightMask = gt(mFrames,2000);
        Refl = zeros(1,size(Frames,3));
        for F = 1:size(Frames,3)
            indFrame = Frames(:,:,F);
            Refl(F) = median(indFrame(and(mask(:,:,UD),BrightMask)));
        end
        clear Frames
        [S,timevec,freqs] = mtspecgramc(diff(Refl-mean(Refl)),SpGramWindow,params);
        NaNpad = NaN*ones(1,(floor(timevec(1)/SpGramWindow(2))-1));
        [~,freqinds] = max(S,[],2);
        HR = freqs(freqinds);
        filtHR = filtfilt(sos,g,HR);
        ProcData.HR = [NaNpad HR NaNpad];
        ProcData.HR_timevec = timevec;
        ProcData.Fs.HR_fs = Fs;
        save([prevdir '\' ProcFileNames(PFI,:)],'ProcData');
        imagesc(timevec,freqs,S'); axis xy;
        hold on; 
        plot(timevec,HR,'r','linewidth',2);
        plot(timevec,filtHR,'k','Linewidth',2); 
        pause(0.001);
        hold off;
    end
    cd(prevdir)
end