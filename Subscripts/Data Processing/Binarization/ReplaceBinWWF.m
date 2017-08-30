function [] = ReplaceBinWWF()
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

filenames = uigetfile('*ProcData.mat','multiselect','on');
[~,~,FileDate,~] = GetFileInfo(filenames);
DateCell = mat2cell(FileDate,ones(1,size(FileDate,1)),size(FileDate,2));
UniqueDates = unique(DateCell);

for UD = 1:length(UniqueDates)
    DateFilter = strcmp(DateCell,UniqueDates{UD});
    DayFiles = filenames(DateFilter);
    for DF = 1:length(DayFiles)
        load(DayFiles{DF})
        [animal,hem,numday,~] = GetFileInfo(DayFiles{DF});
        strday = ConvertDate(numday);
        
        wwf = ProcData.wwf;
        wwf_fs = ProcData.Fs.wwf_fs;
        if DF == 1
            Threshfile = ls('*Thresholds.mat');
            if not(isempty(Threshfile));
                load(Threshfile)
            end
            [wwf_thresh1,wwf_thresh2] = CreateWhiskThreshold(wwf, wwf_fs, strday);
            Thresholds.(['BinWWF_Lower_' strday]) = wwf_thresh1;
            Thresholds.(['BinWWF_Upper_' strday]) = wwf_thresh2;
            save([animal '_' hem '_Thresholds.mat'],'Thresholds');
        end
        load([animal '_' hem '_Thresholds.mat']);
        bin_wwf = binarize_wwf(wwf,wwf_fs,Thresholds.(['BinWWF_Lower_' strday]),...
            Thresholds.(['BinWWF_Upper_' strday]));
        close gcf;
        
        Th1 = find(bin_wwf==1);
        Th2 = find(bin_wwf==0.5);
        
        plot(wwf); axis tight;
        hold on;
        scatter(Th2,15*ones(size(Th2)),'c.');
        scatter(Th1,15*ones(size(Th1)),'k.');
        ylim([-30 30])
        ylabel('Whisker Angle')
        title(strday);
        hold off;
        
             
%         bin_inds = find(bin_wwf);
%         plot(wwf); axis tight;
%         hold on; scatter(bin_inds,45*ones(size(bin_inds)),'r.');
%         hold off;
%         ylim([-50 50])
        pause(0.1)
        ProcData.Bin_wwf = bin_wwf;
        
        save(DayFiles{DF},'ProcData')
    end
end