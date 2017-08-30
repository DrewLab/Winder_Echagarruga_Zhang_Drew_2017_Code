function Muscimol_graph(~)
% close all;
clear all
for abc=1:2
    %open excl file to read in files
    [ndata, texter, alldata]=xlsread('Muscimol_Data.xlsx');
    s=size(alldata,1);
    s=s-1;
    
    %% Set up table for mixed effects ANOVA- ATW
    InfusionTypes = cell((size(alldata,1)-1)*2,1);
    InfusionTypes(1:(size(alldata,1)-1)) = {'aCSF'};
    InfusionTypes(size(alldata,1):end) = {'Muscimol'};
    Vessel = [1:(size(alldata,1)-1), 1:(size(alldata,1)-1)]';
    Animals = [alldata(2:end,3); alldata(2:end,3)]; % start at 2 to omit header
    
    clear resting_raw_diameter
    clear resting_raw_diameter_normalized
    %load each file individually and exclude all running times
    for n=1:s
        filename(n,1)=strcat(alldata(n+1,3),alldata(n+1,abc),alldata(n+1,4));
        str=char(filename(n,1));
        load(str);


        %% mv_mpP.T_stand_seg is filled with rest that is separated based on being greater than 10sec and at least 2 sec after running ends
        if size(mv_mpP.velocity.T_stand_seg,2)>0 
            started1=1;
            for filling_in=1:size(mv_mpP.velocity.T_stand_seg,2)
                start_tr=mv_mpP.velocity.T_stand_seg(1,filling_in);
                end_tr=mv_mpP.velocity.T_stand_seg(2,filling_in);
                lel=end_tr-start_tr;
                resting_raw_diameter(n,started1:started1+lel)=mv_mpP.Vessel.diameter_norm(start_tr:end_tr);   
                started1=started1+lel;
            end
        else
            resting_raw_diameter(n,:)=NaN;
        end

        %% filling in empty spaces with NaN
        for jo=1:size(resting_raw_diameter,1)
            for jo2=1:size(resting_raw_diameter,2)
                if resting_raw_diameter(jo,jo2)==0
                    resting_raw_diameter(jo,jo2)=NaN;
                end
            end
        end
    end

    %% normalizing resting diameter by mean of resting diameter
    for ihm=1:size(resting_raw_diameter,1)
        resting_raw_diameter_normalized(ihm,:)=resting_raw_diameter(ihm,:)./nanmean(resting_raw_diameter(ihm,:));
    end
        %% taking the standard dev of normalized resting diameter
    for ihm=1:size(resting_raw_diameter,1)  
        y_normalized(abc,ihm)=nanstd(sgolayfilt(resting_raw_diameter_normalized(ihm,:),3,9));      
    end
    
    %% filling in mat files to use the normalized resting raw diameter for graph example
    if abc==1
        acsf_ex=resting_raw_diameter_normalized(6,1:160);
    elseif abc==2
        mus_ex=resting_raw_diameter_normalized(6,1:160);
    end
end

%% 2P individual example
plot((1:160)./8,detrend(sgolayfilt(acsf_ex,3,15)),'k','Linewidth',1.5);
hold on
plot((1:160)./8,detrend(sgolayfilt(mus_ex,3,15)),'c','Linewidth',1.5);
ylim([-0.05 0.05])
ylabel('\DeltaD/D');
xlabel('Time (s)')



%% graphing per animal_normalized
figure
bar(1,nanmean(y_normalized(2,1:5)./y_normalized(1,1:5)),'BarWidth',0.3,'FaceColor',[0.2 0.2 0.2])
hold on
scatter(ones(1,5),y_normalized(2,1:5)./y_normalized(1,1:5),'ko','MarkerFaceColor','w')
bar(2,nanmean(y_normalized(2,6:9)./y_normalized(1,6:9)),'BarWidth',0.3,'FaceColor',[0.4 0.4 0.4])
scatter(2*ones(1,4),y_normalized(2,6:9)./y_normalized(1,6:9),'ko','MarkerFaceColor','w')
bar(3,nanmean(y_normalized(2,10:13)./y_normalized(1,10:13)),'BarWidth',0.3,'FaceColor',[0.6 0.6 0.6])
scatter(3*ones(1,4),y_normalized(2,10:13)./y_normalized(1,10:13),'ko','MarkerFaceColor','w')
bar(4,nanmean(y_normalized(2,14:19)./y_normalized(1,14:19)),'BarWidth',0.3,'FaceColor',[0.8 0.8 0.8])
scatter(4*ones(1,6),y_normalized(2,14:19)./y_normalized(1,14:19),'ko','MarkerFaceColor','w')
ylabel('\DeltaD/D \sigma_{muscimol} / \DeltaD/D \sigma_{aCSF}')
xlabel('Mouse #')
axis([0.5 4.5 0 2])
ax = gca;
set(ax,'XTick',1:4);
set(ax,'XTickLabel',1:4);
ylim([0 2])
axis([0.5 4.5 0 2])

reduction = y_normalized(2,:)./y_normalized(1,:);
mean_reduc = mean(reduction);
std_reduc = std(reduction);

% Linear Mixed effects model
t = table(categorical(Animals), categorical(Vessel),...
    categorical(InfusionTypes), [y_normalized(1,:)'; y_normalized(2,:)'],...
    'VariableNames',{'ID','V','IT','D'});

% Model Justification: The diameters of 19 vessels were measured after
% muscimol and aCSF infusion from four animals. This is a crossed design
% where every vessel was imaged under both treatments. Thus the first level
% occurs at the level of measurement:
% Diameter(i,j) = B(0,j) + B(1,j)*Vessel + error(i,j)
%   where i indexes the treatment, and j indexes the animal. The mean
%   diameter among animals [B(0,j)] does not change making B(0,j) fixed.
%   The vessels are chosen randomly, adding a measurement uncertainty
%   U(1,j) to the contributions of the vessels (vessels=random effect). 
% Thus:
% B(0,j) = Y(0,0) + Y(0,1)*Treatment,
% B(1,j) = Y(1,0) + Y(1,1)*Treatment + U(1,j).

lme = fitlme(t,'D ~ 1 + IT + (IT|V)');
title(['p=' num2str(lme.Coefficients.pValue(2)) '; t(' ...
    num2str(lme.Coefficients.DF(2)) ')=' num2str(lme.Coefficients.tStat(2))...
    ' mean\pmst.dev = ' num2str(round(mean_reduc*100)/100) '\pm' num2str(round(std_reduc*100)/100)]);
end
    
