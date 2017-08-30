function results = compute_resting_variance(data)
% compute_resting_variance: compute variance of CBV signal and neural
%                           activity signal for each trial
% Input:
%       data: raw data for each trial
% Output:
%       results: 
% See also:
%       findSpeedSeg.m

animalList = fieldnames(data); % a list of mouse ID
drugList = fieldnames(data.(animalList{1})); % a list of drug names
sessionList = fieldnames(data.(animalList{1}).(drugList{1}));

%% Compute resting variance for CBV signals and neural activity signals
for animalIdx = 1:numel(animalList)
    animal = animalList{animalIdx};
    for drugIdx = 1:numel(drugList)
        drug = drugList{drugIdx};
        for sessionIdx = 1:numel(sessionList)
            session = sessionList{sessionIdx};
            trialList = fieldnames(data.(animal).(drug).(session));
            
            CBV_sd = []; % standard deviation of CBV fluctuations
            GammaPower_sd = []; % standard deviation of Gammaband Power fluctuations
            MUAPower_sd = []; % standard deviation of MUA Power fluctuations
            GammaPower_ave = []; % average of Gammaband Power
            MUAPower_ave = []; % average of MUA power
            
            for trialIdx = 1:numel(trialList)
                trial = trialList{trialIdx};
                % get data for specific trial
                Fr = data.(animal).(drug).(session).(trial).Fr; 
                CBV = data.(animal).(drug).(session).(trial).CBV; 
                GammaPower = data.(animal).(drug).(session).(trial).GammaPower;
                MUAPower = data.(animal).(drug).(session).(trial).MUAPower;
                speed = data.(animal).(drug).(session).(trial).speed;
 
                allRestPeriod = []; % all rest period in one trial
                
                % find resting and running segments
                seg = findSpeedSeg(speed, Fr, 14,0);
                
                if ~isempty(seg.Rest) % resting segment exists
                    for idx = 1:size(seg.Rest,1) % for each resting segment
                        restPeriod = (seg.Rest(idx,1)+4*Fr):seg.Rest(idx,2); % one segment of the rest period
                        allRestPeriod = [allRestPeriod, restPeriod];
                        CBV_var_tmp(idx) = nanstd(detrend(CBV(restPeriod)));
                        GammaPower_var_tmp(idx) = nanstd(detrend(GammaPower(restPeriod)));
                        MUAPower_var_tmp(idx) = nanstd(detrend(MUAPower(restPeriod)));
                    end 
                end
                
                
                CBV_sd = [CBV_sd;CBV_var_tmp(:);];
                GammaPower_sd = [GammaPower_sd; GammaPower_var_tmp(:);];
                MUAPower_sd = [MUAPower_sd; MUAPower_var_tmp(:);];

                GammaPower_ave = [GammaPower_ave;GammaPower(allRestPeriod);];
                MUAPower_ave = [MUAPower_ave; MUAPower(allRestPeriod);];
                
                CBV_var_tmp(idx) = [];
                GammaPower_var_tmp(idx) = [];
                MUAPower_var_tmp(idx) = [];
            end
            results.(animal).(drug).(session).GammaPower.ave = nanmean(GammaPower_ave);
            results.(animal).(drug).(session).MUAPower.ave = nanmean(MUAPower_ave);
            
            results.(animal).(drug).(session).CBV.sd = nanmean(CBV_sd);
            results.(animal).(drug).(session).GammaPower.sd = nanmean(GammaPower_sd);
            results.(animal).(drug).(session).MUAPower.sd = nanmean(MUAPower_sd);
       end
    end
end             

end