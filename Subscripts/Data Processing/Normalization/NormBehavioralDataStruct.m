function [PercNormData,ZNormData] = NormBehavioralDataStruct(DataStruct,dataType)
%   function [PercNormData,ZNormData] = NormBehavioralDataStruct(DataStruct,dataType)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Normalizes measured data by periods of rest.
%   
%_______________________________________________________________
%   PARAMETERS:             
%               DataStruct - [struct] contains behavioral data
%
%               dataType - [string] designates the behavioral measurement
%               to be normalized
%_______________________________________________________________
%   RETURN:                     
%               PercNormData - [matrix,array] the data normalized to the
%                               resting baseline.
%
%               ZNormData - [matrix,array] the data normlized to the
%               resting standard deviation.
%_______________________________________________________________

% Load the Baselines Structure
wraptext(['NormBehavioralDataStruct: Normalizing ' dataType '...']);
BaseFile = dir('*Baselines.mat');
if isempty(BaseFile);
    error('No Baselines found in the current directory...');
elseif size(BaseFile,1)>1;
    error('Multiple Baseline files found in the current directory...');
end
load(BaseFile.name);
if not(isfield(Baselines,dataType))
    error(['No Baseline found for ' dataType])
end
PercNormData = DataStruct.Data;
ZNormData = DataStruct.Data;
Sessions = unique(DataStruct.FileDate);
for s = 1:length(Sessions)
    session = Sessions{s};
    session_str = ConvertDate(session);
    session_inds = strcmp(DataStruct.FileDate,session);
    
    % Calculate the baseline differently depending on data type
    if iscell(DataStruct.Data)
        session_data = DataStruct.Data(session_inds);
        norm_session_data = cell(size(session_data));
        znorm_session_data = cell(size(session_data));
        Zbaseline = mean(Baselines.(dataType).(session_str).StDev,1)';
        baseline = mean(Baselines.(dataType).(session_str).Means,1)';
        for sd = 1:size(session_data,1)
            cell_base = baseline*ones(1,size(session_data{sd},2));
            Zcell_base = Zbaseline*ones(1,size(session_data{sd},2));
            norm_session_data{sd} = session_data{sd}./cell_base;
            znorm_session_data{sd} = (session_data{sd}-mean(session_data{sd}))./Zcell_base;
        end
        PercNormData(session_inds) = norm_session_data;
        ZNormData(session_inds) = znorm_session_data;
    else
        Zbaseline = mean(Baselines.(dataType).(session_str).StDev,1);
        baseline = mean(Baselines.(dataType).(session_str).Means,1);
        % Preallocate array and use for permutation
        norm_session_data = DataStruct.Data(session_inds,:,:);
        znorm_session_data = DataStruct.Data(session_inds,:,:);
        
        % Permute norm_session_data to handle both matrix and array (squeeze
        % causes a matrix dimension error if not permuted)
        session_data = permute(norm_session_data,unique([2,1,ndims(norm_session_data)],...
            'stable'));
        for sd = 1:size(session_data,2)
            norm_session_data(sd,:,:) = squeeze(session_data(:,sd,:))./...
                (ones(size(session_data,1),1)*baseline);
            dc = mean(squeeze(session_data(:,sd,:)));
            znorm_session_data(sd,:,:) = (squeeze(session_data(:,sd,:))-dc)./...
                (ones(size(session_data,1),1)*Zbaseline);
        end
        PercNormData(session_inds,:,:) = norm_session_data;
        ZNormData(session_inds,:,:) = znorm_session_data;
    end
end

    