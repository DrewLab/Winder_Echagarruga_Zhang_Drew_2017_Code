function out = findSpeedSeg(s, Fr, tRest, tRun)
    % findSpeedSeg: find rest and running segments in locomotion data
    % 
    % USAGE:
    %       out = findSpeedSeg(s, Fr);
    %       out = findSpeedSeg(s, Fr, tRest);
    %       out = findSpeedSeg(s, Fr, tRest, tRun);
    %
    % INPUTS:
    %       s: binarized speed
    %       Fr: frame rate, in Hz
    %       tRest: criteria for resting period, default = 3 seconds
    %       tRun: criteria for running period, default = 5 seconds
    %
    % OUTPUTS:
    %       out.tRest: resting period length, in seconds
    %       out.tRun: running period length, in seconds
    %       out.Rest: start and end point of SELECTED resting periods,
    %                 N-by-2 matrix (contains resting period meets criteria)
    %       out.Run: start and end point of SELECTED running periods,
    %                N-by-2 matrix (contains running period meets criteria)
    %       out.RestR: start and end point of SELECTED resting periods (OFFSET)
    %                 N-by-2 matrix (contains resting period meets criteria)
    %       out.RunR: start and end point of SELECTED running periods (OFFSET)
    %                N-by-2 matrix (contains running period meets criteria)
    %       out.durationRest: duration of each resting period, N-by-1 vector 
    %                         (contains all resting period)
    %       out.durationRun: duration of each running period, N-by-1 vector
    %                        (contains all running period)
    %       out.stRest: start point of each resting period, N-by-1 vector
    %                   (contains all resting period)
    %       out.finRest: final point of each resting period, N-by-1 vector
    %                   (contains all resting period)
    %       out.stRun: start point of each running period, N-by-1 vector
    %                  (contains all resting period)
    %       out.finRun: final point of each running period, N-by-1 vector
    %                   (contains all resting period)
    %
    % Last Modified: 10-10-2016
    % Modified By: Qingguang Zhang (qingguang.zhang@gmail.com)
    

    if nargin < 2
        error('Not enough input arguments');
    end
    
    if nargin < 3 
        tRest = 3; % default = 3 seconds
        tRun = 5; % default = 5 seconds
    end
    
    
    out.tRest = tRest; % output resting period length
    out.tRun = tRun; % output running period length
    
    % Make sure there is no NaN in the speed data
    idx = find(isnan(s)); s(idx) = 0; clear idx;
    
    %% First, get the start and end points of ALL resting and running periods
    if s(1) == 0 % resting at the beginning of the trial
        a = 1;
    else % otherwise, running at the beginning of the trial
        a = 0;
    end
    
    if s(end) == 0 % resting at the end of the trial
        b = 1;
    else % otherwise, running at the end of the trial
        b = 0;
    end
    
    s1 = [a; s; b;]; % form a new speed time series
    
    A1 = diff(s1); % first derivative
    A1(end) = []; % remove extra time points we added in the previous step
    
    stRest = find(A1 == -1); % start point of rest
    stRun = find(A1 == 1); % start point of run

    s2 = flipud(s1); % flip the formed speed time series
    
    A2 = diff(s2); % first derivative
    A2 = flipud(A2); % flip again
    A2(1) = [];
    
    finRest = find(A2 == -1); % end point of rest
    finRun = find(A2 == 1); % end point of running
    
    % -------------------------------------------------------- %
    % start and final point of all resting and running periods
    % -------------------------------------------------------- %
    out.allSeg.stRest = stRest; % resting start point
    out.allSeg.finRest = finRest; % resting final point
    out.allSeg.stRun = stRun; % running start point
    out.allSeg.finRun = finRun; % running final point
    
    durationRest = finRest-stRest+1; % duration of each resting segment 
    durationRun = finRun-stRun+1; % duration of each running segment
    
    out.allSeg.durationRest = durationRest(:);
    out.allSeg.durationRun = durationRun(:);

    %% Convert short (< 1 second) resting periods to continuous running period
    % find all rest period less than 1 second
    shortRest = find(durationRest < 1 * Fr); % location of short rest period

    if ~isempty(shortRest) % short resting period exist
        
        if isequal(stRest(shortRest(1)), 1) % first segment is rest &  less than 1 second
            shortRest(1) = NaN; % remove first short period, because we can not convert it to continuous running period
            stRest(1) = NaN; % remove start point of first rest period
            finRest(1) = NaN; % remove end point of first rest period
        end
        
        if isequal(finRest(shortRest(end)),length(s)) % last segment is rest and less than 1 second
            shortRest(end) = NaN; % remove this short period, because we can not convert it to continuous running period
            stRest(end) = NaN; 
            finRest(end) = NaN;          
        end
        
        % Remove NaN elements, so we don't need to index these values
        shortRest = shortRest(~isnan(shortRest));
        % DO NOT remove stRest and finRest, because we need to use the index value again

        s3 = s; % copy a new speed data
        % if less than 1 second, the resting period become running period
        for i = 1:length(shortRest)
            s3(stRest(shortRest(i)):finRest(shortRest(i))) = 1; % generate a new one without < 1 second rest
        end
        
        % Re-search the time series for all running and resting period
        if s3(1) == 0, a = 1; else a = 0; end
        
        if s3(end) == 0, b = 1; else b = 0; end
        
        s4 = [a; s3; b;];
        
        A4 = diff(s4);
        A4(end) = [];
        
        stRest = find(A4 == -1); % start point of rest, onset
        stRun = find(A4 == 1); % start point of run, onset
        
        s5 = flipud(s4);
        
        A5 = diff(s5);
        A5 = flipud(A5);
        A5(1) = [];
        
        finRest = find(A5 == -1);
        finRun = find(A5 == 1);
        
        durationRest = finRest-stRest+1; % duration of each resting segment
        durationRun = finRun-stRun+1; % duration of each running segment
        
        % -------------------------------------------------------- %
        % start and final point of adjusted resting and running periods
        % -------------------------------------------------------- %
        out.adjustSeg.stRest = stRest; % resting start point
        out.adjustSeg.finRest = finRest; % resting final point
        out.adjustSeg.stRun = stRun; % running start point
        out.adjustSeg.finRun = finRun; % running final point        
        out.adjustSeg.durationRest = durationRest(:);
        out.adjustSeg.durationRun = durationRun(:);
    end
    
    %% Select Rest + Run pair for Locomotion Triggered Response
    
    % the first segment should not be running stage, and
    % the last segment should not be resting stage
    if s(1) == 1 % first segment is running period, remove this segment
        stRun(1) = [];
        finRun(1) = [];
        durationRun(1) = [];
    end
    
    if s(end) == 0 % last segment is resting period, remove this segment
        stRest(end) = [];
        finRest(end) = [];
        durationRest(end) = [];
    end

    % now we get the Rest + Run pair, check if each pair meets criteria
    segRest = find(durationRest >= tRest*Fr);
    segRun = find(durationRun >= tRun*Fr);
    
    % find the pair meets both criteria
    seg = intersect(segRest, segRun);
    
    if isempty(seg) % No pair meets criteria
        disp('findSpeedSeg -> No segment meets the ONSET selection criteria');
        out.Rest = [];
        out.Run = [];
    else 
        out.Rest = [stRest(seg), finRest(seg);]; % N-by-2 matrix, [st, fin;] X # of segments
        out.Run = [stRun(seg), finRun(seg);]; % N-by-2 matrix       
    end    
end







