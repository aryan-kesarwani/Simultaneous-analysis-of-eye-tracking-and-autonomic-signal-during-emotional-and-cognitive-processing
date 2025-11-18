function rmSacData = rmSac(sacData, tw, thr, vthr, dur, intv)

if nargin < 2 || isempty(tw),   tw   = model.epoch; end
if nargin < 3 || isempty(thr),  thr  = [0.1 2];     end
if nargin < 4 || isempty(vthr), vthr = [0 200];     end
if nargin < 5 || isempty(dur),  dur  = false;       end
if nargin < 6 || isempty(intv), intv = false;       end
%%
numTrial  = length(sacData);
rmSacData = cell(numTrial,1);
for nTrial = 1:numTrial
    sac = sacData{nTrial};

    if isempty(sac)
        continue
    end

    if exist("tw","var")
        rmIdx = sac(:,1)<tw(1) | sac(:,1)>tw(2);
        sac(rmIdx,:) = [];
    end

    if exist("thr","var")
        amp = sqrt(sac(:,6).^2 + sac(:,7).^2);
        rmIdx = amp < thr(1) | amp > thr(2);
        sac(rmIdx,:) = [];
    end
    
    if exist("vthr","var")
        rmIdx = sac(:,3) < vthr(1) | sac(:,3) > vthr(2);
        sac(rmIdx,:) = [];
    end

    if dur
        rmIdx = sac(:,2)-sac(:,1) > dur(2) | sac(:,2)-sac(:,1) < dur(1);
        sac(rmIdx,:) = [];
    end

    if intv && size(sac,1)>1
        rmIdx = [nan; sac(2:end,1) - sac(1:end-1,2)] < intv;
        sac(rmIdx,:) = [];
    end

    rmSacData{nTrial} = sac;
end
