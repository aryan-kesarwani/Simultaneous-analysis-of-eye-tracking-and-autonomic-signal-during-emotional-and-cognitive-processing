function [nData, rmFlag] = detectSac(Data, evtList)
%%
numTrial = length(Data);
nData    = cell(numTrial,1);
rmFlag = false(numTrial,1);

for nTrial = 1:numTrial
    subData = Data{nTrial};

    try
    l_sac = getMsac(subData(:,2),subData(:,3));
    r_sac = getMsac(subData(:,7),subData(:,8));
    [sac, ~, ~] = binsacc(l_sac,r_sac);
    catch
        rmFlag(nTrial) = true;
        continue
    end

    if isempty(sac)
        continue
    end

    onsetTP  = find(subData(:,1) == evtList(nTrial),1);
    if isempty(onsetTP)
        onsetTP  = find(subData(:,1) == evtList(nTrial)+1,1);
    end

    sac(:,1) = sac(:,1) - onsetTP;
    sac(:,2) = sac(:,2) - onsetTP;
    sac(:,8) = sac(:,8) - onsetTP;
    sac(:,9) = sac(:,9) - onsetTP;
    nData{nTrial} = sac;
end

function sac =getMsac(x,y)
VFAC     = 6;
MINDUR   = 6;
fs = 1000;
% [B,A] = butter(3,200/500,"low");
% x = filtfilt(B,A,double(x));
% y = filtfilt(B,A,double(y));

eye_pos = [x,y];
eye_vel = zeros(size(eye_pos));
eye_vel(2:end-1,:) = (eye_pos(3:end,:)-eye_pos(1:end-2,:))*fs/2;

[sac, ~] = microsacc(eye_pos,eye_vel,VFAC,MINDUR);

