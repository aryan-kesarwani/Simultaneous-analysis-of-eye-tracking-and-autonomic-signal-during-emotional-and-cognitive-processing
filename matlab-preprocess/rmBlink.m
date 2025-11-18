function [nData, rmFlag] = rmBlink(Data,type,eye)
%%
if ~ismember(type, ["interp","nan"])
    error("Type error")
end

if strcmp(eye,"left")
    eyeIdx = 4;
    clmIdx = 2:7;
elseif strcmp(eye,"right")
    eyeIdx = 9;
    clmIdx = 8:13;
end

nData = Data;
rmFlag = false(length(Data),1);

numTrial = length(Data);
for nTrial = 1:numTrial
    subData = Data{nTrial};

    diamData = subData(:,eyeIdx);
    timeData = subData(:,1);
    while any(diff(timeData)<=0)
        timeData = rmDupes(timeData);
    end
    [valOut,~,~] = rawDataFilter(timeData, diamData);
    valOut([1 end]) = 1;

    validTimeData = timeData(valOut);
    validDiamData = diamData(valOut);

    if strcmp(type,"interp")
        interpData = interp1(validTimeData,validDiamData,timeData,'linear');
        nData{nTrial}(:,eyeIdx) = interpData;

    elseif strcmp(type, "nan")
        [~, idx] = ismember(validTimeData, timeData);
        nanData = nan(length(timeData), length(clmIdx));
        nanData(idx,:) = subData(valOut,clmIdx);
        nData{nTrial}(:,clmIdx) = nanData;
    end

    rmFlag(nTrial) = mean(valOut) < 0.3;
end





