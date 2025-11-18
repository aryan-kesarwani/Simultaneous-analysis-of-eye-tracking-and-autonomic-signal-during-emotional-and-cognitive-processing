function nData = upSample(Data)

nData = Data;

numTrial = length(Data);
for nTrial = 1:numTrial
    subData = Data{nTrial};
    tp = subData(:,1);
    while any(diff(tp)<=0)
        tp = rmDupes(tp);
    end
    subData = interp1(tp,subData,tp(1):1:tp(end));
    subData(:,1) = tp(1):1:tp(end);
    nData{nTrial} = subData;
end

