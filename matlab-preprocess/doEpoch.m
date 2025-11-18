function [l_epData, r_epData, rmFlag] = doEpoch(Data,evtList,tw)
%%
numTrial = length(Data);
l_epData = nan(tw(2)-tw(1),numTrial);
r_epData = nan(tw(2)-tw(1),numTrial);
rmFlag   = false(numTrial,1);

for nTrial = 1:numTrial
    subData = Data{nTrial};
    subEvt  = evtList(nTrial);

    tp = find(subData(:,1) == subEvt);

    try
        l_epData(:,nTrial) = subData(tp+tw(1)+1:tp+tw(2),4);
        r_epData(:,nTrial) = subData(tp+tw(1)+1:tp+tw(2),9);
    catch
        rmFlag(nTrial) = true;
    end
end
    

