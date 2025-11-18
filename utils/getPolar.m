function [op,x] = getPolar(allModels,nbins,DV,ep,olVars)
%%
itv = 2*pi/nbins;
x = 0:itv:2*pi;
x = x(1:end-1);

op = cell2mat(arrayfun(@(m) getBinValues(m,DV,ep,nbins,olVars), allModels, ...
    'UniformOutput',false));


end

function temp = getBinValues(model,DV,ep,nbins,olVars)
itv = 2*pi/nbins;
tw = itv/2:itv:2*pi-itv/2;
rmIdx = reduceOR(model.rmFlag, olVars);
allVal = cellfun(@(x) x, ...
    model.rmSacData(~cellfun(@isempty, model.rmSacData)&~rmIdx), ...
    'UniformOutput', false);
allVal = cell2mat(allVal);
temp = nan(nbins,1);
if isempty(allVal)
    return
end
idx = allVal(:,1) >= ep(1) & allVal(:,1) <= ep(2);
allVal = allVal(idx,:);
ang = cart2pol(allVal(:,6),allVal(:,7));
ang(ang<0) = ang(ang<0) + 2*pi;

for k = 1:nbins
    if k == 1
        rng = [tw(end) tw(1)];
        id = ang < rng(2) | ang >= rng(1);
    else
        rng = [tw(k-1) tw(k)];
        id = ang > rng(1) & ang <= rng(2);
    end
    temp(k) = getDV(allVal,DV,id);
end
end

function op = getDV(allVal,DV,id)
switch DV
    case 'Rate'
        op = sum(id)/size(allVal,1);
    case 'Amp'
        sac = allVal(id,:);
        op = mean(sqrt(sac(:,6).^2 + sac(:,7).^2),'omitnan');
    case 'Vel'
        op = mean(allVal(id,3),'omitnan');
end
end


function idx = reduceOR(rmFlagStruct, olVars)
    if isempty(olVars)
        idx = false(size(rmFlagStruct.ACC));
        return;
    end

    idx = false(size(rmFlagStruct.(olVars{1})));
    for i = 1:length(olVars)
        idx = idx | rmFlagStruct.(olVars{i});
    end
end



