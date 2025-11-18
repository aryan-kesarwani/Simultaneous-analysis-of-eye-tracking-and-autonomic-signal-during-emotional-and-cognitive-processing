function op = plotPolar(projectDir,tw,iv,ivVal,nbins,col,lty,rmFlat)
%%
precDir = dir(projectDir + "/derivatives/prec/sub*");
numSubj = length(precDir);
op_sub = cell(numSubj, 1); % Preallocate cell array for efficiency

itv = 2*pi/nbins;
itv_tw = [itv/2:itv:2*pi-itv/2] - pi;

for nSubj = 1:numSubj
    pFlat = false;
    subIdx = sprintf('//sub-%02d',nSubj);
    subjDir = [precDir(nSubj).folder subIdx '//'];

    subFile = dir([subjDir '/*info.mat']);
    load([subFile.folder '//' subFile.name])
    if subInfo.olFlag
        continue
    end

    [labels,rmIdx] = getIV(subjDir,iv,ivVal);

    Ang = getDV(subjDir, 'Ang', []);
    stw = tw - subInfo.epoch(1);
    val = nan(length(rmIdx),nbins);
    num = nan(length(rmIdx),1);
    for i = 1:length(rmIdx)
        temp = Ang(stw(1):stw(2),rmIdx{i});
        temp = temp(~isnan(temp));
        if length(unique(temp)) < 3
            pFlat = true;
        end
        for k = 1:nbins
            if k == 1
                rng = [itv_tw(end) itv_tw(1)];
                id = temp < rng(2) | temp >= rng(1);
            else
                rng = [itv_tw(k-1) itv_tw(k)];
                id = temp > rng(1) & temp <= rng(2);
            end
            val(i,k) = sum(id)/length(temp);
            num(i) = length(temp);
        end
    end
    if pFlat && rmFlat
        continue
    end
    op = [table(num, 'VariableNames', {'n'}), ...
    labels, table(val, 'VariableNames', {'Polar'})];
    sub = table(repmat(nSubj, height(op), 1), 'VariableNames', "Subj");
    op_sub{nSubj} = [sub, op];
end
op = vertcat(op_sub{:});
%%
x = 0:itv:2*pi;
x = x(1:end-1);
Vars      = op.Properties.VariableNames;
groupings = op(:, Vars(3:end-1));      % 所有分組變數
[~, ~, groupIdx] = unique(groupings, 'rows');
numCond   = max(groupIdx);

for nCond = 1:numCond
    condData = op.(Vars{end})(groupIdx == nCond,:);
    errData  = std(condData)/sqrt(size(condData,1)-1);
    avgData  = mean(condData);
    polarplot(x([1:end 1]),avgData([1:end 1]), ...
        'Color',col(nCond,:), 'LineStyle',lty{nCond}); hold on
    for i = 1:length(x)
        polarplot([x(i), x(i)], ...
            [avgData(i)-errData(i), avgData(i)+errData(i)], ...
            'Color', col(nCond,:))
    end
end
