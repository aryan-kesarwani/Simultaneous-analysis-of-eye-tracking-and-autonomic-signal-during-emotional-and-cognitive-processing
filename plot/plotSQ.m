function op = plotSQ(projectDir,tw,iv,ivVal,col)
%%
precDir = dir(projectDir + "/derivatives/prec/sub*");
numSubj = length(precDir);
op_sub = cell(numSubj, 1); % Preallocate cell array for efficiency

for nSubj = 1:numSubj

    subIdx = sprintf('//sub-%02d',nSubj);
    subjDir = [precDir(nSubj).folder subIdx '//'];

    subFile = dir([subjDir '/*info.mat']);
    load([subFile.folder '//' subFile.name])
    if subInfo.olFlag
        continue
    end

    [labels,rmIdx] = getIV(subjDir,iv,ivVal);

    Amp = getDV(subjDir, 'Amp', []);
    Vel = getDV(subjDir, 'Vel', []);

    allVal = cell(length(rmIdx),2);
    num = nan(length(rmIdx),1);
    stw = tw - subInfo.epoch(1);
    for i = 1:length(rmIdx)
        temp = Amp(stw(1):stw(2),rmIdx{i});
        allVal{i,1} = temp(~isnan(temp));
        temp = Vel(stw(1):stw(2),rmIdx{i});
        allVal{i,2} = temp(~isnan(temp));
        num(i) = length(temp(~isnan(temp)));
    end
    op = [table(num, 'VariableNames', {'n'}), ...
    labels, table(allVal, 'VariableNames', {'Seq'})];
    sub = table(repmat(nSubj, height(op), 1), 'VariableNames', "Subj");
    op_sub{nSubj} = [sub, op];
end
op = vertcat(op_sub{:});

%%
figure; hold on
Vars      = op.Properties.VariableNames;
groupings = op(:, Vars(3:end-1));      % 所有分組變數
[~, ~, groupIdx] = unique(groupings, 'rows');
numCond   = max(groupIdx);

for nCond = 1:numCond
    condData = cell2mat(op.(Vars{end})(groupIdx == nCond, :));
    scatter(condData(:,1),condData(:,2), 5, ...
        'MarkerFaceColor', col(nCond,:), ...
        'MarkerEdgeColor','none', ...
        'MarkerFaceAlpha',0.2)
end

ylabel('Peak velocity (deg/s)','fontsize',12)
xlabel('Amplitude (deg)','fontsize',12)


