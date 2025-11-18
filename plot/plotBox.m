function plotBox(Data,varargin)
%%
figure; hold on

Vars      = Data.Properties.VariableNames;
groupings = Data(:, Vars(3:end-1));      % 所有分組變數
[~, ~, groupIdx] = unique(groupings, 'rows');
numCond   = max(groupIdx);

p = inputParser;
addParameter(p, 'col', {[0 .5 .8];[.8 .1 0];[0 .5 0]});
addParameter(p, 'xlab', '');
addParameter(p, 'ylab', '');
addParameter(p, 'dot', false);
parse(p, varargin{:});

col = p.Results.col;
dot = p.Results.dot;

labelStr = cell(numCond,1);
for k = 1:numCond
    thisRow  = groupings(find(groupIdx==k,1), :);     % 抓其中一列模板
    labelStr{k} = strjoin(string(table2cell(thisRow)), ', '); % 例如 "Low, Compatible"
end

for nCond = 1:numCond
    condData = Data.(Vars{end})(groupIdx == nCond, :);
    errData = std(condData)/sqrt(length(condData)-1);
    avgData = mean(condData);

    if dot
        scatter(nCond+(rand(length(condData),1)-.5)*0.5, ...
            condData,10, ...
            'LineWidth',.2,'MarkerFaceColor','none', ...
            'MarkerEdgeColor',col(nCond,:));
    end
    errorbar(nCond,avgData,errData, ...
        '-s','CapSize',10,"MarkerFaceColor",col(nCond,:), ...
        'MarkerSize',6,'Color','k','MarkerEdgeColor','k');
end
    xticks(1:numel(labelStr));
    xticklabels(labelStr);
%%
xlim([.5 numCond+.5])
ylabel(p.Results.ylab,'fontsize',12)
% xlabel(p.Results.xlab,'fontsize',12)



