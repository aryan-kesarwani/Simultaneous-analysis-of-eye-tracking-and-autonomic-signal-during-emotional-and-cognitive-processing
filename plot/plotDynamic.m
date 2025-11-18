function plotDynamic(Data,varargin)
%%
figure; hold on

Vars      = Data.Properties.VariableNames;
groupings = Data(:, Vars(3:end-1));      % 所有分組變數
[~, ~, groupIdx] = unique(groupings, 'rows');
numCond   = max(groupIdx);

p = inputParser;
addParameter(p, 'diff', false);
addParameter(p, 'tw', [1, 5000, 1]);
addParameter(p, 'col', jet(numCond));
addParameter(p, 'xlab', '');
addParameter(p, 'ylab', '');
addParameter(p, 'lnstyle', repmat({'-'},numCond,1))
parse(p, varargin{:});

tw  = p.Results.tw;
col = p.Results.col;
lnstyle = p.Results.lnstyle;

% 先把每一條件的「標籤」轉成字串 cell array
labelStr = cell(numCond,1);
for k = 1:numCond
    thisRow  = groupings(find(groupIdx==k,1), :);     % 抓其中一列模板
    labelStr{k} = strjoin(string(table2cell(thisRow)), ', '); % 例如 "Low, Compatible"
end

tp  = tw(1):tw(2);
for nCond = 1:numCond
    condData = Data.(Vars{end})(groupIdx == nCond, :);
    condData = condData(:, tw(1)-tw(3)+1 : tw(2)-tw(3)+1);
    errData  = std(condData)/sqrt(size(condData,1)-1);
    avgData  = mean(condData);

    %── 畫平均線並指定 DisplayName ──
    plot(tp, avgData, ...
        'LineWidth', 1, ...
        'Color', col(nCond,:), ...
        'LineStyle', lnstyle{nCond}, ...
        'DisplayName', labelStr{nCond});   % ←關鍵

    %── 畫 ±SEM 帶 ──
    hFill = fill([tp, fliplr(tp)], ...
         [avgData-errData, fliplr(avgData+errData)], ...
         col(nCond,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    set(hFill,'HandleVisibility','off');
end

xlim([tp(1) tp(end)]);
legend('show', 'Location', 'best'); % 自動抓 DisplayName
ylabel(p.Results.ylab,'fontsize',12)
xlabel(p.Results.xlab,'fontsize',12)
