function plotBox_app(ax, Data, varargin)
    hold(ax, 'on')

    % Input parser
    p = inputParser;
    addParameter(p, 'col', {[0 .5 .8]; [.8 .1 0]; [0 .5 0]});
    addParameter(p, 'xlab', '');
    addParameter(p, 'ylab', '');
    addParameter(p, 'DV', '');
    addParameter(p, 'BT', '');
    addParameter(p, 'WT', '');
    addParameter(p, 'dot', false);  % 是否畫出所有資料點
    parse(p, varargin{:});

    col = p.Results.col;
    dot = p.Results.dot;
    DV  = p.Results.DV;
    BT  = p.Results.BT;
    WT  = p.Results.WT;

    % 將 table 拆成 cell array（4 條件）
    if istable(Data)
        op = cell(4,1);
        n = 1;
        for i = 1:2
            for j = 1:2
                op{n} = Data.(DV)(Data.(BT) == i & Data.(WT) == j)';
                n = n + 1;
            end
        end
        Data = op;
    end

    numCond = length(Data);

    % 主畫圖迴圈
    for nCond = 1:numCond
        condData = Data{nCond};
        errData = std(condData) / sqrt(length(condData) - 1);
        avgData = mean(condData);

        % 畫 dot scatter（可選）
        if dot
            scatter(ax, nCond + (rand(length(condData),1)-0.5)*0.5, ...
                condData, 10, ...
                'LineWidth', 0.2, ...
                'MarkerFaceColor', 'none', ...
                'MarkerEdgeColor', col{nCond});
        end

        % 畫 error bar + 平均值
        errorbar(ax, nCond, avgData, errData, ...
            '-s', 'CapSize', 10, ...
            'MarkerFaceColor', col{nCond}, ...
            'MarkerSize', 6, ...
            'Color', 'k', ...
            'MarkerEdgeColor', 'k');
    end

    % 座標軸設定
    xlim(ax, [0.5, numCond + 0.5]);
    xlabel(ax, p.Results.xlab, 'FontSize', 12);
    ylabel(ax, p.Results.ylab, 'FontSize', 12);
end
