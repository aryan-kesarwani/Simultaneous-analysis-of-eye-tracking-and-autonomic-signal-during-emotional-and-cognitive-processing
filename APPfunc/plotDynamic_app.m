function h = plotDynamic_app(ax, Data, varargin)
    hold(ax, 'on')

    % 參數解析
    p = inputParser;
    addParameter(p, 'diff', false);  % 是否畫差異波形
    addParameter(p, 'tw', [1, 1, 1]);  % [start, end, baseline]
    addParameter(p, 'col', {[0 .5 .8]; [0.8 .1 0]; [0 .5 0]});
    addParameter(p, 'xlab', '');
    addParameter(p, 'ylab', '');
    addParameter(p, 'lnstyle', {'-', '-'});
    addParameter(p, 'DV', '');
    addParameter(p, 'BT', '');
    addParameter(p, 'WT', '');
    parse(p, varargin{:});

    tw       = p.Results.tw;
    col      = p.Results.col;
    lnstyle  = p.Results.lnstyle;
    DV       = p.Results.DV;
    BT       = p.Results.BT;
    WT       = p.Results.WT;
    doDiff   = p.Results.diff;

    % 若為 table 則轉為 cell array
    if istable(Data)
        op = cell(4,1);
        n = 1;
        for i = 1:2
            for j = 1:2
                op{n} = Data.(DV)(Data.(BT) == i & Data.(WT) == j, :)';
                n = n + 1;
            end
        end
        Data = op;
    end

    tp = tw(1):tw(2);  % 時間軸
    numCond = length(Data);

    for nCond = 1:numCond
        condData = Data{nCond};
        condData = condData(tw(1) - tw(3) + 1 : tw(2) - tw(3) + 1, :)';
        avgData = mean(condData, 1);
        errData = std(condData, 0, 1) / sqrt(size(condData, 1));

        plot(ax, tp, avgData, 'LineWidth', 1, ...
            'Color', col{nCond}, 'LineStyle', lnstyle{nCond});

        fill(ax, [tp, fliplr(tp)], ...
                  [avgData - errData, fliplr(avgData + errData)], ...
                  col{nCond}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end

    % 差異波形（通常用在 condition 2 vs 1）
    if doDiff && numCond >= 2
        d = (Data{2} - Data{1})';
        d = d(:, tw(1) - tw(3) + 1 : tw(2) - tw(3) + 1);
        avgDiff = mean(d, 1);
        errDiff = std(d, 0, 1) / sqrt(size(d, 1));

        plot(ax, tp, avgDiff, 'LineWidth', 1, 'Color', col{3});
        fill(ax, [tp, fliplr(tp)], ...
                 [avgDiff - errDiff, fliplr(avgDiff + errDiff)], ...
                 col{3}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end

    % 座標設定
    xlim(ax, [tp(1), tp(end)]);
    xlabel(ax, p.Results.xlab, 'FontSize', 12);
    ylabel(ax, p.Results.ylab, 'FontSize', 12);

    if nargout > 0
        h = ax;
    end
end
