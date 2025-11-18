function result = mixed_anova(op,DV,BT,WT,plotflag)
% MIXED_ANOVA - 2x2 mixed ANOVA with effect size, post hoc and plot.
%
% Input:
%   DV       - dependent variable (Nx1)
%   SubjID   - subject ID (Nx1)
%   Within   - within-subject factor (Nx1)
%   Between  - between-subject factor (Nx1)
%
% Output:
%   result   - struct with stats, effect size, and plots

if nargin < 5; plotflag = false; end

% Prepare
SubjID = op.Subj;
Within = op.(WT);
Between = op.(BT);
DV = op.(DV);

SubjID   = categorical(SubjID(:));
Within   = categorical(Within(:));
Between  = categorical(Between(:));
DV       = DV(:);

% Long-form table
tbl = table(SubjID, Within, Between, DV);

% Rename Within levels
withinLevels = categories(Within);
condNames = strcat("Cond", string(1:numel(withinLevels)));
map = containers.Map(withinLevels, condNames);
tbl.Within = categorical(values(map, cellstr(tbl.Within)));

% Wide-form for fitrm
wideDV = unstack(tbl(:, {'SubjID','Within','DV'}), 'DV', 'Within');
[~, idx] = unique(tbl.SubjID);
BetweenPerSubject = tbl.Between(idx);
wideTbl = [wideDV table(BetweenPerSubject)];
wideTbl.Properties.VariableNames{end} = 'Between';

% Formula string
dvNames = wideTbl.Properties.VariableNames(2:(1+numel(condNames)));
formulaStr = sprintf('%s ~ Between', strjoin(dvNames, ','));
WithinTbl = table(categorical(dvNames'), 'VariableNames', {'Within'});

% Fit repeated-measures model
rm = fitrm(wideTbl, formulaStr, 'WithinDesign', WithinTbl);

% ANOVA
stats = ranova(rm, 'WithinModel', 'Within');
betweenStats = anova(rm);

% ==== Partial eta squared ====
SS = stats.SumSq;
partial_eta2 = SS ./ (SS + stats.SumSq(end));  % Using Within error term
stats.partial_eta2 = partial_eta2;


% ==== Pairwise comparisons ====
% -- (A) 1) Within 主效應 (已經有了) -------------------------------
mc_within   = multcompare(rm, 'Within',  'ComparisonType','bonferroni');

% -- (B) 2) Between 主效應 ----------------------------------------
mc_between  = multcompare(rm, 'Between', 'ComparisonType','bonferroni');

% -- (C) 3) Simple-effects：在各 Between 群組內比較 Within 水平 ----
simple_within = cell(numel(categories(Between)),1);
for g = 1:numel(categories(Between))
    cats = categories(Between);
    rows = wideTbl.Between == cats{g};
    rm_g = fitrm(wideTbl(rows,:), sprintf('%s ~ 1', strjoin(dvNames, ',')), 'WithinDesign', WithinTbl);
    simple_within{g} = multcompare(rm_g, 'Within', ...
        'ComparisonType','bonferroni');
end

% -- (D) 4) Simple-effects：在各 Within 條件比 Between 群組 ---------
simple_between = cell(numel(condNames),1);
for c = 1:numel(condNames)
    dvCol  = dvNames{c};
    tbl_c  = wideTbl(:, {'Between', dvCol});
    rm_c   = fitrm(tbl_c, dvCol + " ~ Between");
    simple_between{c} = multcompare(rm_c, 'Between', ...
                         'ComparisonType','bonferroni');
end

% ==== Plot interaction ====
m = grpstats(tbl, {'Within','Between'}, {'mean'}, 'DataVars', 'DV');
if plotflag
    group0 = m.mean_DV(m.Between=='1');
    group1 = m.mean_DV(m.Between=='2');

    figure;
    plot(1:numel(condNames), [group0 group1], '-o', 'LineWidth', 2);
    xticks(1:numel(condNames));
    xticklabels(condNames);
    xlabel('Within condition');
    ylabel('Mean DV');
    legend({'Group 1','Group 2'}, 'Location', 'best');
    title(['Interaction Plot: ' BT ' × ' WT]);
    xlim([.5 2.5]);
end
% Output
result.ranova = stats;
result.anova_between = betweenStats;
result.plot_data = m;
result.partial_eta2 = partial_eta2;
result.posthoc.within_main    = mc_within;
result.posthoc.between_main   = mc_between;
result.posthoc.simple_within  = simple_within;   % cell，每格是各群組
result.posthoc.simple_between = simple_between;  % cell，每格是各條件


fmt = ['Btw: F(%.0f,%.0f) = %.02f, p = %.03f, pes = %.02f \n' ...
       'Wth: F(%.0f,%.0f) = %.02f, p = %.03f, pes = %.02f \n' ...
       'Int: F(%.0f,%.0f) = %.02f, p = %.03f, pes = %.02f \n' ...
       'Btw1: p(bonferroni) = %.03f \n' ...
       'Btw2: p(bonferroni) = %.03f \n'];
f = sprintf(fmt, ...
    stats.DF(2), ...
    stats.DF(3), ...
    stats.F(2), ...
    stats.pValue(2), ...
    stats.partial_eta2(2), ...
    ...
    stats.DF(2), ...
    stats.DF(3), ...
    stats.F(4), ...
    stats.pValue(4), ...
    stats.partial_eta2(4), ...
    ...
    stats.DF(2), ...
    stats.DF(3), ...
    stats.F(5), ...
    stats.pValue(5), ...
    stats.partial_eta2(5), ...
    ...
    result.posthoc.simple_within{1}.pValue(1), ...
    result.posthoc.simple_within{2}.pValue(1));

fprintf(f)
result.f = f;
end
