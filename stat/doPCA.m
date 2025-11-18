function doPCA(Data,dv,tw)
%%
Data.(dv) = Data.(dv)(:,tw(1)-tw(3):tw(2)-tw(3));
numComp = 2;
[F] = ep_doPCA('temp', 'Promax', 3, 'SVD', 'COV', numComp, Data.(dv), 'K');

%
% Plot comp * SD
figure; hold on
ls = {'-','--'};
for nComp = 1:numComp
    d = F.FacPat(:,nComp).*F.varSD;
    plot(d,'Color','k','LineStyle',ls{nComp})
end
yline(0)
plotLed({'— Comp 1','- - Comp 2'}, {'k';'k'},'t')
xlabel('Time from face onset (ms)','FontSize',12);
ylabel('Loadings * SD','FontSize',12);

%
figure; hold on
x = tw(1):1:tw(2);

d = Data.(dv);
err = std(d)/sqrt(length(d(:,1))-1);
y = mean(d);
plot(x,y,'linewidth',1,'Color',[0 .5 .8]);
fill([x';flip(x')],[y-err,flip(y+err)], ...
    [0 .5 .8], ...
    'FaceAlpha',0.2, ...
    'EdgeColor','none')

d = d.* F.FacPat(:,1)' + d.* F.FacPat(:,2)';
err = std(d)/sqrt(length(d(:,1))-1);
y = mean(d);
plot(x,y,'linewidth',1,'Color',[.8 .1 0]);
fill([x';flip(x')],[y-err,flip(y+err)], ...
    [.8 .1 0], ...
    'FaceAlpha',0.2, ...
    'EdgeColor','none')
yline(0)
xlabel('Time from face onset (ms)','FontSize',12);
ylabel('Pupil size (mm)','FontSize',12);
plotLed({'Oringinal data','Comp 1+2'},{[0 .5 .8];[.8 .1 0]},'b')
%% Plot comp 1

temp = Data;
temp.(dv) = temp.(dv) .*F.FacPat(:,1)';
plotDynamic(temp, ...
    'tw', [1 4000 1], ...
    'col', [[0 .5 .8];[.8 .1 0];[0 .5 .8];[.8 .1 0]], ...
    'lnstyle', {'-','-','--','--'})
yline(0,'HandleVisibility','off')
title('Comp 1')

%% Plot comp 2

temp = Data;
temp.(dv) = temp.(dv) .*F.FacPat(:,2)';
plotDynamic(temp, ...
    'tw', [1 4000 1], ...
    'col', [[0 .5 .8];[.8 .1 0];[0 .5 .8];[.8 .1 0]], ...
    'lnstyle', {'-','-','--','--'})
yline(0,'HandleVisibility','off')
title('Comp 2')

%% Plot score 1

temp = Data;
temp.(dv) = F.FacScr(:,1);
vName = Data.Properties.VariableNames;
if length(vName) == 5
    plotBox(temp, ...
        'col', [[0 .5 .8];[.8 .1 0];[0 .5 .8];[.8 .1 0]], ...
        'dot', true, ...
        'ylab', 'Score');
    fprintf('Componemnt 1 score: \n')
    result = mixed_anova(temp, ...
        dv, ...
        vName{3}, ...
        vName{4});
elseif length(vName) == 4
    plotBox(temp, ...
    'col', [[0 .5 .8];[.8 .1 0]], ...
    'dot', true, ...
    'ylab', 'Score');
end
title('Comp 1')

%% Plot score 2

temp = Data;
temp.(dv) = F.FacScr(:,2);
vName = Data.Properties.VariableNames;
if length(vName) == 5
    plotBox(temp, ...
        'col', [[0 .5 .8];[.8 .1 0];[0 .5 .8];[.8 .1 0]], ...
        'dot', true, ...
        'ylab', 'Score');
    fprintf('Componemnt 2 score: \n')
    result = mixed_anova(temp, ...
        dv, ...
        vName{3}, ...
        vName{4});
elseif length(vName) == 4
    plotBox(temp, ...
    'col', [[0 .5 .8];[.8 .1 0]], ...
    'dot', true, ...
    'ylab', 'Score');
    % [df,tval,pval,cd,f] = doTTest(op,DV,IV,IVValues,false,'both');
end
title('Comp 2')








