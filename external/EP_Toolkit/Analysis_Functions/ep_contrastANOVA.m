function ep_contrastANOVA
% ep_contrastANOVA -
% Allows contrasts to be specified for the robust ANOVA module.

%History
%  by Joseph Dien (7/23/09)
%  jdien07@mac.com
%
%  bugfix 8/27/09 JD
%  Location of windows appearing partly off screen on some systems.  Fixed.
%
%  bugfix 1/20/11 JD
%  Fixed levels of ANOVA factors being listed in the wrong order when there are more than one ANOVA factor.
%  Fixed crash when no between factor specified.
%  Fixed error message incorrectly rejecting between or within contrast specified as being "1", which should mean
%  no contrast of that type.
%
%  bugfix 3/24/11 JD
%  Fixed crash when performing contrasts with dataset having no between group factors.
%
%  bugfix 11/1/13 JD
%  Fixes font sizes on Windows.
%
%  bugfix 1/12/14 JD
%  Workaround for Matlab bug periodically causing screen size to register as having zero size.
%
% modified 3/18/14 JD
% Changed uses of "temp" as a variable name to "tempVar" due to other Matlab programmers often using it as a function
% name, resulting in collisions.
%
% modified 12/11/15 JD
% Checks if p-value variability exceeds two standard deviations over the threshold.
% Includes warning summary of %age of singular and nearly singular matrices in the output.
%
% bugfix 6/12/17 JD
% Fixed crash when using contrast function.
% Only presents the between contrast controls if there are more than one between levels.
%
% bugfix 7/27/18 JD
% Fixed crash when running contrast and there is no between group variable.
%
% modified 11/29/19 JD
% Enabled reading of text files with a greater range of variations in character encoding, end-of-line markers, and field separation markers.
% Also a ccommodating .csv suffix.
%
% modified & bugfix 1/19/21 JD
% Added option to use conventional ANOVA statistics.
% No longer crashes when only a subset of the columns were included.
%
% bugfix 4/4/25 JD
% Fixed crash when running contrast, added effect size and meanR calculation to the robust statistics, and other fixes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 1999-2025  Joseph Dien
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global EPmain EPContrast

scrsz = EPmain.scrsz;
windowHeight=scrsz(4);

contrastFigure=findobj('Name', 'Contrast Window');

if ~isempty(contrastFigure)
    close(contrastFigure)
end

EPContrast.handles.figure=figure('Name', 'Contrast Window', 'NumberTitle', 'off', 'Position',[201 1 700 windowHeight]);
colormap jet;

%initial error checking of settings
for i=1:6
    if xor(isempty(EPmain.anova.data.betweenName{i}),(EPmain.anova.data.between(i)>length(EPmain.anova.data.columnNames)))
        msg{1}='All between-group independent variables need to have a three letter label specified.';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    if ~isempty(EPmain.anova.data.factor{i}) && isempty(EPmain.anova.data.levels{i})
        msg{1}='There is a within-group factor with a name but no levels.';
        [msg]=ep_errorMsg(msg);
        return
    end
end

set(EPmain.handles.anova.contrast,'enable','off');
set(EPmain.handles.anova.load,'enable','off');
set(EPmain.handles.anova.view,'enable','off');
set(EPmain.handles.anova.run,'enable','off');
set(EPmain.handles.anova.done,'enable','off');
drawnow

EPContrast.handles.run = uicontrol('Style', 'pushbutton', 'String', 'Run','FontSize',EPmain.fontsize,...
    'Position', [20 windowHeight-130 60 30], 'Callback', @ANOVAcontrast);

EPContrast.handles.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
    'Position', [20 windowHeight-160 60 30], 'Callback', ['close(''Contrast Window'');', 'ep(''start'');']);

%within-group levels table
totalFactorLevels=EPmain.anova.data.rightColumn-EPmain.anova.data.leftColumn+1;

factorLevels=cell(0);
factorNames=cell(0);
for theFactor=1:length(EPmain.anova.data.factor)
    if ~isempty(EPmain.anova.data.factor{theFactor}) && ~isempty(EPmain.anova.data.levels{theFactor})
        factorLevels{end+1}=EPmain.anova.data.levels{theFactor};
        factorNames{end+1}=EPmain.anova.data.factor{theFactor};
    end
end

repFactor=1;
repCount=1;
lvlCount=1;
numFactors=length(factorLevels);
for theFactor=numFactors:-1:1
    for theRow=1:totalFactorLevels
        tableDataW{theRow,theFactor+1}=factorLevels{theFactor}(lvlCount);
        repCount=repCount+1;
        if repCount > repFactor
            repCount=1;
            lvlCount=lvlCount+1;
            if lvlCount > length(factorLevels{theFactor})
                lvlCount=1;
            end
        end
    end
    repFactor=repFactor*length(factorLevels{theFactor});
end

for i=1:totalFactorLevels
    tableDataW{i,1}=EPmain.anova.data.columnNames{EPmain.anova.data.leftColumn+i-1};
    tableDataW{i,2+numFactors}=0;
    tableDataW{i,3+numFactors}=0;
    tableDataW{i,4+numFactors}=0;
    tableDataW{i,5+numFactors}=0;
    tableDataW{i,6+numFactors}=0;
end

tableNamesW{1}='Cells';
for i=1:numFactors
    tableNamesW{i+1}=factorNames{i};
end

tableNamesW{2+numFactors}='Con1';
tableNamesW{3+numFactors}='Con2';
tableNamesW{4+numFactors}='Con3';
tableNamesW{5+numFactors}='Con4';
tableNamesW{6+numFactors}='Con5';

columnEditableW =  [repmat(false,1,1+numFactors) repmat(true,1,5)];

columnWidthW{1}=135;
for i=2:length(tableNamesW)
    columnWidthW{i}=50;
end

try
    EPContrast.handles.withinTable = uitable('Data',tableDataW,'ColumnName',tableNamesW,'FontSize',EPmain.fontsize,...
        'ColumnWidth',columnWidthW,...
        'ColumnEditable', columnEditableW, 'Position',[100 windowHeight-400 569 300]);
catch
    uicontrol('Style','text','HorizontalAlignment','left','String', 'This function does not work with this version of Matlab.','FontSize',EPmain.fontsize,...
        'Position',[100 windowHeight-400 569 300]);
end

%between-group levels table

betweenFactors=[];
totBetweenLevels=1;
betweenFactorLevels=cell(0);
for iFactor=1:6
    if ~isempty(EPmain.anova.data.betweenName{iFactor}) && (EPmain.anova.data.covariate(iFactor)~=2)
        theLevels={EPmain.anova.data.betweenLvl{:,EPmain.anova.data.between(iFactor)-EPmain.anova.data.rightColumn}};
        betweenFactors(end+1)=EPmain.anova.data.between(iFactor);
        betweenFactorLevels{end+1}=unique(theLevels);
        totBetweenLevels=totBetweenLevels*length(unique(theLevels));
    end
end

if totBetweenLevels>1
    repFactor=1;
    repCount=1;
    lvlCount=1;
    numFactors=length(betweenFactors);
    for iFactor=1:length(betweenFactorLevels)
        for iRow=1:totBetweenLevels
            tableDataB{iRow,iFactor}=betweenFactorLevels{iFactor}{lvlCount};
            repCount=repCount+1;
            if repCount > repFactor
                repCount=1;
                lvlCount=lvlCount+1;
                if lvlCount > length(betweenFactorLevels{iFactor})
                    lvlCount=1;
                end
            end
        end
        repFactor=repFactor*length(betweenFactorLevels{iFactor});
    end
    
    for iFactor=1:numFactors
        tableNamesB{iFactor}=EPmain.anova.data.betweenName{iFactor};
    end
    
    for i=1:totBetweenLevels
        tableDataB{i,1+numFactors}=0;
        tableDataB{i,2+numFactors}=0;
        tableDataB{i,3+numFactors}=0;
        tableDataB{i,4+numFactors}=0;
        tableDataB{i,5+numFactors}=0;
    end
    
    tableNamesB{1+numFactors}='Con1';
    tableNamesB{2+numFactors}='Con2';
    tableNamesB{3+numFactors}='Con3';
    tableNamesB{4+numFactors}='Con4';
    tableNamesB{5+numFactors}='Con5';
    
    columnEditableB =  [repmat(false,1,numFactors) repmat(true,1,5)];
    
    for i=1:length(tableNamesB)
        columnWidthB{i}=50;
    end
    
    try
        EPContrast.handles.betweenTable = uitable('Data',tableDataB,'ColumnName',tableNamesB,'FontSize',EPmain.fontsize,...
            'ColumnWidth',columnWidthB,...
            'ColumnEditable', columnEditableB, 'Position',[100 windowHeight-720 569 300]);
    catch
        uicontrol('Style','text','HorizontalAlignment','left','String', 'This function does not work with this version of Matlab.','FontSize',EPmain.fontsize,...
            'Position',[100 windowHeight-720 569 300]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ANOVAcontrast(src,eventdata) %run the ANOVA contrasts
global EPmain EPContrast EPtictoc

set(EPContrast.handles.run,'enable','off');
set(EPContrast.handles.done,'enable','off');

WithinTableData=get(EPContrast.handles.withinTable,'Data');
if isfield(EPContrast.handles,'betweenTable')
    BetweenTableData=get(EPContrast.handles.betweenTable,'Data');
else
    BetweenTableData=[];
end

numContrasts=0;
withinContrasts=cell(5,1);
betweenContrasts=cell(5,1);
withinWeights=cell(5,1);
betweenWeights=cell(5,1);
for i=1:5
    withinCon=[WithinTableData{:,size(WithinTableData,2)-5+i}];
    withinWeights{i}=withinCon;
    withinCon=normalize(withinCon);
    if ~isempty(BetweenTableData)
        betweenCon=[BetweenTableData{:,size(BetweenTableData,2)-5+i}];
        betweenWeights{i}=betweenCon;
        betweenCon=normalize(betweenCon);
    else
        betweenCon=[];
    end
    if any(withinCon) || any(betweenCon)
        numContrasts=numContrasts+1;
        if any(withinCon)
            withinContrasts{i}=withinCon';
            if (sum(withinCon) ~= 0) && (length(find(withinCon ==1)) ~= length(withinCon))
                msg{1}='Error: within contrast does not sum to zero.';
                [msg]=ep_errorMsg(msg);
                set(EPContrast.handles.run,'enable','on');
                set(EPContrast.handles.done,'enable','on');
                return
            end
        else
            withinContrasts{i}=ones(size(WithinTableData,1),1);
        end
        if any(betweenCon)
            betweenContrasts{i}=betweenCon;
            if (sum(betweenCon) ~= 0) && (length(find(betweenCon ==1)) ~= length(betweenCon))
                msg{1}='Error: between contrast does not sum to zero.';
                [msg]=ep_errorMsg(msg);
                set(EPContrast.handles.run,'enable','on');
                set(EPContrast.handles.done,'enable','on');
                return
            end
        else
            betweenContrasts{i}=ones(1,size(BetweenTableData,1));
        end
    end
end

[outFileName, pathname] = uiputfile('*.html','ANOVA output:');
if outFileName == 0
    return %user hit cancel on file requestor
end
outFileName=[pathname outFileName];
if exist(outFileName,'file')
    delete(outFileName); %user must have clicked "yes" to whether to replace existing file
end
[pathstr, fileName, ext] = fileparts(outFileName);
if ~strcmp(ext,'.html')
    ext='.html';
end

outfid=fopen([pathstr filesep fileName ext],'w');

[ANOVAfiles, pathname] = uigetfile({'*.txt';'*.csv'},'Open:','MultiSelect','on');
activeDirectory=pathname;
if ~iscell(ANOVAfiles)
    tempVar=ANOVAfiles;
    ANOVAfiles=[];
    ANOVAfiles{1}=tempVar;
end
if ANOVAfiles{1}==0
    msg{1}='No filenames selected. You have to click on a name';
    [msg]=ep_errorMsg(msg);
    set(EPContrast.handles.run,'enable','on');
    set(EPContrast.handles.done,'enable','on');
    return
end
if ~iscell(ANOVAfiles)
    tempVar=ANOVAfiles;
    ANOVAfiles=[];
    ANOVAfiles{1}=tempVar;
end
for theFile=1:size(ANOVAfiles,2)
    ANOVAfiles{theFile}=[activeDirectory ANOVAfiles{theFile}];
end

numComps=EPmain.anova.numComps;
if numComps==0
    numComps=1;
end
alpha.uncorrected=.05; %threshold for declaring statistical significance.
alpha.corrected=alpha.uncorrected/numComps; %with Bonferroni correction

PER=EPmain.preferences.anova.trimming;
NUMSIM=EPmain.preferences.anova.bootstrap;
REPS=EPmain.preferences.anova.reps;
SEED=EPmain.preferences.anova.seed;
MISSING=EPmain.preferences.anova.missing;
statsToolbox=0;
verList=ver;
if ~isempty(find(strcmp('Statistics and Machine Learning Toolbox',{verList.Name})))
    statsToolbox=1;
end

if isempty(MISSING)
    MISSING=Inf;
    disp('MISSING set to Inf');
end
switch EPmain.anova.method
    case 1 %robust
        ANOVAmethod=1;
        trimOption=1;
    case 2 %conventional
        ANOVAmethod=0;
        trimOption=0;
    otherwise
        error('oops programmer error.')
end

fprintf(outfid,'<font face="Courier">');
ep_tictoc('begin');
for theFile=1:length(ANOVAfiles)
    ep_tictoc;if EPtictoc.stop;ep('start');return;end
    [pathname, fileName, ext] = fileparts(ANOVAfiles{theFile});
    fprintf(outfid,'%s </BR>',[fileName ext]);
    
    [ANOVAheader, ANOVAdata] = ep_loadANOVA(ANOVAfiles{theFile},EPmain.anova.data);
    if isempty(ANOVAdata)
        set(EPContrast.handles.run,'enable','on');
        set(EPContrast.handles.done,'enable','on');
        return
    end
    Y=ANOVAdata.data;
    NX=ANOVAdata.subjects(:)';
    Y(Y==MISSING)=NaN; %the conventional code recognizes NaN but not missing data code number.
    Y=Y(:,EPmain.anova.data.leftColumn:EPmain.anova.data.rightColumn);
    for iContrast=1:numContrasts
        C=betweenContrasts{iContrast};
        U=withinContrasts{iContrast};
        if length(U)~=size(Y,2)
            msg{1}='Number of variables in contrast does not match number of variables in the dataset';
            [msg]=ep_errorMsg(msg);
            set(EPContrast.handles.run,'enable','on');
            set(EPContrast.handles.done,'enable','on');
            return
        end
        subList=zeros(size(Y,1),1);
        for iGroup=1:length(C)
            subList(1+sum(NX(1:iGroup-1)):sum(NX(1:iGroup)))=C(iGroup);
        end
        Y2=Y([subList~=0],[U~=0]);
        C=C(C~=0);
        U=U(U~=0);
        if all(C==1)
            C=1;
            NX=size(Y,1);
        end
        fprintf(outfid,'Contrast: %d </BR>',iContrast);
        fprintf(outfid,'Within: %s </BR>',num2str(withinWeights{iContrast}));
        fprintf(outfid,'Between: %s </BR>',num2str(betweenWeights{iContrast}));
        if any(EPmain.anova.data.covariate==2)
            fprintf(outfid,'Covariates not included.</BR>');
        end
        
        testResults=[];
        if ANOVAmethod
            %robust statistics
            resultsTot=cell(REPS,1);
            pList=zeros(REPS,1);
            for iRep=1:REPS
                ep_tictoc;if EPtictoc.stop;return;end
                [~, ~, RESULTS, ~]=runRobust(Y2, NX, C, U, trimOption, PER, ANOVAmethod, NUMSIM, SEED*iRep, MISSING, 1, alpha, 0, 0, 1);
                resultsTot{iRep}=RESULTS;
                pList(iRep)=RESULTS(4);
            end
            if all(pList==0)
                return
            end
            [~, B]=sort(pList);
            RESULTS=resultsTot{B==ceil(length(pList)/2)};
            if isempty(RESULTS)
                return
            end
            alpha.pvar=std(pList)*2;
            testResults.statsType='robust';
            testResults.GPOWERtest='';
            testResults.stat=RESULTS(1);
            testResults.numDF=RESULTS(2);
            testResults.denDF=RESULTS(3);
            testResults.pValue=RESULTS(4);
            testResults.epsilon=RESULTS(9);
            testResults.statsTest='N-way';
            testResults.GPOWER=RESULTS(10);
            %effect size as only the matlab code computes effect sizes and also MSE.
            [~, STDIZER, esRESULTS, MSE]=ep_WJGLMml_mat(Y2, NX, C, U, trimOption, PER, ANOVAmethod, NUMSIM, SEED, MISSING, 1, alpha, 0, 0, 1);
            if isempty(MSE)
                msg='Error: ANOVA run aborted.';
                [msg]=ep_errorMsg(msg);
                set(EPContrast.handles.run,'enable','on');
                set(EPContrast.handles.done,'enable','on');
                return
            end
            testResults.MSE=MSE;
            testResults.effectSize=esRESULTS(5);
            testResults.lowerConfid=esRESULTS(6);
            testResults.upperConfid=esRESULTS(7);
            testResults.effectSizeTest='d<sub>R</sub>';
            Rmatrix=pinv(diag(sqrt(diag(STDIZER))))*STDIZER*pinv(diag(sqrt(diag(STDIZER))));
            testResults.meanR=mean(Rmatrix(~tril(ones(size(Rmatrix)))));
        else
            %conventional statistics
            if statsToolbox
                if any(C) && any(U) %split-plot
                    newY=[];
                    newY(:,1)=sum(Y*diag(abs(U.*(U>0))),2);
                    newY(:,2)=sum(Y*diag(abs(U.*(U<0))),2);
                    
                    newY=diag(abs(subList))*newY;
                    newY=newY(subList~=0,:);
                    dataTable=array2table(newY);
                    withinMeas = cell2table({'A';'B'});
                    withinMeas.Properties.VariableNames{'Var1'}='within';
                    
                    dataTable=[dataTable array2table(sign(subList(subList~=0)))];
                    dataTable.Var1=categorical(dataTable.Var1);
                    dataTable.Properties.VariableNames{'Var1'}='between';
                    rm = fitrm(dataTable,'newY1-newY2 ~ between','WithinDesign',withinMeas);
                    ranovatbl = ranova(rm,'WithinModel','within');
                    theEffect='between:within';
                    theError='Error(within)';
                    testResults.pValue=table2array(ranovatbl(theEffect,'pValue'));
                    testResults.stat=table2array(ranovatbl(theEffect,'F'));
                    testResults.numDF=table2array(ranovatbl(theEffect,'DF'));
                    testResults.denDF=table2array(ranovatbl(theError,'DF'));
                    testResults.MSE=table2array(ranovatbl(theError,'MeanSq'));
                    %partial eta-squared according to Lakens(2013)
                    SSeffect=table2array(ranovatbl(theEffect,'SumSq'));
                    SSerror=table2array(ranovatbl(theError,'SumSq'));
                    partEtaSq=SSeffect/(SSeffect+SSerror);

                    %conversion to SPSS partial eta-squared derived from Lakens(2013) appendix
                    %     covMat=table2array(theRM.Covariance);
                    %     corMat=diag(diag(sqrt(covMat.^-1)))*covMat*diag(diag(sqrt(covMat.^-1)));
                    corMat=corrcoef(newY,'Rows','complete'); %drops any observations with missing data.
                    rho=mean(tril(corMat,-1),'all'); %average of all the pairwise correlations
                    sampSize=size(newY,1);
                    numGroups=2;
                    numReps=2;
                    partEtaSqSPSS=(partEtaSq*sampSize*numReps/(2*(1-partEtaSq)*(sampSize-numGroups)*(numReps-1)*(1-rho)));
                    testResults.GPOWER(1)=sqrt(partEtaSq/(1-partEtaSq));
                    testResults.meanR=rho;
                    testResults.effectSize=partEtaSqSPSS;
                    testResults.effectSizeTest='partEtaSqSPSS';
                    %testResults.GPOWER=sqrt(partEtaSqSPSS/(1+partEtaSqSPSS)); %the SPSS version of the f statistic that G*POWER needs.
                    testResults.GPOWERtest='G*POWER f';
                    testResults.statsTest='N-way';
                    testResults.lowerConfid=NaN;
                    testResults.upperConfid=NaN;
                    testResults.statsType='conventional';
                    testResults.epsilon=1; %epsilon always 1 for one degree contrast
                elseif any(U) %within
                    X1=sum(Y.*diag(abs(U.*(U>0))));
                    X2=sum(Y.*diag(abs(U.*(U<0))));
                    
                    [~,p,ci,stats] = ttest(X1,X2); %paired t-test
                    %[RESULTS, MSE] = runANOVA(mainAnovaTable, mainRM, factorGroupTable, factorTable, betweenGroup, withinGroup, contrastGroupName, contrastName, contrast, epsilonType);
                    %RESULTS(1)=mean(Y*orth(U))^2/MSE;
                    testResults=[];
                    testResults.stat=stats.tstat;
                    testResults.numDF=1; %numerator df always 1 for one degree contrast
                    testResults.denDF=stats.df;
                    testResults.pValue=p;
                    testResults.epsilon=1; %epsilon always 1 for one degree contrast
                    testResults.MSE=stats.sd/stats.df;
                    testResults.statsTest='paired samples t-test';
                    testResults.statsType='conventional';
                    testResults.GPOWER=stats.tstat/sqrt(size(Y,1));
                    testResults.GPOWERtest='dz';
                    testResults.lowerConfid=ci(1);
                    testResults.upperConfid=ci(2);
                    X1m=mean(X1,'omitnan');
                    X2m=mean(X2,'omitnan');
                    X1std=std(X1,'omitnan');
                    X2std=std(X2,'omitnan');
                    dz=mean(X1-X2)/std(X1-X2,'omitnan');
                    rho=corrcoef([X1 X2],'Rows','complete');
                    drm=dz*sqrt(2*(1-rho(1,2)));
                    hgrm=drm*(1-(3/(4*(length(X1)*2)-9)));
                    dav=(X1m-X2m)/((X1std+X2std)/2);
                    hgav=dav*(1-(3/(4*(length(X1)*2)-9)));
                    if (abs((abs((X1m-X2m)/(sqrt((((length(X1)-1)*X1std^2)+((length(X1)-1)*X2std^2))/(length(X1)+length(X1)-2)))))-drm)...
                            <abs((abs((X1m-X2m)/(sqrt((((length(X1)-1)*X1std^2)+((length(X1)-1)*X2std^2))/(length(X1)+length(X1)-2)))))-dav))
                        testResults.effectSize=hgrm;
                        testResults.effectSizeTest='Hedges grm';
                    else
                        testResults.effectSize=hgav;
                        testResults.effectSizeTest='Hedges gav';
                    end
                    testResults.meanR=[];
                elseif any(C) %between
                    Y=diag(abs(subList))*Y;
                    X1=mean(Y(subList>0,:),2);
                    X2=mean(Y(subList<0,:),2);
                    
                    [~,p,ci,stats] = ttest2(X1,X2); %two-sample t-test with homogenous variance (since conventional, not robust)
                    
                    if postHocResults(find(all(levelPairs(iPairs,:)==postHocResults(:,1:2),2)),6) <= postHocAlpha.uncorrected
                        postHocAlpha.corrected=inf; %significant after post-hoc
                    else
                        postHocAlpha.corrected=-inf; %not significant after post-hoc
                    end
                    
                    testResults=[];
                    testResults.stat=stats.tstat;
                    testResults.numDF=1; %numerator df always 1 for one degree contrast
                    testResults.denDF=stats.df;
                    testResults.pValue=p;
                    testResults.epsilon=1; %epsilon always 1 for one degree contrast
                    testResults.MSE=stats.sd/stats.df;
                    testResults.statsTest='independent samples t-test';
                    testResults.statsType='conventional';
                    %ds=abs((mean(X1)-mean(X2))/(sqrt((((length(X1)-1)*var(X1))+((length(X2)-1)*var(X2)))/(length(X1)+length(X2)-2))));
                    ds=stats.tstat*sqrt(2/size(Y,1));
                    testResults.GPOWER=ds;
                    testResults.GPOWERtest='ds';
                    testResults.lowerConfid=ci(1);
                    testResults.upperConfid=ci(2);
                    hgs=ds*(1-(3/(4*(length(X1)+length(X2))-9)));
                    testResults.effectSize=hgs;
                    testResults.effectSizeTest='Hedges gs';
                    testResults.meanR=[];
                end
            else
                %handrolled conventional ANOVAs
                disp('Not functional yet');
                return
                %                     [MUHAT, SIGMA, RESULTS, MSE]=ep_WJGLMml_mat(Y, NX, groupContrast{betweenGroup}, contrast{withinGroup}', trimOption, PER, ANOVAmethod, NUMSIM, SEED, MISSING, OPT3, alpha, SCALE, LOC1, LOC2);
            end
            alpha.pvar=0;
        end
        ep_printADFResults(testResults, alpha, 0, outfid, EPmain.preferences.anova);
    end
end
ep_tictoc('end');
set(EPContrast.handles.run,'enable','on');
set(EPContrast.handles.done,'enable','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MUHAT, SIGMA, RESULTS, MSE]=runRobust(Y, NX, groupContrast, contrast, trimOption, PER, ANOVAmethod, NUMSIM, SEED, MISSING, OPT3, alpha, SCALE, LOC1, LOC2)
Y(isnan(Y))=MISSING; %the robust code recognizes missing data code number but not NaN.
try
    [MUHAT, SIGMA, RESULTS, MSE]=ep_WJGLMml(Y, NX, groupContrast, contrast, trimOption, PER, ANOVAmethod, NUMSIM, SEED, MISSING, 0, alpha, SCALE, LOC1, LOC2);
catch ME
    disp(ME.message);
    disp('Mex function failed - attempting Matlab function.');
    [MUHAT, SIGMA, RESULTS, MSE]=ep_WJGLMml_mat(Y, NX, groupContrast, contrast, trimOption, PER, ANOVAmethod, NUMSIM, SEED, MISSING, 0, alpha, SCALE, LOC1, LOC2);
end
RESULTS(10)=NaN;
