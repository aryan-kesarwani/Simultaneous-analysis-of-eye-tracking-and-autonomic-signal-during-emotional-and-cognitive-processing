function ep_ANOVApane(varargin)
% ep_ANOVApane - function ep_ANOVApane(varargin);
% GUI interface for the ANOVA Pane.
%

%History
%  by Joseph Dien (7/23/19)
%  jdien07@mac.com
%
%
% modified 4/24/20 JD
% Added preference buttons to panes.
%
% modified 8/31/20 JD
% Added support for option to run all posthocs even if the effect is not significant.
%
% modified 1/18/21 JD
% Added option to use conventional ANOVA statistics and to specify covariates.
%
% bugfix 12/5/21 JD
% Fixed crash when trying to handle no-header file with between-group variables.
%
% bugfix 5/5/22 JD
% Fixed crash when setting a Between Group factor to "none."
%
% modified 7/27/22 JD
% Added ability to load and save ANOVA table settings and auto ANOVA table button.
%
% bugfix 8/14/22 JD
% Fixed not updating "have" value after using the Auto button in the ANOVA function.
%
% bugfix 6/27/24 JD
% Fixed auto option not working for multi-factor structures where the first letter represents a factor.
% Fixed crash when saving ANOVA table name.
%
% bugfix 8/12/24 JD
% Fixed crash when loading ANOVA table.
%
% bugfix 9/7/24 JD
% Fixed not updating "have" value after using the Load ANOVAtable button.
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

global EPmain

set(EPmain.handles.hMainWindow,'Name', 'Robust ANOVA');

if isempty(EPmain.anova.data)
    EPmain.anova.data.leftColumn=1;
    EPmain.anova.data.rightColumn=1;
    EPmain.anova.data.columnNames{1}='none';
    EPmain.anova.data.totalLevels=0;
    EPmain.anova.numComps=0;
    EPmain.anova.allPosthoc=0;
    EPmain.anova.data.data=[];
    EPmain.anova.data.areaNames=cell(0);
%     EPmain.anova.error=EPmain.errorList{1};
    
    for iFactor=1:6
        EPmain.anova.data.between(iFactor)=2;
        EPmain.anova.data.factor{iFactor}='';
        EPmain.anova.data.levels{iFactor}='';
        EPmain.anova.data.betweenName{iFactor}='';
        EPmain.anova.data.covariate(iFactor)=1;
    end
end

betweenNames=EPmain.anova.data.columnNames;
betweenNames{end+1}='none';
numCols=size(EPmain.anova.data.data,2);
numSpecs=length(find(strcmp('spec',EPmain.anova.data.areaNames)));


EPmain.handles.anova.prefs = uicontrol('Style', 'pushbutton', 'String', '•','FontSize',EPmain.fontsize,...
    'Position', [5 485 10 10], 'Callback', ['global EPmain;','EPmain.prefReturn=''ANOVA'';','EPmain.mode=''preferenceANOVA'';','ep(''start'');']);

uicontrol('Style','text',...
    'String','Data Columns','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[25 470 100 20]);

if isfield(EPmain.anova.data,'name')
    theColNames=EPmain.anova.data.columnNames(1:numCols);
    EPmain.handles.anova.leftColumn = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',theColNames,...
        'Value',EPmain.anova.data.leftColumn,'Position',[20 450 150 20],...
        'Callback', ['global EPmain;','tempVar=get(EPmain.handles.anova.leftColumn,''Value'');','if tempVar ~=0,EPmain.anova.data.leftColumn=tempVar;end;','if isempty(tempVar),EPmain.anova.data.leftColumn=tempVar;end;','ep(''start'');']);
    
    EPmain.handles.anova.rightColumn = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',theColNames,...
        'Value',EPmain.anova.data.rightColumn,'Position',[20 430 150 20],...
        'Callback', ['global EPmain;','tempVar=get(EPmain.handles.anova.rightColumn,''Value'');','if tempVar ~=0,EPmain.anova.data.rightColumn=tempVar;end;','if isempty(tempVar),EPmain.anova.data.rightColumn=tempVar;end;','ep(''start'');']);
else
    theColNames=EPmain.anova.data.columnNames;
    EPmain.handles.anova.leftColumn = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',theColNames,...
        'Value',EPmain.anova.data.leftColumn,'Position',[20 450 150 20],...
        'Callback', ['global EPmain;','tempVar=get(EPmain.handles.anova.leftColumn,''Value'');','if tempVar ~=0,EPmain.anova.data.leftColumn=tempVar;end;','if isempty(tempVar),EPmain.anova.data.leftColumn=tempVar;end;','ep(''start'');']);
    
    EPmain.handles.anova.rightColumn = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',theColNames,...
        'Value',EPmain.anova.data.rightColumn,'Position',[20 430 150 20],...
        'Callback', ['global EPmain;','tempVar=get(EPmain.handles.anova.rightColumn,''Value'');','if tempVar ~=0,EPmain.anova.data.rightColumn=tempVar;end;','if isempty(tempVar),EPmain.anova.data.rightColumn=tempVar;end;','ep(''start'');']);
    set(EPmain.handles.anova.leftColumn,'enable','off');
    set(EPmain.handles.anova.rightColumn,'enable','off');
end

uicontrol('Style','text','FontSize',EPmain.fontsize,...
    'String',sprintf('have: %d',EPmain.anova.data.totalLevels),'HorizontalAlignment','left',...
    'Position',[25 410 100 20]);

if isfield(EPmain.anova.data,'name')
    needLevels=EPmain.anova.data.rightColumn-EPmain.anova.data.leftColumn+1;
    if needLevels == 1
        needLevels=0; %no within factors
    end
else
    needLevels=0;
end

uicontrol('Style','text','FontSize',EPmain.fontsize,...
    'String',sprintf('need: %d',needLevels),'HorizontalAlignment','left',...
    'Position',[100 410 100 20]);

uicontrol('Style','text','FontSize',EPmain.fontsize,...
    'String','Between Group','HorizontalAlignment','left',...
    'Position',[25 380 100 20]);

for iFactor=1:6
    
    betweenRange=[setdiff([1:length(betweenNames)-1],[1:EPmain.anova.data.rightColumn,EPmain.anova.data.between(setdiff(1:6,iFactor))]) length(betweenNames)];
    %range of values betweenGroup can take.  Ones already taken are not an option.  "none" is always an option.
    
    if ~ismember(EPmain.anova.data.between(iFactor),betweenRange)
        EPmain.anova.data.between(iFactor)=length(betweenNames); %the last in the list of between names is "none"
    end
    
    betweenList=betweenNames(betweenRange);
    EPmain.handles.anova.between(iFactor) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',betweenList,...
        'Value',find(EPmain.anova.data.between(iFactor)==betweenRange),'Position',[5 380-iFactor*20 120 20],...
        'Callback', @betweenGroup);
    
    if ~isfield(EPmain.anova.data,'name')
        set(EPmain.handles.anova.between(iFactor),'enable','off');
    end
    
    if (EPmain.anova.data.between(iFactor) == length(betweenNames))
        EPmain.anova.data.betweenName{iFactor}='';
    end
    
    EPmain.handles.anova.betweenName(iFactor)= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
        'String',EPmain.anova.data.betweenName{iFactor},...
        'Position',[115 380-iFactor*20 40 20],...
        'Callback', @betweenName);
    
    if ~isfield(EPmain.anova.data,'name') || (EPmain.anova.data.between(iFactor) == length(betweenNames))
        set(EPmain.handles.anova.betweenName(iFactor),'enable','off');
    end
    
    if EPmain.anova.method==2
        EPmain.handles.anova.covariate(iFactor) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'F';'C'},...
            'Value',EPmain.anova.data.covariate(iFactor),'Position',[150 380-iFactor*20 55 20],...
            'Callback', ['global EPmain;','tempVar=get(EPmain.handles.anova.covariate(' num2str(iFactor) '),''Value'');','if tempVar ~=0,EPmain.anova.data.covariate(' num2str(iFactor) ')=tempVar;end;','if isempty(tempVar),EPmain.anova.data.covariate(' num2str(iFactor) ')=tempVar;end;','ep(''start'');']);
    end
end

uicontrol('Style','text','FontSize',EPmain.fontsize,...
    'String','Factor','HorizontalAlignment','left',...
    'Position',[20 230 50 20]);

uicontrol('Style','text','FontSize',EPmain.fontsize,...
    'String','Levels','HorizontalAlignment','left',...
    'Position',[80 230 50 20]);

uicontrol('Style','frame',...
    'Position',[15 105 148 130]);

for iFactor=1:6
    EPmain.handles.anova.factor(iFactor)= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
        'String',EPmain.anova.data.factor{iFactor},...
        'Position',[20 230-iFactor*20 50 20],...
        'Callback', @ANOVATable);
    
    if ~isfield(EPmain.anova.data,'name')
        set(EPmain.handles.anova.factor(iFactor),'enable','off');
    end
    
    EPmain.handles.anova.levels(iFactor)= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
        'String',EPmain.anova.data.levels{iFactor},...
        'Position',[80 230-iFactor*20 77 20],...
        'Callback', @ANOVATable);
    
    if ~isfield(EPmain.anova.data,'name')
        set(EPmain.handles.anova.levels(iFactor),'enable','off');
    end    
end

EPmain.handles.anova.loadANOVAtable = uicontrol('Style', 'pushbutton', 'String', 'Load','FontSize',EPmain.fontsize,...
    'Position', [165 210 40 20], 'Callback', @loadANOVAtable);

EPmain.handles.anova.saveANOVAtable = uicontrol('Style', 'pushbutton', 'String', 'Save','FontSize',EPmain.fontsize,...
    'Position', [165 190 40 20], 'Callback', @saveANOVAtable);

EPmain.handles.anova.autoANOVAtable = uicontrol('Style', 'pushbutton', 'String', 'Auto','FontSize',EPmain.fontsize,...
    'Position', [165 170 40 20], 'Callback',@autoANOVAtable);

EPmain.handles.anova.clearANOVAtable = uicontrol('Style', 'pushbutton', 'String', 'Clear','FontSize',EPmain.fontsize,...
    'Position', [165 150 40 20], 'Callback', @clearANOVAtable);

EPmain.handles.anova.undoANOVAtable = uicontrol('Style', 'pushbutton', 'String', 'Undo','FontSize',EPmain.fontsize,...
    'Position', [165 130 40 20], 'Callback', @undoANOVAtable);

if isempty(EPmain.anova.undo)
    set(EPmain.handles.anova.undoANOVAtable,'enable','off');
end

uicontrol('Style','text','FontSize',EPmain.fontsize,...
    'String','MCP','HorizontalAlignment','left',...
    'Position',[5 60 70 20]);

EPmain.handles.anova.numComps= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
    'String',EPmain.anova.numComps,...
    'Position',[40 60 30 20],...
    'Callback', ['global EPmain;','EPmain.anova.numComps=str2num(get(EPmain.handles.anova.numComps,''String''));']);

% EPmain.handles.anova.error = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
%     'String',EPmain.errorList,...
%     'Value',find(strcmp(EPmain.anova.error,EPmain.errorList)),'Position',[5 80 90 20],...
%     'Callback', ['global EPmain;','tempVar=get(EPmain.handles.anova.error,''Value'');','if tempVar ~=0,EPmain.anova.error=tempVar;end;','if isempty(tempVar),EPmain.anova.data.error=EPmain.errorList{tempVar};end;','ep(''start'');']);
% 
verList=ver;
if ~isempty(find(strcmp('Statistics and Machine Learning Toolbox',{verList.Name})))
    EPmain.handles.anova.method = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',{'Robust';'Conventional'},...
        'Value',EPmain.anova.method,'Position',[100 80 100 20],...
        'Callback', ['global EPmain;','tempVar=get(EPmain.handles.anova.method,''Value'');','if tempVar ~=0,EPmain.anova.method=tempVar;end;','if isempty(tempVar),EPmain.anova.data.method=tempVar;end;','ep(''start'');']);
else
    EPmain.handles.anova.method = uicontrol('Style','text','FontSize',EPmain.fontsize,...
        'String','Robust','Position',[100 80 100 20]);
    EPmain.anova.method=1;
end

EPmain.handles.anova.allPosthoc= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
    'String','All Posthoc',...
    'CallBack',['global EPmain;','EPmain.anova.allPosthoc=get(EPmain.handles.anova.allPosthoc,''Value'');','ep(''start'');'],'FontSize',EPmain.fontsize,...
    'Value',EPmain.anova.allPosthoc,'Position',[100 60 100 20],'TooltipString','Check to run all post-hocs even when effect was not significant.');

EPmain.handles.anova.contrast = uicontrol('Style', 'pushbutton', 'String', 'Contrast','FontSize',EPmain.fontsize,...
    'Position', [60 30 60 30], 'Callback', 'ep_contrastANOVA');

EPmain.handles.anova.view = uicontrol('Style', 'pushbutton', 'String', 'View','FontSize',EPmain.fontsize,...
    'Position', [5 30 60 30], 'Callback', 'ep_viewANOVA');

EPmain.handles.anova.load = uicontrol('Style', 'pushbutton', 'String', 'Load','FontSize',EPmain.fontsize,...
    'Position', [5 0 60 30], 'Callback', 'ep(''loadANOVA'')');

EPmain.handles.anova.run = uicontrol('Style', 'pushbutton', 'String', 'Run','FontSize',EPmain.fontsize,...
    'Position', [60 0 60 30], 'Callback', 'ep(''runANOVA'')');

if  (EPmain.anova.data.totalLevels==1) || ((EPmain.anova.data.totalLevels ~= (EPmain.anova.data.rightColumn-EPmain.anova.data.leftColumn+1)) && ~(EPmain.anova.data.totalLevels < 2 && length(char([EPmain.anova.data.betweenName])) > 0))
    set(EPmain.handles.anova.run,'enable','off');
    set(EPmain.handles.anova.view,'enable','off');
end

EPmain.handles.anova.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
    'Position', [120 0 60 30], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ANOVATable(src,eventdata)
%process entry into ANOVA factor table.

global EPmain

theFactor=find(src == EPmain.handles.anova.factor);
theLevel=find(src == EPmain.handles.anova.levels);

tempVar=[];

if theFactor
    theEntry=get(EPmain.handles.anova.factor(theFactor),'string');
    if length(theEntry)==3
        tempVar=EPmain.anova.data;
        EPmain.anova.data.factor{theFactor}=theEntry;
    else
        EPmain.anova.data.factor{theFactor}=[];
    end
elseif theLevel
    theEntry=get(EPmain.handles.anova.levels(theLevel),'string');
    EPmain.anova.data.levels{theLevel}=theEntry;
    tempVar=EPmain.anova.data;
end

ANOVAlevels=0;
for theFactor=1:6
    if ~isempty(EPmain.anova.data.levels{theFactor}) && ~isempty(EPmain.anova.data.factor{theFactor})
        if ANOVAlevels==0
            ANOVAlevels=1;
        end
        ANOVAlevels=ANOVAlevels*length(EPmain.anova.data.levels{theFactor});
    end
end

EPmain.anova.data.totalLevels=ANOVAlevels;

if ~isempty(tempVar)
    EPmain.anova.undo.factor=tempVar.factor;
    EPmain.anova.undo.levels=tempVar.levels;
end

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function betweenName(src,eventdata)
%process entry of between factor name into ANOVA pane

global EPmain

theFactor=find(src == EPmain.handles.anova.betweenName);
if theFactor
    theEntry=get(EPmain.handles.anova.betweenName(theFactor),'string');
    if length(theEntry)==3
        EPmain.anova.data.betweenName{theFactor}=theEntry;
    else
        EPmain.anova.data.betweenName{theFactor}='';
    end
end

ep('start');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function betweenGroup(src,eventdata)
%process selection of between group level column in ANOVA pane

global EPmain

theFactor=find(src == EPmain.handles.anova.between);
if ~isempty(theFactor)
    tempVar=get(EPmain.handles.anova.between(theFactor),'value');
    if tempVar ~=0
        betweenNames=EPmain.anova.data.columnNames;
        betweenNames{end+1}='none';
        betweenRange=[setdiff([1:length(betweenNames)],[1:EPmain.anova.data.rightColumn,EPmain.anova.data.between(setdiff(1:6,theFactor))]) length(betweenNames)];
        EPmain.anova.data.between(theFactor)=betweenRange(tempVar);
        if betweenRange(tempVar) < length(betweenNames) %if not "none"
            for iRow=1:size(EPmain.anova.data.data,1)
                EPmain.anova.data.betweenLvl{iRow,theFactor} = EPmain.anova.data.origData{iRow,betweenRange(tempVar)};
            end
        end
    end
    if isempty(tempVar)
        betweenNames=EPmain.anova.data.columnNames;
        betweenNames{end+1}='none';
        betweenRange=[setdiff([1:length(betweenNames)],[1:EPmain.anova.data.rightColumn,EPmain.anova.data.between(setdiff(1:6,theFactor))]) length(betweenNames)];
        EPmain.anova.data.between(theFactor)=betweenRange(tempVar);
    end
end


ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadANOVAtable(src,eventdata)
%load settings in ANOVA table.

global EPmain

[FileName,PathName,FilterIndex] = uigetfile('','Load ANOVA Table');
if FileName ~= 0
    try
        eval(['tempVar=load( ''' PathName FileName ''');']);
    catch
        tempVar=[];
    end
    if isfield(tempVar,'EPanovaTable')
        EPmain.anova.undo.factor=EPmain.anova.data.factor;
        EPmain.anova.undo.levels=EPmain.anova.data.levels;
        EPmain.anova.data.factor=tempVar.EPanovaTable.factor;
        EPmain.anova.data.levels=tempVar.EPanovaTable.levels;

        ANOVAlevels=0;
        for theFactor=1:6
            if ~isempty(EPmain.anova.data.levels{theFactor}) && ~isempty(EPmain.anova.data.factor{theFactor})
                if ANOVAlevels==0
                    ANOVAlevels=1;
                end
                ANOVAlevels=ANOVAlevels*length(EPmain.anova.data.levels{theFactor});
            end
        end

        EPmain.anova.data.totalLevels=ANOVAlevels;

    else
        msg{1}='Error: not an ANOVA table file.';
        [msg]=ep_errorMsg(msg);
    end
    ep('start');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveANOVAtable(src,eventdata)
%save settings in ANOVA table.

global EPmain

nameStart=strfind(EPmain.anova.data.name,': ');
if isempty(nameStart)
    nameStart=1;
else
    if length(nameStart) > 1
        nameStart=nameStart(1);
    end
    nameStart=nameStart+2;
end
datasetName=EPmain.anova.data.name(nameStart:end);

[fileName,PathName,FilterIndex] = uiputfile('*.mat','Save ANOVA Table',[datasetName '_ANOVAtable']);

if fileName ~= 0
    EPanovaTable.factor=EPmain.anova.data.factor;
    EPanovaTable.levels=EPmain.anova.data.levels;
    eval(['save ''' PathName fileName ''' EPanovaTable']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clearANOVAtable(src,eventdata)
%clears settings in ANOVA table.

global EPmain

EPmain.anova.undo.factor=EPmain.anova.data.factor;
EPmain.anova.undo.levels=EPmain.anova.data.levels;
for iFactor=1:6
    EPmain.anova.data.factor{iFactor}='';
    EPmain.anova.data.levels{iFactor}='';
end

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function undoANOVAtable(src,eventdata)
%undo settings in ANOVA table.

global EPmain

tempVar=EPmain.anova.data;
EPmain.anova.data.factor=EPmain.anova.undo.factor;
EPmain.anova.data.levels=EPmain.anova.undo.levels;
EPmain.anova.undo.factor=tempVar.factor;
EPmain.anova.undo.levels=tempVar.levels;

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function autoANOVAtable(src,eventdata)
%clears settings in ANOVA table.

global EPmain

ep_tictoc('begin');

tempVar=EPmain.anova.data;
condNames=tempVar.columnNames(tempVar.leftColumn:tempVar.rightColumn);
if strcmp(condNames,'none')
    return
end

for iName=1:length(condNames)
    condNames{iName}=strtok(condNames{iName},'-');
end

for iFactor=1:6
    tempVar.factor{iFactor}='';
    tempVar.levels{iFactor}='';
end

if (length(condNames)==1) || any(diff(cellfun(@length,condNames)))
    %names are not all the same length so just one factor
    if length(condNames) <= 26
        tempVar.factor{1}='AF1';
        for iLevel=1:length(condNames)
            tempVar.levels{1}(iLevel)=condNames{iLevel}(1);
        end
        if length(tempVar.levels{1}) ~= length(unique(tempVar.levels{1}))
            for iLevel=1:length(condNames)
                tempVar.levels{1}(iLevel)=96+iLevel;
            end
        end
    else
        ep_tictoc('end');
        return
    end
else
    numLetters=length(condNames{1});
    noDiff=ones(numLetters,1);
    for iLetter=1:numLetters
        for iName=1:length(condNames)
            if ~strcmp(condNames{1}(iLetter),condNames{iName}(iLetter))
                noDiff(iLetter)=0;
                continue
            end
        end
    end
    if all(noDiff)
        disp('All the names are the same.')
        return
    end
    for iName=1:length(condNames)
        condNames{iName}=condNames{iName}(~noDiff); %keep only the part of the names that changes
    end
    
    numFactors=length(condNames{1});
    if numFactors > 6
        disp('Maximum of six ANOVA factors.')
        return
    end
    condNamesChar=char(condNames);
    levelList=cell(numFactors,1);
    for iGrand=1:length(condNames)
        indexList{iGrand}=iGrand;
    end
    comboCount=1;
    for iFac=1:numFactors
        levelList{iFac}=unique(condNamesChar(:,iFac)','stable');
        comboCount=comboCount*length(levelList{iFac});
    end
    if comboCount ==  length(condNames)
        %there are enough conditions that it could be fully crossed
        for iFac=1:numFactors
            tempVar.factor{iFac}=['AF' num2str(iFac)];
            tempVar.levels{iFac}=levelList{iFac};
        end
    else
        %not fully crossed
        if length(condNames) <= 26
            tempVar.factor{1}='AF1';
            for iLevel=1:length(condNames)
                tempVar.levels{1}(iLevel)=condNames{iLevel}(1);
            end
            if length(tempVar.levels{1}) ~= length(unique(tempVar.levels{1}))
                for iLevel=1:length(condNames)
                    tempVar.levels{1}(iLevel)=96+iLevel;
                end
            end
        else
            ep_tictoc('end');
            return
        end
    end
end

ANOVAlevels=0;
for theFactor=1:6
    if ~isempty(tempVar.levels{theFactor}) && ~isempty(tempVar.factor{theFactor})
        if ANOVAlevels==0
            ANOVAlevels=1;
        end
        ANOVAlevels=ANOVAlevels*length(tempVar.levels{theFactor});
    end
end

EPmain.anova.data.totalLevels=ANOVAlevels;

EPmain.anova.undo.factor=EPmain.anova.data.factor;
EPmain.anova.undo.levels=EPmain.anova.data.levels;
EPmain.anova.data.factor=tempVar.factor;
EPmain.anova.data.levels=tempVar.levels;
ep_tictoc('end');

ep('start');

