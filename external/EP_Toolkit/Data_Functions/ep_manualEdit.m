function ep_manualEdit
% ep_manualEdit - ep_manualEdit -
% Provides window for performing manual editing.
%
%Input:
%

%History
%  by Joseph Dien (1/20/13)
%  jdien07@mac.com
%
%
% bugfix 4/1/13 JD
% Fixed crash when using secondary datasets with differing fields, such as PCA and not.
% Fixed crash when right shifting the cell and a secondary dataset is already at the maximum cell.
%
% modified 5/9/13 JD
% Added table which lists %age of bad channels and allows channels to be marked globally good or bad.
%
% modified 10/15/13 JD
% Added table which lists trial specifics for single trial data and boundary lines for continuous data.
%
% bugfix 11/1/13 JD
% Fixes font sizes on Windows.
%
% modified 11/3/13 JD
% For continuous data, baseline correct each channel by entire one second epoch so that waves will be visible.
% + and - change scale buttons added.
%
% modified 11/6/13 JD
% Scan and Waves functions can now present event markings.
%
% modified 11/2/13 JD
% Added toggle to center data and a toggle to switch between clicking to mark bad channels and clicking to expand the
% waveform into a separate window.
%
% bugfix 12/5/13 JD
% Fixed showing bad channel markings for first cell regardless of the currently displayed cell.
%
% bugfix 1/12/14 JD
% Workaround for Matlab bug periodically causing screen size to register as having zero size.
%
% modified 2/26/14 JD
% Added View function option to plot or erpimage all trials and all subjects.
%
% modified 3/19/14 JD
% Eliminated noTable option for old versions of Matlab.
%
% modified 4/24/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure, including support for complex numbers.
%
% bugfix 5/1/14 JD
% Fixed crash when saving edits.
% Fixed all subjects shown as having a bad cell if the first subject in an average dataset has a bad cell.
%
% bugfix 12/2/14 JD
% Fixed crash when the datasets have different numbers of channels
% (including regional ones).
%
% bugfix 12/31/14 JD
% Fixed crash when Scanning single-trial dataset where one condition has
% fewer trials than the others and the current trial goes beyond the number
% that it has available.
% Fixed crash when Scanning single-trial dataset where one condition has
% fewer trials than the others and one switches to that cell while the
% trial is set on a number beyond the number that is available.
%
%  modified 5/25/14 JD
%  Set colormap to jet even for Matlab 2014b onwards.
%
% bugfix 5/29/15 JD
% Fixed crash in Scan function when dataset has no events.
% Fixed navigating around one-second epochs in continuous data not working properly.
%
% bugfix 8/14/15 JD
% Fixed centering of continuous data being based on entire dataset rather
% than just the one-second epoch being viewed, rendering it sometimes not
% useful.
%
% bugfix 8/30/15 JD
% Fixed amplitude spectral density calculated as divided by Hz rather than
% by square root of Hz.
% Fixed dB of amplitude data not converted to power first.
%
% bugfix 1/15/16 JD
% Now allows power scaled data to be displayed as amplitudes.
% Now handles complex FFT numbers.
%
% modified 1/21/16 JD
% Consolidated spectral unit controls so just four options (cmp, asd, psd, dB).
%
% modified 3/8/16 JD
% Eliminated support for .power field.
%
% bugfix 3/10/17 JD
% Wrong channels being marked when clicking on channels using Edit mode for cells beyond the first in single-trial data.
%
% bugfix 6/16/17 JD
% Fixed conversion to spectral density dividing by bin width rather than sqrt(bin width).
%
% modified 6/15/18 JD
% Added option to add marks for RT and selected events to View figures.
% Eliminated the event control.
%
% modified 8/8/18 JD
% Sped up speed of scans after addition of event marking slowed it down unacceptably.
%
% bugfix 3/9/19 JD
% Fixed crash in View>Scan when zooming in (+) or zooming out (-) segmented data that has no events.
%
% bugfix 4/17/19 JD
% Fixed crash and incorrect bad channel and bad trial info in View>Scan for single-trial data where a trial has been deleted via the Edit function.
%
% bugfix 1/12/20 JD
% Fixed crash in Scan when there are no events in the present waveform.
% Fixed crash in Scan when zooming in or out.
% Fixed Scan not showing single event lines at the correct latency.
% Fixed Scan not displaying event line if the initial epoch did not have an event or leaving it unchanged if changing from an epoch with an event to one without.
%
% modified 4/13/20 JD
% ep_saveEPdataset now handles replacing existing EPdataset entries.
%
% modified 4/22/20 JD
% Added support for up to eight colors in waveform figures.
%
% bugfix 8/17/20 JD
% Fixed centering to handle NaN values.
%
% bugfix 3/19/21 JD
% Fixed subsequent crashes because fields are in the wrong order after manual editing.
% Fixed edit fields not appearing until epoch is changed.
%
% modified 3/31/25 JD
% Added support for virtual grand averages.
% Made various fixes to buggy single-trial interface.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 1999-2018  Joseph Dien
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

global EPdataset EPmain EPmanualEdit EPtictoc

scrsz = EPmain.scrsz;

EPmanualEdit.dataType=EPdataset.dataset(EPmain.view.dataset(1)).dataType;
EPmanualEdit.edited=0;
% EPmanualEdit.eventValue=length(EPmain.view.eventList)+1;

EPmain.handles.manualEditing = figure('Name', 'Manual Editing', 'NumberTitle', 'off', 'Position',[scrsz(1) scrsz(4)-1100 200 500], 'MenuBar', 'none', 'CloseRequestFcn',@confirmCloseEdit);
colormap jet;
drawnow

firstColor=0;
for iColor=1:EPmain.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        if firstColor==0
            firstColor=iColor;
        end
        EPdata=ep_loadEPdataset(EPmain.view.dataset(iColor));
        if ~isempty(EPdata.GAVsubs)
            numSubs=length(EPdata.subNames);
            numVsubs=size(EPdata.GAVsubs,1)-1;
            numRsubs=numSubs-numVsubs;
            %convert virtual GAVEs to normal form so waveforms are available.
            EPdata=ep_combineData(EPdata,'convert',{[],[],[],[],[],[]},[],[],[]);
            if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
            if isempty(EPdata)
                msg{1}='Error: Could not convert virtual grand average into real one.';
                [msg]=ep_errorMsg(msg);
                return
            end
        end
        EPmanualEdit.pca{iColor}=EPdata.pca;
        EPdata=rmfield(EPdata,'pca'); %pca field contents too variable to merge into one structure so temporarily remove
        EPmanualEdit.EPdata(iColor)=EPdata;
    end
end
EPmanualEdit.maxSegs=length(EPdataset.dataset(EPmain.view.dataset(1)).trialNames);
EPmanualEdit.maxCells=length(EPdataset.dataset(EPmain.view.dataset(1)).cellNames);
EPmanualEdit.maxSubs=length(EPdataset.dataset(EPmain.view.dataset(1)).subNames);

figure(EPmain.handles.manualEditing)

EPmanualEdit.handles.leftCell = uicontrol('Style', 'pushbutton', 'String', '<--','FontSize',EPmain.fontsize,...
    'Position', [5 460 50 30], 'Callback', @leftCell);
if EPmain.view.cell(1) ==1
    set(EPmanualEdit.handles.leftCell,'enable','off');
end

EPmanualEdit.handles.rightCell = uicontrol('Style', 'pushbutton', 'String', '-->','FontSize',EPmain.fontsize,...
    'Position', [150 460 50 30], 'Callback', @rightCell);
if EPmain.view.cell(1) == EPmanualEdit.maxCells
    set(EPmanualEdit.handles.rightCell,'enable','off');
end

EPmanualEdit.handles.cell = uicontrol('Style','popupmenu',...
    'String',EPmain.view.theCells{EPmain.view.dataset(1)},'ForegroundColor','blue','FontSize',EPmain.fontsize,...
    'Value',EPmain.view.cell(1),'Position',[55 465 95 20],...
    'Callback', @menuCell);

EPmanualEdit.handles.leftSub = uicontrol('Style', 'pushbutton', 'String', '<--','FontSize',EPmain.fontsize,...
    'Position', [5 430 50 30], 'Callback', @leftSub);
if EPmain.view.subject(1) ==1
    set(EPmanualEdit.handles.leftSub,'enable','off');
end

EPmanualEdit.handles.rightSub = uicontrol('Style', 'pushbutton', 'String', '-->','FontSize',EPmain.fontsize,...
    'Position', [150 430 50 30], 'Callback', @rightSub);
if EPmain.view.subject(1) == EPmanualEdit.maxSubs
    set(EPmanualEdit.handles.rightSub,'enable','off');
end

EPmanualEdit.handles.sub = uicontrol('Style','popupmenu',...
    'String',EPdataset.dataset(EPmain.view.dataset(1)).subNames,'ForegroundColor','blue','FontSize',EPmain.fontsize,...
    'Value',EPmain.view.subject(1),'Position',[55 435 95 20],...
    'Callback', @menuSub);

EPmanualEdit.handles.leftTrial = uicontrol('Style', 'pushbutton', 'String', '<--','FontSize',EPmain.fontsize,...
    'Position', [5 400 50 30], 'Callback', @leftTrial);
if (EPmain.view.trial(1) ==1) || strcmp(EPmanualEdit.dataType,'average')
    set(EPmanualEdit.handles.leftTrial,'enable','off');
end

EPmanualEdit.handles.rightTrial = uicontrol('Style', 'pushbutton', 'String', '-->','FontSize',EPmain.fontsize,...
    'Position', [150 400 50 30], 'Callback', @rightTrial);
if (EPmain.view.trial(1) == EPmanualEdit.maxSegs) || strcmp(EPmanualEdit.dataType,'average')
    set(EPmanualEdit.handles.rightTrial,'enable','off');
end

if ~isempty(EPdataset.dataset(EPmain.view.dataset(1)).trialNames)
    EPmanualEdit.handles.trialNum = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',EPdataset.dataset(EPmain.view.dataset(1)).trialNames(strcmp(EPmain.view.theCells{EPmain.view.dataset(1)}(EPmain.view.cell(1)),EPdataset.dataset(EPmain.view.dataset(1)).cellNames)),'ForegroundColor','blue',...
        'Value',EPmain.view.trial(1),'Position',[55 405 95 20],...
        'Callback', [@menuTrial]);
elseif strcmp(EPmanualEdit.dataType,'continuous')
    theEpochs=cellstr(num2str([1:ceil(length(EPdataset.dataset(EPmain.view.dataset(1)).timeNames)/EPdataset.dataset(EPmain.view.dataset(1)).Fs)]'));
    for i=1:length(theEpochs)
        theEpochs{i}=['Sec: ' theEpochs{i}];
    end
    EPmanualEdit.handles.trialNum = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',theEpochs,'ForegroundColor','blue',...
        'Value',EPmain.view.trial(1),'Position',[55 405 95 20],...
        'Callback', [@menuTrial]);
else
    EPmanualEdit.handles.trialNum = uicontrol('Style','text','HorizontalAlignment','left','String', '     No Trials','FontSize',EPmain.fontsize,...
        'ForegroundColor','blue','Position',[55 405 95 20]);
end

if strcmp(EPmanualEdit.dataType,'average')
    set(EPmanualEdit.handles.trialNum,'enable','off');
end

EPmanualEdit.handles.zoomOut = uicontrol('Style', 'pushbutton', 'String', '-','FontSize',EPmain.fontsize*2,...
    'Position', [130 360 30 30], 'Callback', @zoomOut);

EPmanualEdit.handles.zoomIn = uicontrol('Style', 'pushbutton', 'String', '+','FontSize',EPmain.fontsize*2,...
    'Position', [160 360 30 30], 'Callback', @zoomIn);

EPmanualEdit.handles.edit = uicontrol('Style', 'togglebutton', 'String', 'Edit','FontSize',EPmain.fontsize,...
    'Position', [130 330 30 30], 'Callback', @scanEditMode);

if ~isempty(EPmanualEdit.EPdata(1).facNames)
    set(EPmanualEdit.handles.edit,'enable','off');
end

EPmanualEdit.handles.center = uicontrol('Style', 'togglebutton', 'String', 'Ctr','FontSize',EPmain.fontsize,...
    'Position', [160 330 30 30], 'Callback', @centerLines);

if strcmp(EPmanualEdit.dataType,'continuous')
    theLabel='Second';
elseif strcmp(EPmanualEdit.dataType,'single_trial')
    theLabel='Trial';
else
    theLabel='Cell';
end

if strcmp(EPmanualEdit.EPdata(firstColor).dataType,'average')
    theCell=EPmain.view.cell(firstColor);
else  %if single_trial data or continuous
    cellList=find(strcmp(EPmain.view.theCells{EPmain.view.dataset(firstColor)}(EPmain.view.cell(firstColor)),EPdataset.dataset(EPmain.view.dataset(firstColor)).cellNames));
    if EPmain.view.allTrials(firstColor) || strcmp(EPmanualEdit.EPdata(firstColor).dataType,'continuous')
        theCell=cellList;
    else
        if EPmain.view.trial(firstColor) > length(cellList)
            theCell=cellList(end);
        else
            theCell=cellList(EPmain.view.trial(firstColor));
        end
    end
end

EPmanualEdit.handles.blinkLabel = uicontrol('Style','text','HorizontalAlignment','left','String', ['Blink ' theLabel],'FontSize',EPmain.fontsize,...
    'ForegroundColor','black','Position',[25 365 100 20]);

EPmanualEdit.handles.blink = uicontrol('Style','text','FontSize',EPmain.fontsize,...
    'String',num2str(EPmanualEdit.EPdata(1).analysis.blinkTrial(1,theCell)),'Position',[5 365 20 20]);

EPmanualEdit.handles.saccadeLabel = uicontrol('Style','text','HorizontalAlignment','left','String', ['Saccade ' theLabel],'FontSize',EPmain.fontsize,...
    'ForegroundColor','black','Position',[25 345 100 20]);

EPmanualEdit.handles.saccade = uicontrol('Style','text','FontSize',EPmain.fontsize,...
    'String',num2str(EPmanualEdit.EPdata(1).analysis.saccadeTrial(1,theCell)),'Position',[5 345 20 20]);

EPmanualEdit.handles.moveLabel = uicontrol('Style','text','HorizontalAlignment','left','String', ['Movement ' theLabel],'FontSize',EPmain.fontsize,...
    'ForegroundColor','black','Position',[25 325 100 20]);

EPmanualEdit.handles.move = uicontrol('Style','text','FontSize',EPmain.fontsize,...
    'String',num2str(EPmanualEdit.EPdata(1).analysis.moveTrial(1,theCell)),'Position',[5 325 20 20]);

EPmanualEdit.handles.moveLabel = uicontrol('Style','text','HorizontalAlignment','left','String', ['Bad ' theLabel],'FontSize',EPmain.fontsize,...
    'ForegroundColor','black','Position',[25 305 100 20]);

if strcmp(EPmanualEdit.dataType,'average')
    badTrial=0;
    if EPmanualEdit.EPdata(1).avgNum(EPmain.view.subject(1),theCell) == -1
        badTrial=1;
    end
    EPmanualEdit.handles.bad = uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
        'Value',badTrial,'Position',[5 305 20 20],...
        'Callback', ['global EPmanualEdit;','EPmanualEdit.EPdata(1).avgNum(1,' num2str(theCell) ')=-get(EPmanualEdit.handles.bad,''Value'');','EPmanualEdit.edited=1;']);
else
    EPmanualEdit.handles.bad = uicontrol('Style','checkbox',...
        'Value',EPmanualEdit.EPdata(1).analysis.badTrials(1,theCell),'Position',[5 305 20 20],'FontSize',EPmain.fontsize,...
        'Callback', ['global EPmanualEdit;','EPmanualEdit.EPdata(1).analysis.badTrials(1,' num2str(theCell) ')=get(EPmanualEdit.handles.bad,''Value'');','EPmanualEdit.edited=1;']);
end

if ~isempty(EPmanualEdit.EPdata(1).facNames)
    set(EPmanualEdit.handles.bad,'enable','off');    
end

% EPmanualEdit.handles.events = uicontrol('Style', 'popupmenu', 'String', [EPmain.view.eventList; 'All'],'FontSize',EPmain.fontsize,...
%     'Value',EPmanualEdit.eventValue,...
%     'Position', [125 300 80 20], 'Callback', @changeEvent);
% 
% set(EPmain.handles.view.events,'enable','off');
% 
if any(strcmp(EPmanualEdit.dataType,{'average','grand_average'}))
    badChansPercent=squeeze(sum(isnan(EPmanualEdit.EPdata(1).analysis.badChans(EPmain.view.subject(1),:,:)),2))./length(EPmanualEdit.EPdata(1).cellNames);
else
    badChansPercent=squeeze(sum(EPmanualEdit.EPdata(1).analysis.badChans(EPmain.view.subject(1),:,:)==-1,2))./length(EPmanualEdit.EPdata(1).cellNames);
end

EPmanualEdit.tableData(:,1)=EPdataset.dataset(EPmain.view.dataset(1)).chanNames;
EPmanualEdit.tableData(:,2)=num2cell(badChansPercent);
EPmanualEdit.tableData(:,3)=num2cell(badChansPercent==1);

tableNames{1}='chan';
tableNames{2}='%bad';
tableNames{3}='global';

if ~isempty(EPmanualEdit.EPdata(1).facNames)
    columnEditable =  [false false false];
else
    columnEditable =  [false false true];
end
ColumnFormat{1}='char';
ColumnFormat{2}='numeric';
ColumnFormat{3}='logical';

columnWidth=num2cell([60 40 40]);

EPmanualEdit.handles.hTable = uitable('Data',EPmanualEdit.tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
    'ColumnEditable', columnEditable,'columnWidth',columnWidth,...
    'RearrangeableColumns','on',...
    'CellEditCallback',@changeGlobal,'Position',[5 200 200 100]);

if strcmp(EPmanualEdit.dataType,'single_trial') && ~isempty(EPdataset.dataset(EPmain.view.dataset(1)).trialSpecNames)
    EPmanualEdit.specTableData(1:length(EPdataset.dataset(EPmain.view.dataset(1)).trialSpecNames),1)=EPdataset.dataset(EPmain.view.dataset(1)).trialSpecNames;
    EPmanualEdit.specTableData(1:length(EPdataset.dataset(EPmain.view.dataset(1)).trialSpecNames),2)=EPdataset.dataset(EPmain.view.dataset(1)).trialSpecs(EPmain.view.trial(1),:);
    
    columnEditable =  [false false];
    ColumnFormat{1}='char';
    ColumnFormat{2}='char';
    
    columnWidth=num2cell([60 100]);
    
    EPmanualEdit.handles.specTable = uitable('Data',EPmanualEdit.specTableData,'FontSize',EPmain.fontsize,...
        'ColumnEditable', columnEditable,'columnWidth',columnWidth,...
        'Position',[5 40 200 150]);
    
end

EPmanualEdit.handles.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
    'Position', [92 0 50 35], 'Callback', @cancelEdit);

EPmanualEdit.handles.done = uicontrol('Style', 'pushbutton', 'String', 'Keep','FontSize',EPmain.fontsize,...
    'Position', [152 0 50 35], 'Callback', @done);

if ~isempty(EPmanualEdit.EPdata(1).facNames)
    set(EPmanualEdit.handles.done,'enable','off');    
end

% updateDisplay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cancelEdit(~,~)
%quit from manual editing window without saving edits

global EPmain EPwaves EPmanualEdit 

delete(EPmain.handles.manualEditing);
try
    delete(EPwaves.handles.waves.hWaveWindow);
catch
    disp('Wave window is missing, presumably already closed manually by user.');
end

EPmanualEdit=[];

ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function done(~,~)
%quit from manual editing window and save edits

global EPmain EPwaves EPmanualEdit EPdataset

close(EPmain.handles.manualEditing);
try
    close(EPwaves.handles.waves.hWaveWindow);
catch
    disp('Wave window is missing, presumably already closed manually by user.');
end

if ~exist([EPdataset.EPwork filesep 'EPwork'],'dir')
    msg{1}='The work directory cannot be found.';
    [msg]=ep_errorMsg(msg);
    return
end

if EPmanualEdit.edited
    EPmanualEdit.EPdata(1).pca=EPmanualEdit.pca{1}; %put pca field back in, which was temporarily taken out since contents are not standardized.
    EPdata=ep_updateEPfile(EPmanualEdit.EPdata(1)); %ensure that EPdata fields are in the correct order.
    ep_saveEPdataset(EPdata,EPmain.view.dataset(1),'no');
end

EPmanualEdit=[];

ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function leftCell(~,~)
%move to previous cell

global EPmain EPdataset EPmanualEdit

for iColor=1:EPmain.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        if EPmain.view.cell(iColor) > 1
            EPmain.view.cell(iColor)=EPmain.view.cell(iColor)-1;
            set(EPmain.handles.view.cell(iColor),'Value',EPmain.view.cell(iColor));
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames)
                theEpochs=EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames));
            elseif strcmp(EPmanualEdit.dataType,'continuous')
                theEpochs=cellstr(num2str([1:ceil(length(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames)/EPdataset.dataset(EPmain.view.dataset(iColor)).Fs)]'));
                for iEpoch=1:length(theEpochs)
                    theEpochs{iEpoch}=['Sec: ' theEpochs{iEpoch}];
                end
            else
                theEpochs=cell(1);
                theEpochs{1}='     No Trials';
            end
            set(EPmain.handles.view.trial(iColor),'Value',1);
            set(EPmain.handles.view.trial(iColor),'String',theEpochs);
            EPmain.view.trial(iColor)=1;
        end
    end
end

updateManualEdit
updateDisplay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rightCell(~,~)
%move to next cell

global EPmain EPdataset EPmanualEdit

for iColor=1:EPmain.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        if EPmain.view.cell(iColor) < length(EPmain.view.theCells{EPmain.view.dataset(iColor)})
            EPmain.view.cell(iColor)=EPmain.view.cell(iColor)+1;
            set(EPmain.handles.view.cell(iColor),'Value',EPmain.view.cell(iColor));
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames)
                theEpochs=EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames));
            elseif strcmp(EPmanualEdit.dataType,'continuous')
                theEpochs=cellstr(num2str([1:ceil(length(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames)/EPdataset.dataset(EPmain.view.dataset(iColor)).Fs)]'));
                for iEpoch=1:length(theEpochs)
                    theEpochs{iEpoch}=['Sec: ' theEpochs{iEpoch}];
                end
            else
                theEpochs=cell(1);
                theEpochs{1}='     No Trials';
            end
            set(EPmain.handles.view.trial(iColor),'Value',1);
            set(EPmain.handles.view.trial(iColor),'String',theEpochs);
            EPmain.view.trial(iColor)=1;
        end
    end
end

updateManualEdit
updateDisplay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function leftSub(~,~)
%move to previous subject

global EPmain EPdataset EPmanualEdit

for iColor=1:EPmain.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        if EPmain.view.subject(iColor) > 1 && (EPmain.view.subject(iColor) <= length(EPmanualEdit.EPdata(iColor).subNames))
            EPmain.view.subject(iColor)=EPmain.view.subject(iColor)-1;
            set(EPmain.handles.view.subject(iColor),'Value',EPmain.view.subject(iColor));
        end
    end
end

updateTable
updateManualEdit
updateDisplay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rightSub(~,~)
%move to next subject

global EPmain EPdataset

for iColor=1:EPmain.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        if EPmain.view.subject(iColor) < length(EPdataset.dataset(EPmain.view.dataset(iColor)).subNames)
            EPmain.view.subject(iColor)=EPmain.view.subject(iColor)+1;
            set(EPmain.handles.view.subject(iColor),'Value',EPmain.view.subject(iColor));
        end
    end
end

updateTable
updateManualEdit
updateDisplay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function leftTrial(~,~)
%move to previous trial

global EPmain EPdataset EPmanualEdit EPwaves

for iColor=1:EPmain.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        if (EPmain.view.trial(iColor) > 1) && (EPmain.view.trial(iColor) <= length(get(EPmanualEdit.handles.trialNum,'String')))
            if EPmain.view.trial(iColor) > 1
                EPmain.view.trial(iColor)=EPmain.view.trial(iColor)-1;
            end
            set(EPmain.handles.view.trial(iColor),'Value',EPmain.view.trial(iColor));
            if strcmp(EPmanualEdit.dataType,'continuous')
                EPwaves.startSamp(iColor)=min(find(1000*(EPmain.view.trial(iColor)-1) <= EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames));
                EPwaves.lastSamp(iColor)=max(find((1000*(EPmain.view.trial(iColor))-EPwaves.sampleSize) >= EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames));
                EPwaves.numPoints=EPwaves.lastSamp(iColor)-EPwaves.startSamp(iColor)+1;
                EPmain.view.startSamp=1000*(EPmain.view.trial(iColor)-1); %for continuous data, each "trial" is 1000 ms
                EPmain.view.endSamp=1000*(EPmain.view.trial(iColor));
                set(EPmain.handles.view.startSamp,'String',EPmain.view.startSamp);
                set(EPmain.handles.view.endSamp,'String',EPmain.view.endSamp);
            end
        end
    end
end
    
updateManualEdit
updateDisplay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rightTrial(~,~)
%move to next trial

global EPmain EPdataset EPmanualEdit EPwaves

for iColor=1:EPmain.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        if EPmain.view.trial(iColor) < length(get(EPmanualEdit.handles.trialNum,'String'))
            cellList=find(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames));
            if (EPmain.view.trial(iColor) == length(cellList)) && ~strcmp(EPmanualEdit.dataType,'continuous')
                EPmain.view.trial(iColor)=length(cellList);
            else
                EPmain.view.trial(iColor)=EPmain.view.trial(iColor)+1;
            end
            set(EPmain.handles.view.trial(iColor),'Value',EPmain.view.trial(iColor));
            if strcmp(EPmanualEdit.dataType,'continuous')
                EPwaves.startSamp(iColor)=min(find(1000*(EPmain.view.trial(iColor)-1) <= EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames));
                EPwaves.lastSamp(iColor)=max(find((1000*(EPmain.view.trial(iColor))-EPwaves.sampleSize) >= EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames));
                EPwaves.numPoints=EPwaves.lastSamp(iColor)-EPwaves.startSamp(iColor)+1;
                EPmain.view.startSamp=1000*(EPmain.view.trial(iColor)-1); %for continuous data, each "trial" is 1000 ms
                EPmain.view.endSamp=1000*(EPmain.view.trial(iColor));
                set(EPmain.handles.view.startSamp,'String',EPmain.view.startSamp);
                set(EPmain.handles.view.endSamp,'String',EPmain.view.endSamp);
            end
        end
    end
end

updateManualEdit
updateDisplay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateDisplay(~,~)
%update displays after changing the data to be viewed

global EPmain EPdataset EPwaves EPmanualEdit EPtictoc

figure(EPwaves.handles.waves.hWaveWindow)

EPmanualEdit.eventLines=ep_collateEventLines('EPwaves');

%organize the data for plotting
theMax=0;
trialCount=0;
for iColor=1:EPmain.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)

        theSub=EPmain.view.subject(iColor);
        if strcmp(EPmanualEdit.EPdata(iColor).dataType,'average')
            theCell=EPmain.view.cell(iColor);
            numSubs=length(EPdataset.dataset(EPmain.view.dataset(iColor)).subNames);
            if EPmain.view.allTrials(iColor)
                theSub=[1:numSubs];
            end
        else  %if single_trial data or continuous
            cellList=find(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames));
            if EPmain.view.allTrials(iColor) || strcmp(EPmanualEdit.EPdata(iColor).dataType,'continuous')
                theCell=cellList;
            else
                if EPmain.view.trial(iColor) > length(cellList)
                    theCell=cellList(end);
                else
                    theCell=cellList(EPmain.view.trial(iColor));
                end
            end
        end
        
        tempData=ep_expandFacs(EPmanualEdit.EPdata(iColor),[],EPwaves.startSamp(iColor):EPwaves.lastSamp(iColor),theCell,theSub,EPmain.view.factor(iColor),EPwaves.startBins(iColor):EPwaves.lastBins(iColor));
        ep_tictoc;if EPtictoc.stop;return;end
        if get(EPmanualEdit.handles.center,'Value')
            if strcmp('VLT',EPmain.view.dataTransform)
                for iChan=1:size(tempData,1)
                    tempData(iChan,:,:,:,:,:,:)=tempData(iChan,:,:,:,:,:,:)-mean(tempData(iChan,:,:,:,:,:,:),2,'omitnan');
                end
            end
        end
        
        tempEvents=[];
        if strcmp(EPmanualEdit.EPdata(iColor).dataType,'continuous')
            if strcmp('VLT',EPmain.view.dataTransform)
                for iChan=1:size(tempData,1)
                    tempData(iChan,:,:,:,:)=tempData(iChan,:,:,:,:)-mean(tempData(iChan,:,:,:,:),2,'omitnan'); %center the waveforms if continuous
                end
            end
            if ~isempty(EPmanualEdit.EPdata(iColor).events{EPmain.view.subject(iColor),1})
                tempEvents=EPmanualEdit.EPdata(iColor).events{EPmain.view.subject(iColor),1}(([EPmanualEdit.EPdata(iColor).events{1}.sample]>=EPwaves.startSamp(iColor)) & ([EPmanualEdit.EPdata(iColor).events{1}.sample]<=EPwaves.lastSamp(iColor)));
            end
        else
            for iCell=1:length(theCell)
                for iSub=1:length(theSub)
                    tempEvents=[tempEvents EPmanualEdit.EPdata(iColor).events{theSub(iSub),theCell(iCell)}];
                end
            end
        end
        EPwaves.eventWave{iColor}=cell(1);
        EPwaves.eventWave{iColor}{1}=zeros(1,size(tempData,2));
        EPwaves.boundary{iColor}=[];
        if ~isempty(EPmanualEdit.eventLines{iColor})
            EPwaves.eventWave{iColor}{1}=histc(EPmanualEdit.eventLines{iColor},[1:size(tempData,2)]);
            theMax=max([theMax EPwaves.eventWave{iColor}{1}]);
        end
        if ~isempty(tempEvents)
            boundaryEvents=find(strcmp('boundary',{tempEvents.value}));
            if ~isempty(boundaryEvents)
                EPwaves.boundary{iColor}=tempEvents(boundaryEvents).sample;
            end
        end
%         if ~any(any(EPwaves.eventWave{iColor}{1}))
%             EPwaves.eventWave{iColor}{1}=[];
%         end

        if ~isreal(EPmanualEdit.EPdata(iColor).data) %if the data has an imaginary component, as in spectral data
            tempDataImag=imag(tempData);
            tempData=real(tempData);
        else
            tempDataImag=[];
        end
        if ~isempty(EPmanualEdit.EPdata(iColor).relNames) %collapse over the relations dimension if any to provide a mean
            if strcmp(EPmanualEdit.EPdata(iColor).dataType,'average')
                goodRelChans=find(squeeze(any(any(~isnan(EPmanualEdit.EPdata(iColor).analysis.badChans(theSub,theCell,:)),2),1)));
            else
                goodRelChans=find(squeeze(any(any((EPmanualEdit.EPdata(iColor).analysis.badChans(theSub,theCell,:)~=-1),2),1)));
            end
            if isscalar(EPmanualEdit.EPdata(iColor).reference.current)
                goodRelChans=setdiff(goodRelChans,EPmanualEdit.EPdata(iColor).reference.current); %coherence with a single reference channel is NaN.
            end
            tempData=mean(abs(tempData(:,:,:,:,:,:,goodRelChans)),7,'omitnan');
            if ~isempty(tempDataImag)
                tempDataImag=mean(abs(tempDataImag(:,:,:,:,:,:,goodRelChans)),7,'omitnan');
            end
        end
        if strcmp(EPmanualEdit.EPdata(iColor).dataType,'average')
            numWaves=size(tempData,4);
        else
            numWaves=size(tempData,3);
        end
        %the 4 dimensions of totalData are chans, points, waves(trials/cells/subjects), and frequencies
        EPwaves.totalData(EPwaves.chanIX{iColor},:,trialCount+1:trialCount+numWaves,:)=squeeze(tempData); %rearrange order of channels to be consistent with other datasets
        trialCount=trialCount+numWaves;
        if ~isempty(tempDataImag)
            EPwaves.totalData(EPwaves.chanIX{iColor},:,trialCount+1:trialCount+numWaves,:)=squeeze(tempDataImag); %rearrange order of channels to be consistent with other datasets
            trialCount=trialCount+numWaves;
        end
    end
end

if strcmp(EPmanualEdit.dataType,'average') %averaged data
    theCell=EPmain.view.cell(1);
elseif strcmp(EPmanualEdit.dataType,'continuous')
    theCell=EPmain.view.trial(1);
else  %if single_trial data
    theCell=intersect(find(strcmp(EPmain.view.theCells{EPmain.view.dataset(1)}(EPmain.view.cell(1)),EPmanualEdit.EPdata(1).cellNames)),...
        find(EPmanualEdit.EPdata(1).trialNames(EPmain.view.trial(1))==EPmanualEdit.EPdata(1).trialNames));
end

if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
    nonCorrWaves=find(~EPwaves.correlPlotScaleIndex);
    if (EPmain.view.FFTunits > 1)
        EPwaves.totalData(:,:,nonCorrWaves,:)=abs(EPwaves.totalData(:,:,nonCorrWaves,:)); %convert complex number to real number
    end
    EPwaves.totalData(:,:,nonCorrWaves,:)=EPwaves.totalData(:,:,nonCorrWaves,:)/sqrt(mean(diff(EPdata.freqNames))); %convert to spectral density
    if EPmain.view.FFTunits > 2
        EPwaves.totalData(:,:,nonCorrWaves,:)=EPwaves.totalData(:,:,nonCorrWaves,:).^2; %convert amplitude to power
    end
    if (EPmain.view.FFTunits == 4)
        if ~all(EPwaves.totalData(:,:,nonCorrWaves,:) >=0)
            disp('Some values negative so will take absolute value prior to log transform.  This can happen with PCA results due to minor imprecisions in PCA.');
        end
        EPwaves.totalData(:,:,nonCorrWaves,:)=log10(abs(EPwaves.totalData(:,:,nonCorrWaves,:)))*10; %convert to dB log scaling
        tempVar=EPwaves.totalData(:,:,nonCorrWaves,:);
        tempVar(isinf(tempVar))=-flintmax;
        EPwaves.totalData(:,:,nonCorrWaves,:)=tempVar; %log10 of zero is -inf.  Replace with maximum possible double-precision negative number.
    end
end

if ~all(EPwaves.correlPlotScaleIndex) && any(EPwaves.correlPlotScaleIndex) %if not all the waveforms are correlations but some are, then rescale the correlations to match the plots
    EPwaves.totalData(:,:,find(EPwaves.correlPlotScaleIndex),:)=(EPwaves.totalData(:,:,find(EPwaves.correlPlotScaleIndex),:)*(EPwaves.plotMVmax-EPwaves.plotMVmin))+EPwaves.plotMVmin;
end

switch EPmain.view.dataTransform
    case 'VLT'
        for iChan=1:length(EPwaves.handles.waves.hLines)
            for iColor=1:length(EPwaves.handles.waves.hLines{iChan})
                refreshdata(EPwaves.handles.waves.hLines{iChan}(iColor),'caller');
            end
            %refresh event markers
            for iColor=1:length(EPwaves.eventWave) 
                if ~isempty(EPwaves.eventWave{iColor})
                    theColor=EPwaves.plotColors(iColor);
                    if ~isempty(EPwaves.eventWave{iColor}{1})
                        plotPoints=find(EPwaves.eventWave{theColor}{1}>min(EPwaves.eventWave{theColor}{1}));
                        if isempty(plotPoints)
                            EPwaves.eventData{theColor}{1}=(EPwaves.eventWave{theColor}{1}*(abs(EPwaves.plotMVmin/2)))+EPwaves.plotMVmin;
                            EPwaves.eventTimes{theColor}{1}=[EPwaves.firstTime:EPwaves.spacing:EPwaves.firstTime+(EPwaves.spacing*(EPwaves.numPoints-1))];
                            plotWidth=.1;
                        elseif isscalar(plotPoints)
                            EPwaves.eventData{theColor}{1}=[EPwaves.plotMVmin:((EPwaves.plotMVmin/2)-EPwaves.plotMVmin)/9:EPwaves.plotMVmin/2];
                            theTime=EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames(plotPoints);
                            EPwaves.eventTimes{theColor}{1}=[repmat(theTime,10,1)];
                            plotWidth=2;
                        else
                            EPwaves.eventData{theColor}{1}=(EPwaves.eventWave{theColor}{1}*(abs(EPwaves.plotMVmin/2)))+EPwaves.plotMVmin;
                            EPwaves.eventTimes{theColor}{1}=[EPwaves.firstTime:EPwaves.spacing:EPwaves.firstTime+(EPwaves.spacing*(EPwaves.numPoints-1))];
                            plotWidth=5;
                        end
                        %                         for iLine=1:length(EPwaves.handles.waves.eventLines{iChan,theColor})
                        %                             set(EPwaves.handles.waves.eventLines{iChan,theColor}(iLine),'XDataSource',['EPwaves.eventTimes{' num2str(theColor) '}{1}(' num2str(iLine) ')']);
                        %                             set(EPwaves.handles.waves.eventLines{iChan,theColor}(iLine),'YDataSource',['EPwaves.eventData{' num2str(theColor) '}{1}(' num2str(iLine) ')']);
                        %                         end
                        if isfield(EPwaves.handles.waves,'eventLines')
                            refreshdata(EPwaves.handles.waves.eventLines{iChan,theColor}(1),'caller');
                            set(EPwaves.handles.waves.eventLines{iChan,theColor}(1),'LineWidth',plotWidth);
                        end
                    end
                end
            end
            if ~isempty(EPwaves.boundary{iColor})
                hold on
                for iBoundary=1:length(EPwaves.boundary{iColor})
                    theSample=EPwaves.boundary{iColor}(iBoundary);
                    EPwaves.handles.waves.boundary{iChan,iColor} = line([theSample theSample],[EPwaves.plotMVmin EPwaves.plotMVmax],'Color',EPwaves.thePlotColors(iColor,:),'LineWidth',1);
                end
                hold off
            end
        end
        
    case 'FFT'
        for iChan=1:length(EPwaves.handles.waves.hLines)
            for iColor=1:length(EPwaves.handles.waves.hLines{iChan})
                refreshdata(EPwaves.handles.waves.hLines{iChan}(iColor),'caller');
            end
        end
    case 'TFT'
        for iChan=1:length(EPwaves.handles.waves.hLines)
            for iColor=5:length(EPwaves.handles.waves.hLines{iChan})
                refreshdata(EPwaves.handles.waves.hLines{iChan}(iColor),'caller');
            end
        end
end

figure(EPmain.handles.manualEditing)

set(EPmanualEdit.handles.blink,'String',num2str(EPmanualEdit.EPdata(1).analysis.blinkTrial(1,theCell)));

set(EPmanualEdit.handles.saccade,'String',num2str(EPmanualEdit.EPdata(1).analysis.saccadeTrial(1,theCell)));

set(EPmanualEdit.handles.move,'String',num2str(EPmanualEdit.EPdata(1).analysis.moveTrial(1,theCell)));

if strcmp(EPmanualEdit.dataType,'average')
    badTrial=0;
    if EPmanualEdit.EPdata(1).avgNum(EPmain.view.subject(1),theCell) == -1
        badTrial=1;
    end
    theBad=badTrial;
    theCallBack=['global EPmanualEdit;','EPmanualEdit.EPdata(1).avgNum(1,' num2str(theCell) ')=-get(EPmanualEdit.handles.bad,''Value'');','EPmanualEdit.edited=1;'];
else
    theBad=EPmanualEdit.EPdata(1).analysis.badTrials(1,theCell);
    theCallBack=['global EPmanualEdit;','EPmanualEdit.EPdata(1).analysis.badTrials(1,' num2str(theCell) ')=get(EPmanualEdit.handles.bad,''Value'');','EPmanualEdit.edited=1;'];
end
set(EPmanualEdit.handles.bad,'Value',theBad);
set(EPmanualEdit.handles.bad,'Callback',theCallBack);

if ~isempty(EPmanualEdit.EPdata(1).facNames)
    set(EPmanualEdit.handles.bad,'enable','off');    
end

if isempty(EPmanualEdit.EPdata(1).trialNames) %averaged data or continuous data
    theCell=EPmain.view.cell(1);
else  %if single_trial data
    theCell=intersect(find(strcmp(EPmain.view.theCells{EPmain.view.dataset(1)}(EPmain.view.cell(1)),EPmanualEdit.EPdata(1).cellNames)),...
        find(EPmanualEdit.EPdata(1).trialNames(EPmain.view.trial(1))==EPmanualEdit.EPdata(1).trialNames));
end
for iChan=1:length(EPwaves.handles.waves.hWave)
    if ~isempty(intersect(iChan,EPwaves.chanIX{1}))
        if ((EPmanualEdit.EPdata(1).analysis.badChans(EPmain.view.subject(1),theCell,EPwaves.chanIX{1}(iChan))==-1)&&~strcmp(EPmanualEdit.dataType,'average')) || (isnan(EPmanualEdit.EPdata(1).analysis.badChans(EPmain.view.subject(1),theCell,EPwaves.chanIX{1}(iChan)))&&strcmp(EPmanualEdit.dataType,'average'))
            set(EPwaves.handles.waves.hWave(iChan),'Color',[1 .5 .5])
        else
            set(EPwaves.handles.waves.hWave(iChan),'Color',[1 1 1])
        end
    else
        set(EPwaves.handles.waves.hWave(iChan),'Color',[1 1 1])
    end
end

if strcmp(EPmanualEdit.dataType,'single_trial') && ~isempty(EPdataset.dataset(EPmain.view.dataset(1)).trialSpecNames)
    EPmanualEdit.specTableData(1:length(EPdataset.dataset(EPmain.view.dataset(1)).trialSpecNames),2)=EPdataset.dataset(EPmain.view.dataset(1)).trialSpecs(EPmain.view.trial(1),:);
    
    set(EPmanualEdit.handles.specTable,'Data',EPmanualEdit.specTableData);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function menuTrial(~,~)
%update displays after changing the trial setting

global EPmain EPdataset EPmanualEdit EPwaves

newTrial=get(EPmanualEdit.handles.trialNum,'Value');

for iColor=1:EPmain.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset) && (EPmain.view.trial(iColor) <= length(get(EPmanualEdit.handles.trialNum,'String')))
        EPmain.view.trial(iColor)=min(newTrial,length(get(EPmanualEdit.handles.trialNum,'String')));
        cellList=find(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames));
        if (EPmain.view.trial(iColor) > length(cellList)) && ~strcmp(EPmanualEdit.dataType,'continuous')
            EPmain.view.trial(iColor)=length(cellList);
        end
        set(EPmain.handles.view.trial(iColor),'Value',EPmain.view.trial(iColor));
        if strcmp(EPmanualEdit.dataType,'continuous')
            EPwaves.startSamp(iColor)=min(find(1000*(EPmain.view.trial(iColor)-1) <= EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames));
            EPwaves.lastSamp(iColor)=max(find((1000*(EPmain.view.trial(iColor))-EPwaves.sampleSize) >= EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames));
            EPwaves.numPoints=EPwaves.lastSamp(iColor)-EPwaves.startSamp(iColor)+1;
            EPmain.view.startSamp=1000*(EPmain.view.trial(iColor)-1); %for continuous data, each "trial" is 1000 ms
            EPmain.view.endSamp=1000*(EPmain.view.trial(iColor));
            set(EPmain.handles.view.startSamp,'String',EPmain.view.startSamp);
            set(EPmain.handles.view.endSamp,'String',EPmain.view.endSamp);
        end
    end
end

if EPmain.view.trial(1) == length(get(EPmanualEdit.handles.trialNum,'String'))
    set(EPmanualEdit.handles.rightTrial,'enable','off');
else
    set(EPmanualEdit.handles.rightTrial,'enable','on');
end

if EPmain.view.trial(1) ==1
    set(EPmanualEdit.handles.leftTrial,'enable','off');
else
    set(EPmanualEdit.handles.leftTrial,'enable','on');
end

updateManualEdit
updateDisplay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function menuCell(~,~)
%update displays after changing the cell setting

global EPmain EPdataset EPmanualEdit

newCell=get(EPmanualEdit.handles.cell,'Value');

for iColor=1:EPmain.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        EPmain.view.cell(iColor)=min(newCell,length(EPmain.view.theCells{EPmain.view.dataset(iColor)}));
        set(EPmain.handles.view.cell(iColor),'Value',EPmain.view.cell(iColor));
        if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames)
            theEpochs=EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames));
        elseif strcmp(EPmanualEdit.dataType,'continuous')
            theEpochs=cellstr(num2str([1:ceil(length(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames)/EPdataset.dataset(EPmain.view.dataset(iColor)).Fs)]'));
            for iEpoch=1:length(theEpochs)
                theEpochs{iEpoch}=['Sec: ' theEpochs{iEpoch}];
            end
        else
            theEpochs=cell(1);
            theEpochs{1}='     No Trials';
        end
        EPmain.view.trial(iColor)=1;
        set(EPmain.handles.view.trial(iColor),'Value',1);
        set(EPmain.handles.view.trial(iColor),'String',theEpochs);
    end
end

updateManualEdit
updateDisplay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function menuSub(~,~)
%update displays after changing the subject setting

global EPmain EPdataset EPmanualEdit

newSubject=get(EPmanualEdit.handles.sub,'Value');

for iColor=1:EPmain.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset) && (EPmain.view.subject(iColor) <= length(EPmanualEdit.EPdata(iColor).subNames))
        EPmain.view.subject(iColor)=min(newSubject,length(EPdataset.dataset(EPmain.view.dataset(iColor)).subNames));
        set(EPmain.handles.view.subject(iColor),'Value',EPmain.view.subject(iColor));
    end
end

if EPmain.view.subject(1) == length(EPdataset.dataset(EPmain.view.dataset(1)).subNames)
    set(EPmanualEdit.handles.rightSub,'enable','off');
else
    set(EPmanualEdit.handles.rightSub,'enable','on');
end

if EPmain.view.subject(1) ==1
    set(EPmanualEdit.handles.leftSub,'enable','off');
else
    set(EPmanualEdit.handles.leftSub,'enable','on');
end

updateTable
updateManualEdit
updateDisplay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function badChannel(src,eventdata)
global EPwaves EPmain EPdataset EPmanualEdit

chan=find([EPwaves.handles.waves.hWave] == src);

if isempty(chan)
    for iChan=1:length(EPwaves.handles.waves.hLines)
        for theLine=1:length(EPwaves.handles.waves.hLines{iChan})
            if src == EPwaves.handles.waves.hLines{iChan}(theLine)
                chan=iChan;
            end
        end
    end
end

if isempty(chan)
    disp('Could not find the channel for some reason.')
    return
end

if ~any(strcmp(EPmanualEdit.dataType,{'average','grand_average'}))
    if isempty(EPmanualEdit.EPdata(1).trialNames) %continuous data
        theCell=EPmain.view.cell(1);
    else  %if single_trial data
        theCell=intersect(find(strcmp(EPmain.view.theCells{EPmain.view.dataset(1)}(EPmain.view.cell(1)),EPmanualEdit.EPdata(1).cellNames)),...
            find(EPmanualEdit.EPdata(1).trialNames(EPmain.view.trial(1))==EPmanualEdit.EPdata(1).trialNames));
    end
    if EPmanualEdit.EPdata(1).analysis.badChans(EPmain.view.subject(1),theCell,EPwaves.chanIX{1}(chan))==-1
        set(EPwaves.handles.waves.hWave(chan),'Color',[1 1 1])
        EPmanualEdit.EPdata(1).analysis.badChans(EPmain.view.subject(1),theCell,EPwaves.chanIX{1}(chan))=0;
    else
        set(EPwaves.handles.waves.hWave(chan),'Color',[1 .5 .5])
        EPmanualEdit.EPdata(1).analysis.badChans(EPmain.view.subject(1),theCell,EPwaves.chanIX{1}(chan))=-1;
    end
else
    if isnan(EPmanualEdit.EPdata(1).analysis.badChans(EPmain.view.subject(1),EPmain.view.cell(1),EPwaves.chanIX{1}(chan)))
        set(EPwaves.handles.waves.hWave(chan),'Color',[1 1 1])
        EPmanualEdit.EPdata(1).analysis.badChans(EPmain.view.subject(1),EPmain.view.cell(1),EPwaves.chanIX{1}(chan))=0;
    else
        set(EPwaves.handles.waves.hWave(chan),'Color',[1 .5 .5])
        EPmanualEdit.EPdata(1).analysis.badChans(EPmain.view.subject(1),EPmain.view.cell(1),EPwaves.chanIX{1}(chan))=NaN;
    end
end

EPmanualEdit.edited=1;

updateTable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function updateTable
%draw global bad channel table

global EPmain EPdataset EPmanualEdit

if any(strcmp(EPmanualEdit.dataType,{'average','grand_average'}))
    badChansPercent=squeeze(sum(isnan(EPmanualEdit.EPdata(1).analysis.badChans(EPmain.view.subject(1),EPmanualEdit.EPdata(1).avgNum(EPmain.view.subject(1),:)~=-1,:)),2))./length(find(EPmanualEdit.EPdata(1).avgNum(EPmain.view.subject(1),:) ~=-1));
else
    badChansPercent=squeeze(sum(EPmanualEdit.EPdata(1).analysis.badChans(EPmain.view.subject(1),EPmanualEdit.EPdata(1).analysis.badTrials==0,:)==-1,2))./length(find(EPmanualEdit.EPdata(1).analysis.badTrials==0));
end

EPmanualEdit.tableData(:,1)=EPdataset.dataset(EPmain.view.dataset(1)).chanNames;
EPmanualEdit.tableData(:,2)=num2cell(badChansPercent);
EPmanualEdit.tableData(:,3)=num2cell(badChansPercent==1);

set(EPmanualEdit.handles.hTable,'Data',EPmanualEdit.tableData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeGlobal(src,eventdata)
%flip global bad channel checkboxes

global EPmain EPmanualEdit

chan=eventdata.Indices(1);

if any(strcmp(EPmanualEdit.dataType,{'average','grand_average'}))
    if eventdata.PreviousData==0
        newValue=NaN;
        newPercent=1;
    else
        newValue=0;
        newPercent=0;
    end
else
    if eventdata.PreviousData==0
        newValue=-1;
        newPercent=1;
    else
        newValue=0;
        newPercent=0;
    end
end

EPmanualEdit.EPdata(1).analysis.badChans(EPmain.view.subject(1),:,chan)=newValue;
EPmanualEdit.tableData{eventdata.Indices(1),3}=eventdata.NewData;
EPmanualEdit.tableData{eventdata.Indices(1),2}=newPercent;

EPmanualEdit.edited=1;

updateTable

updateDisplay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function centerLines(src,eventdata)
%zoom in on the data

global EPwaves EPmanualEdit EPmain

if get(EPmanualEdit.handles.center,'Value')
    if strcmp('VLT',EPmain.view.dataTransform)
        maxLim=max(abs(EPwaves.plotMVmin),EPwaves.plotMVmax);
        for chan=1:length(EPwaves.handles.waves.hWave)
            set(EPwaves.handles.waves.hWave(chan),'YLim',[-maxLim maxLim])
        end
    end
else
    for chan=1:length(EPwaves.handles.waves.hWave)
        set(EPwaves.handles.waves.hWave(chan),'YLim',[EPwaves.plotMVmin EPwaves.plotMVmax])
    end
end

updateDisplay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zoomIn(src,eventdata)
%zoom in on the data

global EPwaves EPmain EPdataset EPmanualEdit

EPwaves.plotMVmin=EPwaves.plotMVmin/2;
EPwaves.plotMVmax=EPwaves.plotMVmax/2;

for iColor=1:length(EPwaves.eventWave)
    if ~isempty(EPwaves.eventWave{iColor})
        if ~isempty(EPwaves.eventWave{iColor}{1})
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                EPwaves.eventWave{iColor}{1}=EPwaves.eventWave{iColor}{1}/2; %rescale event lines
%                 boundary=find(strcmp(EPwaves.eventType{iColor},'boundary'));
%                 if ~isempty(boundary)
%                     EPwaves.eventWave{iColor}(boundary,(EPwaves.eventWave{iColor}(boundary,:) ~= EPwaves.plotMVmin))=EPwaves.plotMVmax;
%                 end
            end
        end
    end
end

if get(EPmanualEdit.handles.center,'Value') && strcmp('VLT',EPmain.view.dataTransform)
    maxLim=max(abs(EPwaves.plotMVmin),EPwaves.plotMVmax);
    for chan=1:length(EPwaves.handles.waves.hWave)
        set(EPwaves.handles.waves.hWave(chan),'YLim',[-maxLim maxLim])
    end
else
    for chan=1:length(EPwaves.handles.waves.hWave)
        set(EPwaves.handles.waves.hWave(chan),'YLim',[EPwaves.plotMVmin EPwaves.plotMVmax])
    end
end

if strcmp(EPmain.view.dataTransform,'VLT')
    for chan=1:length(EPwaves.handles.waves.hWave)
        if isfield(EPwaves.handles.waves,'eventLines')
            for iColor=1:size(EPwaves.handles.waves.eventLines,2)
                for iLine=1:length(EPwaves.handles.waves.eventLines{chan,iColor})
                    refreshdata(EPwaves.handles.waves.eventLines{chan,iColor}(iLine),'caller');
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function zoomOut(src,eventdata)
%zoom out on the data

global EPwaves EPmain EPdataset EPmanualEdit

EPwaves.plotMVmin=EPwaves.plotMVmin*2;
EPwaves.plotMVmax=EPwaves.plotMVmax*2;

for iColor=1:length(EPwaves.eventWave)
    if ~isempty(EPwaves.eventWave{iColor})
        if ~isempty(EPwaves.eventWave{iColor}{1})
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                EPwaves.eventWave{iColor}{1}=EPwaves.eventWave{iColor}{1}*2; %rescale event lines
%                 boundary=find(strcmp(EPwaves.eventType{iColor},'boundary'));
%                 if ~isempty(boundary)
%                     EPwaves.eventWave{iColor}(boundary,(EPwaves.eventWave{iColor}(boundary,:) ~= EPwaves.plotMVmin))=EPwaves.plotMVmax;
%                 end
            end
        end
    end
end

if get(EPmanualEdit.handles.center,'Value') && strcmp('VLT',EPmain.view.dataTransform)
    maxLim=max(abs(EPwaves.plotMVmin),EPwaves.plotMVmax);
    for iChan=1:length(EPwaves.handles.waves.hWave)
        set(EPwaves.handles.waves.hWave(iChan),'YLim',[-maxLim maxLim])
    end
else 
    for iChan=1:length(EPwaves.handles.waves.hWave)
        set(EPwaves.handles.waves.hWave(iChan),'YLim',[EPwaves.plotMVmin EPwaves.plotMVmax])
    end
end

if strcmp(EPmain.view.dataTransform,'VLT')
    for iChan=1:length(EPwaves.handles.waves.hWave)
        if isfield(EPwaves.handles.waves,'eventLines')
            for iColor=1:size(EPwaves.handles.waves.eventLines,2)
                for iLine=1:length(EPwaves.handles.waves.eventLines{iChan,iColor})
                    refreshdata(EPwaves.handles.waves.eventLines{iChan,iColor}(iLine),'caller');
                end
            end
        end
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function changeEvent(~,~)
% %update displays after changing the events setting
% 
% global EPmanualEdit
% 
% EPmanualEdit.eventValue=get(EPmanualEdit.handles.events,'Value');
% 
% updateDisplay

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function scanEditMode(~,~)
%toggle whether clicking on channel marks it as bad or expands it.

global EPmanualEdit EPwaves

for iChan=1:length(EPwaves.handles.waves.hWave)
    if get(EPmanualEdit.handles.edit,'Value')
        set(EPwaves.handles.waves.hWave(iChan),'ButtonDownFcn',@badChannel);
        for iLine=1:length(EPwaves.handles.waves.hLines{iChan})
            set(EPwaves.handles.waves.hLines{iChan}(iLine),'ButtonDownFcn',@badChannel);
        end
    else
        set(EPwaves.handles.waves.hWave(iChan),'ButtonDownFcn',@ep_expandChan);
        for iLine=1:length(EPwaves.handles.waves.hLines{iChan})
            set(EPwaves.handles.waves.hLines{iChan}(iLine),'ButtonDownFcn',@ep_expandChan);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Confirms closing of the Edit pane.
function confirmCloseEdit(src,event)

global EPmanualEdit

if isempty(EPmanualEdit)
    selection = questdlg('Close The Manual Edit Window Without Saving Changes?', ...
        'Confirm Closure', ...
        'Yes','No','Yes');
    switch selection
        case 'Yes'
            cancelEdit
        case 'No'
            return
    end
else
    cancelEdit
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Confirms closing of the Edit pane.
function updateManualEdit(src,event)

global EPmanualEdit EPwaves EPmain

set(EPmanualEdit.handles.cell,'Value',EPmain.view.cell(1));
set(EPmanualEdit.handles.trialNum,'String',get(EPmain.handles.view.trial(1),'String'));
set(EPmanualEdit.handles.trialNum,'Value',get(EPmain.handles.view.trial(1),'Value'));

if EPmain.view.cell(1) ==1
    set(EPmanualEdit.handles.leftCell,'enable','off');
else
    set(EPmanualEdit.handles.leftCell,'enable','on');
end
if EPmain.view.cell(1) == length(EPmain.view.theCells{EPmain.view.dataset(1)})
    set(EPmanualEdit.handles.rightCell,'enable','off');
else
    set(EPmanualEdit.handles.rightCell,'enable','on');
end
if EPmain.view.trial(1) ==1
    set(EPmanualEdit.handles.leftTrial,'enable','off');
else
    set(EPmanualEdit.handles.leftTrial,'enable','on');
end
if EPmain.view.trial(1) == length(get(EPmanualEdit.handles.trialNum,'String'))
    set(EPmanualEdit.handles.rightTrial,'enable','off');
else
    set(EPmanualEdit.handles.rightTrial,'enable','on');
end

[EPwaves.legendNames, ~]=ep_waveNames;