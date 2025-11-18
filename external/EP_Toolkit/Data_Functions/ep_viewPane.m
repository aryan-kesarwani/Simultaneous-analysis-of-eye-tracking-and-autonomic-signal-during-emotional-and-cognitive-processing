function ep_viewPane
% ep_viewPane - function ep_viewPane;
% GUI interface for the View Pane.
%

%History
%  by Joseph Dien (7/7/19)
%  jdien07@mac.com
%
% bugfix 7/23/19 JD
% Fixed crash when changing datasets and there was no RT data in the prior dataset.
% Fixed manually changing voltage settings not working.
%
% bugfix 8/27/19 JD
% Fixed crash because Scan button was enabled while -all- setting chosen for the main dataset.
%
% bugfix 9/17/19 JD
% Fixed crash when color other than blue had -all- selected for cells or subjects and the data had events and then a dataset was changed.
%
% bugfix 9/23/19 JD
% Fixed not defaulting to correct subject for colors other than blue when changing datasets.
% Fixed defaults for factors and subjects when switching dataset types to a factor dataset.
%
% modified 11/5/19 JD
% Added sessNums sessNames fields.
%
% modified 12/11/19 JD
% Added average-wise std measure as a banding option.
%
% modified 2/25/20 JD
% Added support for viewing BOSC data.
%
% bugfix & modified 4/24/20 JD
% Added support for up to eight colors in waveform figures.
% Fixed View pane voltage values not updating correctly when there are colors that have been set to "none."
% Added preference buttons to panes.
%
% bugfix 7/20/20 JD
% Not computing max sample correctly when data is continuous and overall length is less than one second.
% Fixed not being able to plot multiple continuous datasets.
%
% bugfix 3/3/21 JD
% Fixed problems after using all eight colors to display factor results.
%
% bugfix 12/8/22 JD
% Fixed crash in user interface when the data to be displayed are not voltage measures or continuous data.
%
% bugfix 2/2/25 JD
% Fixed bands options not working for averaged data when anything like noise or SD were dropped during averaging.
%
% modified 4/10/25 JD
% Added support for lastChange field to determine if the pane should be reset.
% Fixed various things not restricting themselves to just the displayed colors.
%
% bugfix 4/30/25 JD
% Fixed Trials popup menu disappearing when data settings changed to one where the current option (e.g., -CI-) is not applicable.
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

global EPdataset EPmain

set(EPmain.handles.hMainWindow,'Name', 'View EEG');
refresh

if ~isempty(EPmain.view)
    if length(EPmain.view.lastChange) ~=length(EPdataset.dataset)
        EPmain.view=[];
    else
        for iDataset=1:length(EPdataset.dataset)
            if ~strcmp(EPdataset.dataset(iDataset).lastChange,EPmain.view.lastChange{iDataset})
                EPmain.view=[];
                break
            end
        end
    end
end

if isempty(EPmain.view)    
    %no option to switch to this panel unless there are data to view.
    for iDataset=1:length(EPdataset.dataset)
        EPmain.view.theFileNames{iDataset,1}=EPdataset.dataset(iDataset).dataName;
        EPmain.view.lastChange{iDataset}=EPdataset.dataset(iDataset).lastChange;
        [uniqueCells waveOrder m]=unique(EPdataset.dataset(iDataset).cellNames,'first'); %waveOrder is the true order of the unique cells (wave numbers of the first occurences)
        [a b cellOrder]=unique(waveOrder,'first'); %cellOrder is the true order of the unique cells (numbered by cells)
        EPmain.view.theCells{iDataset}(cellOrder)=uniqueCells;
    end
    
    EPmain.view.trialList=cell(EPmain.maxColors,1);
    EPmain.view.allTrials=zeros(EPmain.maxColors,1); %1=all trials/subs,2=trial/subs erpimage,3=all facs,4=fac erpimage,5=all cells,6=cells erpimage
    EPmain.view.allCells=cell(EPmain.maxColors,1);
    EPmain.view.allCellsList={'-none-';'-GFP-';'-noise-';'-Trl StDev-';'-CI-';'-Sub StDev-';'-all-';'-erpimage-'};
    EPmain.view.correl=zeros(EPmain.maxColors,1);
    EPmain.view.STS=zeros(EPmain.maxColors,1);
    EPmain.view.BSC=zeros(EPmain.maxColors,1);
    EPmain.view.rel=zeros(EPmain.maxColors,1);
    EPmain.view.subject=zeros(EPmain.maxColors,1);
    EPmain.view.factor=zeros(EPmain.maxColors,1);
    EPmain.view.dataset=zeros(EPmain.maxColors,1);
    EPmain.view.cell=zeros(EPmain.maxColors,1);
    EPmain.view.trial=zeros(EPmain.maxColors,1);
    EPmain.view.type=cell(EPmain.maxColors,1);
    %initially set up the colors
    for iColor=1:EPmain.maxColors
        if iColor <= length(unique(EPdataset.dataset(end).cellNames))
            EPmain.view.dataset(iColor)=length(EPdataset.dataset);
            EPmain.view.cell(iColor)=min(iColor,length(EPmain.view.theCells{EPmain.view.dataset(iColor)}));
            EPmain.view.trial(iColor)=1;
            EPmain.view.allCells{iColor}='-none-';
            GAVlist=find(strcmp('GAV',EPdataset.dataset(EPmain.view.dataset(iColor)).subTypes) & cellfun(@isempty,strfind(EPdataset.dataset(EPmain.view.dataset(iColor)).subNames,'autoPCA')));
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)
                if ~isempty(GAVlist)
                    EPmain.view.subject(iColor)=GAVlist(end);
                else
                    EPmain.view.subject(iColor)=length(EPdataset.dataset(EPmain.view.dataset(iColor)).subNames);
                end
            else
                STSlist=find(strcmp('sampleTest',EPdataset.dataset(EPmain.view.dataset(iColor)).subNames));
                GAVonly=setdiff(GAVlist,STSlist);
                if ~isempty(GAVonly)
                    EPmain.view.subject(iColor)=GAVonly(end);
                else
                    EPmain.view.subject(iColor)=1;
                end
            end
            if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'single_trial')
                EPmain.view.trialList{iColor}=EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames));
            end
            EPmain.view.factor(iColor)=1;
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames) && (length(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)>1)
                EPmain.view.factor(iColor,1)=length(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)+1;
            end
            
            EPmain.view.rel(iColor)=~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).relNames);
            EPmain.view.correl(iColor)=0;
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames)
                if EPmain.view.rel(iColor)
                    EPmain.view.correl(iColor)=1; %coherence measure is correlation and phase-lock is unitless 0 to 1
                end
            end
            if strcmp('STS',EPdataset.dataset(EPmain.view.dataset(iColor)).cellTypes{EPmain.view.cell(iColor)})
                EPmain.view.STS(iColor)=1;
            else
                EPmain.view.STS(iColor)=0;
            end
            if any(strcmp('BSC',EPdataset.dataset(EPmain.view.dataset(iColor)).chanTypes))
                EPmain.view.BSC(iColor)=1;
            else
                EPmain.view.BSC(iColor)=0;
            end
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames)
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames)
                    EPmain.view.type{iColor}='TFT';
                else
                    EPmain.view.type{iColor}='FFT';
                end
            else
                EPmain.view.type{iColor}='VLT'; %the type regardless of whether it is BSC or not
            end
        else
            EPmain.view.dataset(iColor,1)=length(EPdataset.dataset)+1; %all the cells shown so just leave this one blank.
        end
        EPmain.view.changeDatasetFlag(iColor)=0;
        EPmain.view.changeFlag(iColor)=1;
    end
    EPmain.view.FFTunits=4;
    EPmain.view.marker1=[];
    EPmain.view.marker2=[];
    EPmain.view.plotMVmin=zeros(EPmain.maxColors,1);
    EPmain.view.plotMVmax=zeros(EPmain.maxColors,1);
    EPmain.view.plotMVminAbs=zeros(EPmain.maxColors,1);
    EPmain.view.plotMVmaxAbs=zeros(EPmain.maxColors,1);
    EPmain.view.startSamp=0;
    EPmain.view.endSamp=0;
    EPmain.view.startHz=0;
    EPmain.view.endHz=0;
    EPmain.view.binHz=0;
    EPmain.view.edited.bottomVolt=0;
    EPmain.view.edited.topVolt=0;
    EPmain.view.edited.startSamp=0;
    EPmain.view.edited.endSamp=0;
    EPmain.view.edited.startHz=0;
    EPmain.view.edited.endHz=0;
    EPmain.view.manual.bottomVolt=0;
    EPmain.view.manual.topVolt=0;
    EPmain.view.manual.startSamp=0;
    EPmain.view.manual.endSamp=0;
    EPmain.view.manual.startHz=0;
    EPmain.view.manual.endHz=0;
    EPmain.view.flexMode=strcmp(EPdataset.dataset(end).timeUnits,'per');
    
    if ~isempty(EPdataset.dataset(end).freqNames)
        if ~isempty(EPdataset.dataset(end).timeNames)
            EPmain.view.dataTransform='TFT'; %the dataset will count as TFT even if BOSC, although non-BOSC can be displayed at the same time.
        else
            EPmain.view.dataTransform='FFT';
        end
    else
        EPmain.view.dataTransform='VLT';
    end
    
    EPmain.view.dataType=EPdataset.dataset(end).dataType;
    
    EPmain.view.eventList=cell(0);
    EPmain.view.RT=0;
    for iColor=1:EPmain.numColors
        if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
            if isempty(EPmain.view.eventList) || ~any(ismember(EPmain.view.dataset(iColor),EPmain.view.dataset(1:iColor-1)))
%                 if ismember(EPmain.view.allTrials(iColor),[5 6])
                    cellList=[1:length(EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames)];
%                 else
%                     cellList=EPmain.view.cell(iColor);
%                 end
                if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'single_trial')
%                    if ismember(EPmain.view.allTrials(iColor),[1 2])
                        cellList=EPmain.view.trialList{iColor};
%                     else
%                         cellList=EPmain.view.trial(iColor);
%                     end
                    subList=EPmain.view.subject(iColor);
                else
%                     if ismember(EPmain.view.allTrials(iColor),[1 2])
                        subList=[1:length(EPdataset.dataset(EPmain.view.dataset(iColor)).subNames)];
%                     else
%                         subList=EPmain.view.subject(iColor);
%                     end
                end
                typeValues=cell(0);
                for iSub=1:length(subList)
                    theSub=subList(iSub);
                    for iCell=1:length(cellList)
                        theCell=cellList(iCell);
                        if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'single_trial')
                            theEvents=EPdataset.dataset(EPmain.view.dataset(iColor)).events{theSub,theCell};
                        else
                            theEvents=EPdataset.dataset(EPmain.view.dataset(iColor)).events{theSub,theCell};
                        end
                        for iEvent=1:length(theEvents)
                            typeValues{end+1}=[num2str(theEvents(iEvent).type) '-' num2str(theEvents(iEvent).value)];
                        end
                    end
                end
                if any(strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).trialSpecNames,'RT'))
                    EPmain.view.RT=1;
                end
            end
        end
        EPmain.view.eventList(end+1:end+length(typeValues),1)=typeValues;
    end
    EPmain.view.eventList=unique(EPmain.view.eventList);
    EPmain.view.events=length(EPmain.view.eventList)+1; %add -none-
    if EPmain.view.RT
        EPmain.view.events=EPmain.view.events+1; %add -RT-
        EPmain.view.eventList{end+1}='-RT-';
    end
    if EPmain.view.events > 2
        EPmain.view.events=EPmain.view.events+1;  %add -all-
        EPmain.view.eventList{end+1}='-all-';
    end
    EPmain.view.eventList{end+1}='-none-';
end

activeList=find(EPmain.view.dataset(1:EPmain.numColors) <= length(EPdataset.dataset));

if any(EPmain.view.changeDatasetFlag) %if the dataset was just changed
    theColor=min(find(EPmain.view.changeDatasetFlag));
    if EPmain.view.dataset(theColor) <= length(EPdataset.dataset)
        if strcmp(EPdataset.dataset(EPmain.view.dataset(theColor)).dataType,'single_trial')
            EPmain.view.cell(theColor)=1;
            EPmain.view.trialList{theColor}=EPdataset.dataset(EPmain.view.dataset(theColor)).trialNames(strcmp(EPmain.view.theCells{EPmain.view.dataset(theColor)}(EPmain.view.cell(theColor)),EPdataset.dataset(EPmain.view.dataset(theColor)).cellNames));
        elseif (EPmain.view.changeDatasetFlag(theColor) <= length(EPdataset.dataset)) && any(strcmp(EPmain.view.theCells{EPmain.view.changeDatasetFlag(theColor)}(EPmain.view.cell(theColor)),EPmain.view.theCells{EPmain.view.dataset(theColor)}))
            EPmain.view.cell(theColor)=find(strcmp(EPmain.view.theCells{EPmain.view.changeDatasetFlag(theColor)}(EPmain.view.cell(theColor)),EPmain.view.theCells{EPmain.view.dataset(theColor)}));
        else
            EPmain.view.cell(theColor)=min(theColor,length(EPdataset.dataset(EPmain.view.dataset(theColor)).cellNames));
        end
        EPmain.view.subject(theColor)=1;
        GAVlist=find(strcmp('GAV',EPdataset.dataset(EPmain.view.dataset(theColor)).subTypes) & cellfun(@isempty,strfind(EPdataset.dataset(EPmain.view.dataset(theColor)).subNames,'autoPCA')));
        if ~isempty(EPdataset.dataset(EPmain.view.dataset(theColor)).facNames)
            if ~isempty(GAVlist)
                EPmain.view.subject(theColor)=GAVlist(end);
            else
                EPmain.view.subject(theColor)=length(EPdataset.dataset(EPmain.view.dataset(theColor)).subNames);
            end
        else
            
            STSlist=find(strcmp('sampleTest',EPdataset.dataset(EPmain.view.dataset(theColor)).subNames));
            GAVonly=setdiff(GAVlist,STSlist);
            if ~isempty(GAVonly)
                EPmain.view.subject(theColor)=GAVonly(end);
            end
        end
        %gavSubs=find(strcmp(EPdataset.dataset(EPmain.view.dataset(theColor)).subTypes,'GAV'));
        %                 if strcmp(EPdataset.dataset(EPmain.view.dataset(theColor)).dataType,'average') && (length(gavSubs)==1)
        %                     EPmain.view.subject(theColor)=gavSubs;
        %                 end
        EPmain.view.factor(theColor)=1;
        if ~isempty(EPdataset.dataset(EPmain.view.dataset(theColor)).facNames) && (length(EPdataset.dataset(EPmain.view.dataset(theColor)).facNames)>1)
            theFac=length(EPdataset.dataset(EPmain.view.dataset(theColor)).facNames);
            if any(EPmain.view.allTrials(1:EPmain.numColors)==4)
                theFac=theFac+2;
            elseif any(EPmain.view.allTrials(1:EPmain.numColors)==3) || all(EPmain.view.allTrials(1:EPmain.numColors)==0)
                    theFac=theFac+1;
            end
            EPmain.view.factor(theColor)=theFac;
        end
        
        EPmain.view.trial(theColor)=1;
        
        EPmain.view.allTrials(theColor)=0;
        
        EPmain.view.rel(theColor)=~isempty(EPdataset.dataset(EPmain.view.dataset(theColor)).relNames);
        
        %if changed transform type (voltage, FFT, or TFT), then change over all the other colors too.
        EPmain.view.correl(theColor)=0;
        if EPmain.view.rel(theColor)
            EPmain.view.correl(theColor)=1; %coherence measure is correlation and phase-lock is unitless 0 to 1
        end
        if strcmp('STS',EPdataset.dataset(EPmain.view.dataset(theColor)).cellTypes{EPmain.view.cell(theColor)})
            EPmain.view.STS(theColor)=1;
        else
            EPmain.view.STS(theColor)=0;
        end
        if any(strcmp('BSC',EPdataset.dataset(EPmain.view.dataset(theColor)).chanTypes))
            EPmain.view.BSC(theColor)=1;
        else
            EPmain.view.BSC(theColor)=0;
        end
        if ~isempty(EPdataset.dataset(EPmain.view.dataset(theColor)).freqNames)
            if ~isempty(EPdataset.dataset(EPmain.view.dataset(theColor)).timeNames)
                newType='TFT';
            else
                newType='FFT';
            end
        else
            newType='VLT';
        end
        EPmain.view.type{theColor}=newType;
        
        if (~strcmp(newType,EPmain.view.dataTransform) && ~(EPmain.view.BSC(theColor) || all(EPmain.view.BSC(setdiff(activeList,theColor))))) || (xor(strcmp(EPmain.view.dataType,'continuous'),strcmp(EPdataset.dataset(EPmain.view.dataset(theColor)).dataType,'continuous'))) || (xor(strcmp(EPdataset.dataset(EPmain.view.dataset(theColor)).timeUnits,'per'),EPmain.view.flexMode))
            %if new dataset is different from existing ones then change them too to keep consistent
            %BSC can be displayed with all three types
            disp('The new dataset is different in type than current datasets so all four colors will all be replaced with the new dataset.');
            EPmain.view.dataTransform=newType;
            EPmain.view.dataType=EPdataset.dataset(EPmain.view.dataset(theColor)).dataType;
            EPmain.view.flexMode=strcmp(EPdataset.dataset(EPmain.view.dataset(theColor)).timeUnits,'per');
            for iColor=1:EPmain.maxColors
                if iColor ~= theColor
                    if iColor <= length(EPmain.view.theCells{EPmain.view.dataset(theColor)})
                        EPmain.view.dataset(iColor)=EPmain.view.dataset(theColor);
                        EPmain.view.cell(iColor)=min(iColor,length(EPmain.view.theCells{EPmain.view.dataset(iColor)}));
                        if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'single_trial')
                            EPmain.view.trialList{iColor}=EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames));
                        end
                        EPmain.view.trial(iColor)=1;
                        EPmain.view.allCells{iColor}='-none-';
                        EPmain.view.subject(iColor)=1;
                        GAVlist=find(strcmp('GAV',EPdataset.dataset(EPmain.view.dataset(iColor)).subTypes) & cellfun(@isempty,strfind(EPdataset.dataset(EPmain.view.dataset(iColor)).subNames,'autoPCA')));
                        if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)
                            if ~isempty(GAVlist)
                                EPmain.view.subject(iColor)=GAVlist(end);
                            else
                                EPmain.view.subject(iColor)=length(EPdataset.dataset(EPmain.view.dataset(iColor)).subNames);
                            end
                        else
                            
                            STSlist=find(strcmp('sampleTest',EPdataset.dataset(EPmain.view.dataset(iColor)).subNames));
                            GAVonly=setdiff(GAVlist,STSlist);
                            if ~isempty(GAVonly)
                                EPmain.view.subject(iColor)=GAVonly(end);
                            end
                        end
                        EPmain.view.factor(iColor)=1;
                        if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames) && (length(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)>1)
                            theFac=length(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames);
                            if any(EPmain.view.allTrials(1:EPmain.numColors)==4)
                                theFac=theFac+2;
                            elseif any(EPmain.view.allTrials(1:EPmain.numColors)==3) || all(EPmain.view.allTrials(1:EPmain.numColors)==0)
                                theFac=theFac+1;
                            end
                            EPmain.view.factor(iColor)=theFac;
                        end
                        
                        EPmain.view.changeFlag(iColor)=1;
                        EPmain.view.correl(iColor)=EPmain.view.correl(theColor);
                        EPmain.view.STS(iColor)=EPmain.view.STS(theColor);
                        EPmain.view.rel(iColor)=EPmain.view.rel(theColor);
                    else
                        EPmain.view.dataset(iColor)=length(EPdataset.dataset)+1; %all the cells shown so just leave this one blank.
                    end
                end
            end
        else
            if ~strcmp('TFT',EPmain.view.dataTransform) && strcmp('TFT',newType)
                EPmain.view.dataTransform='TFT'; %if the new one is BSC, then change the type to TFT
            end
        end
    else
        EPmain.view.correl(theColor)=0;
        EPmain.view.rel(theColor)=0;
        EPmain.view.STS(theColor)=0;
        EPmain.view.BSC(theColor)=0;
    end
    
    previousEvent=EPmain.view.eventList{EPmain.view.events};
    EPmain.view.eventList=cell(0);
    EPmain.view.RT=0;
    typeValues=cell(0);
    %set up the list of events
    for iColor=1:EPmain.numColors
        if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
            if isempty(EPmain.view.eventList) || ~any(ismember(EPmain.view.dataset(iColor),EPmain.view.dataset(1:iColor-1)))
%                 if ismember(EPmain.view.allTrials(iColor),[5 6])
                    cellList=[1:length(EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames)];
%                 else
%                     cellList=EPmain.view.cell(iColor);
%                 end
                if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'single_trial')
%                    if ismember(EPmain.view.allTrials(iColor),[1 2])
                        cellList=EPmain.view.trialList{iColor};
%                     else
%                         cellList=EPmain.view.trial(iColor);
%                     end
                    subList=EPmain.view.subject(iColor);
                else
%                     if ismember(EPmain.view.allTrials(iColor),[1 2])
                        subList=[1:length(EPdataset.dataset(EPmain.view.dataset(iColor)).subNames)];
%                     else
%                         subList=EPmain.view.subject(iColor);
%                     end
                end
                typeValues=cell(0);
                for iSub=1:length(subList)
                    theSub=subList(iSub);
                    for iCell=1:length(cellList)
                        theCell=cellList(iCell);
                        if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'single_trial')
                            theEvents=EPdataset.dataset(EPmain.view.dataset(iColor)).events{theSub,theCell};
                        else
                            theEvents=EPdataset.dataset(EPmain.view.dataset(iColor)).events{theSub,theCell};
                        end
                        for iEvent=1:length(theEvents)
                            typeValues{end+1}=[num2str(theEvents(iEvent).type) '-' num2str(theEvents(iEvent).value)];
                        end
                    end
                end
                if any(strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).trialSpecNames,'RT'))
                    EPmain.view.RT=1;
                end
            end
        end
        EPmain.view.eventList(end+1:end+length(typeValues),1)=typeValues;
    end
    EPmain.view.eventList=unique(EPmain.view.eventList);
    EPmain.view.events=length(EPmain.view.eventList)+1;
    if EPmain.view.RT
        EPmain.view.events=EPmain.view.events+1; %add -RT-
        EPmain.view.eventList{end+1}='-RT-';
    end
    if EPmain.view.events > 2
        EPmain.view.events=EPmain.view.events+1;  %add -all-
        EPmain.view.eventList{end+1}='-all-';
    end
    EPmain.view.eventList{end+1}='-none-';
    if any(strcmp(previousEvent,EPmain.view.eventList))
        EPmain.view.events=find(strcmp(EPmain.view.eventList,previousEvent));
    end
    
    EPmain.view.plotMVmin(theColor)=0;
    EPmain.view.plotMVmax(theColor)=0;
    EPmain.view.plotMVminAbs(theColor)=0;
    EPmain.view.plotMVmaxAbs(theColor)=0;
    EPmain.view.changeDatasetFlag(theColor)=0;
    EPmain.view.changeFlag(theColor)=1;
    EPmain.view.edited.bottomVolt=0;
    EPmain.view.edited.topVolt=0;
    EPmain.view.edited.startSamp=0;
    EPmain.view.edited.endSamp=0;
    EPmain.view.edited.startHz=0;
    EPmain.view.edited.endHz=0;
    EPmain.view.manual.bottomVolt=0;
    EPmain.view.manual.topVolt=0;
    EPmain.view.manual.startSamp=0;
    EPmain.view.manual.endSamp=0;
    EPmain.view.manual.startHz=0;
    EPmain.view.manual.endHz=0;
end

marker1=num2str(EPmain.view.marker1);
marker2=num2str(EPmain.view.marker2);

%compute voltages and epoch for the changed colors
if any(EPmain.view.changeFlag)
    for iColor=1:EPmain.maxColors
        if EPmain.view.changeFlag(iColor)
            EPmain.view.changeFlag(iColor)=0;
            EPmain.view.allTrials(iColor)=0;
            
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                
                if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'continuous') %continuous data
                    theCell=EPmain.view.trial(iColor);
                    if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVmin(theCell,1,EPmain.view.factor(iColor)))
                        EPmain.view.plotMVmin(iColor)=EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVmin(theCell,1,EPmain.view.factor(iColor));
                    else
                        EPmain.view.plotMVmin(iColor)=[];
                    end
                    if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVmax(theCell,1,EPmain.view.factor(iColor)))
                        EPmain.view.plotMVmax(iColor)=EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVmax(theCell,1,EPmain.view.factor(iColor));
                    else
                        EPmain.view.plotMVmax(iColor)=[];
                    end
                    if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVminAbs(theCell,1,EPmain.view.factor(iColor)))
                        EPmain.view.plotMVminAbs(iColor)=EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVminAbs(theCell,1,EPmain.view.factor(iColor));
                    else
                        EPmain.view.plotMVminAbs(iColor)=[];
                    end
                    if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVmaxAbs(theCell,1,EPmain.view.factor(iColor)))
                        EPmain.view.plotMVmaxAbs(iColor)=EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVmaxAbs(theCell,1,EPmain.view.factor(iColor));
                    else
                        EPmain.view.plotMVmaxAbs(iColor)=[];
                    end
                    
                elseif strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'average') %averaged data
                    numSubs=length(EPdataset.dataset(EPmain.view.dataset(iColor)).subNames);
                    if EPmain.view.subject(iColor) > numSubs
                        theSubject=[1:numSubs];
                        if EPmain.view.subject(iColor) == numSubs+1
                            EPmain.view.allTrials(iColor)=1; %'all'
                        end
                        if EPmain.view.subject(iColor) == numSubs+2
                            EPmain.view.allTrials(iColor)=2; %'erpimage'
                        end
                    else
                        theSubject=EPmain.view.subject(iColor);
                    end
                    numFacs=length(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames);
                    if (EPmain.view.factor(iColor) > numFacs) && (numFacs > 0)
                        theFactor=[1:numFacs];
                        if EPmain.view.factor(iColor) == numFacs+1
                            EPmain.view.allTrials(iColor)=3; %'all'
                        end
                        if EPmain.view.factor(iColor) == numFacs+2
                            EPmain.view.allTrials(iColor)=4; %'erpimage'
                        end
                    else
                        theFactor=EPmain.view.factor(iColor);
                    end
                    numCells=length(EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames);
                    theCell=EPmain.view.cell(iColor);
                    if theCell > numCells
                        theCell=[1:numCells];
                        if EPmain.view.cell(iColor) == numCells+1
                            EPmain.view.allTrials(iColor)=5; %'all'
                        end
                        if EPmain.view.cell(iColor) == numCells+2
                            EPmain.view.allTrials(iColor)=6; %'erpimage'
                        end
                    end
                    EPmain.view.plotMVmin(iColor)=min(EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVmin(theCell,theSubject,theFactor));
                    EPmain.view.plotMVmax(iColor)=max(EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVmax(theCell,theSubject,theFactor));
                    EPmain.view.plotMVminAbs(iColor)=min(EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVminAbs(theCell,theSubject,theFactor));
                    EPmain.view.plotMVmaxAbs(iColor)=max(EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVmaxAbs(theCell,theSubject,theFactor));
                    
                else  %if single_trial data
                    EPmain.view.trialList{iColor}=EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames));
                    if EPmain.view.trial(iColor) > length(EPmain.view.trialList{iColor})
                        trialList=EPmain.view.trialList{iColor};
                        if EPmain.view.trial(iColor) == length(EPmain.view.trialList{iColor})+1
                            EPmain.view.allTrials(iColor)=1; %'all'
                        end
                        if EPmain.view.trial(iColor) == length(EPmain.view.trialList{iColor})+2
                            EPmain.view.allTrials(iColor)=2; %'erpimage'
                        end
                    else
                        trialList=EPmain.view.trial(iColor);
                    end
                    numFacs=length(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames);
                    if (EPmain.view.factor(iColor) > numFacs) && (numFacs > 0)
                        theFactor=[1:numFacs];
                        if EPmain.view.factor(iColor) == numFacs+1
                            EPmain.view.allTrials(iColor)=3; %'all'
                        end
                        if EPmain.view.factor(iColor) == numFacs+2
                            EPmain.view.allTrials(iColor)=4; %'erpimage'
                        end
                    else
                        theFactor=EPmain.view.factor(iColor);
                    end
                    trialsInCell=find(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames));
                    EPmain.view.plotMVmin(iColor)=min(EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVmin(trialsInCell(trialList),1,theFactor));
                    EPmain.view.plotMVmax(iColor)=max(EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVmax(trialsInCell(trialList),1,theFactor));
                    EPmain.view.plotMVminAbs(iColor)=min(EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVminAbs(trialsInCell(trialList),1,theFactor));
                    EPmain.view.plotMVmaxAbs(iColor)=max(EPdataset.dataset(EPmain.view.dataset(iColor)).plotMVmaxAbs(trialsInCell(trialList),1,theFactor));
                end
            end
        end
    end
    
    
    if strcmp('FFT',EPmain.view.dataTransform)
        EPmain.view.startSamp=NaN;
        EPmain.view.endSamp=NaN;
    else
        EPmain.view.startSamp=0;
        EPmain.view.endSamp=0;
        % Time is counted for an epoch as, say, -200 to 800 ms.  In this case the samples have some length to them depending on
        % the digitization rate, such as 4ms for a 250Hz sampling rate.  The points prior to the baseline are counted from the
        % left side of the baseline sample and the points following the baseline are counted from the right side of the baseline
        % sample.  Thus, when the actual samples are calculated, one first adds the baseline to the numbers, yielding 0 to
        % 1000ms.  One then must add the length of the baseline to the first number to obtain 4 to 1000ms.  Then one would
        % divide the numbers by the sample size, resulting in 1 to 250 samples, which are the correct samples to use.
        % For continuous data or other such data where there is no baseline, one must by this convention start with 0 ms,
        % as in 0-1000ms and 1000-2000ms.
        % Allows for different sized epochs.  View Waves will just drop the excess time points based on the time names.
        for iColor=1:EPmain.numColors
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                if ~isnan(EPmain.view.startSamp) && ~isnan(EPmain.view.endSamp)
                    if ~EPmain.view.startSamp && ~EPmain.view.endSamp %if neither has been set yet
                        if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'continuous') %continuous data
                            EPmain.view.startSamp=1000*(EPmain.view.trial(iColor)-1); %for continuous data, each "trial" is 1000 ms
                            if length(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames)<EPdataset.dataset(EPmain.view.dataset(iColor)).Fs
                                EPmain.view.endSamp=(1000/EPdataset.dataset(EPmain.view.dataset(iColor)).Fs)*length(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames);
                            else
                                EPmain.view.endSamp=1000*(EPmain.view.trial(iColor));
                            end
                        else
                            EPmain.view.startSamp=min(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames);
                            EPmain.view.endSamp=max(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames)+(1000/EPdataset.dataset(EPmain.view.dataset(iColor)).Fs);
                        end
                    else
                        if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'continuous') %continuous data
                            newStartSamp=1000*(EPmain.view.trial(iColor)-1); %for continuous data, each "trial" is 1000 ms
                            if length(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames)<EPdataset.dataset(EPmain.view.dataset(iColor)).Fs
                                newEndSamp=(1000/EPdataset.dataset(EPmain.view.dataset(iColor)).Fs)*length(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames);
                            else
                                newEndSamp=1000*(EPmain.view.trial(iColor));
                            end

                            if EPmain.view.startSamp ~= newStartSamp
                                EPmain.view.startSamp=NaN; %set to NaN if the samples are different for the different datasets
                            end
                            if length(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames)<EPdataset.dataset(EPmain.view.dataset(iColor)).Fs
                                if EPmain.view.endSamp ~= (newEndSamp/EPdataset.dataset(EPmain.view.dataset(iColor)).Fs)*length(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames)
                                    EPmain.view.endSamp=NaN; %set to NaN if the samples are different for the different datasets
                                end
                            else
                                if EPmain.view.endSamp ~= newEndSamp
                                    EPmain.view.endSamp=NaN; %set to NaN if the samples are different for the different datasets
                                end
                            end
                        else
                            EPmain.view.startSamp = min([EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames; EPmain.view.startSamp]);
                            EPmain.view.endSamp= max([EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames+(1000/EPdataset.dataset(EPmain.view.dataset(iColor)).Fs); EPmain.view.endSamp]);
                        end
                    end
                end
            end
        end
    end
    
    if strcmp('VLT',EPmain.view.dataTransform)
        EPmain.view.startHz=NaN;
        EPmain.view.endHz=NaN;
    else
        EPmain.view.startHz=0;
        EPmain.view.endHz=0;
        for iColor=1:EPmain.numColors
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames) %it is possible to be in TFT mode and yet for a dataset to not have frequency data if it is due to a BOSC dataset
                    if ~isnan(EPmain.view.startHz) && ~isnan(EPmain.view.endHz)
                        if ~EPmain.view.startHz && ~EPmain.view.endHz %if neither has been set yet
                            EPmain.view.startHz=min(EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames);
                            EPmain.view.endHz=max(EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames);
                            EPmain.view.binHz=EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames(2)-EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames(1); %making assumption all have the same bin resolution
                        else
                            if EPmain.view.startHz ~= min(EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames)
                                EPmain.view.startHz=NaN;
                            end
                            if EPmain.view.endHz ~= max(EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames)
                                EPmain.view.endHz=NaN;
                            end
                        end
                    end
                end
            end
        end
    end
end
if ~isempty(activeList)
    if all(EPmain.view.correl(activeList))
        plotMVmin=min(EPmain.view.plotMVmin(activeList & ~EPmain.view.STS(activeList) & ~EPmain.view.BSC(activeList)));
        plotMVmax=max(EPmain.view.plotMVmax(activeList & ~EPmain.view.STS(activeList)) & ~EPmain.view.BSC(activeList));
    elseif all(EPmain.view.BSC(activeList))
        plotMVmin=min(EPmain.view.plotMVmin(activeList));
        plotMVmax=max(EPmain.view.plotMVmax(activeList));
    else
        %if any datasets are not correlations, then ignore correlations when calculating min and max values
        if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'})) && (EPmain.view.FFTunits > 1) && any(~strcmp('VLT',EPmain.view.type(activeList)) & ~EPmain.view.BSC(activeList))
            plotMVmin=min(EPmain.view.plotMVminAbs(activeList(~EPmain.view.correl(activeList) & ~EPmain.view.STS(activeList) & ~EPmain.view.BSC(activeList))));
            plotMVmax=max(EPmain.view.plotMVmaxAbs(activeList(~EPmain.view.correl(activeList) & ~EPmain.view.STS(activeList) & ~EPmain.view.BSC(activeList))));
        else
            plotMVmin=min(EPmain.view.plotMVmin(activeList(~EPmain.view.correl(activeList) & ~EPmain.view.STS(activeList) & ~EPmain.view.BSC(activeList))));
            plotMVmax=max(EPmain.view.plotMVmax(activeList(~EPmain.view.correl(activeList) & ~EPmain.view.STS(activeList) & ~EPmain.view.BSC(activeList))));
        end
    end
else
    plotMVmin=0;
    plotMVmax=0;
end

EPmain.handles.view.prefs = uicontrol('Style', 'pushbutton', 'String', '•','FontSize',EPmain.fontsize,...
    'Position', [5 485 10 10], 'Callback', ['global EPmain;','EPmain.prefReturn=''view'';','EPmain.mode=''preferenceView'';','ep(''start'');']);

menuStart=20;
menuX=150;
menuY=18;

for iColor=1:EPmain.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        numSubs=length(EPdataset.dataset(EPmain.view.dataset(iColor)).subNames);
        numVsubs=max(0,size(EPdataset.dataset(EPmain.view.dataset(iColor)).GAVsubs,1)-1);
        numRsubs=numSubs-numVsubs;
        numCells=length(EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames);
        numVcells=max(0,size(EPdataset.dataset(EPmain.view.dataset(iColor)).GAVsubs,2)-1);
        numRcells=numCells-numVcells;
    else
        numSubs=0;
        numVsubs=0;
        numRsubs=0;
        numCells=0;
        numVcells=0;
        numRcells=0;
    end
    if iColor <=4 %only four fit on the pane
        leftStart=menuStart;
        topStart=iColor;
    else
        leftStart=EPmain.panesize(1)+1+menuStart;
        topStart=iColor-4;
    end
    theDatasets=[EPmain.view.theFileNames; {'none'}];
    EPmain.handles.view.dataset(iColor) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',theDatasets,'ForegroundColor',EPmain.preferences.view.color(iColor).RGB,...
        'Value',EPmain.view.dataset(iColor),'Position',[leftStart 500-menuY*((topStart-1)*5+1) menuX menuY],...
        'Callback', ['pause(.2);','global EPmain;','EPmain.view.changeDatasetFlag(' num2str(iColor) ')=EPmain.view.dataset(' num2str(iColor) ');','tempVar=get(EPmain.handles.view.dataset(' num2str(iColor) '),''value'');','if tempVar ~=0,EPmain.view.dataset(' num2str(iColor) ')=tempVar;end;','if isempty(tempVar),EPmain.view.dataset(' num2str(iColor) ')=tempVar;end;','ep(''start'')']);
    
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        theCells=EPmain.view.theCells{EPmain.view.dataset(iColor)};
        if ~strcmp('TFT',EPmain.view.dataTransform) && length(theCells) > 1 && strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'average') && ~ismember(EPmain.view.allTrials(iColor),[1 2 3 4])
            theCells{end+1}='-all-';
            theCells{end+1}='-erpimage-';
        end
        EPmain.handles.view.cell(iColor) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',theCells,'ForegroundColor',EPmain.preferences.view.color(iColor).RGB,...
            'Value',EPmain.view.cell(iColor),'Position',[leftStart 500-menuY*((topStart-1)*5+2) menuX menuY],...
            'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.cell(' num2str(iColor) '),''value'');','if tempVar ~=0,EPmain.view.cell(' num2str(iColor) ')=tempVar;end;','if isempty(tempVar),EPmain.view.cell(' num2str(iColor) ')=tempVar;end;','EPmain.view.changeFlag(' num2str(iColor) ')=1;','ep(''start'')']);
        
        theSubs=EPdataset.dataset(EPmain.view.dataset(iColor)).subNames;
        if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).sessNums)
            for iSub=1:length(theSubs)
                if EPdataset.dataset(EPmain.view.dataset(iColor)).sessNums(iSub)
                    theSubs{iSub}=[theSubs{iSub} '-' EPdataset.dataset(EPmain.view.dataset(iColor)).sessNames{EPdataset.dataset(EPmain.view.dataset(iColor)).sessNums(iSub)}];
                end
            end
        end
        if ~strcmp('TFT',EPmain.view.dataTransform) && length(theSubs) > 1 && ~ismember(EPmain.view.allTrials(iColor),[3 4 5 6])
            theSubs{end+1}='-all-';
            theSubs{end+1}='-erpimage-';
        end
        EPmain.handles.view.subject(iColor) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',theSubs,'ForegroundColor',EPmain.preferences.view.color(iColor).RGB,...
            'Value',EPmain.view.subject(iColor),'Position',[leftStart 500-menuY*((topStart-1)*5+3) menuX menuY],...
            'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.subject(' num2str(iColor) '),''value'');','if tempVar ~=0,EPmain.view.subject(' num2str(iColor) ')=tempVar;end;','if isempty(tempVar),EPmain.view.subject(' num2str(iColor) ')=tempVar;end;','EPmain.view.changeFlag(' num2str(iColor) ')=1;','ep(''start'')']);
        
        if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames)
            theTrials=num2cell(EPmain.view.trialList{iColor});
            if ~strcmp('TFT',EPmain.view.dataTransform) && ~ismember(EPmain.view.allTrials(iColor),[3 4 5 6])
                theTrials{end+1}='-all-';
                theTrials{end+1}='-erpimage-';
            end
            EPmain.handles.view.trial(iColor) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',theTrials,'ForegroundColor',EPmain.preferences.view.color(iColor).RGB,...
                'Value',EPmain.view.trial(iColor),'Position',[leftStart 500-menuY*((topStart-1)*5+4) menuX menuY],...
                'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.trial(' num2str(iColor) '),''value'');','if tempVar ~=0,EPmain.view.trial(' num2str(iColor) ')=tempVar;end;','if isempty(tempVar),EPmain.view.trial(' num2str(iColor) ')=tempVar;end;','EPmain.view.changeFlag(' num2str(iColor) ')=1;','ep(''start'')']);
        elseif strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'continuous')
            theEpochs=cellstr(num2str([1:ceil(length(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames)/ceil(EPdataset.dataset(EPmain.view.dataset(iColor)).Fs))]'));
            for i=1:length(theEpochs)
                theEpochs{i}=['Sec: ' theEpochs{i}];
            end
            EPmain.handles.view.trial(iColor) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',theEpochs,'ForegroundColor',EPmain.preferences.view.color(iColor).RGB,...
                'Value',EPmain.view.trial(iColor),'Position',[leftStart 500-menuY*((topStart-1)*5+4) menuX menuY],...
                'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.trial(' num2str(iColor) '),''value'');','if tempVar ~=0,EPmain.view.trial(' num2str(iColor) ')=tempVar;end;','if isempty(tempVar),EPmain.view.trial(' num2str(iColor) ')=tempVar;end;','EPmain.view.changeFlag(' num2str(iColor) ')=1;','ep(''start'')']);
        elseif strcmp('VLT',EPmain.view.dataTransform)
            theTrials=cell(0);
            theTrials{end+1}='-none-';
            theTrials{end+1}='-GFP-';
            if ~EPmain.view.allTrials(iColor)
                if EPdataset.dataset(EPmain.view.dataset(iColor)).isNoise
                    theTrials{end+1}='-noise-';
                elseif strcmp(EPmain.view.allCells{iColor},'-noise-')
                    EPmain.view.allCells{iColor}='-none-';
                end
                if EPdataset.dataset(EPmain.view.dataset(iColor)).isCovAVE
                    theTrials{end+1}='-Trl StDev-';
                elseif strcmp(EPmain.view.allCells{iColor},'-Trl StDev-')
                    EPmain.view.allCells{iColor}='-none-';
                end
                if  (strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).subTypes{EPmain.view.subject(iColor)},'GAV') && ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).GAVsubs) && (EPmain.view.subject(iColor) > numRsubs) && ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).GAVsubs{EPmain.view.subject(iColor)-numRsubs+1,max(1,EPmain.view.cell(iColor)-numRcells+1),EPmain.view.factor(iColor)})) ||...
                        (strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).subTypes{EPmain.view.subject(iColor)},'AVG') && EPdataset.dataset(EPmain.view.dataset(iColor)).isCovAVE)
                    theTrials{end+1}='-CI-';
                elseif strcmp(EPmain.view.allCells{iColor},'-CI-')
                    EPmain.view.allCells{iColor}='-none-';
                end
                if (strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).subTypes{EPmain.view.subject(iColor)},'GAV') && ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).GAVsubs) && (EPmain.view.subject(iColor) > numRsubs) && ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).GAVsubs{EPmain.view.subject(iColor)-numRsubs+1,max(1,EPmain.view.cell(iColor)-numRcells+1),EPmain.view.factor(iColor)}))
                    theTrials{end+1}='-Sub StDev-';
                elseif strcmp(EPmain.view.allCells{iColor},'-Sub StDev-')
                    EPmain.view.allCells{iColor}='-none-';
                end
            end
            EPmain.handles.view.trial(iColor) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',theTrials,'ForegroundColor',EPmain.preferences.view.color(iColor).RGB,...
                'Value',find(strcmp(EPmain.view.allCells{iColor},theTrials)),'Position',[leftStart 500-menuY*((topStart-1)*5+4) menuX menuY],...
                'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.trial(' num2str(iColor) '),''value'');','tempVar2=get(EPmain.handles.view.trial(' num2str(iColor) '),''String'');','if tempVar ~=0,EPmain.view.allCells{' num2str(iColor) '}=tempVar2{tempVar};end;','if isempty(tempVar),EPmain.view.allCells(' num2str(iColor) ')=tempVar;end;','EPmain.view.changeFlag(' num2str(iColor) ')=1;','ep(''start'')']);
        else
            EPmain.handles.view.trial(iColor) = uicontrol('Style','text','HorizontalAlignment','left','String', '     No Trials','FontSize',EPmain.fontsize,...
                'ForegroundColor','blue','Position',[leftStart 500-menuY*((topStart-1)*5+4) menuX menuY]);
        end
        
        if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)
            theFacs=EPdataset.dataset(EPmain.view.dataset(iColor)).facNames;
            if length(theFacs) > 1 && ~ismember(EPmain.view.allTrials(iColor),[1 2 5 6])
                theFacs{end+1}='-all-';
                theFacs{end+1}='-erpimage-';
            end
            EPmain.handles.view.factor(iColor) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',theFacs,'ForegroundColor',EPmain.preferences.view.color(iColor).RGB,...
                'Value',EPmain.view.factor(iColor),'Position',[leftStart 500-menuY*((topStart-1)*5+5) menuX menuY],...
                'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.factor(' num2str(iColor) '),''value'');','if tempVar ~=0,EPmain.view.factor(' num2str(iColor) ')=tempVar;end;','if isempty(tempVar),EPmain.view.factor(' num2str(iColor) ')=tempVar;end;','EPmain.view.changeFlag(' num2str(iColor) ')=1;','ep(''start'')']);
        else
            EPmain.handles.view.factor(iColor) = uicontrol('Style','text','HorizontalAlignment','left','String', '     No Factors','FontSize',EPmain.fontsize,...
                'ForegroundColor',EPmain.preferences.view.color(iColor).RGB,'Position',[leftStart 500-menuY*((topStart-1)*5+5) menuX menuY]);
        end
    end
end

if length(EPmain.view.eventList) > 1
    EPmain.handles.view.events = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',EPmain.view.eventList,'ForegroundColor','magenta',...
        'Value',EPmain.view.events,'Position',[menuStart 500-menuY*21 menuX menuY],...
        'Callback', ['global EPmain;','tempVar=get(EPmain.handles.view.events,''value'');','if tempVar ~=0,EPmain.view.events=tempVar;end;','if isempty(tempVar),EPmain.view.events=tempVar;end;','ep(''start'')']);
else
    EPmain.handles.view.events = uicontrol('Style','text','HorizontalAlignment','left','String', '     No Events','FontSize',EPmain.fontsize,...
        'ForegroundColor','magenta','Position',[menuStart 500-menuY*21 menuX menuY]);
end

uicontrol('Style','frame',...
    'Position',[5 36 155 62]);

theMax=plotMVmax;
theMin=plotMVmin;

theTextColor.bottomVolt='black';
if EPmain.view.edited.bottomVolt
    theMin=EPmain.view.manual.bottomVolt;
    theTextColor.bottomVolt='blue';
end
theTextColor.topVolt='black';
if EPmain.view.edited.topVolt
    theMax=EPmain.view.manual.topVolt;
    theTextColor.topVolt='blue';
end

if all(EPmain.view.correl(activeList))
    EPmain.handles.view.FFTunits = uicontrol('Style','text','HorizontalAlignment','left','String', 'Corr','FontSize',EPmain.fontsize,...
        'ForegroundColor','black','Position',[10 77 50 20]);
elseif all(EPmain.view.BSC(activeList))
    EPmain.handles.view.FFTunits = uicontrol('Style','text','HorizontalAlignment','left','String', 'Bosc','FontSize',EPmain.fontsize,...
        'ForegroundColor','black','Position',[10 77 50 20]);
elseif any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'})) && any(~strcmp('VLT',EPmain.view.type(activeList)) & ~EPmain.view.BSC(activeList))
    %the cached min and max numbers are based on the greater of the min and max of the real and imaginary parts
    theMax=theMax/sqrt(EPmain.view.binHz); %convert to spectral density
    theMin=theMin/sqrt(EPmain.view.binHz); %convert to spectral density
    if EPmain.view.FFTunits > 2
        theMax=theMax.^2; %convert amplitude to power
        theMin=theMin.^2; %convert amplitude to power
    end
    if (EPmain.view.FFTunits == 4)
        theMax=log10(abs(theMax))*10; %convert to dB log scaling
        theMin=log10(abs(theMin))*10; %convert to dB log scaling
        if isinf(theMin) && (theMax > -100)
            theMin=-100; %log10 of zero is -inf.  Replace with -2 to maintain useful range.
            theTextColor.bottomVolt='blue';
        end
    end
    EPmain.handles.view.FFTunits = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'Value',EPmain.view.FFTunits,'Position',[1 77 65 20],...
        'String',{'cm','am','pw','dB'},...
        'Callback', ['global EPmain;','EPmain.view.FFTunits=get(EPmain.handles.view.FFTunits,''value'');','ep(''start'')'],...
        'TooltipString','Units for spectral data.');
    
else
    EPmain.handles.view.FFTunits = uicontrol('Style','text','HorizontalAlignment','left','String', 'Voltage','FontSize',EPmain.fontsize,...
        'ForegroundColor','black','Position',[10 77 50 20]);
end

if EPmain.view.edited.startSamp
    startSamp=EPmain.view.manual.startSamp;
    theTextColor.startSamp='blue';
else
    startSamp=EPmain.view.startSamp;
    theTextColor.startSamp='black';
end
if EPmain.view.edited.endSamp
    endSamp=EPmain.view.manual.endSamp;
    theTextColor.endSamp='blue';
else
    endSamp=EPmain.view.endSamp;
    theTextColor.endSamp='black';
end

if EPmain.view.edited.startHz
    startHz=EPmain.view.manual.startHz;
    theTextColor.startHz='blue';
else
    startHz=EPmain.view.startHz;
    theTextColor.startHz='black';
end
if EPmain.view.edited.endHz
    endHz=EPmain.view.manual.endHz;
    theTextColor.endHz='blue';
else
    endHz=EPmain.view.endHz;
    theTextColor.endHz='black';
end

EPmain.handles.view.topVolt = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%4.4f', theMax),'FontSize',EPmain.fontsize,...
    'Position',[10 58 40 20],'ForegroundColor',theTextColor.topVolt,'Callback',@checkViewSettings);
EPmain.handles.view.bottomVolt = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%4.4f', theMin),'FontSize',EPmain.fontsize,...
    'Position',[10 38 40 20],'ForegroundColor',theTextColor.bottomVolt,'Callback',@checkViewSettings);

if EPmain.view.flexMode
    theUnit='%';
else
    theUnit='Ms';
end
h = uicontrol('Style','text','HorizontalAlignment','left','String', theUnit,'FontSize',EPmain.fontsize,...
    'ForegroundColor','black','Position',[55 77 50 20]);
EPmain.handles.view.startSamp = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', round(startSamp)),'FontSize',EPmain.fontsize,...
    'Position',[55 58 35 20],'ForegroundColor',theTextColor.startSamp,'Callback',@checkViewSettings);
EPmain.handles.view.endSamp = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', round(endSamp)),'FontSize',EPmain.fontsize,...
    'Position',[55 38 35 20],'ForegroundColor',theTextColor.endSamp,'Callback',@checkViewSettings);
if strcmp('FFT',EPmain.view.dataTransform)
    set(EPmain.handles.view.startSamp ,'enable','off');
    set(EPmain.handles.view.endSamp ,'enable','off');
end

h = uicontrol('Style','text','HorizontalAlignment','left','String', 'Hz','FontSize',EPmain.fontsize,...
    'ForegroundColor','black','Position',[90 77 50 20]);
EPmain.handles.view.startHz = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', round(startHz)),'FontSize',EPmain.fontsize,...
    'Position',[90 58 35 20],'ForegroundColor',theTextColor.startHz,'Callback',@checkViewSettings);
EPmain.handles.view.endHz = uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', round(endHz)),'FontSize',EPmain.fontsize,...
    'Position',[90 38 35 20],'ForegroundColor',theTextColor.endHz,'Callback',@checkViewSettings);
if strcmp('VLT',EPmain.view.dataTransform)
    set(EPmain.handles.view.startHz ,'enable','off');
    set(EPmain.handles.view.endHz ,'enable','off');
end

h = uicontrol('Style','text','HorizontalAlignment','left','String', 'Mark','FontSize',EPmain.fontsize,...
    'ForegroundColor','black','Position',[125 77 30 20]);

EPmain.handles.view.marker1 = uicontrol('Style','edit','HorizontalAlignment','left','String', marker1,'FontSize',EPmain.fontsize,...
    'Position',[125 58 30 20],...
    'Callback',@checkViewSettings);

EPmain.handles.view.marker2 = uicontrol('Style','edit','HorizontalAlignment','left','String', marker2,'FontSize',EPmain.fontsize,...
    'Position',[125 38 30 20],...
    'Callback',@checkViewSettings);

%          h = uicontrol('Style','text','HorizontalAlignment','left','String', 'evt','FontSize',EPmain.fontsize,...
%             'ForegroundColor','black','Position',[160 57 25 20]);
%         EPmain.handles.view.evt = uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
%             'Value',EPmain.view.evt,'Position',[180 80 50 20],...
%             'Callback', ['global EPmain;','EPmain.view.evt=get(EPmain.handles.view.evt,''value'');','ep(''start'')'],...
%             'TooltipString','Include event lines in graphs.');

%         if ~strcmp('VLT',EPmain.view.dataTransform)
%             set(EPmain.handles.view.evt,'enable','off');
%         end

if EPmain.numColors<EPmain.maxColors
    EPmain.handles.view.expand = uicontrol('Style', 'pushbutton', 'String', '>','FontSize',EPmain.fontsize,...
        'Position', [175 250 20 20], 'Callback', ['global EPmain;','EPmain.numColors=EPmain.maxColors;','temp=get(EPmain.handles.hMainWindow,''Position'');','temp(3)=temp(3)*2;','set(EPmain.handles.hMainWindow,''Position'',temp);','ep(''start'');']);
else
    EPmain.handles.view.expand = uicontrol('Style', 'pushbutton', 'String', '<','FontSize',EPmain.fontsize,...
        'Position', [175 250 20 20], 'Callback', ['global EPmain;','EPmain.numColors=4;','temp=get(EPmain.handles.hMainWindow,''Position'');','temp(3)=temp(3)/2;','set(EPmain.handles.hMainWindow,''Position'',temp);','ep(''start'');']);
end

EPmain.handles.view.waves= uicontrol('Style', 'pushbutton', 'String', 'Waves','FontSize',EPmain.fontsize,...
    'Position', [2 0 50 35], 'Callback', 'ep(''viewWaves'')');

EPmain.handles.view.topos= uicontrol('Style', 'pushbutton', 'String', 'Topos','FontSize',EPmain.fontsize,...
    'Position', [52 0 50 35], 'Callback', 'ep(''viewTopos'')');

EPmain.handles.view.scan = uicontrol('Style', 'pushbutton', 'String', 'Scan','FontSize',EPmain.fontsize,...
    'Position', [102 0 50 35], 'Callback', 'ep(''viewEdit'')');

if strcmp('FFT',EPmain.view.dataTransform) || any(EPmain.view.allTrials(1:EPmain.numColors) == 2) || (EPmain.view.allTrials(1) > 0)
    set(EPmain.handles.view.scan ,'enable','off'); %can't scan with erpimages or FFT or if main is 'all'
end

if any(ismember(EPmain.view.allTrials(1:EPmain.numColors),[2 4 6]))
    set(EPmain.handles.view.topos ,'enable','off'); %can't view topos with erpimage options
end

for iColor=1:EPmain.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        if (strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'average') && (EPmain.view.trial(iColor) > 1))
            %can't view topos or scan with any band settings
            set(EPmain.handles.view.topos ,'enable','off');
            set(EPmain.handles.view.scan ,'enable','off');
        end
    end
end

EPmain.handles.view.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
    'Position', [152 0 50 35], 'Callback', ['global EPmain;','EPmain.mode=''main'';','temp=get(EPmain.handles.hMainWindow,''Position'');','temp(3)=EPmain.panesize(1);','set(EPmain.handles.hMainWindow,''Position'',temp);','ep(''start'');']);

if (theMin >= theMax) || (isnan(theMin) && isnan(theMax))
    set(EPmain.handles.view.waves,'enable','off');
    set(EPmain.handles.view.topos,'enable','off');
    disp('There is no data to view.');
end

if startSamp > endSamp
    set(EPmain.handles.view.waves,'enable','off');
    set(EPmain.handles.view.topos,'enable','off');
    disp('Start of the window needs to be before the end of the window.');
end

if (startHz > endHz) || (strcmp('FFT',EPmain.view.dataTransform) && (startHz == endHz))
    set(EPmain.handles.view.waves,'enable','off');
    set(EPmain.handles.view.topos,'enable','off');
    disp('Start of the window needs to be before the end of the window.');
end

if (startSamp == endSamp) && (startHz == endHz)
    set(EPmain.handles.view.waves,'enable','off');
    set(EPmain.handles.view.topos,'enable','off');
    disp('Start of the window needs to be before the end of the window.');
end

for iColor=1:EPmain.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        if ((endSamp-startSamp)==(1000/EPdataset.dataset(EPmain.view.dataset(iColor)).Fs)) && any(strcmp(EPmain.view.dataTransform,{'VLT','TFT'}))
            set(EPmain.handles.view.waves,'enable','off');
            set(EPmain.handles.view.scan,'enable','off');
            disp('Needs at least two time points.');
        end
    end
end

if ~isempty(marker1)
    if any(strcmp(EPmain.view.dataTransform,{'VLT','TFT'}))
        if (str2double(marker1) < startSamp) || (str2double(marker1) > endSamp)
            set(EPmain.handles.view.waves,'enable','off');
            set(EPmain.handles.view.topos,'enable','off');
            disp('Marker 1 needs to be within the window.');
        end
    else
        if (str2double(marker1) < startHz) || (str2double(marker1) > endHz)
            set(EPmain.handles.view.waves,'enable','off');
            set(EPmain.handles.view.topos,'enable','off');
            disp('Marker 1 needs to be within the window.');
        end
    end
end

if ~isempty(marker2)
    if any(strcmp(EPmain.view.dataTransform,{'VLT','TFT'}))
        if (str2double(marker2) < startSamp) || (str2double(marker2) > endSamp)
            set(EPmain.handles.view.waves,'enable','off');
            set(EPmain.handles.view.topos,'enable','off');
            disp('Marker 2 needs to be within the window.');
        end
    else
        if (str2double(marker2) < startHz) || (str2double(marker2) > endHz)
            set(EPmain.handles.view.waves,'enable','off');
            set(EPmain.handles.view.topos,'enable','off');
            disp('Marker 2 needs to be within the window.');
        end
    end
end

if EPmain.view.flexMode
    if (startSamp < 0) || (startSamp > 100)
        set(EPmain.handles.view.waves,'enable','off');
        set(EPmain.handles.view.topos,'enable','off');
        disp('Flexmode window boundaries must be within the range of zero to 100.');
    end
end

% if strcmp(EPmain.view.dataType,'continuous')
%     if (startSamp < 0) || (startSamp > 1000) || isnan(startSamp)
%         set(EPmain.handles.view.waves,'enable','off');
%         set(EPmain.handles.view.topos,'enable','off');
%         disp('Continuous data window boundaries must be within the range of zero to 1000.');
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function checkViewSettings(src,eventdata)
%check the settings in the view pane to see if the Waves control should be enabled.

global EPmain EPdataset

set(EPmain.handles.view.waves,'enable','on');

bottomVolt=str2double(get(EPmain.handles.view.bottomVolt,'string'));
topVolt=str2double(get(EPmain.handles.view.topVolt,'string'));
startSamp=str2double(get(EPmain.handles.view.startSamp,'string'));
endSamp=str2double(get(EPmain.handles.view.endSamp,'string'));
startHz=str2double(get(EPmain.handles.view.startHz,'string'));
endHz=str2double(get(EPmain.handles.view.endHz,'string'));
marker1=str2double(get(EPmain.handles.view.marker1,'string'));
marker2=str2double(get(EPmain.handles.view.marker2,'string'));

if isnan(bottomVolt)
    bottomVolt=[];
end
if isnan(topVolt)
    topVolt=[];
end
if isnan(startSamp)
    startSamp=[];
end
if isnan(endSamp)
    endSamp=[];
end
if isnan(startHz)
    startHz=[];
end
if isnan(endHz)
    endHz=[];
end
if isnan(marker1)
    marker1=[];
end
if isnan(marker2)
    marker2=[];
end

activeList=find(EPmain.view.dataset(1:EPmain.numColors) <= length(EPdataset.dataset));

%convert Hz values to non-dB amplitude measures
if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'})) && ~all(EPmain.view.correl(activeList)) && any(~strcmp('VLT',EPmain.view.type(activeList)) & ~EPmain.view.BSC(activeList))
    if (EPmain.view.FFTunits == 4)
        if bottomVolt == -flintmax
            bottomVolt=-inf; %log10 of zero is -inf.  Reverse replacing with maximum possible double-precision negative number.
        end
        bottomVolt=10^(bottomVolt/10);
        if topVolt == -flintmax
            topVolt=-inf; %log10 of zero is -inf.  Reverse replacing with maximum possible double-precision negative number.
        end
        topVolt=10^(topVolt/10);
    end
    if (EPmain.view.FFTunits > 2)
        bottomVolt=sqrt(bottomVolt); %convert power to amplitude
        topVolt=sqrt(topVolt); %convert power to amplitude
    end
    bottomVolt=bottomVolt*sqrt(EPmain.view.binHz); %convert from spectral density
    topVolt=topVolt*sqrt(EPmain.view.binHz); %convert from spectral density
end

if isempty(startHz) || isempty(endHz) || isempty(startSamp) || isempty(endSamp) || isempty(bottomVolt) || isempty(topVolt)
    EPmain.view.changeFlag(1:4)=1;
end

if src==EPmain.handles.view.bottomVolt
    if isempty(bottomVolt) || isnan(bottomVolt)
        EPmain.view.edited.bottomVolt=0;
    else
        EPmain.view.edited.bottomVolt=1;
        EPmain.view.manual.bottomVolt=bottomVolt;
    end
end

if src==EPmain.handles.view.topVolt
    if isempty(topVolt) || isnan(topVolt)
        EPmain.view.edited.topVolt=0;
    else
        EPmain.view.edited.topVolt=1;
        EPmain.view.manual.topVolt=topVolt;
    end
end

if src==EPmain.handles.view.startSamp
    if isempty(startSamp) || isnan(startSamp)
        EPmain.view.edited.startSamp=0;
    else
        EPmain.view.edited.startSamp=1;
        EPmain.view.manual.startSamp=startSamp;
    end
end

if src==EPmain.handles.view.endSamp
    if isempty(endSamp) || isnan(endSamp)
        EPmain.view.edited.endSamp=0;
    else
        EPmain.view.edited.endSamp=1;
        EPmain.view.manual.endSamp=endSamp;
    end
end

if src==EPmain.handles.view.startHz
    if isempty(startHz) || isnan(startHz)
        EPmain.view.edited.startHz=0;
    else
        EPmain.view.edited.startHz=1;
        EPmain.view.manual.startHz=startHz;
    end
end

if src==EPmain.handles.view.endHz
    if isempty(endHz) || isnan(endHz)
        EPmain.view.edited.endHz=0;
    else
        EPmain.view.edited.endHz=1;
        EPmain.view.manual.endHz=endHz;
    end
end

EPmain.view.startSamp=startSamp;
EPmain.view.endSamp=endSamp;
EPmain.view.startHz=startHz;
EPmain.view.endHz=endHz;
EPmain.view.marker1=marker1;
EPmain.view.marker2=marker2;

ep('startView');