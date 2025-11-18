function [EPdataOut]=ep_selectData(EPdataIn,keep)
%  [EPdataOut]=ep_selectData(EPdataIn,keep)
%       Selects out specified members from each dimension of the data.  When selecting out trials from a continuous
%       file, selects out one second segments.  Will still be continuous data, so normally just specify a single epoch. 
%       If continuous, one can also specify both points and cells.  In this case,
%       the cells will represent the one second epochs and the points will indicate how many points is in this
%       one-second epoch.  If points is not specified for a continuous file, then EPdata.Fs will be used to determine how many points in a second.
%       If continuous and no cells are specified, just points, then will just trim the data.
%       If the trimmed starting samples are an even multiple of one second, will just drop those epochs.
%       Otherwise, the .analysis fields will be reinitialized starting with the new starting sample.
%
%Inputs:
%  EPdata         : Structured array with the data and accompanying information in EP file format.  See readData.
%  keep           : Which entries to keep from each dimension (cell array with six cells, one per dimension).  Empty
%                   means keep everything.
%
%Outputs:
%  EPdata        : Structured array with the data and accompanying information in EP file format.  See readData.

%History:
%  by Joseph Dien (5/25/11)
%  jdien07@mac.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% modified 1/28/12 JD
% Eliminated REF channel type and added reference field to keep track of original and current reference scheme.
% 
% modified 2/22/12 JD
% Noise and std fields no longer optional.  Now set to empty if not used.
%
% modified 5/23/12 JD
% Added support for 6th dimension of frequency for time-frequency analyses
%
% modified 6/14/12 JD
% Added support for selecting epochs out of continuous data.
%
% bugfix 9/6/12 JD
% Fixed selection of points from a continuous file.
%
% bugfix 9/15/12 JD
% Fixed crash when deleting factors from data with a combined factor add.  Fixed crash when deleting factors.
%
% bugfix 2/6/13 JD
% Fixed crashes after Edit function used to trim the range of frequencies due to not applying trimming to std or FacVecF
% fields.
%
% modified 10/9/13 JD
% Added recTime field.
%
% bugfix 11/14/13 JD
% .analysis edit fields updated when points dropped from a continuous file.  Boundary events added as needed.
%
% bufix 3/12/14 JD
% Handles decimal sampling rates gracefully.
%
% modified 3/24/14 JD
% Added .cov field.
%
% bufix 3/27/14 JD
% Fixed crash when selecting time points and there are empty event cells.
%
% modified 4/9/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure.
%
% modified 4/14/14 JD
% Added .covNum field.
%
% bufix 6/12/14 JD
% Fixed blank keys field of events being produced without .key (e.g., .keys.keyCod instead of .keys.key.keyCode)
%
% bufix 6/29/14 JD
% Fixed subject selection not being applied to .cov.Nq field, resulting in crashes down the line for combined subject average files.
%
% modified 7/16/14 JD
% Simplified keys code field structure.
%
% bufix 7/28/14 JD
% Fixed crash when selecting points from a continuous file such that the points align with the one-second epochs in the
% .analysis fields.
% Fixed only first event sample being updated when points selected from continuous data, as in trimming data.
% Fixed boundary events falling on edge of selected time range not being deleted.
%
% bugfix 10/27/14 JD
% Fixed crash when selecting points in continuous data where the number of
% points is not an even multiple of one second.
%
% bugfix 3/20/15 JD
% Fixed crash when channels selected in the Chans tab and then the working
% copy was saved.
%
% bugfix 7/4/15 JD
% Fixed crash when selecting less than a second of data.
%
% modified 9/4/15 JD
% Added trial specs for average files.
%
% modified 10/16/16 JD
% Added .stims field.
%
% modified 5/4/17 JD
% Optimized .stims field processing.
%
% modified 11/22/17 JD
% Added support for impedances field.
%
% modified 2/11/18 JD
% Added support for stdCM field.
%
% bufix & modified 1/9/19 JD
% Now allows CMB and SGL factors to be intermixed.
% FacVar and FacVarQ now include CMB factors.
%
% modified 4/9/19 JD
% Added support for task level performance measures.
%
% bugfix 5/3/19 JD
% Fixed crash when selecting subjects from a dataset with no task specs.
%
% modified 9/10/19 JD
% Added AOI subfield to stims field.
%
% modified 11/4/19 JD
% Added sessNums sessNames fields.
%
% modified 12/24/19 JD
% Upgraded support of std information by adding .covAVE and .GAVsubs fields and eliminating .std and .stdCM fields.
%
% bugfix 3/10/20 JD
% Fixed adding redundant copies to the stims entries.
%
% modified 1/27/25 JD
% Added video field.
%
% modified 4/5/25 JD
% Added support for virtual grand averages.
% Fixed crash when data has relationship channels.
%
% bugfix 6/14/25 JD
% Fixed crash when dropping constituents of retained virtual cells or subjects, by converting to real data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 1999-2016  Joseph Dien
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

global EPtictoc

EPdataOut=EPdataIn;

if length(keep) ~= 6
    disp('Cell array for data selection needs to be six long.');
    EPdataOut=[];
    return
end

chans=keep{1};
points=keep{2};
cells=keep{3};
subs=keep{4};
facs=keep{5};
if isempty(EPdataIn.facNames)
    facs=1;
end
freqs=keep{6};
if isempty(EPdataIn.relNames)
    relChans=1;
else
    relChans=[1:length(EPdataIn.relNames)];
end

numSubs=length(EPdataIn.subNames);
numVsubs=max(0,size(EPdataIn.GAVsubs,1)-1);
numRsubs=numSubs-numVsubs;
numCells=length(EPdataIn.cellNames);
numVcells=max(0,size(EPdataIn.GAVsubs,2)-1);
numRcells=numCells-numVcells;

%if the constituents of any retained virtual cells or subjects are dropped, then will need
%to convert the virtual data to real.
if (~isempty(cells) && any(cells>numRcells))
    cellConvert=zeros(numRcells,1);
    cellCounter=0;
    for iCell=1:numRcells
        if any(cells==iCell)
            cellCounter=cellCounter+1;
            cellConvert(iCell)=cellCounter;
        else
            cellConvert(iCell)=NaN;
        end
    end
    for iCell=1:length(cells)
        theCell=cells(iCell);
        if theCell>numRcells
            theGAVcell=theCell-numRcells;
            if ~any(ismember(EPdataOut.GAVsubs{1,theGAVcell+1,1}(:,1),cells(cells<=numRcells)))
                %if any of the constituents of a virtual GAVE are deleted, then convert to a normal GAVE
                EPdataOut=ep_combineData(EPdataOut,'cells',{[],[],[],[],[],[]},[],EPdataOut.cellNames{theCell},[]);
                numVcells=numVcells-1;
                numRcells=numRcells+1;
            else
                %convert subject indices to the new ones accounting for deletions
                for iFac=1:size(EPdataOut.GAVsubs,3)
                    tempGAVsubs=[];
                    for iRow=1:size(EPdataOut.GAVsubs{1,theGAVcell+1,iFac},1)
                        if ~isnan(cellConvert(EPdataOut.GAVsubs{1,theGAVcell+1,iFac}(iRow,1)))
                            tempGAVsubs(end+1,:)=[cellConvert(EPdataOut.GAVsubs{1,theGAVcell+1,iFac}(iRow,1)) EPdataOut.GAVsubs{1,theGAVcell+1,iFac}(iRow,2)];
                        end
                    end
                    EPdataOut.GAVsubs{1,theGAVcell+1,iFac}=tempGAVsubs;
                end
            end
        end
    end
end
if (~isempty(subs) && any(subs>numRsubs))
    subConvert=zeros(numRsubs,1);
    subCounter=0;
    for iSub=1:numRsubs
        if any(subs==iSub)
            subCounter=subCounter+1;
            subConvert(iSub)=subCounter;
        else
            subConvert(iSub)=NaN;
        end
    end
    for iSub=1:length(subs)
        theSub=subs(iSub);
        if theSub>numRsubs
            theGAV=theSub-numRsubs;
            if ~any(ismember(EPdataOut.GAVsubs{theGAV+1,1,1}(:,1),subs(subs<=numRsubs)))
                %if any of the constituents of a virtual GAVE are deleted, then convert to a normal GAVE
                EPdataOut=ep_combineData(EPdataOut,'subjects',{[],[],[],[],[],[]},[],EPdataOut.subNames{theSub},[]);
                numVsubs=numVsubs-1;
                numRsubs=numRsubs+1;
            else
                %convert subject indices to the new ones accounting for deletions
                for iCell=1:size(EPdataOut.GAVsubs,2) %including the first column
                    for iFac=1:size(EPdataOut.GAVsubs,3)
                        tempGAVsubs=[];
                        for iRow=1:size(EPdataOut.GAVsubs{theGAV+1,iCell,iFac},1)
                            if ~isnan(subConvert(EPdataOut.GAVsubs{theGAV+1,iCell,iFac}(iRow,1)))
                                tempGAVsubs(end+1,:)=[subConvert(EPdataOut.GAVsubs{theGAV+1,iCell,iFac}(iRow,1)) EPdataOut.GAVsubs{theGAV+1,iCell,iFac}(iRow,2)];
                            end
                        end
                        EPdataOut.GAVsubs{theGAV+1,iCell,iFac}=tempGAVsubs;
                    end
                end
            end
        end
    end
end

if ~isempty(chans)
    if ~isempty(EPdataOut.facVecS)
        EPdataOut.facVecS=EPdataOut.facVecS(chans,:);
    else
        EPdataOut.data=EPdataOut.data(chans,:,:,:,:,:,relChans);
    end
    if ~isempty(EPdataOut.facData)
        EPdataOut.facData=EPdataOut.facData(chans,:,:,:,:,:,relChans);
    end
    if ~isempty(EPdataOut.noise)
        EPdataOut.noise=EPdataOut.noise(chans,:,:,:,:);
    end
    if ~isempty(EPdataOut.covAVE)
        if size(EPdataOut.covAVE,7)==1
            EPdataOut.covAVE=EPdataOut.covAVE(chans,:,:,:,:,:,1);
        else
            EPdataOut.covAVE=EPdataOut.covAVE(chans,:,:,:,:,:,chans);
        end
    end
    if ~isempty(EPdataOut.cov)
        EPdataOut.cov.covMatrix=EPdataOut.cov.covMatrix(:,chans,chans);
    end
    EPdataOut.chanNames=EPdataOut.chanNames(chans);
    EPdataOut.chanTypes=EPdataOut.chanTypes(chans);
    if ~isempty(EPdataOut.analysis.badChans)
        EPdataOut.analysis.badChans=EPdataOut.analysis.badChans(:,:,chans);
    end
    if ~isempty(EPdataOut.eloc)
        EPdataOut.eloc=EPdataOut.eloc(chans);
    end
    
    reference.type=EPdataIn.reference.type;
    reference.original=[];
    reference.current=[];
    for i=1:length(EPdataIn.reference.original)
        if any(EPdataIn.reference.original(i) == chans)
            reference.original(end+1)=sum(chans <= EPdataIn.reference.original(i));
        end
    end
    for i=1:length(EPdataIn.reference.current)
        if any(EPdataIn.reference.current(i) == chans)
            reference.current(end+1)=sum(chans <= EPdataIn.reference.current(i));
        end
    end
    EPdataOut.reference=reference;
    if ~isempty(EPdataIn.impedances.channels)
        EPdataOut.impedances.channels=EPdataOut.impedances.channels(chans,:);
    end
end

%if points are specified and not the case of continuous data being chopped up into one-second epochs.
if ~isempty(points) && ~(strcmp(EPdataOut.dataType,'continuous') && ~isempty(cells))
    if ~ismember(1,points)
        EPdataOut.recTime=EPdataOut.recTime+min(points)-1;
    end
    if ~isempty(EPdataOut.facVecT)
        EPdataOut.facVecT=EPdataOut.facVecT(points,:);
    else
        EPdataOut.data=EPdataOut.data(:,points,:,:,:,:,:);
    end
    if ~isempty(EPdataOut.facData)
        EPdataOut.facData=EPdataOut.facData(:,points,:,:,:,:,:);
    end
    if ~isempty(EPdataOut.noise)
        EPdataOut.noise=EPdataOut.noise(:,points,:,:,:,:);
    end
    if ~isempty(EPdataOut.covAVE)
        EPdataOut.covAVE=EPdataOut.covAVE(:,points,:,:,:,:,:);
    end
    EPdataOut.timeNames=EPdataOut.timeNames(points);
    if ~isempty(EPdataOut.video)
        for iTrial=1:length(EPdataOut.video)
            numLess=length(find(EPdataOut.video.times<EPdataOut.timeNames(points(1))));
            if numLess > 1
                EPdataOut.video.frames(1:numLess-1)=[];
                EPdataOut.video.times(1:numLess-1)=[];
            end
            numMore=length(find(EPdataOut.video.times>EPdataOut.timeNames(points(end))));
            if numMore > 1
                EPdataOut.video.frames(end-numMore+2:end)=[];
                EPdataOut.video.times(end-numMore+2:end)=[];
            end
        end
    end
    %delete events
    oldEvents=EPdataOut.events;
    numPoints=length(EPdataIn.timeNames);
    dropPoints=find(~ismember([1:numPoints],points));
    for subject = 1:length(EPdataOut.subNames)
        for wave = 1:length(EPdataOut.cellNames)
            if ~isempty(EPdataOut.events{subject,wave})
                numEvents=length(EPdataOut.events{subject,wave});
                keepEvents=true(numEvents,1);
                for iEvent=1:numEvents
                    eventSample=EPdataOut.events{subject,wave}(iEvent).sample;
                    if any(ismember(dropPoints,eventSample)) || (any(ismember(dropPoints,floor(eventSample))) && any(ismember(dropPoints,ceil(eventSample)))) || (any(ismember(dropPoints,1)) && eventSample<1) || (any(ismember(dropPoints,numPoints)) && eventSample>numPoints) 
                        keepEvents(iEvent)=false;
                    end
                end
                EPdataOut.events{subject,wave}=EPdataOut.events{subject,wave}(keepEvents);
            end
        end
    end    
    if min(points) > 1
        EPdataOut.baseline = EPdataOut.baseline - min(points) +1;
        for subject = 1:length(EPdataOut.subNames)
            for wave = 1:length(EPdataOut.cellNames)
                deleteEvents=[];
                for event = 1:length(EPdataOut.events{subject,wave})
                    if ~isempty(EPdataOut.events{subject,wave})
                        EPdataOut.events{subject,wave}(event).sample=EPdataOut.events{subject,wave}(event).sample-min(points)+1;
                        if strcmp('boundary',EPdataOut.events{subject,wave}(event).type)
                            if (EPdataOut.events{subject,wave}(event).sample==1) || (EPdataOut.events{subject,wave}(event).sample==length(EPdataOut.timeNames))
                                deleteEvents(end+1)=event; %delete boundary if it falls on outer edges of selected range
                            end
                        end
                    end
                end
                EPdataOut.events{subject,wave}(deleteEvents)=[];
            end
        end
    end
    %insert boundary event to mark the deleted gap
    midGaps=find((diff(points)~=1));
    if ~isempty(midGaps)
        theGaps=diff(points);
        for iGap=1:length(midGaps)
            for subject = 1:length(EPdataOut.subNames)
                for wave = 1:length(EPdataOut.cellNames)
                    deletedEvents=find(ismember([oldEvents{subject,wave}.sample],[points(midGaps(iGap))+1:points(midGaps(iGap)+1)]));
                    boundaryList=find(strcmp('boundary',[oldEvents{subject,wave}(deletedEvents).type]));
                    gapTime=0;
                    if ~isempty(boundaryList)
                        for iBoundary=1:length(boundaryList)
                            gapTime=gapTime+oldEvents{subject,wave}(deletedEvents(boundaryList(iBoundary))).sample; %add existing boundary event durations to total deletion time
                        end
                        %need to delete boundary just after new deletion zone
                        if oldEvents{subject,wave}(deletedEvents(boundaryList(end))).sample == points(midGaps(iGap)+1)
                            EPdataOut.events{subject,wave}(find([EPdataOut.events{subject,wave}(find(strcmp('boundary',[EPdataOut.events{subject,wave}.type]))).sample] == points(midGaps(iGap)+1)))=[];
                        end
                    end
                    EPdataOut.events{subject,wave}(end+1).type='boundary';
                    EPdataOut.events{subject,wave}(end).value='boundary';
                    EPdataOut.events{subject,wave}(end).sample=midGaps(iGap)+1;
                    EPdataOut.events{subject,wave}(end).duration=theGaps(midGaps(iGap))-1+gapTime;
                    EPdataOut.events{subject,wave}(end).keys=struct('code','','data','','datatype','','description','');
                end
            end
        end
    end
    if strcmp(EPdataOut.dataType,'continuous')
        epochKeepPoints=floor((points-1)/EPdataIn.Fs)+1;
        epochKeep=unique(epochKeepPoints);
        epochKeepTot=histc(epochKeepPoints,epochKeep);
        epochAllPoints=floor(([0:length(EPdataIn.timeNames)-1])/EPdataIn.Fs)+1;
        epochAll=unique(epochAllPoints);
        epochAllTot=histc(epochAllPoints,unique(epochAllPoints));
        if epochAllTot(end)<EPdataIn.Fs
            epochAll=setdiff(epochKeep,epochAll(end)); %left over samples are tacked onto the final epoch
        end
        
        if ~isempty(epochKeepTot(1:end-1)) && all(epochKeepTot(1:end-1)==EPdataIn.Fs)
            %if points dropped from start are round seconds, then keep edit codes and just drop the epochs
            if epochKeepTot(end) < EPdataIn.Fs
                epochKeep=epochKeep(1:end-1);
            end
            
            EPdataOut.analysis.blinkTrial=EPdataOut.analysis.blinkTrial(1,epochKeep);
            EPdataOut.analysis.saccadeTrial=EPdataOut.analysis.saccadeTrial(1,epochKeep);
            EPdataOut.analysis.saccadeOnset=EPdataOut.analysis.saccadeOnset(1,epochKeep);
            EPdataOut.analysis.moveTrial=EPdataOut.analysis.moveTrial(1,epochKeep);
            EPdataOut.analysis.badTrials=EPdataOut.analysis.badTrials(1,epochKeep);
            EPdataOut.analysis.badChans=EPdataOut.analysis.badChans(1,epochKeep,:);
            
            if (epochAllTot(end) < EPdataIn.Fs) && (epochAllTot(end) ~= epochKeepTot(end)) && (points(end) ~= length(EPdataIn.timeNames))
                % if there were left over samples tacked on but some were trimmed off, then reset analysis codes to zero
                EPdataOut.analysis.blinkTrial(end)=0;
                EPdataOut.analysis.saccadeTrial(end)=0;
                EPdataOut.analysis.saccadeOnset(end)=0;
                EPdataOut.analysis.moveTrial(end)=0;
                EPdataOut.analysis.badTrials(end)=0;
                EPdataOut.analysis.badChans(1,end,:)=0;
            end
        else
            %else reinitialize the edit codes
            numEpochs=floor(size(EPdataOut.data,2)/ceil(EPdataOut.Fs)); %excess time points are tacked onto final epoch
            if numEpochs == 0
                numEpochs =1;
            end
            EPdataOut.analysis.blinkTrial=zeros(1,numEpochs);
            EPdataOut.analysis.saccadeTrial=zeros(1,numEpochs);
            EPdataOut.analysis.saccadeOnset=zeros(1,numEpochs);
            EPdataOut.analysis.moveTrial=zeros(1,numEpochs);
            EPdataOut.analysis.badTrials=zeros(1,numEpochs);
            EPdataOut.analysis.badChans=zeros(1,numEpochs,size(EPdataOut.data,1));
        end
    end
end

if ~isempty(cells)
    if ~isempty(EPdataOut.GAVsubs)
        EPdataOut.GAVsubs=EPdataOut.GAVsubs(:,[1 cells(cells>numRcells)+1-numRcells],facs);
        rCells=cells(cells<=numRcells);
    else
        rCells=cells;
    end
    if strcmp(EPdataOut.dataType,'continuous') %if continuous and cells have been specified, chop up data into one-second epochs
        if ~isempty(points)
            epochPoints=points;
        else
            epochPoints=[1:ceil(EPdataOut.Fs)];
        end
        thePoints=kron((rCells-1)*length(epochPoints),ones(1,length(epochPoints)))+kron(ones(1,length(rCells)),epochPoints);
        EPdataOut.timeNames=EPdataOut.timeNames(thePoints);
        EPdataOut.data=EPdataOut.data(:,thePoints,1,:,:,:,:);
        if ~isempty(EPdataOut.facData)
            EPdataOut.facData=EPdataOut.facData(:,thePoints,1,:,:,:,:);
        end
        if ~isempty(EPdataOut.noise)
            EPdataOut.noise=EPdataOut.noise(:,thePoints,1,:,:,:);
        end
        if ~isempty(EPdataOut.covAVE)
            EPdataOut.covAVE=EPdataOut.covAVE(:,thePoints,1,:,:,:,:);
        end
    else
        EPdataOut.data=EPdataOut.data(:,:,rCells,:,:,:,:);
        if ~isempty(EPdataOut.facData)
            EPdataOut.facData=EPdataOut.facData(:,:,rCells,:,:,:,:);
        end
        if ~isempty(EPdataOut.noise)
            EPdataOut.noise=EPdataOut.noise(:,:,rCells,:,:,:);
        end
        if ~isempty(EPdataOut.covAVE)
            EPdataOut.covAVE=EPdataOut.covAVE(:,:,rCells,:,:,:,:);
        end
        EPdataOut.cellNames=EPdataOut.cellNames(cells);
        EPdataOut.cellTypes=EPdataOut.cellTypes(cells);
        EPdataOut.events=EPdataOut.events(:,rCells);
        EPdataOut.avgNum=EPdataOut.avgNum(:,rCells);
        EPdataOut.covNum=EPdataOut.covNum(:,rCells);
        EPdataOut.subNum=EPdataOut.subNum(:,rCells);
        EPdataOut.recTime=EPdataOut.recTime(rCells);
    end
    EPdataOut.analysis.badChans=EPdataOut.analysis.badChans(:,rCells,:);
    EPdataOut.analysis.moveTrial=EPdataOut.analysis.moveTrial(:,rCells);
    EPdataOut.analysis.blinkTrial=EPdataOut.analysis.blinkTrial(:,rCells);
    EPdataOut.analysis.saccadeTrial=EPdataOut.analysis.saccadeTrial(:,rCells);
    EPdataOut.analysis.saccadeOnset=EPdataOut.analysis.saccadeOnset(:,rCells);
    EPdataOut.analysis.badTrials=EPdataOut.analysis.badTrials(:,rCells);
    if ~isempty(EPdataOut.trialSpecs)
        EPdataOut.trialSpecs=EPdataOut.trialSpecs(rCells,:,:);
    end
    if strcmp(EPdataOut.dataType,'single_trial')
        if ~isempty(EPdataOut.trialNames)
            EPdataOut.trialNames=EPdataOut.trialNames(rCells);
        end
    end
    if ~isempty(EPdataOut.video)
        for iTrial=1:length(EPdataOut.video)
            numLess=length(find(EPdataOut.video.times<EPdataOut.timeNames(points(1))));
            if numLess > 1
                EPdataOut.video.frames(1:numLess-1)=[];
                EPdataOut.video.times(1:numLess-1)=[];
            end
            numMore=length(find(EPdataOut.video.times>EPdataOut.timeNames(points(end))));
            if numMore > 1
                EPdataOut.video.frames(end-numMore+2:end)=[];
                EPdataOut.video.times(end-numMore+2:end)=[];
            end
        end
    end
end

if ~isempty(subs)
    if ~isempty(EPdataOut.GAVsubs)
        EPdataOut.GAVsubs=EPdataOut.GAVsubs([1 subs(subs>numRsubs)+1-numRsubs],:,:);
        rSubs=subs(subs<=numRsubs);
    else
        rSubs=subs;
    end
    EPdataOut.data=EPdataOut.data(:,:,:,rSubs,:,:,:);
    if ~isempty(EPdataOut.facData)
        EPdataOut.facData=EPdataOut.facData(:,:,:,rSubs,:,:,:);
    end
    if ~isempty(EPdataOut.noise)
        EPdataOut.noise=EPdataOut.noise(:,:,:,rSubs,:,:);
    end
    if ~isempty(EPdataOut.covAVE)
        EPdataOut.covAVE=EPdataOut.covAVE(:,:,:,rSubs,:,:,:);
    end
    if ~isempty(EPdataOut.cov)
        EPdataOut.cov.covMatrix=EPdataOut.cov.covMatrix(rSubs,:,:);
        EPdataOut.cov.Nq=EPdataOut.cov.Nq(rSubs);
    end

    EPdataOut.subNames=EPdataOut.subNames(subs);
    EPdataOut.subTypes=EPdataOut.subTypes(subs);
    if ~isempty(EPdataOut.sessNums)
        EPdataOut.sessNums=EPdataOut.sessNums(rSubs);
    end
    if ~isempty(EPdataOut.subjectSpecs)
        EPdataOut.subjectSpecs=EPdataOut.subjectSpecs(rSubs,:);
    end
    if ~isempty(EPdataOut.taskSpecs)
        EPdataOut.taskSpecs=EPdataOut.taskSpecs(rSubs,:,:);
    end
    EPdataOut.avgNum=EPdataOut.avgNum(rSubs,:);
    EPdataOut.covNum=EPdataOut.covNum(rSubs,:);
    EPdataOut.subNum=EPdataOut.subNum(rSubs,:);
    EPdataOut.events=EPdataOut.events(rSubs,:);
    EPdataOut.analysis.badChans=EPdataOut.analysis.badChans(rSubs,:,:);
    EPdataOut.analysis.moveTrial=EPdataOut.analysis.moveTrial(rSubs,:);
    EPdataOut.analysis.blinkTrial=EPdataOut.analysis.blinkTrial(rSubs,:);
    EPdataOut.analysis.saccadeTrial=EPdataOut.analysis.saccadeTrial(rSubs,:);
    EPdataOut.analysis.saccadeOnset=EPdataOut.analysis.saccadeOnset(rSubs,:);
    EPdataOut.analysis.badTrials=EPdataOut.analysis.badTrials(rSubs,:);
    if ~isempty(EPdataOut.trialSpecs)
        EPdataOut.trialSpecs=EPdataOut.trialSpecs(:,:,rSubs);
    end
    if ~isempty(EPdataIn.impedances.channels)
        EPdataOut.impedances.channels=EPdataOut.impedances.channels(:,rSubs);
    end
    if ~isempty(EPdataIn.impedances.ground)
        EPdataOut.impedances.ground=EPdataOut.impedances.ground(rSubs);
    end
end

if ~isempty(facs) && ~isempty(EPdataOut.facNames)
    SGLfacs=find(strcmp('SGL',EPdataOut.facTypes));
    CMBfacs=find(strcmp('CMB',EPdataOut.facTypes));
    SGLkeepFacs=[];
    CMBkeepFacs=[];
    for iFactor=1:length(facs)
        if strcmp('SGL',EPdataOut.facTypes{facs(iFactor)})
            SGLkeepFacs(end+1)=find(SGLfacs==facs(iFactor));
        elseif strcmp('CMB',EPdataOut.facTypes{facs(iFactor)})
            CMBkeepFacs(end+1)=find(CMBfacs==facs(iFactor));
        end
    end

    if ~isempty(EPdataOut.facVecS)
        EPdataOut.facVecS=EPdataOut.facVecS(:,SGLkeepFacs);
    end
    if ~isempty(EPdataOut.facVecT)
        EPdataOut.facVecT=EPdataOut.facVecT(:,SGLkeepFacs);
    end
    if ~isempty(EPdataOut.facVar)
        EPdataOut.facVar=EPdataOut.facVar(:,facs);
    end
    if ~isempty(EPdataOut.facVarQ)
        EPdataOut.facVarQ=EPdataOut.facVarQ(:,facs);
    end
    if ~isempty(EPdataOut.facData)
        EPdataOut.facData=EPdataOut.facData(:,:,:,:,CMBkeepFacs,:,:);
    end
    EPdataOut.data=EPdataOut.data(:,:,:,:,SGLkeepFacs,:,:);
    if ~isempty(EPdataOut.facNames)
        EPdataOut.facNames=EPdataOut.facNames(facs);
    end
    if ~isempty(EPdataOut.facTypes)
        EPdataOut.facTypes=EPdataOut.facTypes(facs);
    end
    if ~isempty(EPdataOut.covAVE)
        EPdataOut.covAVE=EPdataOut.covAVE(:,:,:,:,facs,:,:);
    end
end

if ~isempty(freqs)
    if ~isempty(EPdataOut.facData)
        EPdataOut.facData=EPdataOut.facData(:,:,:,:,:,freqs,:);
    end
    if ~isempty(EPdataOut.facVecF)
        EPdataOut.facVecT=EPdataOut.facVecF(freqs,:);
    else
        EPdataOut.data=EPdataOut.data(:,:,:,:,:,freqs,:);
    end
    if ~isempty(EPdataOut.covAVE)
        EPdataOut.covAVE=EPdataOut.covAVE(:,:,:,:,:,freqs,:);
    end
    EPdataOut.freqNames=EPdataOut.freqNames(freqs);
end

if ~isempty(cells) || ~isempty(subs)
    EPdataOut.stims=struct('name',{},'image',{},'AOI',{});
    for iStim=1:length(EPdataIn.stims)
        for iSub=1:numSubs
            for iCell=1:numRcells
                for iEvent=1:length(EPdataOut.events{iSub,iCell})
                    if any(strcmp(EPdataIn.stims(iStim).name,{EPdataOut.events{iSub,iCell}(iEvent).keys.data})) && ~any(strcmp(EPdataIn.stims(iStim).name,{EPdataOut.stims.name}))
                        EPdataOut.stims(end+1)=EPdataIn.stims(iStim); %keep only stim images whose events are still in the selected data.
                    end
                end
            end
        end
    end
end

[err]=ep_checkEPfile(EPdataOut);
if err
    disp('EP data format defective after data selection performed.');
    EPdataOut=[];
    return
end


