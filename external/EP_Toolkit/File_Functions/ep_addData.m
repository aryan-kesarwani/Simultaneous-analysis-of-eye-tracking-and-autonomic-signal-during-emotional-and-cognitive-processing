function [EPdataOut]=ep_addData(EPdataIn,EPadd,dataDimension,checkFlag)
%  [EPdataOut]=ep_addData(EPdataIn,EPadd,dataDimension,checkFlag);
%       Adds a level of data to one of the data dimensions, using a partial EPdata variable as the source.
%       Defaults added where not present in the partial EPdata variable.
%       The corresponding names field is required.
%
%Inputs:
%  EPdataIn       : Structured array with the input data and accompanying information in EP file format.  See readData.
%  EPadd          : Structured array with the information to be added.
%                   When adding GAVsubs, need to include the first row and the first column.
%                   Will need to add it twice, once for 'subjects' and once for 'cells'.
%                   Include the full GAVsubs matrix each time.
%  dataDimension  : Which data dimension to add to (i.e., 'channels', 'cells', 'subjects', 'factors', 'points')
%                   'cells' includes trials if single-trial data.
%  checkFlag      : EPcheck the final resulting data structure (default=1).
%
%Outputs:
%  EPdataOut      : Structured array with the output data and accompanying information in EP file format.  See readData.

%History:
%  by Joseph Dien (6/4/14)
%  jdien07@mac.com
%
% modified 10/2/14 JD
% Can add datasets with different trial specs.
%
% bugfix 12/2/14 JD
% If adding two sets of data with cells with the same trialNames, then make
% sure they no longer overlap.
%
% bugfix 5/19/15 JD
% Fixed crash when adding non single-trial datasets.
%
% modified 9/4/15 JD
% Added trial specs for average files.
%
% modified 10/12/15 JD
% Fixed crash when adding an average file with trial specs.
%
% bugfix 10/26/15 JD
% Fixed crash when adding to a data file with an empty trialspecs field.
%
% modified 10/16/16 JD
% Added .stims field.
%
% bugfix 6/1/17 JD
% Fixed trialnames not being correctly numbered when adding cells to a single-trial file.
%
% modified 12/22/17 JD
% Added support for impedances field.
%
% bugfix 12/6/17 JD
% Fixed crash when input file has impedances field and channels are being added but none are being added to the impedances field.
%
% modified & bugfix 2/11/18 JD
% Added support for stdCM field.
% No longer tries to combine std values.  Instead sets to zero.
%
% bugfix 10/18/18 JD
% Fixed crash when there are multiple specs with the same name.
%
% bugfix 12/17/18 JD
% Fixed crash when adding cells or trials with trial specs to an existing file that already has trial specs, as in averaging.
%
% modified 4/9/19 JD
% Added support for task level performance measures.
%
% bugfix 9/5/19 JD
% Fixed crash when adding cells or subjects where the .cov field is present but empty.
%
% modified 10/18/19 JD
% Added support for adding points.
%
% modified 11/4/19 JD
% Added sessNums sessNames fields.
%
% modified 1/13/20 JD
% Upgraded support of std information by adding .covAVE and .GAVsubs fields and eliminating .std and .stdCM fields.
%
% bugfix 4/8/20 JD
% Fixed crash when adding subjects.
%
% bugfix 5/28/20 JD
% Fixed error when adding channels due to bug in covMatrix code.
%
% bugfix 7/11/20 JD
% Fixed error when adding two continuous segments together and the first has a number of points that does not divided evenly into one-second epochs.
% Fixed crash when adding continuous segment with events to one without events.
%
% bugfix 2/4/21 JD
% Fixed not adding covAVE field correctly when adding subjects.
% Fixed when adding channels, the default type was SGL rather than EEG.
%
% modified 3/9/21 JD
% Timeshifts the data via interpolation if the time points of the added data are not the same.
%
% bugfix & modified 11/4/21 JD
% Added checkFlag.
% Updated which fields are required to be fully dimensional.
%
% bugfix & modified 4/2/25 JD
% Added video field.
% Fixed incorrect timeName labels when data appended via timepoint dimension.
% Added support for virtual grand averages.
% Fixed incorrectly initializing subNum, aveNum, and covNum fields with -1 (bad) instead of 0 (unknown).
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

global EPtictoc

EPdataOut=EPdataIn;

if ~exist('checkFlag','var')
    checkFlag=1;
end

numChans=length(EPdataIn.chanNames);
numPoints=length(EPdataIn.timeNames);
numCells=length(EPdataIn.cellNames);
numSubs=length(EPdataIn.subNames);
numVsubs=max(0,size(EPdataIn.GAVsubs,1)-1);
numRsubs=numSubs-numVsubs;
numVcells=max(0,size(EPdataIn.GAVsubs,2)-1);
numRcells=numCells-numVcells;
numFacs=length(EPdataIn.facNames);
numFreqs=length(EPdataIn.freqNames);
numRels=length(EPdataIn.relNames);
if numFacs==0
    numFacs=1;
end
if ~isempty(EPdataIn.facData)
    numCMBfacs=size(EPdataIn.facData,5);
else
    numCMBfacs=0;
end
numSGLfacs=numFacs-numCMBfacs;

if isfield(EPadd,'timeNames')
    if ~isempty(EPdataIn.timeNames) && ~isempty(setxor(EPdataIn.timeNames,EPadd.timeNames))
        disp(['Adjusting timing  to match that of the original data.']);
        [EPadd]=ep_interpTime(EPadd,EPdataIn.timeNames);
        ep_tictoc;if EPtictoc.stop;return;end
        if isempty(EPadd)
            disp('Error: The file was unable to be time shifted.');
            return
        end
    end
end

switch dataDimension
    case 'channels'
        if ~isfield(EPadd,'chanNames')
            disp('Error: No channel names specified.');
            return
        end
        numAdded=length(EPadd.chanNames);
        if isempty(EPdataIn.facVecS)
            if isfield(EPadd,'data')
                EPdataOut.data(end+1:end+numAdded,:,:,:,:,:,:)=EPadd.data;
            else
                EPdataOut.data(end+1:end+numAdded,:,:,:,:,:,:)=0;
            end
        end
        if ~isempty(EPdataIn.noise)
            if isfield(EPadd,'noise')
                EPdataOut.noise(end+1:end+numAdded,:,:,:,:,:,:)=EPadd.noise;
            else
                EPdataOut.noise(end+1:end+numAdded,:,:,:,:,:,:)=0;
            end
        end
        if ~isempty(EPdataIn.covAVE)
            EPdataOut.covAVE(end+1:end+numAdded,:,:,:,:,:,:)=NaN;
            if size(EPdataIn.covAVE,7)==1
                if isfield(EPadd,'covAVE')
                    EPdataOut.covAVE(numChans+1:numChans+numAdded,:,:,:,:,:,1)=EPadd.covAVE;
                end
            else
                EPdataOut.covAVE(:,:,:,:,:,:,end+1:end+numAdded)=NaN;
                if isfield(EPadd,'covAVE')
                    EPdataOut.covAVE(numChans+1:numChans+numAdded,:,:,:,:,:,numChans+1:numChans+numAdded)=EPadd.covAVE;
                end
            end
        end
        if ~isempty(EPdataIn.facVecS)
            if isfield(EPadd,'facVecS')
                EPdataOut.facVecS(end+1:end+numAdded,:)=EPadd.facVecS;
            else
                EPdataOut.facVecS(end+1:end+numAdded,:)=0;
            end
        end
        if ~isempty(EPdataIn.facData)
            if isfield(EPadd,'facData')
                EPdataOut.facData(end+1:end+numAdded,:,:,:,:,:,:)=EPadd.facData;
            else
                EPdataOut.facData(end+1:end+numAdded,:,:,:,:,:,:)=0;
            end
        end
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'badChans')
                EPdataOut.analysis.badChans(:,:,end+1:end+numAdded)=EPadd.analysis.badChans;
            else
                EPdataOut.analysis.badChans(:,:,end+1:end+numAdded)=0;
            end
        else
            EPdataOut.analysis.badChans(:,:,end+1:end+numAdded)=0;
        end
        if ~isempty(EPdataIn.cov)
            EPdataOut.cov.covMatrix(:,1:end+numAdded,:)=NaN;
            EPdataOut.cov.covMatrix(:,:,1:end+numAdded)=NaN;
            if isfield(EPadd,'cov')
                EPdataOut.cov.covMatrix(:,numChans+1:numChans+numAdded,numChans+1:numChans+numAdded)=EPadd.cov.covMatrix;
            end
        end
        EPdataOut.chanNames(end+1:end+numAdded,1)=EPadd.chanNames;
        if isfield(EPadd,'chanTypes')
            EPdataOut.chanTypes(end+1:end+numAdded,1)=EPadd.chanTypes;
        else
            [EPdataOut.chanTypes{end+1:end+numAdded,1}]=deal('EEG');
        end
        
        if ~isempty(EPdataIn.eloc)
            if isfield(EPadd,'eloc')
                EPdataOut.eloc(end+1:end+numAdded)=EPadd.eloc;
            else
                EPdataOut.eloc(end+1:end+numAdded)=cell2struct(cell(numAdded,length(fieldnames(EPdataIn.eloc))),fieldnames(EPdataIn.eloc),2);
            end
        end
        
        if ~isempty(EPdataIn.relNames)
            EPdataOut.relNames{end+1:end+numAdded,1}=EPadd.relNames;
            EPdataOut.data(end+1:end+numAdded,:,:,:,:,:,:)=NaN;
            EPdataOut.data(:,:,:,:,:,:,end+1:end+numAdded)=NaN;
            if ~isfield(EPadd,'data')
                for iChan=1:numAdded
                    EPdataOut.data(numChans+iChan,:,:,:,:,:,numChans+iChan)=1; %coherence of channel with itself is real number one.
                end
            end
        end
        
        if ~isempty(EPdataIn.impedances.channels)
            if isfield(EPadd,'impedances') && isfield(EPadd.impedances,'channels') && ~isempty(EPadd.impedances.channels)
                EPdataOut.impedances.channels(end+1:end+numAdded,:)=EPadd.impedances.channels;
            else
                EPdataOut.impedances.channels(end+1:end+numAdded,:)=NaN;
            end
        end
        
    case 'cells'
        
        if ~isfield(EPadd,'cellNames')
            disp('Error: No cell names specified.');
            return
        end
        if ~isfield(EPadd,'trialNames') && strcmp(EPdataIn.dataType,'single_trial')
            disp('Error: No trial names specified.');
            return
        end
        if ~isempty(EPdataIn.trialNames)
            if ~isfield(EPadd,'trialNames')
                disp('Error: No trial names specified.');
                return
            end
        end

        numAdded=length(EPadd.cellNames);
        numVadded=0;
        if isfield(EPadd,'GAVsubs')
            numVadded=size(EPadd.GAVsubs,2);
        end
        numRadded=numAdded-numVadded;
        if isfield(EPadd,'avgNum')
            EPdataOut.avgNum(:,end+1:end+numRadded)=EPadd.avgNum;
        else
            EPdataOut.avgNum(:,end+1:end+numRadded)=0;
        end
        if isfield(EPadd,'covNum')
            EPdataOut.covNum(:,end+1:end+numRadded)=EPadd.covNum;
        else
            EPdataOut.covNum(:,end+1:end+numRadded)=0;
        end
        if isfield(EPadd,'subNum')
            EPdataOut.subNum(:,end+1:end+numRadded)=EPadd.subNum;
        else
            EPdataOut.subNum(:,end+1:end+numRadded)=0;
        end
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'blinkTrial')
                EPdataOut.analysis.blinkTrial(:,end+1:end+numRadded)=EPadd.analysis.blinkTrial;
            else
                EPdataOut.analysis.blinkTrial(:,end+1:end+numRadded)=0;
            end
        else
            EPdataOut.analysis.blinkTrial(:,end+1:end+numRadded)=0;
        end
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'saccadeTrial')
                EPdataOut.analysis.saccadeTrial(:,end+1:end+numRadded)=EPadd.analysis.saccadeTrial;
            else
                EPdataOut.analysis.saccadeTrial(:,end+1:end+numRadded)=0;
            end
        else
            EPdataOut.analysis.saccadeTrial(:,end+1:end+numRadded)=0;
        end
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'saccadeOnset')
                EPdataOut.analysis.saccadeOnset(:,end+1:end+numRadded)=EPadd.analysis.saccadeOnset;
            else
                EPdataOut.analysis.saccadeOnset(:,end+1:end+numRadded)=0;
            end
        else
            EPdataOut.analysis.saccadeOnset(:,end+1:end+numRadded)=0;
        end
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'moveTrial')
                EPdataOut.analysis.moveTrial(:,end+1:end+numRadded)=EPadd.analysis.moveTrial;
            else
                EPdataOut.analysis.moveTrial(:,end+1:end+numRadded)=0;
            end
        else
            EPdataOut.analysis.moveTrial(:,end+1:end+numRadded)=0;
        end
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'badTrials')
                EPdataOut.analysis.badTrials(:,end+1:end+numRadded)=EPadd.analysis.badTrials;
            else
                EPdataOut.analysis.badTrials(:,end+1:end+numRadded)=0;
            end
        else
            EPdataOut.analysis.badTrials(:,end+1:end+numRadded)=0;
        end
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'badChans')
                EPdataOut.analysis.badChans(:,end+1:end+numRadded,:)=EPadd.analysis.badChans;
            else
                EPdataOut.analysis.badChans(:,end+1:end+numRadded,:)=0;
            end
        else
            EPdataOut.analysis.badChans(:,end+1:end+numRadded,:)=0;
        end
        if isfield(EPadd,'recTime')
            EPdataOut.recTime(end+1:end+numRadded)=EPadd.recTime;
        else
            EPdataOut.recTime(end+1:end+numRadded)=1;
        end
        
        if isfield(EPadd,'data')
            EPdataOut.data(:,:,end+1:end+numRadded,:,:,:,:)=EPadd.data;
        else
            EPdataOut.data(:,:,end+1:end+numRadded,:,:,:,:)=0;
        end
        
        if ~isempty(EPdataIn.noise)
            if isfield(EPadd,'noise') && ~isempty(EPadd.noise)
                EPdataOut.noise(:,:,end+1:end+numRadded,:,:,:,:)=EPadd.noise;
            else
                EPdataOut.noise(:,:,end+1:end+numRadded,:,:,:,:)=0;
            end
        end
        if ~isempty(EPdataIn.covAVE)
            EPdataOut.covAVE(:,:,end+1:end+numRadded,:,:,:,:)=NaN;
            if isfield(EPadd,'covAVE')
                EPdataOut.covAVE(:,:,end+1:end+numRadded,:,:,:,:)=EPadd.covAVE;
            end
        end
        if ~isempty(EPdataIn.facData)
            if isfield(EPadd,'facData') && ~isempty(EPadd.facData)
                EPdataOut.facData(:,:,end+1:end+numRadded,:,:,:,:)=EPadd.facData;
            else
                EPdataOut.facData(:,:,end+1:end+numRadded,:,:,:,:)=0;
            end
        end
        EPdataOut.cellNames=[EPdataOut.cellNames(1:numRcells);EPadd.cellNames(1:numRadded);EPdataOut.cellNames(end-numVcells+1:end);EPadd.cellNames(end-numVadded+1:end)];
        if isfield(EPadd,'cellTypes')
            EPdataOut.cellTypes=[EPdataOut.cellTypes(1:numRcells);EPadd.cellTypes(1:numRadded);EPdataOut.cellTypes(end-numVcells+1:end);EPadd.cellTypes(end-numVadded+1:end)];
        else
            tempCellTypes=repmat({'SGL'},numAdded,1);
            EPdataOut.cellTypes=[EPdataOut.cellTypes(1:numRcells);repmat({'SGL'},numRadded,1);EPdataOut.cellTypes(end-numVcells+1:end);repmat({'CMB'},numVadded,1)];
        end
        if isfield(EPadd,'trialNames') && ~isempty(EPadd.trialNames)
            cellNameList=unique(EPadd.cellNames);
            newTrials=zeros(numRadded,1);
            for iCell=1:length(cellNameList)
                theCell=cellNameList{iCell};
                whichCells=find(strcmp(theCell,EPadd.cellNames));
                theAddTrials=EPadd.trialNames(whichCells);
                theInTrials=EPdataIn.trialNames(find(strcmp(theCell,EPdataIn.cellNames)));
                if ~isempty(intersect(theAddTrials,theInTrials))
                    theAddTrials=theAddTrials+max(theInTrials);
                end
                newTrials(whichCells)=theAddTrials;
            end
            EPdataOut.trialNames(end+1:end+numRadded,1)=newTrials;
        end
        
        if isfield(EPadd,'trialSpecs') && ~isempty(EPadd.trialSpecs)
            if ~iscell(EPdataOut.trialSpecs) && isempty(EPdataOut.trialSpecs)
                EPdataOut.trialSpecs=cell(0);
            end
            EPdataOut.trialSpecs(end+1:end+numRadded,:,:)=cell(numRadded,size(EPdataIn.trialSpecs,2),size(EPdataIn.trialSpecs,3));
            for iSpec=1:length(EPadd.trialSpecNames)
                oldSpec=max(find(strcmp(EPadd.trialSpecNames{iSpec},EPdataOut.trialSpecNames)));
                if ~isempty(oldSpec)
                    EPdataOut.trialSpecs(end-numRadded+1:end,oldSpec,:)=EPadd.trialSpecs(:,iSpec,:);
                else
                    EPdataOut.trialSpecNames(end+1)=EPadd.trialSpecNames(iSpec);
                    EPdataOut.trialSpecs(end-numRadded+1:end,end+1,:)=EPadd.trialSpecs(:,iSpec,:);
                end
            end
        else
            if ~isempty(EPdataIn.trialSpecs)
                EPdataOut.trialSpecs(end+1:end+numRadded,:,:)=cell(numRadded,size(EPdataIn.trialSpecs,2),size(EPdataIn.trialSpecs,3));
            end
        end
        
%         if ~isempty(EPdataIn.events)
            if isfield(EPadd,'events') && ~isempty(EPadd.events)
                EPdataOut.events(:,end+1:end+numRadded)=EPadd.events;
            else
                EPdataOut.events(:,end+1:end+numRadded)=cell(size(EPdataIn.events,1),numRadded);
            end
            
            if ~isempty(EPdataIn.stims)
                for iStim=1:length(EPdataIn.stims)
                    if ~any(strcmp(EPdataIn.stims(iStim).name,{EPdataOut.stims.name}))
                        EPdataOut.stims(end+1)=EPdataIn.stims(iStim); %keep only stim images whose events are still in the segmented data.
                    end
                end
            end
            %         end

            if isfield(EPadd,'video')
                for iTrial=1:length(EPadd.video)
                    EPdataOut.video(numCells+iTrial).frames=EPadd.video(iTrial).frames;
                    EPdataOut.video(numCells+iTrial).times=EPadd.video(iTrial).times;
                end
            end

            if isfield(EPadd,'GAVsubs') && ~isempty(EPadd.GAVsubs)
                if isempty(EPdataOut.GAVsubs)
                    %initialize a new GAVsubs structure
                    EPdataOut.GAVsubs=cell(1,1,numFacs);
                end
                EPdataOut.GAVsubs(1,end+1:end+size(EPadd.GAVsubs,2),size(EPadd.GAVsubs,3))=EPadd.GAVsubs(1,:,:);
            end

    case 'subjects'

        if ~isfield(EPadd,'subNames')
            disp('Error: No subject names specified.');
            return
        end
        numAdded=length(EPadd.subNames);
        numVadded=0;
        if isfield(EPadd,'GAVsubs')
            numVadded=size(EPadd.GAVsubs,1);
        end
        numRadded=numAdded-numVadded;

        numSubs=length(EPdataIn.subNames);
        numVsubs=0;
        if ~isempty(EPdataIn.GAVsubs)
            numVsubs=max(0,size(EPdataIn.GAVsubs,1)-1);
        end
        numRsubs=numSubs-numVsubs;
        
        if ~isempty(EPdataIn.subjectSpecs)
            if isfield(EPadd,'subjectSpecs')
                EPdataOut.subjectSpecs(end+1:end+numRadded,:)=EPadd.subjectSpecs;
            else
                EPdataOut.subjectSpecs(end+1:end+numRadded,:)=cell(numRadded,size(EPdataIn.subjectSpecs,2));
            end
        end
        
        if ~isempty(EPdataIn.taskSpecs)
            if isfield(EPadd,'taskSpecs')
                EPdataOut.taskSpecs(end+1:end+numRadded,:,:)=EPadd.taskSpecs;
            else
                EPdataOut.taskSpecs(end+1:end+numRadded,:,:)=zeros(numRadded,size(EPdataIn.taskSpecs,2),size(EPdataIn.taskSpecs,3));
            end
        end
        
%         if ~isempty(EPdataIn.avgNum)
            if isfield(EPadd,'avgNum')
                EPdataOut.avgNum(end+1:end+numRadded,:)=EPadd.avgNum;
            else
                EPdataOut.avgNum(end+1:end+numRadded,:)=0;
            end
%         end
%         if ~isempty(EPdataIn.covNum)
            if isfield(EPadd,'covNum')
                EPdataOut.covNum(end+1:end+numRadded,:)=EPadd.covNum;
            else
                EPdataOut.covNum(end+1:end+numRadded,:)=0;
            end
%         end
%         if ~isempty(EPdataIn.subNum)
            if isfield(EPadd,'subNum')
                EPdataOut.subNum(end+1:end+numRadded,:)=EPadd.subNum;
            else
                EPdataOut.subNum(end+1:end+numRadded,:)=0;
            end
%         end
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'blinkTrial')
                EPdataOut.analysis.blinkTrial(end+1:end+numRadded,:)=EPadd.analysis.blinkTrial;
            else
                EPdataOut.analysis.blinkTrial(end+1:end+numRadded,:)=0;
            end
        else
            EPdataOut.analysis.blinkTrial(end+1:end+numRadded,:)=0;
        end
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'saccadeTrial')
                EPdataOut.analysis.saccadeTrial(end+1:end+numRadded,:)=EPadd.analysis.saccadeTrial;
            else
                EPdataOut.analysis.saccadeTrial(end+1:end+numRadded,:)=0;
            end
        else
            EPdataOut.analysis.saccadeTrial(end+1:end+numRadded,:)=0;
        end
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'saccadeOnset')
                EPdataOut.analysis.saccadeOnset(end+1:end+numRadded,:)=EPadd.analysis.saccadeOnset;
            else
                EPdataOut.analysis.saccadeOnset(end+1:end+numRadded,:)=0;
            end
        else
            EPdataOut.analysis.saccadeOnset(end+1:end+numRadded,:)=0;
        end
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'moveTrial')
                EPdataOut.analysis.moveTrial(end+1:end+numRadded,:)=EPadd.analysis.moveTrial;
            else
                EPdataOut.analysis.moveTrial(end+1:end+numRadded,:)=0;
            end
        else
            EPdataOut.analysis.moveTrial(end+1:end+numRadded,:)=0;
        end
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'badTrials')
                EPdataOut.analysis.badTrials(end+1:end+numRadded,:)=EPadd.analysis.badTrials;
            else
                EPdataOut.analysis.badTrials(end+1:end+numRadded,:)=0;
            end
        else
            EPdataOut.analysis.badTrials(end+1:end+numRadded,:)=0;
        end
        if isfield(EPadd,'analysis')
            if isfield(EPadd.analysis,'badChans')
                EPdataOut.analysis.badChans(end+1:end+numRadded,:,:)=EPadd.analysis.badChans;
            else
                EPdataOut.analysis.badChans(end+1:end+numRadded,:,:)=0;
            end
        else
            EPdataOut.analysis.badChans(end+1:end+numRadded,:,:)=0;
        end
        
        if isfield(EPadd,'data')
            EPdataOut.data(:,:,:,end+1:end+numRadded,:,:,:)=EPadd.data;
        else
            EPdataOut.data(:,:,:,end+1:end+numRadded,:,:,:)=0;
        end
        
        if ~isempty(EPdataIn.noise)
            if isfield(EPadd,'noise') && ~isempty(EPadd.noise)
                EPdataOut.noise(:,:,:,end+1:end+numRadded,:,:,:)=EPadd.noise;
            else
                EPdataOut.noise(:,:,:,end+1:end+numRadded,:,:,:)=0;
            end
        end
        if ~isempty(EPdataIn.covAVE)
            if isfield(EPadd,'covAVE')
                EPdataOut.covAVE(:,:,:,end+1:end+numRadded,:,:,:)=EPadd.covAVE;
            else
                EPdataOut.covAVE(:,:,:,end+1:end+numRadded,:,:,:)=NaN;
            end
        end
        if ~isempty(EPdataIn.facData)
            if isfield(EPadd,'facData') && ~isempty(EPadd.facData)
                EPdataOut.facData(:,:,:,end+1:end+numRadded,:,:,:)=EPadd.facData;
            else
                EPdataOut.facData(:,:,:,end+1:end+numRadded,:,:,:)=0;
            end
        end
        EPdataOut.subNames=[EPdataOut.subNames(1:numRsubs);EPadd.subNames(1:numRadded);EPdataOut.subNames(end-numVsubs+1:end);EPadd.subNames(end-numVadded+1:end)];
        if isfield(EPadd,'subTypes')
            EPdataOut.subTypes=[EPdataOut.subTypes(1:numRsubs);EPadd.subTypes(1:numRadded);EPdataOut.subTypes(end-numVsubs+1:end);EPadd.subTypes(end-numVadded+1:end)];
        else
            if strcmp(EPdataIn.dataType,'average')
                [EPdataOut.subTypes{end+1:end+numAdded,1}]=deal('AVG');
            else
                [EPdataOut.subTypes{end+1:end+numAdded,1}]=deal('RAW');
            end
        end
        if ~isempty(EPdataIn.cov)
            if isfield(EPadd,'cov') && ~isempty(EPadd.cov)
                EPdataOut.cov.covMatrix(end+1:end+numRadded,:,:)=EPadd.cov.covMatrix;
                EPdataOut.cov.Nq(end+1:end+numRadded)=EPadd.cov.Nq;
            else
                EPdataOut.cov.covMatrix(end+1:end+numRadded,:,:)=NaN;
                EPdataOut.cov.Nq(end+1:end+numRadded)=NaN;
            end
        end
        if isfield(EPadd,'events') && ~isempty(EPadd.events)
            EPdataOut.events(end+1:end+numRadded,:)=EPadd.events;
        else
            EPdataOut.events(end+1:end+numRadded,:)=cell(numRadded,size(EPdataIn.events,2));
        end

        if ~isempty(EPdataIn.stims)
            for iStim=1:length(EPdataIn.stims)
                if ~any(strcmp(EPdataIn.stims(iStim).name,{EPdataOut.stims.name}))
                    EPdataOut.stims(end+1)=EPdataIn.stims(iStim); %keep only stim images whose events are still in the segmented data.
                end
            end
        end
        
        if isfield(EPadd,'trialSpecs') && ~isempty(EPadd.trialSpecs)
            if ~iscell(EPdataOut.trialSpecs) && isempty(EPdataOut.trialSpecs)
                EPdataOut.trialSpecs=cell(0);
            end
            EPdataOut.trialSpecs(:,:,end+1:end+numRadded)=cell(size(EPdataIn.trialSpecs,1),size(EPdataIn.trialSpecs,2),numRadded);
            for iSpec=1:length(EPadd.trialSpecNames)
                oldSpec=find(strcmp(EPadd.trialSpecNames{iSpec},EPdataOut.trialSpecNames));
                if ~isempty(oldSpec)
                    EPdataOut.trialSpecs(:,oldSpec,end-numRadded+1:end)=EPadd.trialSpecs(:,iSpec,:);
                else
                    EPdataOut.trialSpecNames(end+1)=EPadd.trialSpecNames(iSpec);
                    EPdataOut.trialSpecs(:,end+1,end-numRadded+1:end)=EPadd.trialSpecs(:,iSpec,:);
                end
            end
        else
            if ~isempty(EPdataIn.trialSpecs)
                EPdataOut.trialSpecs(:,:,end+1:end+numRadded)=cell(size(EPdataIn.trialSpecs,1),size(EPdataIn.trialSpecs,2),numRadded);
            end
        end
        
        if ~isempty(EPdataIn.impedances.channels)
            if isfield(EPadd,'impedances') && isfield(EPadd.impedances,'channels') && ~isempty(EPadd.impedances.channels)
                EPdataOut.impedances.channels(:,end+1:end+numRadded)=EPadd.impedances.channels;
            else
                EPdataOut.impedances.channels(:,end+1:end+numRadded)=NaN;
            end
        end
        
        if ~isempty(EPdataIn.impedances.ground)
            if isfield(EPadd,'impedances') && isfield(EPadd.impedances,'ground') && ~isempty(EPadd.impedances.ground)
                EPdataOut.impedances.ground(end+1:end+numRadded)=EPadd.impedances.ground;
            else
                EPdataOut.impedances.ground(end+1:end+numRadded)=NaN;
            end
        end
        
        if ~isempty(EPdataOut.sessNums) || (isfield(EPadd,'sessNums') && ~isempty(EPadd.sessNums))
            EPdataOut.sessNums(end+1:end+numRadded)=0;
        end
        
        if isfield(EPadd,'sessNames')
            for iName=1:length(EPadd.sessNames)
                if ~isempty(find(strcmp(EPadd.sessNames{iName},EPdataOut.sessNames)))
                    EPdataOut.sessNames{end+1}=EPadd.sessNames{iName};
                end
            end
            if isfield(EPadd,'sessNums') && ~isempty(EPadd.sessNums)
                for iSess=1:length(EPadd.sessNums)
                    EPdataOut.sessNums(end-numRadded+iSess)=find(strcmp(EPadd.sessNames{EPadd.sessNum(iSess)},EPdataOut.sessNames));
                end
            end
        end

        if isfield(EPadd,'GAVsubs') && ~isempty(EPadd.GAVsubs)
            if isempty(EPdataOut.GAVsubs)
                %initialize a new GAVsubs structure
                EPdataOut.GAVsubs=cell(1,1,numFacs);
            end
            EPdataOut.GAVsubs(end+1:end+size(EPadd.GAVsubs,1),1,size(EPadd.GAVsubs,3))=EPadd.GAVsubs;
        end

    case 'factors'

        if ~isfield(EPadd,'facNames')
            disp('Error: No factor names specified.');
            return
        end
        numAdded=length(EPadd.facNames);
        
        if ~isempty(EPdataIn.facData)
            if isfield(EPadd,'facData')
                EPdataOut.facData(:,:,:,:,end+1:end+numAdded,:,:)=EPadd.facData;
            else
                EPdataOut.facData(:,:,:,:,end+1:end+numAdded,:,:)=0;
            end
        end
        
        EPdataOut.facNames(end+1:end+numAdded,1)=EPadd.facNames;
        if isfield(EPadd,'facTypes')
            EPdataOut.facTypes(end+1:end+numAdded,1)=EPadd.facTypes;
        else
            [EPdataOut.facTypes{end+1:end+numAdded,1}]=deal('SGL');
        end
        
    case 'points'
        if ~isfield(EPadd,'timeNames')
            disp('Error: No time names specified.');
            return
        end
        numAdded=length(EPadd.timeNames);
        if isempty(EPdataIn.facVecT)
            if isfield(EPadd,'data')
                EPdataOut.data(:,end+1:end+numAdded,:,:,:,:,:)=EPadd.data;
            else
                EPdataOut.data(:,end+1:end+numAdded,:,:,:,:,:)=0;
            end
        end
        if ~isempty(EPdataIn.noise)
            if isfield(EPadd,'noise')
                EPdataOut.noise(:,end+1:end+numAdded,:,:,:,:,:)=EPadd.noise;
            else
                EPdataOut.noise(:,end+1:end+numAdded,:,:,:,:,:)=0;
            end
        end
        if ~isempty(EPdataIn.covAVE)
            EPdataOut.covAVE(:,end+1:end+numAdded,:,:,:,:,:)=NaN;
            if isfield(EPadd,'covAVE')
                EPdataOut.covAVE(:,end+1:end+numAdded,:,:,:,:,:)=EPadd.covAVE;
            end
        end
        if ~isempty(EPdataIn.facVecT)
            if isfield(EPadd,'facVecT')
                EPdataOut.facVecS(end+1:end+numAdded,:)=EPadd.facVecT;
            else
                EPdataOut.facVecS(end+1:end+numAdded,:)=0;
            end
        end
        if ~isempty(EPdataIn.facData)
            if isfield(EPadd,'facData')
                EPdataOut.facData(:,end+1:end+numAdded,:,:,:,:,:)=EPadd.facData;
            else
                EPdataOut.facData(:,end+1:end+numAdded,:,:,:,:,:)=0;
            end
        end
        if strcmp(EPdataIn.dataType,'continuous')
            %adding two continuous segments together is funky because if the number of points do not divide evenly into one second epochs,
            %then the epochs in the new one can become disaligned.  They will therefore be reset to zero in this scenario.
            if rem(numPoints,ceil(EPdataIn.Fs))
                disp('Resetting edit fields to zero.')
                numEpochs=max(1,floor((numPoints+numAdded)/ceil(EPdataIn.Fs)));
                EPdataOut.analysis.blinkTrial=zeros(1,numEpochs);
                EPdataOut.analysis.saccadeTrial=zeros(1,numEpochs);
                EPdataOut.analysis.saccadeOnset=zeros(1,numEpochs);
                EPdataOut.analysis.moveTrial=zeros(1,numEpochs);
                EPdataOut.analysis.badTrials=zeros(1,numEpochs);
                EPdataOut.analysis.badChans=zeros(1,numEpochs,numChans);
            else
                numEpochAdded=floor(numAdded/ceil(EPdataIn.Fs));
                if isfield(EPadd,'analysis')
                    if isfield(EPadd.analysis,'badChans')
                        EPdataOut.analysis.badChans(:,end+1:end+numEpochAdded,:)=EPadd.analysis.badChans;
                    else
                        EPdataOut.analysis.badChans(:,end+1:end+numEpochAdded,:)=0;
                    end
                else
                    EPdataOut.analysis.badChans(:,end+1:end+numEpochAdded,:)=0;
                end
                if isfield(EPadd,'analysis')
                    if isfield(EPadd.analysis,'blinkTrial')
                        EPdataOut.analysis.blinkTrial(:,end+1:end+numEpochAdded)=EPadd.analysis.blinkTrial;
                    else
                        EPdataOut.analysis.blinkTrial(:,end+1:end+numEpochAdded)=0;
                    end
                else
                    EPdataOut.analysis.blinkTrial(:,end+1:end+numEpochAdded)=0;
                end
                if isfield(EPadd,'analysis')
                    if isfield(EPadd.analysis,'saccadeTrial')
                        EPdataOut.analysis.saccadeTrial(:,end+1:end+numEpochAdded)=EPadd.analysis.saccadeTrial;
                    else
                        EPdataOut.analysis.saccadeTrial(:,end+1:end+numEpochAdded)=0;
                    end
                else
                    EPdataOut.analysis.saccadeTrial(:,end+1:end+numEpochAdded)=0;
                end
                if isfield(EPadd,'analysis')
                    if isfield(EPadd.analysis,'saccadeOnset')
                        EPdataOut.analysis.saccadeOnset(:,end+1:end+numEpochAdded)=EPadd.analysis.saccadeOnset;
                    else
                        EPdataOut.analysis.saccadeOnset(:,end+1:end+numEpochAdded)=0;
                    end
                else
                    EPdataOut.analysis.saccadeOnset(:,end+1:end+numEpochAdded)=0;
                end
                if isfield(EPadd,'analysis')
                    if isfield(EPadd.analysis,'moveTrial')
                        EPdataOut.analysis.moveTrial(:,end+1:end+numEpochAdded)=EPadd.analysis.moveTrial;
                    else
                        EPdataOut.analysis.moveTrial(:,end+1:end+numEpochAdded)=0;
                    end
                else
                    EPdataOut.analysis.moveTrial(:,end+1:end+numEpochAdded)=0;
                end
                if isfield(EPadd,'analysis')
                    if isfield(EPadd.analysis,'badTrials')
                        EPdataOut.analysis.badTrials(:,end+1:end+numEpochAdded)=EPadd.analysis.badTrials;
                    else
                        EPdataOut.analysis.badTrials(:,end+1:end+numEpochAdded)=0;
                    end
                else
                    EPdataOut.analysis.badTrials(:,end+1:end+numEpochAdded)=0;
                end
            end
        else
            %if there is new information on bad data, then it replaces the original.
            if isfield(EPadd,'analysis')
                if isfield(EPadd.analysis,'badChans')
                    EPdataOut.analysis.badChans=EPadd.analysis.badChans;
                else
                    EPdataOut.analysis.badChans=0;
                end
            else
                EPdataOut.analysis.badChans=0;
            end
            if isfield(EPadd,'analysis')
                if isfield(EPadd.analysis,'blinkTrial')
                    EPdataOut.analysis.blinkTrial=EPadd.analysis.blinkTrial;
                else
                    EPdataOut.analysis.blinkTrial=0;
                end
            else
                EPdataOut.analysis.blinkTrial=0;
            end
            if isfield(EPadd,'analysis')
                if isfield(EPadd.analysis,'saccadeTrial')
                    EPdataOut.analysis.saccadeTrial=EPadd.analysis.saccadeTrial;
                else
                    EPdataOut.analysis.saccadeTrial=0;
                end
            else
                EPdataOut.analysis.saccadeTrial=0;
            end
            if isfield(EPadd,'analysis')
                if isfield(EPadd.analysis,'saccadeOnset')
                    EPdataOut.analysis.saccadeOnset=EPadd.analysis.saccadeOnset;
                else
                    EPdataOut.analysis.saccadeOnset=0;
                end
            else
                EPdataOut.analysis.saccadeOnset=0;
            end
            if isfield(EPadd,'analysis')
                if isfield(EPadd.analysis,'moveTrial')
                    EPdataOut.analysis.moveTrial=EPadd.analysis.moveTrial;
                else
                    EPdataOut.analysis.moveTrial=0;
                end
            else
                EPdataOut.analysis.moveTrial=0;
            end
            if isfield(EPadd,'analysis')
                if isfield(EPadd.analysis,'badTrials')
                    EPdataOut.analysis.badTrials=EPadd.analysis.badTrials;
                else
                    EPdataOut.analysis.badTrials=0;
                end
            else
                EPdataOut.analysis.badTrials=0;
            end
        end
        
        % if isfield(EPadd,'recTime')
        %     addTime=EPadd.recTime-EPdataIn.recTime;
        % else
        %     addTime=0;
        % end
        % addTime=addTime*(1000/EPdataIn.Fs);
        % if (EPadd.timeNames(1)+addTime) < EPdataIn.timeNames(end)
        %     disp('Time ranges overlap so appended times being renumbered.');
        %     addTime=EPdataIn.timeNames(end);
        % end
        % EPdataOut.timeNames(end+1:end+numAdded,1)=EPadd.timeNames+addTime;

        EPdataOut.timeNames(end+1:end+numAdded,1)=EPadd.timeNames+EPdataOut.timeNames(end)+median(diff(EPdataOut.timeNames));
        
        if ~isempty(EPdataIn.events)
            if isfield(EPadd,'events') && ~isempty(EPadd.events)
                for iSub=1:size(EPadd.events,1)
                    for iWave=1:size(EPadd.events,2)
                        newEvents=EPadd.events{iSub,iWave};
                        for iEvent=1:length(newEvents)
                            newEvents(iEvent).sample=newEvents(iEvent).sample+length(EPdataIn.timeNames);
                        end
                        if isempty(EPdataOut.events{iSub,iWave})
                            EPdataOut.events{iSub,iWave}=newEvents;
                        elseif ~isempty(newEvents)
                            EPdataOut.events{iSub,iWave}(end+1:end+length(newEvents))=newEvents;
                        end
                    end
                end
            end
        end

        if isfield(EPadd,'video')
            for iTrial=1:numCells
                if ~isempty(EPadd.video)
                    if ~isempty(EPdataOut.video(iTrial).frames)
                        numFrames=size(EPdataOut.video(iTrial).frames);
                    else
                        numFrames=0;
                    end
                    EPdataOut.video(iTrial).frames(numFrames+1:numFrames+length(EPadd.video(iTrial).frames))=EPadd.video(iTrial).frames;
                    EPdataOut.video(iTrial).times(numFrames+1:numFrames+length(EPadd.video(iTrial).frames))=EPadd.video(iTrial).times;
                end
            end
        end

    otherwise
        disp('Data dimension not recognized.');
        EPdataOut=[];
        return
end

if checkFlag
    [err]=ep_checkEPfile(EPdataOut);
    if err
        EPdataOut=[];
    end
end
