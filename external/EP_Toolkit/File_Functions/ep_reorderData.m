function [EPdataOut]=ep_reorderData(EPdataIn,dataDimension,sortOrder)
%  [EPdataOut]=ep_reorderData(EPdataIn,dataDimension,sortOrder);
%       Reorders data according to specified order and dimension. 
%
%Inputs:
%  EPdataIn       : Structured array with the input data and accompanying information in EP file format.  See readData.
%  dataDimension  : The data dimension ('cells', 'subjects', 'channels', 'factors').
%  sortOrder       : The new order of the selected data dimension.
%
%Outputs:
%  EPdataOut      : Structured array with the output data and accompanying information in EP file format.  See readData.

%History:
%  by Joseph Dien (6/12/18)
%  jdien07@mac.com
%
% bufix 12/2/18 JD
% Fixed not reordering impedance values when reordering channels using Edit function.
%
% bufix 12/13/18 JD
% Fixed crash when reordering by cells data with no trial specs.
%
% bufix & modified 1/8/19 JD
% Fixed not handling reordering of combined factors correctly.
% Now allows CMB and SGL factors to be intermixed.
% FacVar and FacVarQ now include CMB factors.
%
% bugfix 2/12/19 JD
% Fixed crash when sorting data by subject and there are no subject specs.
%
% modified 4/9/19 JD
% Added support for task level performance measures.
%
% modified 11/4/19 JD
% Added sessNums sessNames fields.
%
% modified 12/24/19 JD
% Upgraded support of std information by adding .covAVE and .GAVsubs fields and eliminating .std and .stdCM fields.
%
% bugfix 1/20/20 JD
% Fixed crash when reordering subject and there are no session numbers.
%
% bugfix 5/18/24 JD
% Fixed crash when reordering channels and there are no eloc info.
%
% bugfix & modified 3/21/25 JD
% Added video field.
% Fixed .GAVsubs not being updated when subjects reordered.
% Added support for virtual grand averages.
%
% bugfix 6/11/25 JD
% Fixed not supporting virtual cells.
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

EPdataOut=EPdataIn;

numSubs=length(EPdataIn.subNames);
numVsubs=max(0,size(EPdataIn.GAVsubs,1)-1);
numRsubs=numSubs-numVsubs;
numCells=length(EPdataIn.cellNames);
numVcells=max(0,size(EPdataIn.GAVsubs,2)-1);
numRcells=numCells-numVcells;

%[temp2,sortOrder] = sortrows(sortOrder(:),1);

switch dataDimension
    case 'cells'
        if ~isempty(EPdataIn.GAVsubs)
            if any(sortOrder(1:numRcells)>numRcells) || any(sortOrder(numRcells+1:end)<=numRcells)
                %if any of the virtual cells are moved from the end or vice versa
                disp('Error: Virtual cells must be at the end of the cells slots.')
                EPdataOut=[];
                return
            end
            EPdataOut.GAVsubs(:,2:numVcells+1,:)=EPdataIn.GAVsubs(:,sortOrder(numRcells+1:end)-numRcells+1,:);
            for iCell=1:numVcells
                %convert the list of cells to the new order numbers.
                for iFac=1:length(EPdataOut.facNames)
                    for iRow=1:length(EPdataOut.GAVsubs{1,iCell+1,iFac})
                        EPdataOut.GAVsubs{1,iCell+1,iFac}(iRow,1)=sortOrder(EPdataOut.GAVsubs{1,iCell+1,iFac}(iRow,1));
                    end
                end
            end
        end
        if ~isempty(EPdataIn.trialSpecs)
            EPdataOut.trialSpecs=EPdataIn.trialSpecs(sortOrder(1:numRcells),:,:);
        end
        EPdataOut.avgNum=EPdataIn.avgNum(:,sortOrder(1:numRcells));
        EPdataOut.covNum=EPdataIn.covNum(:,sortOrder(1:numRcells));
        EPdataOut.subNum=EPdataIn.subNum(:,sortOrder(1:numRcells));
        EPdataOut.analysis.blinkTrial=EPdataIn.analysis.blinkTrial(:,sortOrder(1:numRcells));
        EPdataOut.analysis.saccadeTrial=EPdataIn.analysis.saccadeTrial(:,sortOrder(1:numRcells));
        EPdataOut.analysis.saccadeOnset=EPdataIn.analysis.saccadeOnset(:,sortOrder(1:numRcells));
        EPdataOut.analysis.moveTrial=EPdataIn.analysis.moveTrial(:,sortOrder(1:numRcells));
        EPdataOut.analysis.badTrials=EPdataIn.analysis.badTrials(:,sortOrder(1:numRcells));
        EPdataOut.analysis.badChans=EPdataIn.analysis.badChans(:,sortOrder(1:numRcells),:);
        EPdataOut.recTime=EPdataIn.recTime(sortOrder(1:numRcells));
        EPdataOut.data=EPdataIn.data(:,:,sortOrder(1:numRcells),:,:,:,:);
        if ~isempty(EPdataIn.facData)
            EPdataOut.facData=EPdataIn.facData(:,:,sortOrder(1:numRcells),:,:,:,:);
        end
        if ~isempty(EPdataIn.noise)
            EPdataOut.noise=EPdataIn.noise(:,:,sortOrder(1:numRcells),:,:,:,:);
        end
        if ~isempty(EPdataIn.covAVE)
            EPdataOut.covAVE=EPdataIn.covAVE(:,:,sortOrder(1:numRcells),:,:,:,:);
        end
        if ~isempty(EPdataIn.GAVsubs) && any(sortOrder>numRcells)
            EPdataOut.GAVsubs(:,2:end,:)=EPdataIn.GAVsubs(:,sortOrder(sortOrder>numRcells)-numRcells+1,:);
        end
        EPdataOut.cellNames=EPdataIn.cellNames(sortOrder);
        EPdataOut.cellTypes=EPdataIn.cellTypes(sortOrder); 
        if ~isempty(EPdataIn.trialNames)
            EPdataOut.trialNames=EPdataIn.trialNames(sortOrder(1:numRcells));
        end
        EPdataOut.events=EPdataIn.events(:,sortOrder(1:numRcells));
        if ~isempty(EPdataIn.video)
            EPdataOut.video=EPdataIn.video(sortOrder(1:numRcells));
        end

    case 'subjects'
        if ~isempty(EPdataIn.GAVsubs)
            if any(sortOrder(1:numRsubs)>numRsubs) || any(sortOrder(numRsubs+1:end)<=numRsubs)
                %if any of the virtual grand averages are moved from the end or vice versa
                disp('Error: Virtual grand averages must be at the end of the subjects slots.')
                EPdataOut=[];
                return
            end
            EPdataOut.GAVsubs(2:numVsubs+1,:,:)=EPdataIn.GAVsubs(sortOrder(numRsubs+1:end)-numRsubs+1,:,:);
            for iGAV=1:numVsubs
                %convert the list of subjects to the new order numbers.
                for iCell=1:length(EPdataOut.cellNames)
                    for iFac=1:length(EPdataOut.facNames)
                        for iRow=1:length(EPdataOut.GAVsubs{iGAV+1,iCell+1,iFac})
                            EPdataOut.GAVsubs{iGAV+1,iCell+1,iFac}(iRow,1)=sortOrder(EPdataOut.GAVsubs{iGAV+1,iCell+1,iFac}(iRow,1));
                        end
                    end
                end
            end
        end
        if ~isempty(EPdataIn.subjectSpecs)
            EPdataOut.subjectSpecs=EPdataIn.subjectSpecs(sortOrder(1:numRsubs),:);
        end
        if ~isempty(EPdataIn.taskSpecs)
            EPdataOut.taskSpecs=EPdataIn.taskSpecs(sortOrder(1:numRsubs),:,:);
        end
        EPdataOut.avgNum=EPdataIn.avgNum(sortOrder(1:numRsubs),:);
        EPdataOut.covNum=EPdataIn.covNum(sortOrder(1:numRsubs),:);
        EPdataOut.subNum=EPdataIn.subNum(sortOrder(1:numRsubs),:);
        EPdataOut.analysis.blinkTrial=EPdataIn.analysis.blinkTrial(sortOrder(1:numRsubs),:);
        EPdataOut.analysis.saccadeTrial=EPdataIn.analysis.saccadeTrial(sortOrder(1:numRsubs),:);
        EPdataOut.analysis.saccadeOnset=EPdataIn.analysis.saccadeOnset(sortOrder(1:numRsubs),:);
        EPdataOut.analysis.moveTrial=EPdataIn.analysis.moveTrial(sortOrder(1:numRsubs),:);
        EPdataOut.analysis.badTrials=EPdataIn.analysis.badTrials(sortOrder(1:numRsubs),:);
        EPdataOut.analysis.badChans=EPdataIn.analysis.badChans(sortOrder(1:numRsubs),:,:);
        EPdataOut.data=EPdataIn.data(:,:,:,sortOrder(1:numRsubs),:,:,:);
        if ~isempty(EPdataIn.facData)
            EPdataOut.facData=EPdataIn.facData(:,:,:,sortOrder(1:numRsubs),:,:,:);
        end
        if ~isempty(EPdataIn.noise)
            EPdataOut.noise=EPdataIn.noise(:,:,:,sortOrder(1:numRsubs),:,:,:);
        end
        if ~isempty(EPdataIn.covAVE)
            EPdataOut.covAVE=EPdataIn.covAVE(:,:,:,sortOrder(1:numRsubs),:,:,:);
        end
        EPdataOut.subNames=EPdataIn.subNames(sortOrder);
        EPdataOut.subTypes=EPdataIn.subTypes(sortOrder);
        if ~isempty(EPdataIn.sessNums)
            EPdataOut.sessNums=EPdataIn.sessNums(sortOrder(1:numRsubs));
        end
        EPdataOut.events=EPdataIn.events(sortOrder(1:numRsubs),:);
        if ~isempty(EPdataIn.impedances.channels)
            EPdataOut.impedances.channels=EPdataIn.impedances.channels(:,sortOrder(1:numRsubs));
        end
        if ~isempty(EPdataIn.impedances.ground)
            EPdataOut.impedances.ground=EPdataIn.impedances.ground(sortOrder(1:numRsubs));
        end
        
    case 'channels'
        
        if ~isempty(EPdataIn.facVecS)
            EPdataOut.facVecS=EPdataIn.facVecS(sortOrder,:);
        else
            EPdataOut.data=EPdataIn.data(sortOrder,:,:,:,:,:,:);
        end
        if ~isempty(EPdataIn.noise)
            EPdataOut.noise=EPdataIn.noise(sortOrder,:,:,:,:);
        end
        if ~isempty(EPdataIn.cov)
            EPdataOut.cov.covMatrix=EPdataIn.cov.covMatrix(:,sortOrder,sortOrder);
        end
        if ~isempty(EPdataIn.covAVE)
            if size(EPdataOut.covAVE,7)==1
                EPdataOut.covAVE=EPdataIn.covAVE(sortOrder,:,:,:,:,:,1);
            else
                EPdataOut.covAVE=EPdataIn.covAVE(sortOrder,:,:,:,:,:,sortOrder);
            end
        end
        if ~isempty(EPdataIn.facData)
            EPdataOut.facData=EPdataIn.facData(sortOrder,:,:,:,:,:,:);
        end
        
        reference.type=EPdataIn.reference.type;
        reference.original=[];
        reference.current=[];
        
        for i=1:length(EPdataIn.reference.original)
            reference.original(i)=find(EPdataIn.reference.original(i) == sortOrder);
        end
        for i=1:length(EPdataIn.reference.current)
            reference.current(i)=find(EPdataIn.reference.current(i) == sortOrder);
        end
        EPdataOut.reference=reference;
        
        EPdataOut.chanNames=EPdataIn.chanNames(sortOrder);
        EPdataOut.chanTypes=EPdataIn.chanTypes(sortOrder);
        if ~isempty(EPdataIn.eloc)
            EPdataOut.eloc=EPdataIn.eloc(sortOrder);
        end
        EPdataOut.analysis.badChans=EPdataIn.analysis.badChans(:,:,sortOrder);
        
        if ~isempty(EPdataIn.impedances.channels)
            EPdataOut.impedances.channels=EPdataIn.impedances.channels(sortOrder,:);
        end
        
        
    case 'factors'
        EPdataOut.facNames=EPdataIn.facNames(sortOrder);
        EPdataOut.facTypes=EPdataIn.facTypes(sortOrder);
        
        SGLfacs=find(strcmp('SGL',EPdataIn.facTypes));
        CMBfacs=find(strcmp('CMB',EPdataIn.facTypes));
        sortOrderSGL=[];
        sortOrderCMB=[];
        for iFactor=1:length(sortOrder)
            theFactor=sortOrder(iFactor);
            if strcmp('SGL',EPdataIn.facTypes{sortOrder(iFactor)})
                sortOrderSGL(end+1)=find(SGLfacs==theFactor);
            else
                sortOrderCMB(end+1)=find(CMBfacs==theFactor);
            end
        end        
        
        EPdataOut.data=EPdataIn.data(:,:,:,:,sortOrderSGL,:,:);
        if ~isempty(EPdataIn.facData)
            EPdataOut.facData=EPdataIn.facData(:,:,:,:,sortOrderCMB,:,:);
        end
        if ~isempty(EPdataIn.facVecT)
            EPdataOut.facVecT=EPdataIn.facVecT(:,sortOrderSGL);
        end
        if ~isempty(EPdataIn.facVecS)
            EPdataOut.facVecS=EPdataIn.facVecS(:,sortOrderSGL);
        end
        if ~isempty(EPdataIn.facVar)
            EPdataOut.facVar=EPdataIn.facVar(:,sortOrder);
        end
        if ~isempty(EPdataIn.facVarQ)
            EPdataOut.facVarQ=EPdataIn.facVarQ(:,sortOrder);
        end
        
    otherwise
        disp('Data dimension not recognized.');
        EPdataOut=[];
        return
end

[err]=ep_checkEPfile(EPdataOut);
if err
    EPdataOut=[];
end
