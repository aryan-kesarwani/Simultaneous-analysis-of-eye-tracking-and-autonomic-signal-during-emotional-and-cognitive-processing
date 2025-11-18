function [EPdataOut]=ep_interpTime(EPdataIn,newTimeNames)
%  [EPdataOut]=ep_interpTime(EPdataIn,newTimeNames)
%       Interpolates new voltage time series based on the newTimes vector.
%
%Inputs:
%  EPdataIn       : Structured array with the input data and accompanying information in EP file format.  See readData.
%  newTimeNames   : Vector with new times, on the same scale as the old times.
%
%Outputs:
%  EPdataOut      : Structured array with the output data and accompanying information in EP file format.  See readData.
%
% When the interpolations go up against the edge of the time series, the first and last time points are used as the estimates of the t0-1 and tEnd+1 voltages.
% Handles two scenarios.  In the time-shift scenario, the difference in the first sample time is taken as the degree of the shift to be applied across the events.
% In the sampling rate change scenario, the first sample is always taken as being sample "1" regardless of how the rate changed.
%History:
%  by Joseph Dien (3/8/21)
%  jdien07@mac.com
%
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

EPdataOut=[];
oldSampDur = median(diff(EPdataIn.timeNames));

if length(EPdataIn.timeNames) < 2
    disp('Error: Cannot interpolate time without at least two timepoints.');
    return;
end

numPreOut=sum(newTimeNames < (EPdataIn.timeNames(1)-oldSampDur));
numPostOut=sum(newTimeNames > (EPdataIn.timeNames(end)+oldSampDur));
numPoints=length(EPdataIn.timeNames);
EPdataOut=EPdataIn;
EPdataOut.facVecT=nan(size(EPdataIn.facVecT,1),length(newTimeNames),size(EPdataIn.facVecT,3),size(EPdataIn.facVecT,4),size(EPdataIn.facVecT,5),size(EPdataIn.facVecT,6),size(EPdataIn.facVecT,7));
EPdataOut.data=nan(size(EPdataIn.data,1),length(newTimeNames),size(EPdataIn.data,3),size(EPdataIn.data,4),size(EPdataIn.data,5),size(EPdataIn.data,6),size(EPdataIn.data,7));
EPdataOut.noise=nan(size(EPdataIn.noise,1),length(newTimeNames),size(EPdataIn.noise,3),size(EPdataIn.noise,4),size(EPdataIn.noise,5),size(EPdataIn.noise,6),size(EPdataIn.noise,7));
EPdataOut.covAVE=nan(size(EPdataIn.covAVE,1),length(newTimeNames),size(EPdataIn.covAVE,3),size(EPdataIn.covAVE,4),size(EPdataIn.covAVE,5),size(EPdataIn.covAVE,6),size(EPdataIn.covAVE,7));

if (numPreOut+numPostOut) < numPoints %if any time points fall within the range of the data
    newTimeRange=newTimeNames(1+numPreOut:end-numPostOut);
    
    beforeTime=EPdataIn.timeNames(1)-oldSampDur;
    afterTime=EPdataIn.timeNames(end)+oldSampDur;
    oldtimeNames=[beforeTime; EPdataIn.timeNames; afterTime];
    
    if ~isempty(EPdataIn.facVecT)
        EPdataIn.facVecT=interp1(oldtimeNames,[EPdataIn.facVecT(1) EPdataIn.facVecT EPdataIn.facVecT(end)],newTimeNames,'linear','extrap');
    else
        newData=zeros(size(EPdataIn.data,1),length(newTimeRange),size(EPdataIn.data,3),size(EPdataIn.data,4),size(EPdataIn.data,5),size(EPdataIn.data,6),size(EPdataIn.data,7));
        for iCell=1:size(EPdataIn.data,3)
            ep_tictoc;if EPtictoc.stop;return;end
            for iChan=1:size(EPdataIn.data,1)
                for iSub=1:size(EPdataIn.data,4)
                    for iFac=1:size(EPdataIn.data,5)
                        for iFreq=1:size(EPdataIn.data,6)
                            for iRel=1:size(EPdataIn.data,7)
                                newData(iChan,:,iCell,iSub,iFac,iFreq,iRel)=interp1(oldtimeNames,[EPdataIn.data(iChan,1,iCell,iSub,iFac,iFreq,iRel) squeeze(EPdataIn.data(iChan,:,iCell,iSub,iFac,iFreq,iRel)) EPdataIn.data(iChan,end,iCell,iSub,iFac,iFreq,iRel)],newTimeNames,'linear','extrap')';
                            end
                        end
                    end
                end
            end
        end
        EPdataOut.data(:,1+numPreOut:length(newTimeRange),:,:,:,:,:)=newData;
    end
    
    if ~isempty(EPdataIn.noise)
        newData=zeros(size(EPdataIn.noise,1),length(newTimeRange),size(EPdataIn.noise,3),size(EPdataIn.noise,4),size(EPdataIn.noise,5),size(EPdataIn.noise,6),size(EPdataIn.noise,7));
        for iCell=1:size(EPdataIn.noise,3)
            ep_tictoc;if EPtictoc.stop;return;end
            for iChan=1:size(EPdataIn.noise,1)
                for iSub=1:size(EPdataIn.noise,4)
                    for iFac=1:size(EPdataIn.noise,5)
                        for iFreq=1:size(EPdataIn.noise,6)
                            for iRel=1:size(EPdataIn.noise,7)
                                newData(iChan,:,iCell,iSub,iFac,iFreq,iRel)=interp1(oldtimeNames,[EPdataIn.noise(iChan,1,iCell,iSub,iFac,iFreq,iRel) squeeze(EPdataIn.noise(iChan,:,iCell,iSub,iFac,iFreq,iRel)) EPdataIn.data(iChan,end,iCell,iSub,iFac,iFreq,iRel)],newTimeNames,'linear','extrap')';
                            end
                        end
                    end
                end
            end
        end
        EPdataOut.noise(:,1+numPreOut:length(newTimeRange),:,:,:,:,:)=newData;
    end
    
    if ~isempty(EPdataIn.covAVE)
        newData=zeros(size(EPdataIn.covAVE,1),length(newTimeRange),size(EPdataIn.covAVE,3),size(EPdataIn.covAVE,4),size(EPdataIn.covAVE,5),size(EPdataIn.covAVE,6),size(EPdataIn.covAVE,7));
        for iCell=1:size(EPdataIn.covAVE,3)
            ep_tictoc;if EPtictoc.stop;return;end
            for iChan=1:size(EPdataIn.covAVE,1)
                for iSub=1:size(EPdataIn.covAVE,4)
                    for iFac=1:size(EPdataIn.covAVE,5)
                        for iFreq=1:size(EPdataIn.covAVE,6)
                            for iRel=1:size(EPdataIn.covAVE,7)
                                newData(iChan,:,iCell,iSub,iFac,iFreq,iRel)=interp1(oldtimeNames,sqrt([EPdataIn.covAVE(iChan,1,iCell,iSub,iFac,iFreq,iRel) squeeze(EPdataIn.covAVE(iChan,:,iCell,iSub,iFac,iFreq,iRel)) EPdataIn.data(iChan,end,iCell,iSub,iFac,iFreq,iRel)]),newTimeNames,'linear','extrap')'.^2;
                            end
                        end
                    end
                end
            end
        end
        EPdataOut.covAVE(:,1+numPreOut:length(newTimeRange),:,:,:,:,:)=newData;
    end
    
    if ~isempty(EPdataIn.facData)
        newData=zeros(size(EPdataIn.facData,1),length(newTimeRange),size(EPdataIn.facData,3),size(EPdataIn.facData,4),size(EPdataIn.facData,5),size(EPdataIn.facData,6),size(EPdataIn.facData,7));
        for iCell=1:size(EPdataIn.facData,3)
            ep_tictoc;if EPtictoc.stop;return;end
            for iChan=1:size(EPdataIn.facData,1)
                for iSub=1:size(EPdataIn.facData,4)
                    for iFac=1:size(EPdataIn.facData,5)
                        for iFreq=1:size(EPdataIn.facData,6)
                            for iRel=1:size(EPdataIn.facData,7)
                                newData(iChan,:,iCell,iSub,iFac,iFreq,iRel)=interp1(oldtimeNames,[EPdataIn.facData(iChan,1,iCell,iSub,iFac,iFreq,iRel) squeeze(EPdataIn.facData(iChan,:,iCell,iSub,iFac,iFreq,iRel)) EPdataIn.data(iChan,end,iCell,iSub,iFac,iFreq,iRel)],newTimeNames,'linear','extrap')';
                            end
                        end
                    end
                end
            end
        end
        EPdataOut.facData(:,1+numPreOut:length(newTimeRange),:,:,:,:,:)=newData;
    end
end

EPdataOut.timeNames=newTimeNames;
newSampDur=median(diff(newTimeNames));
rateChange=oldSampDur/newSampDur;
timeShift=(newTimeNames(1)-EPdataIn.timeNames(1))/newSampDur; %in new samples

for iSub = 1:length(EPdataIn.subNames)
    for iWave = 1:length(EPdataIn.cellNames)
        if ~isempty(EPdataIn.events{iSub,iWave})
            for iEvent = 1:length(EPdataIn.events{iSub,iWave})
                EPdataOut.events{iSub,iWave}(iEvent).sample=((EPdataIn.events{iSub,iWave}(iEvent).sample)*rateChange)-timeShift-(rateChange-1);
            end
            outRange=([EPdataIn.events{iSub,iWave}.sample]<1) | ([EPdataIn.events{iSub,iWave}.sample]>=numPoints+1);
            EPdataOut.events{iSub,iWave}(outRange)=[];
        end
    end
end

EPdataOut.baseline=EPdataIn.baseline*rateChange;
EPdataOut.Fs=1000/newSampDur;
EPdataOut.recTime=EPdataIn.recTime*rateChange-(rateChange-1);
