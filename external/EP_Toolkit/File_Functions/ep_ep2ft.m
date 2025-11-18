function [FTdata eventHdr]=ep_ep2ft(EPdata)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [FTdata eventHdr]=ep_ep2ft(EPdata)
% Converts EPdata data format to FieldTrip data format.  Assumes for FT that trials are in order and were contiguous.
%
%Inputs
%   EPdata: EP data format structured variable (see ep_readData)
%
%Outputs
%	FTdata: FieldTrip data format structured variable.
%   eventHdr: FieldTrip event structured variable.

% History:
%
% by Joseph Dien (5/25/11)
% jdien07@mac.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% bugfix 3/17/14 JD
% Fixed formation of sampleinfo field.
%
% bugfix 8/5/14 JD
% Changed output to be appropriate for both raw and average formats with the time field in the correction orientation.
% Added dimord field.
%
% modified 3/16/20 JD
% Added support for frequency-domain files.
%
% bugfix 5/28/20 JD
% Fixed crash when converting data with multiple events per cell, such as continuous data.
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

FTdata=[];
eventHdr=[];

if isempty(EPdata)
    return
end

if length(EPdata.subNames) > 1
    msg{1}=['Error: Cannot convert multi-subject EPdata data into FieldTrip data format.'];
    [msg]=ep_errorMsg(msg);
    return
end

if ~isempty(EPdata.facNames)
    msg{1}=['Error: Cannot convert EPdata factor data into FieldTrip data format.'];
    [msg]=ep_errorMsg(msg);
    return
end

if ~isempty(EPdata.relNames)
    msg{1}=['Error: Cannot convert EPdata coherence data into FieldTrip data format.'];
    [msg]=ep_errorMsg(msg);
    return
end

numTrials=length(EPdata.cellNames);
numPoints=length(EPdata.timeNames);
numChans=length(EPdata.chanNames);
numFreqs=length(EPdata.freqNames);

FTdata.hdr.Fs=EPdata.Fs;
FTdata.hdr.nChans=numChans;
FTdata.hdr.label=EPdata.chanNames;
FTdata.hdr.nTrials=numTrials;
FTdata.hdr.nSamplesPre=EPdata.baseline/(1000/EPdata.Fs);
FTdata.hdr.nSamples=numPoints;
FTdata.label=EPdata.chanNames;

if numFreqs
    FTdata.fourierspctrm=permute(EPdata.data,[3 1 6 2 4 5]);
    FTdata.dimord='rpt_chan_freq_time';
    FTdata.freq=EPdata.freqNames;
    FTdata.time=EPdata.timeNames;
else
    if strcmp(EPdata.dataType,'average')
        FTdata.time=zeros(numTrials,numPoints);
        FTdata.avg=zeros(numTrials,numChans,numPoints);
        FTdata.sampleinfo=zeros(numTrials,2);
        for theTrial=1:numTrials
            FTdata.time(theTrial,:)=(EPdata.timeNames')/1000; %time in seconds
            FTdata.avg(theTrial,:,:)=EPdata.data(:,:,theTrial);
            FTdata.sampleinfo(theTrial,1)=((theTrial-1)*numPoints)+1;
            FTdata.sampleinfo(theTrial,2)=theTrial*numPoints;
        end
        FTdata.dimord='rpt_chan_time';
    else
        FTdata.time=cell(numTrials,1);
        FTdata.trial=cell(numTrials,1);
        FTdata.sampleinfo=zeros(numTrials,2);
        for theTrial=1:numTrials
            FTdata.time{theTrial,1}=(EPdata.timeNames')/1000; %time in seconds
            FTdata.trial{theTrial,1}=EPdata.data(:,:,theTrial);
            FTdata.sampleinfo(theTrial,1)=((theTrial-1)*numPoints)+1;
            FTdata.sampleinfo(theTrial,2)=theTrial*numPoints;
        end
        FTdata.dimord='{rpt}_chan_time';
    end
    FTdata.fsample=EPdata.Fs;
end

for iSub=1:length(EPdata.subNames)
    for iCell=1:length(EPdata.cellNames)
        %add trial events that mark the start of each epoch
        eventHdr(end+1).type='trial';
        eventHdr(end).value=EPdata.cellNames{iCell};
        eventHdr(end).sample=EPdata.recTime(iCell);
        eventHdr(end).duration=0;
        eventHdr(end).offset=-EPdata.baseline;
        for iEvent=1:length(EPdata.events{iSub,iCell})
            eventHdr(end+1).type=EPdata.events{iSub,iCell}(iEvent).type;
            eventHdr(end).value=EPdata.events{iSub,iCell}(iEvent).value;
            eventHdr(end).sample=EPdata.events{iSub,iCell}(iEvent).sample+EPdata.recTime(iCell);
            eventHdr(end).duration=EPdata.events{iSub,iCell}(iEvent).duration;
            eventHdr(end).offset=[];
        end
    end
end