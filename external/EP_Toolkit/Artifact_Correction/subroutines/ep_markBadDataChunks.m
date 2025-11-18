function [outputLog] = ep_markBadDataChunks(inFile, startChunk, endChunk, badChans, theSubject);
% ep_markBadDataChunks(inFile, startChunk, endChunk, badChans, theSubject);
%	Marks bad channels and trials with a flat line interrupted by a huge spike for another program like NetStation to fix.
%
%Inputs
%	inFile:     filename (not including the .mat suffix or chunk number.  e.g., "NT5") and sourcepath.
%	startChunk: starting chunk (usually 1)
%   endChunk:   ending chunk
%   badChans:   list of globally bad channels.  Will be set to a flat line.
%   theSubject: which subject of the file is being processed.
%
%   The input chunks include: dataChunk in EP data format.
%
%Outputs
%   outputLog: output messages from bad data process.
%
%   Updated output chunks: dataChunk in EP data format.
%
% History:
%
% by Joseph Dien (2/09)
% jdien07@mac.com
%
% modified 3/14/09 JD
% Changed to use EP format data to provide more flexibility with I/O functions.
%
% modified 4/17/09 JD
% Dropped eloc as separate input parameter (now part of data).
%
% bugfix 12/8/09 JD
% Fixed crash when there are bad channels.  Thanks to Alex Lamey.
%
% modified 2/11/10 JD
% Will now work with subject average files with multiple subjects.
% Gets bad channel and bad trial info from the data chunk rather than from the function call.
%
% bufix 3/11/14 JD
% Handles decimal sampling rates gracefully.
%
% bugfix 12/5/18 JD
% Fixed good epochs in average files being treated as bad (and being zeroed) erroneously while bad epochs not zeroed.
% Fixed good channels in average files being treated as bad (and being replaced) erroneously while bad channels not replaced.
% Fixed bad EOG channels that were temporarily interpolated not being fixed by the replace bad channel option.
%
% bugfix 2/20/19 JD
% Fixed crash when data are chunked.
% 
% modified 6/23/19 JD
% Accelerates artifact correction by adding option to keep chunks in RAM.
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

global EPchunk

msg='Marking bad channels and trials.';
disp(msg);
outputLog{1}=msg;

for iChunk = startChunk:endChunk
    if exist('EPchunk','var') && ~isempty(EPchunk)
        dataChunk=EPchunk{iChunk};
    else
        ep_tictoc('ioStart');
        eval(['load ''' deblank(inFile) '-' num2str(iChunk) '.mat''']);
        ep_tictoc('ioFinish');
    end;
    if strcmp(dataChunk.dataType,'continuous')
        numTrials=floor(size(dataChunk.data,2)/ceil(dataChunk.Fs)); %excess time points are tacked onto final epoch
        numSamples = min(ceil(dataChunk.Fs),size(dataChunk.data,2)); %one second epochs
    else
        numTrials = length(dataChunk.cellNames);
        numSamples=length(dataChunk.timeNames);
    end;
    numChans=length(dataChunk.chanNames);
    spike=[1000 zeros(1,numSamples-1)];
        
    if strcmp(dataChunk.dataType,'average')
        badTrialNum=(dataChunk.avgNum==-1);
        badChanNum=isnan(dataChunk.analysis.badChans);
    else
        badTrialNum=dataChunk.analysis.badTrials;
        badChanNum=(dataChunk.analysis.badChans==-1);
    end;
    
    badChanNum(:,:,badChans)=1; %set global bad channels to all be replaced.
    
    for iTrial=1:numTrials
        if badTrialNum(iTrial)
            if strcmp(dataChunk.dataType,'continuous')
                if iTrial == numTrials %excess time points are tacked onto final epoch
                    dataChunk.data(:,(iTrial-1)*numSamples+1:end,1,theSubject)=repmat([1000 zeros(1,size(dataChunk.data,2)-(numSamples*(numTrials-1))-1)],numChans,1);
                else
                    dataChunk.data(:,(iTrial-1)*numSamples+1:iTrial*numSamples,1,theSubject)=repmat(spike,numChans,1);
                end;
            else
                dataChunk.data(:,:,iTrial,theSubject)=repmat(spike,numChans,1);
            end;
        else
            if strcmp(dataChunk.dataType,'continuous')
                if iTrial == numTrials %excess time points are tacked onto final epoch
                    theData=dataChunk.data(:,(iTrial-1)*numSamples+1:end,1,theSubject);
                else
                    theData=dataChunk.data(:,(iTrial-1)*numSamples+1:iTrial*numSamples,1,theSubject);
                end;
            else
                theData=dataChunk.data(:,:,iTrial,theSubject);
            end;
            trialBadChans=find(badChanNum(1,iTrial,:));
            trialGoodChans=setdiff([1:numChans],trialBadChans);
            theData(trialBadChans,:)=repmat(spike,length(trialBadChans),1);
            if strcmp(dataChunk.dataType,'continuous')
                if iTrial == numTrials %excess time points are tacked onto final epoch
                    dataChunk.data(:,(iTrial-1)*numSamples+1:end,1,theSubject)=theData;
                else
                    dataChunk.data(:,(iTrial-1)*numSamples+1:iTrial*numSamples,1,theSubject)=theData;
                end;
            else
                dataChunk.data(:,:,iTrial,theSubject)=theData;
            end;
        end;
    end;
    if exist('EPchunk','var') && ~isempty(EPchunk)
        EPchunk{iChunk}=dataChunk;
    else
        ep_tictoc('ioStart');
        eval (['save ''' inFile '-' num2str(iChunk) '.mat'' dataChunk']);
        ep_tictoc('ioFinish');
    end;
end;
