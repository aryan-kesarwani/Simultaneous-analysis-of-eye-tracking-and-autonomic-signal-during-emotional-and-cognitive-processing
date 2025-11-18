function [outputLog, graphCounter] = ep_fixMains(inFile, startChunk, endChunk, theSubject, mainsFreq, butterflyFig, graphCounter, numGraphs)
% [outputLog, graphCounter] = ep_fixMains(inFile, startChunk, endChunk, theSubject, mainsFreq, butterflyFig, graphCounter, numGraphs)
%	Performs fix on mains noise.
%
%Inputs
%	inFile:     filename (not including the .mat suffix or chunk number.  e.g., "NT5") and sourcepath.
%	startChunk: starting chunk (usually 1)
%   endChunk:   ending chunk
%   theSubject: which subject of the file is being processed.
%   mainsFreq:  the frequency for the mains noise (typically 50 or 60).
%   butterflyFig:  the handle for the output figure.  Otherwise, will open a new figure.
%   graphCounter: the current subplot for the summary figure.
%   numGraphs: the total number of subgraphs in the summary figure.
%
%   The input chunks include: dataChunk in EP data format.
%
%Outputs
%   outputLog: output messages from bad data process.
%   graphCounter: the current subplot for the summary figure.
%
%   Updated output chunks: dataChunk in EP data format.
%
% History:
%
% by Joseph Dien (5/9/20)
% jdien07@mac.com
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

global EPchunk EPtictoc

msg=['Fixing mains noise at ' num2str(mainsFreq) ' Hz.'];
disp(msg);
outputLog{1}=msg;

if ~exist('butterflyFig','var')
    butterflyFig=figure('Name','Artifact Correction','NumberTitle','off');
    colormap jet;
    standAlone=1;
else
    standAlone=0;
end

for iChunk = startChunk:endChunk
    if exist('EPchunk','var') && ~isempty(EPchunk)
        dataChunk=EPchunk{iChunk};
    else
        ep_tictoc('ioStart');
        eval(['load ''' deblank(inFile) '-' num2str(iChunk) '.mat''']);
        ep_tictoc('ioFinish');
    end
    numChans=length(dataChunk.chanNames);
    numPoints=length(dataChunk.timeNames);
    EEGchans=find(strcmp('EEG',dataChunk.chanTypes));
    lineNoiseIn.fPassBand=[0, dataChunk.Fs/2];
    lineNoiseIn.Fs=dataChunk.Fs;
    lineNoiseIn.fScanBandWidth=2;
    lineNoiseIn.lineFrequencies=[mainsFreq:mainsFreq:floor(dataChunk.Fs/2)]; %include harmonics up to the Nyquist.
    lineNoiseIn.lineNoiseChannels=EEGchans';
    lineNoiseIn.maximumIterations=10;
    lineNoiseIn.p = 0.01;
    lineNoiseIn.pad=0;
    lineNoiseIn.taperBandWidth=2;
    lineNoiseIn.tau=100;
    
    originalData=dataChunk.data;
    
    for iCell=1:length(dataChunk.cellNames)
        ep_tictoc;if EPtictoc.stop;return;end
        if strcmp(dataChunk.dataType,'continuous')
            [ALLEEG]=ep_ep2alleeg(ep_selectData(dataChunk,{EEGchans,[],[],[],[],[]}));
            SlidingWinLength=round(min(numPoints,dataChunk.Fs*3)/dataChunk.Fs);
            SlidingWinStep=round(SlidingWinLength/2);
        else
            [ALLEEG]=ep_ep2alleeg(ep_selectData(dataChunk,{EEGchans,[],iCell,theSubject,[],[]}));
            SlidingWinLength=round(numPoints/dataChunk.Fs);
            SlidingWinStep=SlidingWinLength;
        end
        lineNoiseIn.taperWindowSize=SlidingWinLength;
        lineNoiseIn.taperWindowStep=SlidingWinStep;
        [ALLEEG(1), ~] = cleanLineNoise(ALLEEG(1), lineNoiseIn);
        dataChunk.data(EEGchans,:,iCell,theSubject,:,:,:)=ALLEEG(1).data;
    end
    
    if exist('EPchunk','var') && ~isempty(EPchunk)
        EPchunk{iChunk}=dataChunk;
    else
        ep_tictoc('ioStart');
        eval (['save ''' inFile '-' num2str(iChunk) '.mat'' dataChunk']);
        ep_tictoc('ioFinish');
    end
    if ~isempty(butterflyFig) && ishandle(butterflyFig{iChunk})
        trialdata=reshape(dataChunk.data(:,:,:,theSubject),numChans,[]);
        originalData=reshape(originalData(:,:,:,theSubject),numChans,[]);
        displayPeriod=size(trialdata,2);    %Number of timepoints to graph in display.
        totalDisplayPeriod=displayPeriod*size(dataChunk.data,4);
        decimateSamples=ceil(max(1,totalDisplayPeriod/10000));
        figure(butterflyFig{iChunk});
        
        theTitle='Mains Noise';
        plotData=ep_makePlotData(butterflyFig{iChunk},displayPeriod,totalDisplayPeriod,decimateSamples,theTitle,originalData-trialdata,EEGchans,theSubject);
        subplot(numGraphs,1,graphCounter), plot([1:decimateSamples:totalDisplayPeriod],plotData);
        title(theTitle,'Interpreter','none');
        axis([1 totalDisplayPeriod -200 200])
        set(gca,'XTickLabel','','XTick',[]);
        
        theTitle='Without Mains Noise';
        plotData=ep_makePlotData(butterflyFig{iChunk},displayPeriod,totalDisplayPeriod,decimateSamples,theTitle,trialdata,EEGchans,theSubject);
        subplot(numGraphs,1,graphCounter+1), plot([1:decimateSamples:totalDisplayPeriod],plotData);
        title(theTitle,'Interpreter','none');
        axis([1 totalDisplayPeriod -200 200])
        set(gca,'XTickLabel','','XTick',[]);
        
        drawnow
    end
    if standAlone
        try
            MATLABver=ver('MATLAB');
            [a b]=strtok(MATLABver.Version,'.');
            b=b(2:end);
            if ~isprop(butterflyFig,'Number')
                eval (['print -f' num2str(butterflyFig{iChunk}) ' -djpeg ''' inFile '''-' num2str(iChunk) 'mains.jpg']);
            else
                eval (['print -f' num2str(butterflyFig{iChunk}.Number) ' -djpeg ''' inFile '''-' num2str(iChunk) 'mains.jpg']);
            end
        catch
            disp('Couldn''t save a copy of the artifact correction figure.  Perhaps your version of Matlab is not current.');
        end
        close(butterflyFig{iChunk});
    end
end
if ~isempty(butterflyFig)
    graphCounter=graphCounter+2;
end