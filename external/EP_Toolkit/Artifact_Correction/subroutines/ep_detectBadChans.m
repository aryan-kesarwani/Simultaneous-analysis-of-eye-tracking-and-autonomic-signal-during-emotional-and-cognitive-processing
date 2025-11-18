function [badChans shortChans outputLog, graphCounter]=ep_detectBadChans(inFile, startChunk, endChunk, badDataCriteria, butterflyFig, graphCounter, numGraphs, theSubject, chunkSizes);
%  [badChans shortChans outputLog, graphCounter]=ep_detectBadChans(inFile, startChunk, endChunk, badDataCriteria, butterflyFig, graphCounter, numGraphs, theSubject, chunkSizes);
%       Detects global bad channels by identifying ones that cannot be readily predicted by the neighboring channels and by detecting flat data channels.
%       Also notes shorted channels.  Bases routine based on the full dataset.  If chunked, will first reassemble it.
%       After identifying bad channels, for every time point where there at least one good channel is outside the "saturation" range,
%       all the channels for that timepoint are set to NaN.
%
% Cohen, J., & Cohen, P. (1983). Applied multiple regression/correlation analysis for the behavioral sciences. Hillsdale, NJ: Lawrence Erlbaum Associates.
%
%Inputs:
%	inFile:     filename (not including the .mat suffix or chunk number.  e.g., "NT5") and sourcepath.
%               or else EPdata file structure.
%	startChunk: starting chunk (usually 1)
%   endChunk:   ending chunk
%   badDataCriteria:  Criteria for detecting bad data.
%       .neighbors: number of electrodes considered to be neighbors
%       .badchan:   minimum predictability (multiple R) from neighbors to not be considered globally bad
%       .badtrials: percentage of good trials chan is bad to declare a channel globally bad
%       .saturation: followed by range of acceptable data values.  Time points with a channel outside this range will be excluded.
%   butterflyFig:  the handle for the output figure.  Otherwise, will open a new figure.
%   graphCounter: the current subplot for the summary figure.
%   numGraphs: the total number of subgraphs in the summary figure.
%   theSubject: which subject of the file is being processed.
%   chunkSizes: array of lengths of the chunks.
%
%Outputs:
%  badChans : List of bad channels.
%  shortChans: List of shorted channels.
%  outputLog: output messages from bad channel detection process
%  graphCounter: the current subplot for the summary figure.

%History:
%  by Joseph Dien (2/8/09)
%  jdien07@mac.com
%
% modified 3/14/09 JD
% Changed to use EP format data to provide more flexibility with I/O functions.
%
% modified 5/15/09 JD
% Treats flat channels as bad data and excludes from bad channel detection routine.
% Dropped eloc as separate input parameter (now part of data).
% Flat channel not bad if it is the reference channel.
%
% modified 9/4/09 JD
% Added support for multiple refChans to deal with mean mastoid data where the presence of the two reference channels (correlated -1)
% was causing ICA problems.
%
% modified 10/28/09 JD
% Added detection of channels perfectly correlated with a reference channel and which were therefore flat prior to rereferencing.
%
% modified 11/12/09 JD
% Correlated channels can be only nearly perfect (e.g., .9999) and still trigger bad channel code, to account for rounding errors etc.
%
% bugfix 11/20/09 JD
% Replaced "union" commands with "unique" commands because certain situations caused the "union" command to crash in
% Matlab 2007.
%
% modified & bugfix 12/3/09 JD
% Detects non-reference channels that are perfectly correlated and identifies them as bad channels as they must be
% shorted together.  Fixed bug where test of correlation with reference only detecting +1 correlation, not -1 correlation.
% Fixed bug where if there is an explicit reference channel and it is flat, then all reference channels marked bad and
% real bad channels are no longer marked bad.
% Don't apply correlated neighbors test to the reference channels as distant reference channels will always be labeled
% bad.
%
% modified & bugfix 2/24/10 JD
% Now works on average files.
% Fixed bug where neighboring channels for determing whether a channel is not correlating with its neighbors sometimes not chosen
% correctly, which could lead to too many channels being dubbed globally bad.
% No longer treating shorted channels as being bad (too conservative).  Instead just displaying a warning message.
% If there are two reference channels (as in mean mastoids), then no longer require that they have a -1 correlation as one may just be bad.
% If there are two reference channels (as in mean mastoids), then they are still marked as bad channels if they are
% flat.
% Added log output.
% When there are shorted channels, prints out the channel pairs.
%
% bugfix 6/7/11 JD
% Fixed crash when bad channels detection in artifact correction routine generated error message.
%
% modified 1/25/12 JD
% Eliminated REF channel type and added reference field to keep track of original and current reference scheme.
%
% modified 9/22/13 JD
% Restricted bad channel detection to EEG channels.
%
% modified 5/23/16 JD
% Detects channels that went flat partway through a session (more than 10%) and labels them as global bad channels.
%
% modified 4/26/17 JD
% Excludes time points beyond a certain range from global bad channel detection, blink, and saccade routines.
% Divides the data into sections according to the badtrials parameter and checks the badchan multiple R predictability threshold for each section to detect channels that went bad partway through the session.
%
% bugfix 9/30/17 JD
% Median corrects data prior to saturation check to ensure channels with merely high offsets are not treated as bad data.
% Median based on the section rather than the entire dataset.
%
% bugfix 11/26/17 JD
% Fixed calculation of section median bad channel.
%
% bugfix 1/13/19 JD
% Fixed not taking previously edited bad trials into account when identifying channels bad due to being flat.
%
% Bugfix & modified 7/23/19 JD
% Added summary subplots for global bad channels.
% Global bad channel detection now based on baselined/detrended copy of data.
% Fixed crash when chunking single-trial or averaged data.
% Save change out-of-range points to NaN and save the EPchunk for the next step.
%
% bugfix 9/21/20 JD
% Fixed crash when chunking single-trial data.
%
% bugfix 10/10/20 JD
% Fixed crash when chunking continuous data.
%
% bugfix 12/23/20 JD
% Fixed crash when preprocessing average data.
%
% bugfix 9/3/22 JD
% Fixed crash when preprocessing data with NaN values.
%
% modified 9/25/24 JD
% Tweaked bad channel detection algorithm so less sensitive to solitary bad channels with high amplitude.
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

global EPchunk;

badChans{1}=-1;
outputLog={};

if ~exist('theSubject','var')
    theSubject =1;
end

if isstruct(inFile)
    EPdata=inFile;
else
    %reassemble chunked datasets to apply detection routines to the full dataset since this step isn't memory-intensive.
    for iChunk = startChunk:endChunk
        disp([deblank(inFile) '-' num2str(iChunk)]);
        if exist('EPchunk','var') && ~isempty(EPchunk)
            dataChunk=EPchunk{iChunk};
        else
            ep_tictoc('ioStart');
            eval(['load ''' deblank(inFile) '-' num2str(iChunk) '.mat''']);
            ep_tictoc('ioFinish');
        end
        if strcmp('continuous',dataChunk.dataType)
            newNumPoints=length(dataChunk.timeNames);
            numPointsList(iChunk)=newNumPoints;
            numChunkTrials=floor(size(dataChunk.data,2)/ceil(dataChunk.Fs)); %for last chunk, extra points will be attached to final segment
            numTrialsList(iChunk)=numChunkTrials;
            if iChunk==1
                EPdata=dataChunk;
            else
                EPdata.data(:,end+1:end+newNumPoints,:,:,:,:)=dataChunk.data;
                EPdata.timeNames(end+1:end+newNumPoints)=dataChunk.timeNames;
            end
        else
            numChunkTrials=length(dataChunk.cellNames);
            numTrialsList(iChunk)=numChunkTrials;
            if iChunk==1
                EPdata=dataChunk;
            else
                EPdata.data(:,:,end+1:end+numChunkTrials,:,:,:)=dataChunk.data;
            end
        end
        EPdata.analysis.badTrials(:,1+sum(numTrialsList(1:iChunk-1)):sum(numTrialsList(1:iChunk)))=dataChunk.analysis.badTrials;
    end
end

if strcmp(EPdata.dataType,'factors')
    msg='This function does not support factor files.';
    disp(msg);
    outputLog{end+1}=msg;
    badChans=-1;
    return;
end

elecDistances=ep_closestChans(EPdata.eloc);

numSubs=length(EPdata.subNames);
shortChans=[];
numChans=length(EPdata.chanNames);
if length(EPdata.eloc) ~= numChans
    msg='Error: The number of channels in the electrode coordinates file is different from that of the data.';
    disp(msg);
    outputLog{end+1}=msg;
    badChans=-1;
    return;
end

EEGchans=find(strcmp('EEG',EPdata.chanTypes));
testData=EPdata.data(EEGchans,:,:,theSubject);
testData=reshape(testData,length(EEGchans),[]);
numPoints=size(testData,2);
%goodPoints = find((max(testData-repmat(median(testData','omitnan')',1,size(testData,2))) < badDataCriteria.saturation(2)) & (min(testData-repmat(median(testData','omitnan')',1,size(testData,2))) > badDataCriteria.saturation(1)));

%exclude bad points that were subject to a major movement artifact or whatever if it affected more than 10% of the good channels
bcTestData=testData-repmat(median(testData,2,'omitnan'),1,numPoints);
goodPoints = find(~(sum((bcTestData<badDataCriteria.saturation(1))|(bcTestData> badDataCriteria.saturation(2)),1)>(numChans/10)));
testData=testData(:,goodPoints);
allData=reshape(EPdata.data(EEGchans,:,:,theSubject),length(EEGchans),[]); %includes "bad" points

%flat channels are bad.
badChans=find(~std(testData,1,'omitnan')); 
badChans=intersect(badChans,EEGchans);
goodChans=setdiff(EEGchans,badChans);

normData=pinv(diag(std(testData,0,2,'omitnan')))*testData; %standardized data
normData=normData-diag(mean(normData,2))*ones(size(normData)); %normalized data (standardized and mean-corrected)
elecDistances(badChans,:)=inf;
elecDistances(:,badChans)=inf;
origRefChan=EPdata.reference.original;
currRefChan=EPdata.reference.current;

if ~exist('butterflyFig','var') && ~isstruct(inFile)
    butterflyFig=figure('Name','Global Bad Channel Detection','NumberTitle','off');
    colormap jet;
    standAlone=1;
else
    standAlone=0;
end

%global bad channel is flat for more than 10% of the segments
if strcmp(EPdata.dataType,'continuous')
    numGoodTrials=sum(~EPdata.analysis.badTrials); %ignore excess
    flatTrialsChans=zeros(length(goodChans),numGoodTrials);
    trialList=find(~EPdata.analysis.badTrials);
    for iTrial = 1:length(trialList)
        theTrial=trialList(iTrial);
        flatTrialsChans(:,iTrial)=var(squeeze(EPdata.data(goodChans,(theTrial-1)*EPdata.Fs+1:theTrial*EPdata.Fs,1,theSubject))');
    end
else
    numGoodTrials=sum(~EPdata.analysis.badTrials(theSubject,:));
    flatTrialsChans=zeros(length(goodChans),numGoodTrials);
    trialList=find(~EPdata.analysis.badTrials(theSubject,:));
    for iTrial = 1:length(trialList)
        theTrial=trialList(iTrial);
        flatTrialsChans(:,iTrial)=var(squeeze(EPdata.data(goodChans,:,theTrial,theSubject))');
    end
end
badChans=[badChans; goodChans(find(sum(flatTrialsChans'==0)>(numGoodTrials/10)))];
goodChans=setdiff(goodChans,badChans);

if isempty(goodChans)
    msg='Error: No good EEG channels left.';
    disp(msg);
    outputLog{end+1}=msg;
    badChans=-1;
    return;
end

corrs=corrcoef([testData']);
corrSigns=sign(corrs);
corrs=corrSigns.*ceil(abs(corrs)*1000)/1000;

%look for bad channels that are perfectly correlated with each other and must therefore be shorted together.

goodNonRefChans=setdiff(goodChans,currRefChan);

shortPairs=[];
for chan1= 1:length(goodNonRefChans)
    theChan1= goodNonRefChans (chan1);
    for chan2=chan1+1:length(goodNonRefChans)
        theChan2= goodNonRefChans (chan2);
        if abs(corrs(theChan1, theChan2))== 1
            shortChans=unique([shortChans theChan1, theChan2]);
            shortPairs=[shortPairs EPdata.chanNames{theChan1} '-' EPdata.chanNames{theChan2} '; '];
        end
    end
end
if ~isempty(shortPairs)
    msg=['Warning: shorted channels: ' shortPairs];
    disp(msg);
    outputLog{end+1}=msg;
end

%look for bad channels that are not correlated with other channels
numRegressors=min(badDataCriteria.neighbors,length(goodChans)-1);
neighbors=zeros(length(EEGchans),numRegressors);
for iChan=1:length(goodChans)
    chan=goodChans(iChan);
    [E IX]=sort(elecDistances(chan,goodChans));
    neighbors(chan,:)=goodChans(IX(2:numRegressors+1));
end

secLength=floor((badDataCriteria.badtrials/100)*length(goodPoints));
secAllLength=floor((badDataCriteria.badtrials/100)*size(EPdata.data,2));
if secLength >= EPdata.Fs %if section is at least a second long
    numSecs=floor(100/badDataCriteria.badtrials);
else
    numSecs=1;
    secLength=length(goodPoints);
    secAllLength=size(EPdata.data,2);
end

chanPred=zeros(length(EEGchans),numSecs);
badPred=zeros(length(EEGchans),1);
for iSec=1:numSecs
    secPoints=((iSec-1)*secLength)+1:iSec*secLength;
    secAllPoints=((iSec-1)*secAllLength)+1:iSec*secAllLength;
    for iChan=1:length(goodChans)
        theChan=goodChans(iChan);
        Y=normData(theChan,secPoints)';
        X=normData(neighbors(theChan,:),secPoints)';
        R=corrcoef([Y X]);
        R=R(2:end,1);
        B=[ones(secLength,1) X]\Y;
        chanPred(theChan,iSec)=sqrt(sum(B(2:end).*R)); %Multiple R Cohen & Cohen (1983), p. 86
        if chanPred(theChan,iSec) < badDataCriteria.badchan
            badPred(theChan)=1;
        end
        if (length(find(abs(allData(theChan,secAllPoints)-median(allData(theChan,secAllPoints)))>500))/secLength) > (badDataCriteria.badtrials/100)
            badPred(theChan)=1;
        end
    end
end

nonRefGoodChans=setdiff(goodChans,currRefChan);
badChans=[badChans; nonRefGoodChans(find(badPred(nonRefGoodChans)))]; %don't check ref chans for being bad by local channel predictability

if length(currRefChan) == 1
    if std(normData(currRefChan,:)','omitnan')==0
        badChans=setdiff(badChans,currRefChan); %flat channel is not bad if it is the reference channel
    end
end

goodChans=setdiff(goodChans,badChans);

if isempty(goodChans)
    msg='Error: No good EEG channels left.';
    disp(msg);
    outputLog{end+1}=msg;
    badChans=-1;
    return;
end

%look for bad channels that are perfectly correlated with the reference channel(s) and were therefore flat prior to
%rereferencing.

if length(origRefChan) > 2
    msg='There are more than two original reference channels indicated.';
    disp(msg);
    outputLog{end+1}=msg;
    badChans=-1;
    return;
elseif length(origRefChan) == 2
    refCorrs=abs(corrs(origRefChan(1),:))';
    badChans=unique([badChans; setdiff(find(refCorrs == 1),origRefChan)]);
    badChans=intersect(badChans,EEGchans);
elseif length(origRefChan) == 1
    refCorrs=abs(corrs(origRefChan(1),:))';
    badChans=unique([badChans; setdiff(find(refCorrs == 1),origRefChan)]);
    badChans=intersect(badChans,EEGchans);
end

if ~strcmp(EPdata.reference.type,'AVG')
    %If there are two current reference channels and they are not perfectly inversely correlated, there is a problem
    if length(origRefChan) == 2
        if (corrs(currRefChan(1),currRefChan(2)) ~= -1) && isempty(intersect(currRefChan,badChans))
            msg=['The two current reference channels should have a perfect inverse correlation and do not (' num2str(corrs(currRefChan(1),currRefChan(2))) ') so something is wrong.'];
            disp(msg);
            outputLog{end+1}=msg;
            badChans=-1;
            return;
        end
    end
end

goodChans=setdiff(goodChans,badChans);

if isempty(goodChans)
    msg='Error: No good channels left.';
    disp(msg);
    outputLog{end+1}=msg;
    badChans=-1;
    return;
end

%Set points to NaN if they are outside the saturation preferences range for any of the good channels.
% msg='Identifying data points outside the acceptable voltage range and setting them to NaN.';
% disp(msg);
% outputLog{end+1}=msg;

%temporarily rereference data to Cz so that this operation is independent of the data reference.
[czChan, theOrder] = ep_findChan(EPdata.eloc, EPdata.implicit, EPdata.chanNames, EPdata.ced, badChans, 'Cz', EPdata.montage);

% if theOrder > 1
%     msg=['Cz is a bad channel so instead using channel ' EPdata.eloc(czChan).labels '.'];
% elseif strcmpi(EPdata.eloc(czChan).labels,'cz')
%     msg=['Channel Cz identified.'];
% else
%     msg=['Channel ' EPdata.eloc(czChan).labels ' is assumed to be Cz.'];
% end
% outputLog{end+1}=msg;
% disp(msg);
% 
allData=allData-repmat(allData(EEGchans(EEGchans==czChan),:),length(EEGchans),1); %rereference to Cz
badPoints = find((max(allData(goodChans,:)-repmat(median(allData(goodChans,:)','omitnan')',1,size(allData,2))) >= badDataCriteria.saturation(2)) | (min(allData(goodChans,:)-repmat(median(allData(goodChans,:)','omitnan')',1,size(allData,2))) <= badDataCriteria.saturation(1)));
fixData=reshape(EPdata.data(EEGchans,:,:,theSubject),length(EEGchans),[]); %includes "bad" points
fixData(:,badPoints)=NaN;
EPdata.data(EEGchans,:,:,theSubject)=reshape(fixData,length(EEGchans),size(EPdata.data,2),size(EPdata.data,3),1);
pointsDropped=length(badPoints);
if pointsDropped > 0
    msg=[num2str(pointsDropped) ' time points (' num2str(pointsDropped/size(allData,2)*100) '%) were set to NaN as being out-of-range, which has been set at ' num2str(badDataCriteria.saturation(1)) ' and ' num2str(badDataCriteria.saturation(2)) '.'];
else
    msg='No out-of-range time points set to NaN.';
end
outputLog{end+1}=msg;
disp(msg);

if strcmp(EPdata.dataType,'continuous')
    numTrials=1;
    trialSize = length(EPdata.timeNames);
else
    trialSize = length(EPdata.timeNames);
    numTrials = length(EPdata.cellNames);
end

if isempty(EPdata) || isempty(EPdata.data)
    disp('Warning: No file saved due to program error.');
else
    for iChunk = startChunk:endChunk
        if length(startChunk:endChunk)==1 %if no chunking
            dataChunk=EPdata;
        else
            if exist('EPchunk','var') && ~isempty(EPchunk)
                dataChunk=EPchunk{iChunk};
            else
                ep_tictoc('ioStart');
                eval(['load ''' deblank(inFile) '-' num2str(iChunk) '.mat''']);
                ep_tictoc('ioFinish');
            end
            if strcmp('continuous',dataChunk.dataType)
                dataChunk.data=EPdata.data(:,1+sum(numPointsList(1:iChunk-1)):sum(numPointsList(1:iChunk)),:,:,:,:); %don't need to extract rest of dataChunk for continuous data as not changed.
                dataChunk.timeNames=EPdata.timeNames(1+sum(numPointsList(1:iChunk-1)):sum(numPointsList(1:iChunk)));
            else
                dataChunk.data=EPdata.data(:,:,1+sum(numTrialsList(1:iChunk-1)):sum(numTrialsList(1:iChunk)),:,:,:);
            end
            dataChunk.analysis.badTrials=EPdata.analysis.badTrials(:,1+sum(numTrialsList(1:iChunk-1)):sum(numTrialsList(1:iChunk)));
        end
        if exist('EPchunk','var') && ~isempty(EPchunk)
            EPchunk{iChunk}=dataChunk;
        else
            ep_tictoc('ioStart');
            eval (['save ''' outFile '-' num2str(iChunk) '.mat'' dataChunk']);
            ep_tictoc('ioFinish');
        end
    end
end

if ~isstruct(inFile)
    for iChunk=1:length(chunkSizes)
        totalDisplayPeriod=chunkSizes(iChunk); %Number of timepoints to graph in display.
        decimateSamples=ceil(max(1,totalDisplayPeriod/10000));
        displayPeriod=totalDisplayPeriod/size(EPdata.data,4);
        
        if ~isempty(badChans)
            if ~isempty(butterflyFig)
                if standAlone
                    figure(butterflyFig{iChunk});
                    subplot(3,1,2), plot([1:decimateSamples:displayPeriod],reshape(allData(badChans,sum(chunkSizes(1:iChunk-1))+1:decimateSamples:sum(chunkSizes(1:iChunk-1))+displayPeriod),length(badChans),[]));
                    title('global bad channels','Interpreter','none');
                    axis([1 totalDisplayPeriod -200 200])
                    set(gca,'XTickLabel','','XTick',[]);
                    
                    subplot(3,1,3), plot([1:decimateSamples:displayPeriod],reshape(allData(goodChans,sum(chunkSizes(1:iChunk-1))+1:decimateSamples:sum(chunkSizes(1:iChunk-1))+displayPeriod),length(goodChans),[]));
                    title('without global bad channels','Interpreter','none');
                    axis([1 totalDisplayPeriod -200 200])
                    set(gca,'XTickLabel','','XTick',[]);
                else
                    %workaround for weird Matlab 2017b causing it to not recognize the reshape command.
                    temp=zeros(size(EPdata.data,1),size(EPdata.data,2)*size(EPdata.data,3));
                    for iTrial=1:numTrials
                        temp(:,(iTrial-1)*trialSize+1:iTrial*trialSize)=EPdata.data(:,:,iTrial);
                    end
                    figure(butterflyFig{iChunk});
                    theTitle='global bad channels';
                    plotData=ep_makePlotData(butterflyFig{iChunk},displayPeriod,totalDisplayPeriod,decimateSamples,theTitle,temp,badChans,theSubject);
                    subplot(numGraphs,1,graphCounter), plot([1:decimateSamples:totalDisplayPeriod],plotData);
                    title(theTitle,'Interpreter','none');
                    axis([1 totalDisplayPeriod -200 200])
                    set(gca,'XTickLabel','','XTick',[]);
                    
                    %workaround for weird Matlab 2017b causing it to not recognize the reshape command.
                    theTitle='without global bad channels';
                    plotData=ep_makePlotData(butterflyFig{iChunk},displayPeriod,totalDisplayPeriod,decimateSamples,theTitle,temp,goodChans,theSubject);
                    subplot(numGraphs,1,graphCounter+1), plot([1:decimateSamples:totalDisplayPeriod],plotData);
                    title(theTitle,'Interpreter','none');
                    axis([1 totalDisplayPeriod -200 200])
                    set(gca,'XTickLabel','','XTick',[]);
                end
            end
        else
            if ~isempty(butterflyFig)
                if standAlone
                    figure(butterflyFig{iChunk});
                    subplot(3,1,2), plot([1:decimateSamples:displayPeriod],ones(1,ceil(displayPeriod/decimateSamples)));
                    axis([1 totalDisplayPeriod -200 200])
                    set(gca,'XTickLabel','','XTick',[]);
                    subplot(3,1,3), plot([1:decimateSamples:displayPeriod],allData(EEGchans,1:displayPeriod));
                    title('no global bad channels','Interpreter','none');
                    axis([1 totalDisplayPeriod -200 200])
                    set(gca,'XTickLabel','','XTick',[]);
                else
                    figure(butterflyFig{iChunk});
                    if numSubs > 1
                        theTitle='global bad channels';
                    else
                        theTitle='no global bad channels';
                    end
                    plotData=ep_makePlotData(butterflyFig{iChunk},displayPeriod,totalDisplayPeriod,decimateSamples,theTitle,zeros(size(allData,1),size(allData,2)),EEGchans,theSubject);
                    subplot(numGraphs,1,graphCounter), plot([1:decimateSamples:totalDisplayPeriod],plotData);
                    title(theTitle,'Interpreter','none');
                    axis([1 totalDisplayPeriod -200 200])
                    set(gca,'XTickLabel','','XTick',[]);
                    
                    if numSubs > 1
                        theTitle='without global bad channels';
                    else
                        theTitle='with all channels';
                    end
                    temp=zeros(size(EPdata.data,1),size(EPdata.data,2)*size(EPdata.data,3));
                    for iTrial=1:numTrials
                        temp(:,(iTrial-1)*trialSize+1:iTrial*trialSize)=EPdata.data(:,:,iTrial);
                    end
                    plotData=ep_makePlotData(butterflyFig{iChunk},displayPeriod,totalDisplayPeriod,decimateSamples,theTitle,temp,EEGchans,theSubject);
                    subplot(numGraphs,1,graphCounter+1), plot([1:decimateSamples:totalDisplayPeriod],plotData);
                    title(theTitle,'Interpreter','none');
                    axis([1 totalDisplayPeriod -200 200])
                    set(gca,'XTickLabel','','XTick',[]);
                end
            end
        end
        
        if standAlone
            try
                MATLABver=ver('MATLAB');
                [a b]=strtok(MATLABver.Version,'.');
                b=b(2:end);
                if ~isprop(butterflyFig,'Number')
                    eval (['print -f' num2str(butterflyFig{iChunk}) ' -djpeg ''' inFile '''-' num2str(iChunk) 'globalBadChannels.jpg']);
                else
                    eval (['print -f' num2str(butterflyFig{iChunk}.Number) ' -djpeg ''' inFile '''-' num2str(iChunk) 'globalBadChannels.jpg']);
                end
            catch
                disp('Couldn''t save a copy of the alpha correction figure.  Perhaps your version of Matlab is not current.');
            end
        end
    end
    
    if standAlone
        close(butterflyFig);
    end
    if ~isempty(butterflyFig)
        graphCounter=graphCounter+2;
    end
end
