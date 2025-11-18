function cellNums=ep_segmentData(cellTable,inputFiles,importFormat,outputFormat,preview,specData,segSuffix)
% cellNums=ep_segmentData(cellTable,inputFiles,preview,specData,segSuffix) -
% segments continuous data files into single-trial segmented data.
%
%Input:
%  cellTable        : Cell array with segmenting parameters (name, prestim ms, poststim ms, delay, list of five specs with
%                       event spec (type, value, or key), relation, and value.  Empty when not used.
%  inputFiles       : List of session files including path.  If preview mode, then the EPdataset structured variable.
%  importFormat     : file format code for input data
%  outputFormat     : file format code for output data
%  preview          : indicates files should only be assessed for number of resulting epochs without
%                       generating output files.
%  specData         : Parameters for separate trials specs text file.  Empty if none to be used.
%     .specFormat   : format of the text file (EPM=E-Prime _EPM.txt, TSP=tab-delimited with first row containing the field names _TSP.txt)
%     .excludeSpec  : The spec to be used for excluding unwanted trials.  Empty if not to be used.
%     .excludeValue : The value of the spec to be excluded.
%     .excludeSpec2 : The spec to be used for excluding unwanted trials.  Empty if not to be used.
%     .excludeValue2: The value of the spec to be excluded.
%     .excludeSpec3 : The spec to be used for excluding unwanted trials.  Empty if not to be used.
%     .excludeValue3: The value of the spec to be excluded.
%     .specField    : cell array of specs to be included in the segment's trial specs.  Empty if not to be used.
%     .specLabels   : cell array of the labels to use for the trial specs instead of the original field name.
%     .specTable    : Cell array of the template spec data (trials,spec fields)
%     .specNames    : Cell array of the field names for the template spec data for use with the preview mode.
%  segSuffix        : Suffix for the output segment files.
%Output:
%   cellNums       : Array of number of epochs resulting for each cell (subject,cell).
%
% Will include specs from both TRSPs and spec files.  If both have crits with the same name, will retain TS- prefix for the TRSP field.
% In the case of subject specs where both the file and the spec file have fields with the same name, the spec file info will overwrite the existing info.
%
% When -follows- is specified, it will first be assumed that the given parameter is the .type field, usually with its .value field needing to meet the next spec in order to be valid.
% Intervening events will then be ignored.  If no such .type field exists, then the given parameter will be assumed to be the .value field and that the immediately preceding event is being referred to.
% For ~=, in the former case it will look for the preceding .type field and then compare the .value field.  If there is no spec for the .value field, it is automatically not valid.
% In the latter case, it will compare against the .value field of the prior event.  In this case, there does not need to be a further spec.
% Ditto for < and > even though in principle one could legitimately wish to specify a .type that is not present in the dataset.
% Ditto "starts" and "ends" and "contains".
%
% Non-integer sample sizes continue to be a problem.  The approach here will be to strictly adhere to the boundaries provided (i.e., floor and ceil) rather than to round.
%
% Specification of ACC and RT data can be complicated since there are many variations in datasets.  The EP Toolkit assumes that trial specs named ACC and RT refer to these two types of data.
% If the specs tab is being used, columns in specs files that are identified as being ACC or RT data will be renamed accordingly.
% If the -responseCorr- and -responseErr- options are used, then the ACC and RT values will be generated and inserted into the ACC and RT trial specs (and they will be created if not already existing).
% If the ACC and RT fields are not blank, then the trial spec they specify will have their contents copied to the ACC and RT trial specs, overwriting the original contents if they already exist
% or creating them first if they do not already exist.  If there is no trial spec that matches the contents of these fields, 
% then they will be converted to numbers and directly inserted into the ACC and/or RT fields.  If they are not numbers, then they will end up being NaN values.
% This allows a criteria to flexibly define ACC values.  0 is error and 1 is correct and 2 is timeout.

%History
%  by Joseph Dien (11/16/13)
%  jdien07@mac.com
%
% bugfix 3/19/14 JD
% Fixed crash when segmenting or previewing.  Not sure why the syntax was working and suddenly was not.
%
% modified 6/18/14 JD
% Added starts, ends, contains, and -follows- keywords.
%
% modified 7/16/14 JD
% Simplified keys code field structure.
%
% modified 8/31/14 JD
% Added support for adding additional keys information to events in continuous data and trial specs in single-trial
% data.
%
% modified 9/4/14 JD
% Added delay field to segment function.
%
% modified 9/15/14 JD
% Added contents of type field to time-lock events list.
%
% modified 9/17/14 JD
% Adding event output for files that don't support events.
%
% modified 10/16/14 JD
% Passes screen size to ep_readData call.
%
% modified 10/30/14 JD
% Eliminated requirement that segmented epochs not overlap with each other.
%
% bugfix 8/14/15 JD
% Fixed criterion comparisons not being evaluated correctly for < and >.
%
% modified 8/17/15 JD
% Added ability to resegment (and thus reassign cells of) single-trial data.
%
% modified & bugfix 9/3/15 JD
% Adds next TRSP event after segmentation event without regard to whether
% the TRSP event fell within the segmentation period.
% The sectionTrial and SectionNumbers changed to start with 1 rather than 0.
% Added capability to resegment single-trial data to reassign cell
% conditions.
% Added capability to specify OR criteria for a condition by just giving
% them the same name.
% Trial spec names no longer continue to have TS- prefix when segmented
% data saved, which also resulted in trial specs not being recognized
% during segmentation.
%
% modified 12/15/15 JD
% When it runs into multiple matching specs, takes first one rather than dropping the segment.
%
% bugfix 12/18/15 JD
% Fixed crash when the criterion string is longer than the stimulus string for the "starts" and "ends" relations.
%
% modified 10/16/16 JD
% Added .stims field.
%
% modified 11/5/16 JD
% Added support for writing out subject spec text files.
%
% bugfix 11/8/16 JD
% Fixed segmentation ignoring the 6th criterion.
% Fixed section crit < and > relations being interpreted as the reverse direction.
%
% bugfix 3/16/17 JD
% Fixed error when batching multiple subject files and there are more than one row of criteria for the same cell name.
%
% bugfix 5/21/17 JD
% Fixed crash when resegmenting single-trial data.
% Fixed adding sec TrialSpec fields when already present, as in resegmenting data, resulting in crashes during averaging.
%
% modified & bugfix 6/21/17 JD
% Bad data edits passed on to segmented files.
% Flexible segments implemented.
% Fixed crash when segmenting time-frequency data.
% switch to amplitude scaling when adding freq data together.
% Fixed crash that can happen with continuous data.
%
% modified & bugfix 11/25/17 JD
% Added -precedes- crit option to the Segment function.
% Made fixes to -follows- function for relations other than = and ~=
% Added support for TS-EPoffset field.
%
% bugfix 3/1/18 JD
% Fixed unable to batch multiple files when the ced file contains BAD channels.
%
% bugfix 3/16/18 JD
% Now handles situation where the recording was aborted before the last trial spec was recorded.
%
% bugfix 4/29/18 JD
% No longer misses events when the .value field is a number rather than a string.
%
% bugfix 5/13/18 JD
% Now includes 'ep_segmentData' in .history field record.
%
% bugfix 6/14/18 JD
% Fixed using next available TRSP when a TRSP is missing instead of just dropping the trial.
%
% bugfix 7/23/18 JD
% Fixed recTime values not correct.
%
% modified 10/18/18 JD
% Added support for reading E-Prime text output to add trial spec information during segmentation.
% Changed "TS-EPoffset" to "EPoffset"
%
% bugfix 12/5/18 JD
% Fixed not using third exclusion criteria for specs text file.
% Fixed e-prime specs not being merged into data in case where more than one code is used as the time-lock event (e.g., where the event label codes the experimental condition).
%
% bugfix 3/8/19 JD
% Fixed -follows- and -precedes- options not working as expected.
%
% bugfix & modified 4/10/19 JD
% Fixed preview using contents of specs table even when from a different file.  Instead, user must now choose the text file when running Preview.
% Added support for task level performance measures.
%
% bugfix 6/25/19 JD
% Fixed not using trial specs when resegmenting single-trial data.
%
% bugfix & modified 9/13/19 JD
% Fixed crash when using preview mode with specs table.
% Example specs table is now linked to the example data so if the latter is changed then the former is cleared but otherwise it will be used in preview mode.
%
% bugfix & modified 11/3/19 JD
% Fixed criteria not working when .value or .type fields of the events structure are numeric rather than strings.
% Added -rangesamps- functionality.
% eventlock field is now always only based on the .type field, with additional specs needed to utilize the .value field.
% Fixed crash when segmenting continuous data with at least one recording interruption (a boundary event).
% Fixed to correctly compute number of samples when sample sizes are not integers.
%
% modified 12/23/19 JD
% Upgraded support of std information by adding .covAVE and .GAVsubs fields and eliminating .std and .stdCM fields.
%
% bugfix 3/23/20 JD
% Fixed not segmenting based on information in the value field of events.
%
% bugfix & modified 7/7/20 JD
% Fixed crash when segmenting data with 512Hz sampling rate.
% Fixed -precedes- option accepting all values of events as meeting criteria.
% Added -responseCorr- and -responseErr- options.
% Removed preference option to rotate electrode coordinates when importing mff or eeglab files that have internal coordinates as no longer needed.
% Added control to Segment function to specify the ACC and RT trial specs.
%
% bugfix 12/22/20 JD
% If there are no trial specs, will instead offer the key names, if any, of the event lock event for the ACC and RT popup menus.
% Resegmenting single-trial data no longer results in zero segments.
% Fixed crash when there are no trial specs and the Task field is being used.
%
% modified 3/10/21 JD
% For fractional event times, timeshifts the segments via interpolation so that the start of the epoch exactly corresponds to the event time.
% 
% modified 10/3/21 JD
% Upgraded history field to provide more information on changes.
%
% modified 7/24/22 JD
% Added support for reading Matlab .mat files.
%
% bugfix 4/28/23 JD
% Fixed crash when sampling rate is not a multiple of 500.
%
% modified 11/15/24 JD
% May now specify suffix for output segmented files via Preferences setting.
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

global EPmain EPtictoc

cellNums=[];
numOutCells=size(cellTable,1);
cellCount=zeros(length(inputFiles),numOutCells);
flexMode=strcmpi(cellTable{1,4}(1),'F');

%assume the first session file is representative of the rest.
for iFile=1:length(inputFiles)
    ep_tictoc;if EPtictoc.stop;return;end
    if ~preview
        disp(['Working on #' num2str(iFile) ': ' inputFiles{iFile} '.']);
        if iFile==1
            readArg{1}='format';
            readArg{2}=importFormat;
            readArg{3}='type';
            readArg{4}='continuous';
            readArg{7}='screenSize';
            readArg{8}=EPmain.scrsz;
            readArg{9}='FontSize';
            readArg{10}=EPmain.fontsize;
            readArg{11}='silent';
            readArg{12}='on';
            
            SMIsuffix=EPmain.preferences.general.SMIsuffix;
            if ~isempty(SMIsuffix)
                readArg{end+1}='SMIsuffix';
                readArg{end+1}=SMIsuffix;
            end
            specSuffix=EPmain.preferences.general.specSuffix;
            if ~isempty(specSuffix)
                readArg{end+1}='specSuffix';
                readArg{end+1}=specSuffix;
            end
            subjectSpecSuffix=EPmain.preferences.general.subjectSpecSuffix;
            if ~isempty(subjectSpecSuffix)
                readArg{end+1}='subjectSpecSuffix';
                readArg{end+1}=subjectSpecSuffix;
            end
            BVheader=EPmain.preferences.general.BVheader;
            if ~isempty(BVheader)
                readArg{end+1}='BVheader';
                readArg{end+1}=BVheader;
            end
            noInternal=EPmain.preferences.general.noInternal;
            if ~isempty(noInternal)
                readArg{end+1}='noInternal';
                readArg{end+1}=noInternal;
            end

            Name=deblank(inputFiles{1});
            thisReadArg=readArg;
            thisReadArg{end+1}='file';
            [pathstr, fileName, ext]=fileparts(Name);
            thisReadArg{end+1}=Name;
            [inputData, origEloc, outInfo]=ep_readData(thisReadArg);
            if isempty(inputData) || isempty(inputData.data)
                return;
            end
            readArg{end+1}='silent';
            readArg{end+1}='on';
            readArg{end+1}='ced';
            readArg{end+1}=inputData.ced; %assume all the files to be segmented will be using the same ced file.
%             readArg{end+1}='eloc';
%             readArg{end+1}=inputData.eloc; %assume all the files to be segmented will be using the same eloc info.
            readArg{end+1}='montage';
            readArg{end+1}=inputData.montage; %assume all the files to be segmented will be using the same montage.
        else
            thisReadArg=readArg;
            thisReadArg{end+1}='file';
            if strcmp(importFormat,'matlab_mat')
                thisReadArg{end+1}='matlabDims';
                thisReadArg{end+1}=outInfo.matlabDims;
            end
            Name=deblank(inputFiles{iFile});
            thisReadArg{end+1}=Name;
            [pathstr, fileName, ext]=fileparts(Name);
            [inputData, origEloc, outInfo]=ep_readData(thisReadArg);
        end
        if isempty(inputData) || isempty(inputData.data)
            msg{1}=['Error: The file ' fileName ' failed to be read.'];
            [msg]=ep_errorMsg(msg);
            continue
        end
        if ~any(strcmp(inputData.dataType,{'continuous','single_trial'}))
            msg{1}=['Error: The file ' fileName ' is not a continuous or single-trial file.'];
            [msg]=ep_errorMsg(msg);
            continue
        end
        if all(cellfun(@isempty,inputData.events))
            msg{1}=['Error: The file ' fileName ' has no events.'];
            [msg]=ep_errorMsg(msg);
            continue
        end
        numChans=length(inputData.chanNames);
        samplesize=1000/inputData.Fs;
        outputData=inputData;
        outputData.cellNames=cell(0);
        outputData.cellTypes=cell(0);
        outputData.trialSpecs=cell(0);
        outputData.events=cell(0);
        outputData.trialNames=[];
        if strcmp(inputData.dataType,'continuous')
            outputData.dataType='single_trial';
            if flexMode
                numInterPoints=str2num(cellTable{1,4}(2:end));
                outputData.Fs=numInterPoints*10;
                outputData.timeNames=[0:numInterPoints-1]'*samplesize;
                outputData.baseline=0;
                outputData.timeUnits='per';
            else
                outputData.baseline=-str2double(cellTable{1,3})/samplesize;
                outputData.timeNames=[str2double(cellTable{1,3}):samplesize:str2double(cellTable{1,4})-samplesize]';
                numPoints=length(outputData.timeNames);
            end
            outputData.recTime=[];
        end
        outputData.data=zeros(numChans,length(outputData.timeNames),0,1,1,max(1,length(outputData.freqNames)));
        numEpochPoints=length(outputData.timeNames);
    else
        inputData=inputFiles;
%         if ~isempty(specData)
%             if ismac
%                 disp('Load E-Prime text file') %workaround for bug in Mac version of Matlab.
%             end
%             [FileName,PathName,FilterIndex] = uigetfile('*.txt','Load E-Prime text file');
%             if ~isempty(FileName) && ~isnumeric(FileName)
%                 Name = [PathName filesep FileName];
%                 [pathstr, fileName, ext]=fileparts(Name);
%             else
%                 disp('No E-Prime text file chosen so specs table will be ignored.')
%                 Name='preview';
%                 specData=[];
%             end
%         end
        numChans=length(inputData.chanNames);
        samplesize=1000/inputData.Fs;
    end
    
    numEpochs=1;
    numPoints=length(inputData.timeNames);
    if strcmp(inputData.dataType,'single_trial')
        %convert single_trial data so it is organized like continuous   
        tempEvents=cell(1);
        if ~isempty(inputData.recTime)
            [~, trialOrder]=sort(inputData.recTime); %the trials are not necessarily organized chronologically
        else
            trialOrder=[1:length(inputData.cellNames)];
        end
        numEpochs=length(inputData.cellNames);
        if ~preview
            tempData=zeros(size(inputData.data,1),size(inputData.data,2)*size(inputData.data,3),1,size(inputData.data,4),size(inputData.data,5),size(inputData.data,6),size(inputData.data,7));
            for iTrial=1:length(trialOrder)
                theTrial=trialOrder(iTrial);
                tempData(:,(iTrial-1)*numPoints+1:iTrial*numPoints,1,:,:,:,:)=inputData.data(:,:,theTrial,:,:,:,:);
            end
            inputData.data=tempData;
            inputData.timeNames=[1:size(tempData,2)];
        end
        
        recTime=zeros(length(trialOrder),1);
        for iTrial=1:length(trialOrder)
            theTrial=trialOrder(iTrial);
            trialEvents=inputData.events{1,theTrial};
            for iTrialEvent=1:length(trialEvents)
                trialEvents(iTrialEvent).sample=trialEvents(iTrialEvent).sample+(iTrial-1)*numPoints;
            end
            tempEvents{1}=[tempEvents{1} trialEvents];
            recTime(iTrial)=inputData.recTime(theTrial);
        end
        inputData.events=tempEvents;
        outputData.recTime=[];
        outputData.baseline=-str2num(cellTable{1,3})/samplesize;
        if flexMode
            numInterPoints=str2num(cellTable{1,4}(2:end));
            outputData.Fs=numInterPoints*10;
            outputData.timeNames=[0:numInterPoints-1]'*(1000/outputData.Fs);
            outputData.baseline=0;
            outputData.timeUnits='per';
        else
            outputData.timeNames=[str2double(cellTable{1,3}):samplesize:str2double(cellTable{1,4})-samplesize]';
        end
        outputData.data=zeros(numChans,length(outputData.timeNames),0);
    end
    
    trialData=cell(0);
    if ~isempty(specData)
        if ~preview
            specFileName=[pathstr filesep fileName '.txt'];
            if exist(specFileName,'file')
                [experimentFieldNames, experimentData, theHeaders, trialData] = ep_loadSpecFile(specData.specFormat, specFileName);
%                 if ~preview
                    for iSpec=1:length(experimentFieldNames)
                        if ~any(strcmp(experimentFieldNames{iSpec},outputData.subjectSpecNames))
                            outputData.subjectSpecNames(end+1)=experimentFieldNames(iSpec);
                            outputData.subjectSpecs(end+1)=experimentData(iSpec);
                        else
                            outputData.subjectSpecs(max(find(strcmp(experimentFieldNames{iSpec},outputData.subjectSpecNames))))=experimentData(iSpec);
                        end
                    end
%                 end
            else
                theHeaders=cell(0);
                trialData=cell(0);
                experimentFieldNames=cell(0);
                experimentData=cell(0);
                disp(['Error: could not find ' specFileName '.  Perhaps the name of the E-Prime text file does not match up with the name of the EEG file?']);
            end
        else
            %if preview, then just use the example specs
            theHeaders=specData.specNames(1:end-1);
            trialData=specData.specTable;
            experimentFieldNames=cell(0);
            experimentData=cell(0);
        end
        %eliminate practice trials and other such events that do not correspond to even potential segments
        if any(strcmp(specData.excludeSpec,theHeaders))
            trialData(strcmp(specData.excludeValue,trialData(:,strcmp(specData.excludeSpec,theHeaders))),:)=[];
        end
        if any(strcmp(specData.excludeSpec2,theHeaders))
            trialData(strcmp(specData.excludeValue2,trialData(:,strcmp(specData.excludeSpec2,theHeaders))),:)=[];
        end
        if any(strcmp(specData.excludeSpec3,theHeaders))
            trialData(strcmp(specData.excludeValue3,trialData(:,strcmp(specData.excludeSpec3,theHeaders))),:)=[];
        end
        outSpecData=cell(0);
        outSpecLabels=cell(0);
        for iSpec=1:length(specData.specField)
            if ~isempty(specData.specField{iSpec}) && ~isempty(trialData)
                outSpecData(:,end+1)=trialData(:,strcmp(specData.specField{iSpec},theHeaders),:);
                outSpecLabels{end+1}=[specData.specFormat '-' specData.specLabels{iSpec}];
            end
        end
    end
    
    %need to convert .type and .value fields into strings as otherwise criteria might not work
    %since the previous step reorganized the data as if it were continuous, all events are in the same {1} cell.
    eventsString=inputData.events{1};
    for iEvent=1:length(eventsString)
        eventsString(iEvent).value=num2str(eventsString(iEvent).value);
        eventsString(iEvent).type=num2str(eventsString(iEvent).type);
    end
    
    allEvents={eventsString.value};
    [allEvents{find(cellfun(@isempty,allEvents))}]=deal(' ');
    eventSamples=[eventsString.sample];
    [orderedEventSamples, eventOrder]=sort(eventSamples); %the events are not necessarily organized chronologically
    orderedEvents=eventsString(eventOrder);
    badTrials=cell(numOutCells,1);
    badChans=cell(numOutCells,1);
    outTRSPnames=cell(0);
    TRSPindex=find(strcmp('TRSP',{orderedEvents.value}));
    TRSPsamples=[orderedEvents(TRSPindex).sample];
        
    %match up events of interest with the e-prime specs if any
    trialSpecFlag=0;
    if ~isempty(specData)
        typeIndex=find(ismember({orderedEvents.type},cellTable(:,2)));
        valueIndex=find(ismember(cellfun(@num2str,{orderedEvents.value},'UniformOutput',false),cellTable(:,2)));
        specIndex=union(typeIndex,valueIndex); %list of events that correspond to e-prime specs based on temporally ordered events
        
        if size(outSpecData,1) ~=length(specIndex)
            disp(['Warning: Unable to match trial specs file information (' num2str(size(outSpecData,1)) ') with potential time-locking events (' num2str(length(specIndex)) ').'])
        else
            trialSpecFlag=1;
        end
        
    else
        specIndex=[];
    end
    
    if strcmp(inputData.dataType,'single_trial')
        totalNumPoints=numPoints*numEpochs;
    else
        totalNumPoints=numPoints;
    end
    for iCell=1:numOutCells
        ep_tictoc;if EPtictoc.stop;return;end
        %generate list of potential epochs
        typeIndex=find(strcmp(cellTable{iCell,2},{orderedEvents.type}));
        valueIndex=find(strcmp(cellTable{iCell,2},cellfun(@num2str,{orderedEvents.value},'UniformOutput',false)));
        eventIndex=union(valueIndex, typeIndex);
        eventSamples=[orderedEvents(eventIndex).sample];
        delayTime=str2double(cellTable{iCell,5})/(1000/inputData.Fs);
        taskLabel=cellTable{iCell,6};
        
        if ~flexMode
%             if strcmp(inputData.dataType,'continuous')
                temp=str2double(cellTable{iCell,3})/(1000/inputData.Fs);
                prestim=ceil(abs(temp))*sign(temp);
                temp=str2double(cellTable{iCell,4})/(1000/inputData.Fs);
                poststim=floor(abs(temp))*sign(temp);
                epochLength=poststim-prestim;
                epochStart=eventSamples+prestim+delayTime;
                epochEnd=eventSamples+poststim+delayTime-1;
                %check to see if full epoch within recording session
                goodEpochs=intersect(find(epochStart>0),find((epochStart+epochLength-1)<totalNumPoints+1));
                epochRecTime=epochStart(goodEpochs);
%             elseif strcmp(inputData.dataType,'single_trial')
%                 %check to see if new epochs fall within bounds of original epochs
%                 prestim=-length(find(inputData.timeNames<0));
%                 poststim=length(find(inputData.timeNames>=0));
%                 epochLength=poststim-prestim;
%                 epochStart=eventSamples+prestim+delayTime;
%                 epochEnd=eventSamples+poststim+delayTime-1;
%                 eventEpoch=floor((eventSamples-1)/numPoints)+1;
%                 startEpoch=floor((epochStart-1)/numPoints)+1;
%                 endEpoch=floor((epochEnd-1)/numPoints)+1;
%                 goodEpochs=intersect(find(eventEpoch==startEpoch),find(eventEpoch==endEpoch));
%                 epochRecTime=recTime(eventEpoch);
%             end
        else
            %see if flex start event is followed by a flex end event without another intervening start event
            epochStart=eventSamples+delayTime;
            flexEndIndex=union(find(strcmp(cellTable{iCell,3},{orderedEvents.type})),find(strcmp(cellTable{iCell,3},{orderedEvents.value})));
            flexEndSamples=round([orderedEvents(flexEndIndex).sample]);
            epochEnd=[];
            goodEpochs=[];
            for iEvent=1:length(eventIndex)
                endSample=min(flexEndSamples(flexEndSamples>eventSamples(iEvent)));
                epochEnd(end+1)=endSample+delayTime;
                nextEventSample=min(eventSamples(eventSamples>eventSamples(iEvent)));
                if isempty(nextEventSample) || (endSample <= nextEventSample)
                    goodEpochs(end+1)=iEvent;
                end
            end
            epochRecTime=epochStart(goodEpochs);
        end
        eventIndex=eventIndex(goodEpochs);
        epochStart=epochStart(goodEpochs);
        epochEnd=epochEnd(goodEpochs);
                
        if strcmp(inputData.dataType,'continuous')
            %check to see if boundary event falls within epoch
            boundaryIndex=find(strcmp('boundary',{orderedEvents.value}));
            boundarySamples=[orderedEvents(boundaryIndex).sample];
            badEpochs=[];
            for iBoundary=1:length(boundaryIndex)
                badEpochs=[badEpochs find((epochStart<=boundarySamples(iBoundary)) & ((epochStart+epochLength-1)>=boundarySamples(iBoundary)))];
            end
            badEpochs=unique(badEpochs);
            eventIndex(badEpochs)=[];
            epochStart(badEpochs)=[];
            epochEnd(badEpochs)=[];
            epochRecTime(badEpochs)=[];
        end
        
        %determine which potential epochs meets trial spec criteria
        
        goodSpecs=ones(EPmain.segment.numSpecs,1);
        for iSpec=1:EPmain.segment.numSpecs
            if isempty(deblank(cellTable{iCell,EPmain.segment.numFixed+(iSpec-1)*3}))
                goodSpecs(iSpec)=0;
            end
            if isempty(deblank(cellTable{iCell,EPmain.segment.numFixed+1+(iSpec-1)*3}))
                goodSpecs(iSpec)=0;
            elseif any(strcmp(deblank(cellTable{EPmain.segment.numFixed+1+(iSpec-1)*3}),{'=','~='}))
                if isempty(deblank(num2str(cellTable{iCell,EPmain.segment.numFixed+2+(iSpec-1)*3})))
                    goodSpecs(iSpec)=0;
                end
            elseif isempty(deblank(cellTable{iCell,EPmain.segment.numFixed+1+(iSpec-1)*3}))
                goodSpecs(iSpec)=0;
            end
        end
        goodSpecs=find(goodSpecs);
        badEpochs=[];
        allTRSPspecValues=cell(length(eventIndex),1);
        TRSPspecNamesList=cell(0);
        sectionNumber=1;
        sectionTrial=1;
        sectionStart=0;
        sectionNumberList=[];
        sectionTrialList=[];
        for iEvent=1:length(eventIndex) %looping through events potentially defining trials
            theEvent=eventIndex(iEvent);
            stimSpecNames=cell(0); %list of spec names for this event
            stimSpecValues=cell(0);%list of the values for these specs
            stimSpecNames{end+1}='value';
            stimSpecValues{end+1}=orderedEvents(theEvent).value;
            for iKey=1:length(orderedEvents(theEvent).keys)
                stimSpecNames{end+1}=orderedEvents(theEvent).keys(iKey).code;
                stimSpecValues{end+1}=orderedEvents(theEvent).keys(iKey).data;
            end
            if strcmp(inputData.dataType,'single_trial')
                theTrial=ceil(epochStart(iEvent)/numPoints);
                for iTrialSpec=1:length(inputData.trialSpecNames)
                    stimSpecNames{end+1}=inputData.trialSpecNames{iTrialSpec};
                    stimSpecValues{end+1}=inputData.trialSpecs{theTrial,iTrialSpec};
                end
            end
            
            TRSPspecNames=cell(0);
            TRSPspecValues=cell(0);
            if strcmp(inputData.dataType,'continuous')
                epochTRSP=TRSPindex(find((TRSPsamples>=epochStart(iEvent))));
                if ~isempty(epochTRSP)
                    epochTRSP=epochTRSP(1);
                    %if (iEvent==length(eventIndex)) || (orderedEvents(epochTRSP).sample < orderedEvents(eventIndex(iEvent+1)).sample)
                        for iKey=1:length(orderedEvents(epochTRSP).keys)
                            TRSPspecNames{end+1,1}=['TS-' orderedEvents(epochTRSP).keys(iKey).code];
                            TRSPspecValues{end+1,1}=orderedEvents(epochTRSP).keys(iKey).data;
                        end
                    %end
                end
                
%             elseif strcmp(inputData.dataType,'single_trial')
%                 %theEpoch=floor((eventSamples(eventIndex(iEvent))-1)/length(inputData.timeNames))+1;
%                 for iKey=1:length(inputData.trialSpecNames)
%                     TRSPspecNames{end+1,1}=['TS-' inputData.trialSpecNames{iKey}];
%                     TRSPspecValues{end+1,1}=inputData.trialSpecs{iEvent,iKey};
%                 end              
%             else
%                 msg{1}=['Error: The file ' fileName ' is not a continuous or single-trial file.'];
%                 [msg]=ep_errorMsg(msg);
%                 continue
            end
            
            for iSpec=1:length(goodSpecs)
                if any(strcmp(cellTable{iCell,EPmain.segment.numFixed+(iSpec-1)*3},{'-responseCorr-','-responseErr-'}))
                    if ~any(strcmp('TS-RT',TRSPspecNames))
                        TRSPspecNames{end+1,1}='TS-RT';
                        TRSPspecValues{end+1,1}=[];
                    end
                    if ~any(strcmp('TS-ACC',TRSPspecNames))
                        TRSPspecNames{end+1,1}='TS-ACC';
                        TRSPspecValues{end+1,1}=[];
                    end
                end
            end
            
            if trialSpecFlag
                TRSPspecNames(end+1:end+length(outSpecLabels),1)=outSpecLabels;
                TRSPspecValues(end+1:end+size(outSpecData,2),1)=outSpecData(find(specIndex==eventIndex(iEvent)),:);
            end
            
            if ~isempty(TRSPspecNames)
                TRSPspecNamesList=TRSPspecNames;
            end
            if isempty(outTRSPnames)
                outTRSPnames=TRSPspecNames;
                outTRSPnames{end+1,1}='secNum';
                outTRSPnames{end+1,1}='secTrial';
                outTRSPnames{end+1,1}='secTrialFromEnd';
            elseif ~all(ismember(TRSPspecNames,outTRSPnames))
                for iSpec=1:length(TRSPspecNames)
                    if ~any(strcmp(TRSPspecNames{iSpec},outTRSPnames))
                        outTRSPnames(end+1,1)=TRSPspecNames(iSpec);
                    end
                end
            end
            for iSpec=1:length(TRSPspecNames)
                theSpec=find(strcmp(TRSPspecNames{iSpec},outTRSPnames));
                if length(theSpec)>1
                    theSpec=theSpec(1);
                end
                allTRSPspecValues{iEvent}{theSpec}=TRSPspecValues{iSpec};
            end
                        
            keepTrial=1;
            afterPrecede=0;
            firstPrecede=0; %event number (sorted) of the first follow criterion event
            precedeNeg=0; %next spec is following a ~= -precedes- criterion            
            afterFollow=0;
            firstFollow=0; %event number (sorted) of the first follow criterion event
            followNeg=0; %next spec is following a ~= -follows- criterion
            stimTime=[]; %stim sample for -responseCorr- and -responseErr
            theAcc=[]; %accuracy for -responseCorr- and -responseErr
            for iSpec=1:length(goodSpecs)
                theSpec=goodSpecs(iSpec);
                stimSpec=find(strcmp(cellTable{iCell,EPmain.segment.numFixed+(iSpec-1)*3},stimSpecNames));
                if length(stimSpec) > 1
                    disp(['Error: multiple matching stimulus specs (' stimSpecNames{stimSpec(1)} ') for ' Name '.  Will use the first instance.']);
                    stimSpec=stimSpec(1);
                    %keepTrial=0;
                end
                noStim=0;
                badStim=0;
                if afterFollow && ~any(ismember(theSpec-1,goodSpecs))
                    afterFollow=0;
                end
                if afterPrecede && ~any(ismember(theSpec-1,goodSpecs))
                    afterPrecede=0;
                end
                
                %if it is a '-precedes-' spec
                if any(strcmp(cellTable{iCell,EPmain.segment.numFixed+(iSpec-1)*3},{'-precedes-','-responseCorr-','-responseErr-'}))
                    afterPrecede=1;
                    thePrecedeName=cellTable{iCell,EPmain.segment.numFixed+2+(iSpec-1)*3};
                    if firstPrecede %if there is already a precede in effect
                        precede1=min(find(orderedEventSamples==eventsString(theEvent).sample))+1; %ranked order
                        precede2=firstPrecede-1; %ranked order
                    else
                        precede1=min(find(orderedEventSamples==eventsString(theEvent).sample))+1; %ranked order
                        precede2=length(orderedEvents); %ranked order
                        sampleRange=[eventsString(theEvent).sample+1:length(inputData.timeNames)];
                        if theSpec<EPmain.segment.numSpecs && strcmp('-rangeSamps-',cellTable{iCell,EPmain.segment.numFixed+(theSpec-1+2)*3})
                            theTime=cellTable{iCell,EPmain.segment.numFixed+2+(theSpec-1+2)*3};
                            if isnumeric(theTime)
                                theCritTime=eventsString(theEvent).sample+theTime;
                            else
                                theCritTime=eventsString(theEvent).sample+str2double(theTime);
                            end
                            switch deblank(cellTable{iCell,EPmain.segment.numFixed+1+(theSpec-1+2)*3})
                                case '='
                                    sampleRange=theCritTime;
                                case '~='
                                    sampleRange=setdiff(sampleRange,theCritTime);
                                case '<'
                                    sampleRange=sampleRange(sampleRange<theCritTime);
                                case '>'
                                    sampleRange=sampleRange(sampleRange>theCritTime);
                                case '<='
                                    sampleRange=sampleRange(sampleRange<=theCritTime);
                                case '>='
                                    sampleRange=sampleRange(sampleRange>=theCritTime);
                                otherwise
                                    disp(['Error: cell table spec relationship (' cellTable{iCell,EPmain.segment.numFixed+1+(theSpec-1+2)*3} ') not valid.']);
                                    return
                            end
                            if str2double(cellTable{iCell,EPmain.segment.numFixed+2+(theSpec-1+2)*3}) < 1
                                disp(['Error: cell table spec value (' cellTable{iCell,EPmain.segment.numFixed+2+(theSpec-1+2)*3} ') not valid.']);
                                return
                            end
                        end
                    end
                    
                    switch deblank(cellTable{iCell,EPmain.segment.numFixed+1+(theSpec-1)*3})
                        case '='
                            precedeList=find(strcmp(thePrecedeName,{orderedEvents.type}));
                            precedeList=intersect(precedeList,[precede1:precede2]);
                            precedeList=precedeList(ismember([orderedEvents(precedeList).sample],sampleRange));
                            if isempty(precedeList)
                                badStim=1;
                            end
                        case '~='
                            precedeList=find(~strcmp(thePrecedeName,{orderedEvents.type}));
                            precedeList=intersect(precedeList,[precede1:precede2]);
                            precedeList=precedeList(ismember([orderedEvents(precedeList).sample],sampleRange));
                            if isempty(precedeList)
                                badStim=1;
                            else
                                precedeNeg=1;
                            end
                        case '<'
                            precedeList=[];
                            for iPrecEvent=1:length(orderedEvents)
                                if ~isnan(str2double(thePrecedeName))
                                    if str2double(orderedEvents(iPrecEvent).type)<str2double(thePrecedeName)
                                        precedeList(end+1)=iPrecEvent;
                                    end
                                else
                                    if orderedEvents(iPrecEvent).type<thePrecedeName
                                        precedeList(end+1)=iPrecEvent;
                                    end
                                end
                            end
                            precedeList=intersect(precedeList,[precede1:precede2]);
                            precedeList=precedeList(ismember([orderedEvents(precedeList).sample],sampleRange));
                            if isempty(precedeList)
                                badStim=1;
                            end
                        case '>'
                            precedeList=[];
                            for iPrecEvent=1:length(orderedEvents)
                                if ~isnan(str2double(thePrecedeName))
                                    if str2double(orderedEvents(iPrecEvent).type)>str2double(thePrecedeName)
                                        precedeList(end+1)=iPrecEvent;
                                    end
                                else
                                    if orderedEvents(iPrecEvent).type>thePrecedeName
                                        precedeList(end+1)=iPrecEvent;
                                    end
                                end
                            end
                            precedeList=intersect(precedeList,[precede1:precede2]);
                            precedeList=precedeList(ismember([orderedEvents(precedeList).sample],sampleRange));
                            if isempty(precedeList)
                                badStim=1;
                            end
                        case '<='
                            precedeList=[];
                            for iPrecEvent=1:length(orderedEvents)
                                if ~isnan(str2double(thePrecedeName))
                                    if str2double(orderedEvents(iPrecEvent).type)<=str2double(thePrecedeName)
                                        precedeList(end+1)=iPrecEvent;
                                    end
                                else
                                    if orderedEvents(iPrecEvent).type<=thePrecedeName
                                        precedeList(end+1)=iPrecEvent;
                                    end
                                end
                            end
                            precedeList=intersect(precedeList,[precede1:precede2]);
                            precedeList=precedeList(ismember([orderedEvents(precedeList).sample],sampleRange));
                            if isempty(precedeList)
                                badStim=1;
                            end
                        case '>='
                            precedeList=[];
                            for iPrecEvent=1:length(orderedEvents)
                                if ~isnan(str2double(thePrecedeName))
                                    if str2double(orderedEvents(iPrecEvent).type)>=str2double(thePrecedeName)
                                        precedeList(end+1)=iPrecEvent;
                                    end
                                else
                                    if orderedEvents(iPrecEvent).type>=thePrecedeName
                                        precedeList(end+1)=iPrecEvent;
                                    end
                                end
                            end
                            precedeList=intersect(precedeList,[precede1:precede2]);
                            precedeList=precedeList(ismember([orderedEvents(precedeList).sample],sampleRange));
                            if isempty(precedeList)
                                badStim=1;
                            end
                        case 'starts'
                            precedeList=[];
                            for iPrecEvent=1:length(orderedEvents)
                                if strfind(orderedEvents(iPrecEvent).type,thePrecedeName)==1
                                    precedeList(end+1)=iPrecEvent;
                                end
                            end
                            precedeList=intersect(precedeList,[precede1:precede2]);
                            precedeList=precedeList(ismember([orderedEvents(precedeList).sample],sampleRange));
                            if isempty(precedeList)
                                badStim=1;
                            end
                        case 'ends'
                            precedeList=[];
                            for iPrecEvent=1:length(orderedEvents)
                                if strfind(orderedEvents(iPrecEvent).type,thePrecedeName)==(length(allEvents{iPrecEvent})-length(thePrecedeName)+1)
                                    precedeList(end+1)=iPrecEvent;
                                end
                            end
                            precedeList=intersect(precedeList,[precede1:precede2]);
                            precedeList=precedeList(ismember([orderedEvents(precedeList).sample],sampleRange));
                            if isempty(precedeList)
                                badStim=1;
                            end
                        case 'contains'
                            precedeList=[];
                            for iPrecEvent=1:length(orderedEvents)
                                if ~isempty(strfind(theEventsList{iList},thePrecedeName))
                                    precedeList(end+1)=iPrecEvent;
                                end
                            end
                            precedeList=intersect(precedeList,[precede1:precede2]);
                            precedeList=precedeList(ismember([orderedEvents(precedeList).sample],sampleRange));
                            if isempty(precedeList)
                                badStim=1;
                            end
                        otherwise
                            disp(['Error: cell table spec relationship (' cellTable{iCell,EPmain.segment.numFixed+1+(iSpec-1)*3} ') not valid.']);
                            return
                    end
                    if badStim
                        keepTrial=0;
                    else
                        if strcmp(cellTable{iCell,EPmain.segment.numFixed+(iSpec-1)*3},'-responseCorr-')
%                             TRSPspecValues{find(strcmp('TS-ACC',TRSPspecNames))}=1;
%                             TRSPspecValues{find(strcmp('TS-RT',TRSPspecNames))}=0;
                              stimTime=eventsString(theEvent).sample;
                              theAcc=1;
                        end
                        if strcmp(cellTable{iCell,EPmain.segment.numFixed+(iSpec-1)*3},'-responseErr-')
                              stimTime=eventsString(theEvent).sample;
                              theAcc=0;
                        end
                    end
                    
                %if it is a '-follows-' spec
                elseif strcmp('-follows-',cellTable{iCell,EPmain.segment.numFixed+(theSpec-1)*3})
                    afterFollow=1;
                    theFollowName=cellTable{iCell,EPmain.segment.numFixed+2+(theSpec-1)*3};
                    if firstFollow %if there is already a follow in effect
                        follow1=firstFollow+1; %ranked order
                        follow2=min(find(orderedEventSamples==eventsString(theEvent).sample))-1; %ranked order
                    else
                        follow1=1; %ranked order
                        follow2=min(find(orderedEventSamples==eventsString(theEvent).sample))-1; %ranked order
                        sampleRange=[1:eventsString(theEvent).sample-1];
                        if theSpec<EPmain.segment.numSpecs && strcmp('-rangeSamps-',cellTable{iCell,EPmain.segment.numFixed+(theSpec-1+2)*3})
                            theTime=cellTable{iCell,EPmain.segment.numFixed+2+(theSpec-1+2)*3};
                            if isnumeric(theTime)
                                theCritTime=eventsString(theEvent).sample-theTime;
                            else
                                theCritTime=eventsString(theEvent).sample-str2double(theTime);
                            end
                            switch deblank(cellTable{iCell,EPmain.segment.numFixed+1+(theSpec-1+2)*3})
                                case '='
                                    sampleRange=theCritTime;
                                case '~='
                                    sampleRange=setdiff(sampleRange,theCritTime);
                                case '<'
                                    sampleRange=sampleRange(sampleRange>theCritTime);
                                case '>'
                                    sampleRange=sampleRange(sampleRange<theCritTime);
                                case '<='
                                    sampleRange=sampleRange(sampleRange>=theCritTime);
                                case '>='
                                    sampleRange=sampleRange(sampleRange<=theCritTime);
                                otherwise
                                    disp(['Error: cell table spec relationship (' cellTable{iCell,EPmain.segment.numFixed+1+(theSpec-1+2)*3} ') not valid.']);
                                    return
                            end
                            if str2double(cellTable{iCell,EPmain.segment.numFixed+2+(theSpec-1+2)*3}) < 1
                                disp(['Error: cell table spec value (' cellTable{iCell,EPmain.segment.numFixed+2+(theSpec-1+1)*3} ') not valid.']);
                                return
                            end
                        end
                    end
                    
                    switch deblank(cellTable{iCell,EPmain.segment.numFixed+1+(theSpec-1)*3})
                        case '='
                            followList=find(strcmp(theFollowName,{orderedEvents.type}));
                            followList=intersect(followList,[follow1:follow2]);
                            followList=followList(ismember([orderedEvents(followList).sample],sampleRange));
                            if isempty(followList)
                                badStim=1;
                            end
                        case '~='
                            followList=find(~strcmp(theFollowName,{orderedEvents.type}));
                            followList=intersect(followList,[follow1:follow2]);
                            followList=followList(ismember([orderedEvents(followList).sample],sampleRange));
                            if isempty(followList)
                                badStim=1;
                            else
                                followNeg=1;
                            end
                        case '<'
                            followList=[];
                            for iFollEvent=1:length(orderedEvents)
                                if ~isnan(str2double(theFollowName))
                                    if str2double(orderedEvents(iFollEvent).type)<str2double(theFollowName)
                                        followList(end+1)=iFollEvent;
                                    end
                                else
                                    if orderedEvents(iFollEvent).type<theFollowName
                                        followList(end+1)=iFollEvent;
                                    end
                                end
                            end
                            followList=intersect(followList,[follow1:follow2]);
                            followList=followList(ismember([orderedEvents(followList).sample],sampleRange));
                            if isempty(followList)
                                badStim=1;
                            end
                        case '>'
                            followList=[];
                            for iFollEvent=1:length(orderedEvents)
                                if ~isnan(str2double(theFollowName))
                                    if str2double(orderedEvents(iFollEvent).type)>str2double(theFollowName)
                                        followList(end+1)=iFollEvent;
                                    end
                                else
                                    if orderedEvents(iFollEvent).type>theFollowName
                                        followList(end+1)=iFollEvent;
                                    end
                                end
                            end
                            followList=intersect(followList,[follow1:follow2]);
                            followList=followList(ismember([orderedEvents(followList).sample],sampleRange));
                            if isempty(followList)
                                badStim=1;
                            end
                        case '<='
                            followList=[];
                            for iFollEvent=1:length(orderedEvents)
                                if ~isnan(str2double(theFollowName))
                                    if str2double(orderedEvents(iFollEvent).type)<=str2double(theFollowName)
                                        followList(end+1)=iFollEvent;
                                    end
                                else
                                    if orderedEvents(iFollEvent).type<=theFollowName
                                        followList(end+1)=iFollEvent;
                                    end
                                end
                            end
                            followList=intersect(followList,[follow1:follow2]);
                            followList=followList(ismember([orderedEvents(followList).sample],sampleRange));
                            if isempty(followList)
                                badStim=1;
                            end
                        case '>='
                            followList=[];
                            for iFollEvent=1:length(orderedEvents)
                                if ~isnan(str2double(theFollowName))
                                    if str2double(orderedEvents(iFollEvent).type)>=str2double(theFollowName)
                                        followList(end+1)=iFollEvent;
                                    end
                                else
                                    if orderedEvents(iFollEvent).type>=theFollowName
                                        followList(end+1)=iFollEvent;
                                    end
                                end
                            end
                            followList=intersect(followList,[follow1:follow2]);
                            followList=followList(ismember([orderedEvents(followList).sample],sampleRange));
                            if isempty(followList)
                                badStim=1;
                            end
                        case 'starts'
                            followList=[];
                            for iFollEvent=1:length(orderedEvents)
                                if strfind(orderedEvents(iFollEvent).type,theFollowName)==1
                                    followList(end+1)=iFollEvent;
                                end
                            end
                            followList=intersect(followList,[follow1:follow2]);
                            followList=followList(ismember([orderedEvents(followList).sample],sampleRange));
                            if isempty(followList)
                                badStim=1;
                            end
                        case 'ends'
                            followList=[];
                            for iFollEvent=1:length(orderedEvents)
                                if strfind(orderedEvents(iFollEvent).type,theFollowName)==(length(allEvents{iFollEvent})-length(theFollowName)+1)
                                    followList(end+1)=iFollEvent;
                                end
                            end
                            followList=intersect(followList,[follow1:follow2]);
                            followList=followList(ismember([orderedEvents(followList).sample],sampleRange));
                            if isempty(followList)
                                badStim=1;
                            end
                        case 'contains'
                            followList=[];
                            for iFollEvent=1:length(orderedEvents)
                                if ~isempty(strfind(theEventsList{iList},theFollowName))
                                    followList(end+1)=iFollEvent;
                                end
                            end
                            followList=intersect(followList,[follow1:follow2]);
                            followList=followList(ismember([orderedEvents(followList).sample],sampleRange));
                            if isempty(followList)
                                badStim=1;
                            end
                        otherwise
                            disp(['Error: cell table spec relationship (' cellTable{iCell,EPmain.segment.numFixed+1+(theSpec-1)*3} ') not valid.']);
                            return
                    end
                    if badStim
                        keepTrial=0;
                    else
                        if sectionStart ~= firstFollow
                            sectionStart=firstFollow;
                            sectionNumber=sectionNumber+1;
                            sectionTrial=1;
                        end
                    end
                    
                elseif strcmp(cellTable{iCell,EPmain.segment.numFixed+(theSpec-1)*3},'-rangeSamps-')
                    if (theSpec<3) || ~any(strcmp(cellTable{iCell,EPmain.segment.numFixed+(theSpec-1-2)*3},{'-precedes-','-follows-','-responseCorr-','-responseErr-'}))
                        disp(['Error: cell table spec criteria (' cellTable{iCell,EPmain.segment.numFixed+(theSpec-1)*3} ') not valid.']);
                        return
                    end
                    %if it is a regular spec
                elseif ~any(strcmp(cellTable{iCell,EPmain.segment.numFixed+(theSpec-1)*3},{'-secNum-','-secTrial-','-secTrialFromEnd-'}))
                    theCriterion=cellTable{iCell,EPmain.segment.numFixed+2+(theSpec-1)*3};
                    stimList=cell(0);
                    stimEventList=[];
                    if afterPrecede %match to -precedes- event specs rather than the lock event specs
                        for iPrecede=1:length(precedeList)
                            thePrecedeEvent=precedeList(iPrecede);
                            precedeStimSpecNames=cell(0);
                            precedeSpecSpecValues=cell(0);
                            precedeStimSpecNames{1}='value';
                            precedeSpecSpecValues{1}=eventsString(eventOrder(thePrecedeEvent)).value;
                            for iKey=1:length(eventsString(eventOrder(thePrecedeEvent)).keys)
                                theCode=eventsString(eventOrder(thePrecedeEvent)).keys(iKey).code;
                                theData=eventsString(eventOrder(thePrecedeEvent)).keys(iKey).data;
                                if ~isempty(theCode) && ~isempty(theData)
                                    precedeStimSpecNames{end+1}=theCode;
                                    precedeSpecSpecValues{end+1}=theData;
                                end
                            end
                            %if the -precedes- event has no specs and yet a spec condition was specified, then treated as no precede spec was present, resulting in a bad epoch.
                            stimSpec=find(strcmp(cellTable{iCell,EPmain.segment.numFixed+(theSpec-1)*3},precedeStimSpecNames));
                            if ~isempty(stimSpec)
                                stimList{end+1}=num2str(precedeSpecSpecValues{stimSpec});
                                stimEventList(end+1)=thePrecedeEvent;
                                %firstPrecede=thePrecedeEvent;
                            end
                        end
                    elseif afterFollow %match to -follows- event specs rather than the lock event specs
                        for iFollow=length(followList):-1:1
                            theFollowEvent=followList(iFollow);
                            followStimSpecNames=cell(0);
                            followSpecSpecValues=cell(0);
                            followStimSpecNames{1}='value';
                            followSpecSpecValues{1}=eventsString(eventOrder(theFollowEvent)).value;
                            for iKey=1:length(eventsString(eventOrder(theFollowEvent)).keys)
                                theCode=eventsString(eventOrder(theFollowEvent)).keys(iKey).code;
                                theData=eventsString(eventOrder(theFollowEvent)).keys(iKey).data;
                                if ~isempty(theCode) && ~isempty(theData)
                                    followStimSpecNames{end+1}=theCode;
                                    followSpecSpecValues{end+1}=theData;
                                end
                            end
                            %if the -follows- event has no specs and yet a spec condition was specified, then treated as no follow spec was present, resulting in a bad epoch.
                            stimSpec=find(strcmp(cellTable{iCell,EPmain.segment.numFixed+(theSpec-1)*3},followStimSpecNames));
                            if ~isempty(stimSpec)
                                stimList{end+1}=num2str(followSpecSpecValues{stimSpec});
                                stimEventList(end+1)=theFollowEvent;
                                %firstFollow=theFollowEvent;
                            end
                        end
                    else
                        if ~isempty(stimSpec)
                            theStim=num2str(stimSpecValues{stimSpec});
                            stimList{1}=theStim;
                        else
                            theStim=[];
                        end
                        stimEventList=theEvent;
                    end
                    
                    if isempty(stimList)
                        noStim=1;
                    else
                        for iStim=1:length(stimList)
                            badStim=0;
                            theStim=stimList{iStim};
                            cmp=ep_compareStrings(theStim,theCriterion);
                            if isempty(cmp)
                                badStim=1;
                            else
                                switch deblank(cellTable{iCell,EPmain.segment.numFixed+1+(theSpec-1)*3})
                                    case '='
                                        if cmp ~=0
                                            badStim=1;
                                        end
                                    case '~='
                                        if cmp ==0
                                            badStim=1;
                                        end
                                    case '<'
                                        if cmp > -1
                                            badStim=1;
                                        end
                                    case '>'
                                        if cmp < 1
                                            badStim=1;
                                        end
                                    case '<='
                                        if cmp > 0
                                            badStim=1;
                                        end
                                    case '>='
                                        if cmp < 0
                                            badStim=1;
                                        end
                                    case 'starts'
                                        if (length(theCriterion)>length(theStim)) || ~strcmp(theCriterion,theStim(1:length(theCriterion)))
                                            badStim=1;
                                        end
                                    case 'ends'
                                        if (length(theCriterion)>length(theStim)) || ~strcmp(theCriterion,theStim(end-length(theCriterion)+1:end))
                                            badStim=1;
                                        end
                                    case 'contains'
                                        if isempty(findstr(theCriterion,theStim))
                                            badStim=1;
                                        end
                                    otherwise
                                        disp(['Error: cell table spec relationship (' cellTable{iCell,EPmain.segment.numFixed+1+(theSpec-1)*3} ') not valid.']);
                                        return
                                end
                            end
                            if badStim==0
                                if afterPrecede
                                    firstPrecede=stimEventList(iStim);
                                elseif afterFollow
                                    firstFollow=stimEventList(iStim);
                                end
                                break %just need one good one
                            end
                        end
                        if followNeg || precedeNeg 
                            if badStim==0
                                badStim=1; %if any at all match then bad
                            else
                                badStim=0; %if none matched then good
                            end
                        end
                        followNeg=0;
                        precedeNeg=0;
                    end
                    if badStim
                        keepTrial=0;
                    end
                    
                    noTRSP=0;
                    if ~afterFollow && ~afterPrecede %TRSP crits do not succeed a -follows- or -precedes- criteria
                        TRSPspec=find(strcmp(cellTable{iCell,EPmain.segment.numFixed+(theSpec-1)*3},TRSPspecNames));
                        if length(TRSPspec) > 1
                            disp(['Error: multiple matching TRSP specs (' TRSPspecNames{TRSPspec(1)} ') for ' Name '.']);
                            keepTrial=0;
                        end
                        badTRSP=0;
                        if isempty(TRSPspec)
                            noTRSP=1;
                        else
                            cmp=ep_compareStrings(TRSPspecValues{TRSPspec},cellTable{iCell,EPmain.segment.numFixed+2+(theSpec-1)*3});
                            if isempty(cmp)
                                badTRSP=1;
                            else
                                switch deblank(cellTable{iCell,EPmain.segment.numFixed+1+(theSpec-1)*3})
                                    case '='
                                        if cmp ~=0
                                            badTRSP=1;
                                        end
                                    case '~='
                                        if cmp ==0
                                            badTRSP=1;
                                        end
                                    case '<'
                                        if cmp > -1
                                            badTRSP=1;
                                        end
                                    case '>'
                                        if cmp < 1
                                            badTRSP=1;
                                        end
                                    case '<='
                                        if cmp > 0
                                            badTRSP=1;
                                        end
                                    case '>='
                                        if cmp < 0
                                            badTRSP=1;
                                        end
                                    case 'starts'
                                        theStim=num2str(TRSPspecValues{TRSPspec});
                                        theCriterion=cellTable{iCell,EPmain.segment.numFixed+2+(theSpec-1)*3};
                                        if (length(theCriterion)>length(theStim)) || ~strcmp(theCriterion,theStim(1:length(theCriterion)))
                                            badTRSP=1;
                                        end
                                    case 'ends'
                                        theStim=num2str(TRSPspecValues{TRSPspec});
                                        theCriterion=cellTable{iCell,EPmain.segment.numFixed+2+(theSpec-1)*3};
                                        if (length(theCriterion)>length(theStim)) || ~strcmp(theCriterion,theStim(end-length(theCriterion)+1:end))
                                            badTRSP=1;
                                        end
                                    case 'contains'
                                        theStim=num2str(TRSPspecValues{TRSPspec});
                                        theCriterion=cellTable{iCell,EPmain.segment.numFixed+2+(theSpec-1)*3};
                                        if isempty(findstr(theCriterion,theStim))
                                            badTRSP=1;
                                        end
                                    otherwise
                                end
                            end
                            if badTRSP
                                keepTrial=0;
                            end
                        end
                    else
                        afterFollow=0;
                        afterPrecede=0;
                    end
                    if noStim && noTRSP
                        keepTrial=0;
                    end
                end%-follows- or -precedes- versus not
                if ~keepTrial
                    break %no need to check rest of specs
                end
            end%iSpec
            if ~keepTrial
                badEpochs=[badEpochs iEvent];
            else
                sectionNumberList=[sectionNumberList sectionNumber];
                sectionTrialList=[sectionTrialList sectionTrial];
                sectionTrial=sectionTrial+1;
                if ~isempty(stimTime) %add RT and ACC data for -responseCorr- and -responseErr-
                    allTRSPspecValues{iEvent}{find(strcmp('TS-ACC',TRSPspecNames))}=theAcc;
                    allTRSPspecValues{iEvent}{find(strcmp('TS-RT',TRSPspecNames))}=round((eventsString(stimEventList(iStim)).sample-stimTime)*(1000/inputData.Fs));
                end
            end
        end%iEvent
        eventIndex(badEpochs)=[];
        epochStart(badEpochs)=[];
        epochEnd(badEpochs)=[];
        epochRecTime(badEpochs)=[];
        allTRSPspecValues(badEpochs)=[];
        
        sectionNumHist=histc(sectionNumberList,unique(sectionNumberList));
        sectionNumTrials=[];
        for iSection=1:length(sectionNumHist)
            sectionNumTrials=[sectionNumTrials repmat(sectionNumHist(iSection),1,sectionNumHist(iSection))];
        end
        sectionTrialFromEndList=sectionNumTrials-sectionTrialList+1;
        
        %Check for section crits now that the sections have been fully defined
        badEpochs=[];
        for iEvent=1:length(eventIndex) %looping through events potentially defining trials
            ep_tictoc;if EPtictoc.stop;return;end
            badStim=0;
            for iSpec=1:length(goodSpecs)
                theSpec=goodSpecs(iSpec);
                theCrit=cellTable{iCell,EPmain.segment.numFixed+(iSpec-1)*3};
                theStim=cellTable{iCell,EPmain.segment.numFixed+2+(iSpec-1)*3};
                if ischar(theStim)
                    theStim=str2num(theStim);
                end
                if any(strcmp(theCrit,{'-secNum-','-secTrial-','-secTrialFromEnd-'}))
                    switch theCrit
                        case '-secNum-'
                            theValue=sectionNumberList(iEvent);
                        case '-secTrial-'
                            theValue=sectionTrialList(iEvent);
                        case '-secTrialFromEnd-'
                            theValue=sectionTrialFromEndList(iEvent);
                    end
                    switch deblank(cellTable{iCell,EPmain.segment.numFixed+1+(iSpec-1)*3})
                        case '='
                            if theValue ~= theStim
                                badStim=1;
                            end
                        case '~='
                            if theValue == theStim
                                badStim=1;
                            end
                        case '<'
                            if theValue >= theStim
                                badStim=1;
                            end
                        case '>'
                            if theValue <= theStim
                                badStim=1;
                            end
                        case '<='
                            if theValue > theStim
                                badStim=1;
                            end
                        case '>='
                            if theValue < theStim
                                badStim=1;
                            end
                    end
                end
            end
            if badStim
                badEpochs=[badEpochs iEvent];
            end
        end
        eventIndex(badEpochs)=[];
        epochStart(badEpochs)=[];
        epochEnd(badEpochs)=[];
        epochRecTime(badEpochs)=[];
        allTRSPspecValues(badEpochs)=[];
        sectionNumberList(badEpochs)=[];
        sectionTrialList(badEpochs)=[];
        sectionTrialFromEndList(badEpochs)=[];
        
        %if there is a TS-EPoffset field, then adjust the epochs accordingly.
        if any(strcmp('EPoffset',TRSPspecNamesList))
            offsetTRSP=find(strcmp('EPoffset',TRSPspecNamesList));
            disp('Adjusting the epoch timepoints according to the EPoffset values')
            badEpochs=[];
            for iEvent=1:length(eventIndex) %looping through events potentially defining trials
                theOffset=round(str2num(allTRSPspecValues{iEvent}{offsetTRSP})/(1000/inputData.Fs));
                if ~isempty(theOffset)
                    epochStart(iEvent)=epochStart(iEvent)+theOffset;
                    epochEnd(iEvent)=epochEnd(iEvent)+theOffset;
                    if (epochStart(iEvent) < 1) || (epochEnd(iEvent) > length(inputData.timeNames))
                        badEpochs(end+1)=iEvent;
                    end
                end
            end
            eventIndex(badEpochs)=[];
            epochStart(badEpochs)=[];
            epochEnd(badEpochs)=[];
            epochRecTime(badEpochs)=[];
            allTRSPspecValues(badEpochs)=[];
            sectionNumberList(badEpochs)=[];
            sectionTrialList(badEpochs)=[];
            sectionTrialFromEndList(badEpochs)=[];
        end
        
        cellCount(iFile,iCell)=length(eventIndex);
        if ~preview
            numCellTrials=length(eventIndex);
            if numCellTrials > 0
                if ~any(strcmp('secNum',TRSPspecNamesList))
                    TRSPspecNamesList{end+1}='secNum';
                end
                if ~any(strcmp('secTrial',TRSPspecNamesList))
                    TRSPspecNamesList{end+1}='secTrial';
                end
                if ~any(strcmp('secTrialFromEnd',TRSPspecNamesList))
                    TRSPspecNamesList{end+1}='secTrialFromEnd';
                end
                secNumTRSP=find(strcmp('secNum',TRSPspecNamesList));
                secTrialTRSP=find(strcmp('secTrial',TRSPspecNamesList));
                secTrialFromEndTRSP=find(strcmp('secTrialFromEnd',TRSPspecNamesList));
                sectionNumTrials=histc(sectionNumberList,unique(sectionNumberList));
                for iEpoch=1:numCellTrials
                    allTRSPspecValues{iEpoch}{secNumTRSP}=sectionNumberList(iEpoch);
                    allTRSPspecValues{iEpoch}{secTrialTRSP}=sectionTrialList(iEpoch);
                    allTRSPspecValues{iEpoch}{secTrialFromEndTRSP}=sectionTrialFromEndList(iEpoch);
                    if ~isempty(cellTable{iCell,7}) && ~any(strcmp('TS-ACC',cellTable(iCell,7)))
                        newACC=find(strcmp(cellTable{iCell,7},outTRSPnames));
                        oldACC=find(strcmp('TS-ACC',outTRSPnames));
                        if isempty(oldACC)
                            outTRSPnames{end+1,1}='TS-ACC';
                            allTRSPspecValues{iEpoch}{end+1}=[];
                            oldACC=length(outTRSPnames);
                        end
                        if ~isempty(newACC)
                            allTRSPspecValues{iEpoch}{oldACC}=allTRSPspecValues{iEpoch}{newACC};
                        else
                            theEventKeys=orderedEvents(eventIndex(iEpoch)).keys;
                            theKey=find(strcmp(cellTable{iCell,7},{theEventKeys.code}));
                            if ~isempty(theEventKeys) && ~isempty(theKey)
                                allTRSPspecValues{iEpoch}{oldACC}=str2double(theEventKeys(theKey).data);
                            else
                                allTRSPspecValues{iEpoch}{oldACC}=str2double(cellTable{iCell,7}); %enter in ACC directly as 0,1,2
                            end
                        end
                    end
                    if ~isempty(cellTable{iCell,8}) && ~any(strcmp('TS-RT',cellTable(iCell,8)))
                        newRT=find(strcmp(cellTable{iCell,8},outTRSPnames));
                        oldRT=find(strcmp('TS-RT',outTRSPnames));
                        if isempty(oldRT)
                            outTRSPnames{end+1,1}='TS-RT';
                            allTRSPspecValues{iEpoch}{end+1}=[];
                            oldRT=length(outTRSPnames);
                        end
                        if ~isempty(newRT)
                            allTRSPspecValues{iEpoch}{oldRT}=allTRSPspecValues{iEpoch}{newRT};
                        else
                            theEventKeys=orderedEvents(eventIndex(iEpoch)).keys;
                            theKey=find(strcmp(cellTable{iCell,8},{theEventKeys.code}));
                            if ~isempty(theEventKeys) && ~isempty(theKey)
                                allTRSPspecValues{iEpoch}{oldRT}=str2double(theEventKeys(theKey).data);
                            else
                                allTRSPspecValues{iEpoch}{oldRT}=str2double(cellTable{iCell,8}); %enter in RT directly
                            end
                        end
                    end
                end
                
                allEventSamples=round([orderedEvents.sample]);
                trialCountStart=0;
                sameCellName=find(strcmp(cellTable{iCell,1},cellTable(1:iCell-1,1)));
                if ~isempty(sameCellName)
                    trialCountStart=sum(cellCount(iFile,sameCellName));
                end
                if flexMode
                    badPoints=find(kron([inputData.analysis.badTrials(1,:)],ones(1,round(inputData.Fs))));
                    badChanPoints=cell(0);
                    if strcmp(inputData.dataType,'continuous')
                        for iChan=1:length(inputData.chanNames)
                            badChanPoints{iChan}=find(kron([inputData.analysis.badChans(1,:,iChan)==-1],ones(1,round(inputData.Fs))));
                        end
                    else
                        for iChan=1:length(inputData.chanNames)
                            badChanPoints{iChan}=find(kron([inputData.analysis.badChans(1,:,iChan)==-1],ones(1,numPoints)));
                        end
                    end
                end
                badTrials{iCell}=zeros(1,numCellTrials);
                badChans{iCell}=zeros(1,numCellTrials,length(inputData.chanNames));
                for iEpoch=1:numCellTrials
                    if flexMode
                        deltaT=floor((epochEnd(iEpoch)-epochStart(iEpoch)+1)/numInterPoints);
                        if strcmp(inputData.dataType,'continuous')
                            goodPoints=epochStart(iEpoch):epochEnd(iEpoch);
                            goodPoints=setdiff(goodPoints,badPoints);
                        else
                            theInTrial=(floor(epochStart(iEpoch)-1)/numPoints)+1;
                            if inputData.analysis.badTrials(1,theInTrial)
                                goodPoints=[];
                            else
                                goodPoints=epochStart(iEpoch):epochEnd(iEpoch);
                            end
                        end
                        outputData.data(:,:,end+1,1,1,1)=0;
                        if isempty(goodPoints)
                            badTrials{iCell}(1,iEpoch)=1;
                        else
                            for iChan=1:length(inputData.chanNames)
                                goodChanPoints=setdiff(goodPoints,badChanPoints{iChan});
                                for iPoint=1:numInterPoints
                                    epochPoints=intersect([(iPoint-1)*deltaT+epochStart(iEpoch):iPoint*deltaT+epochStart(iEpoch)-1],goodChanPoints);
                                    if ~isempty(epochPoints)
                                        if ~isempty(inputData.freqNames) && any(strcmp(inputData.chanTypes{iChan},{'EEG','REG'}))
                                            %switch to amplitude scaling when adding freq data together.
                                            outputData.data(iChan,iPoint,end,1,:,:)=squeeze(mean(abs(inputData.data(iChan,epochPoints,1,1,:,:)),2));
                                        else
                                            outputData.data(iChan,iPoint,end,1,:,:)=squeeze(mean(inputData.data(iChan,epochPoints,1,1,:,:),2));
                                        end
                                    else
                                        outputData.data(:,iPoint,end,1,1,1)=0;
                                        badChans{iCell}(1,iEpoch,iChan)=-1;
                                    end
                                end
                            end
                        end
                        outputData.events{1,end+1}=orderedEvents(find((allEventSamples>=epochStart(iEpoch)) & (allEventSamples <= epochEnd(iEpoch))));
                        for iEvent=1:length(outputData.events{1,end})
                            outputData.events{1,end}(iEvent).sample=floor((outputData.events{1,end}(iEvent).sample-epochStart(iEpoch)+1-1)/numInterPoints)+1; %change event time to be relative to epoch
                        end
                    else
                        startSample=round(epochStart(iEpoch));
                        endSample=round(epochEnd(iEpoch));
                        epochPoints=[startSample:endSample]';
                        if startSample==1
                            epochPoints=[startSample; epochPoints];
                        else
                            epochPoints=[startSample-1; epochPoints];
                        end
                        if endSample==numPoints
                            epochPoints=[epochPoints; endSample];
                        else
                            epochPoints=[epochPoints; endSample+1];
                        end
                        outputData.data(:,:,end+1,1,:,:)=interp1([startSample-1 startSample:endSample endSample+1],inputData.data(:,epochPoints,1,1,:,:)',[epochStart(iEpoch):epochEnd(iEpoch)])';
                        if ~isempty(inputData.covAVE)
                            outputData.covAVE(:,:,end+1,1,:,:,:,:)=interp1([startSample-1 startSample:endSample endSample+1],sqrt(inputData.covAVE(:,epochPoints,1,1,:,:)'),[epochStart(iEpoch):epochEnd(iEpoch)])'.^2;
                        end
                        outputData.events{1,end+1}=orderedEvents(find((allEventSamples>=epochStart(iEpoch)) & (allEventSamples <= epochStart(iEpoch)+epochLength-1)));
                        for iEvent=1:length(outputData.events{1,end})
                            outputData.events{1,end}(iEvent).sample=outputData.events{1,end}(iEvent).sample-epochStart(iEpoch)+1; %change event time to be relative to epoch
                        end
                        if strcmp(inputData.dataType,'continuous')
                            theInTrial=floor((epochStart(iEpoch)-1)/inputData.Fs)+1;
                        else
                            theInTrial=floor((epochStart(iEpoch)-1)/numPoints)+1;
                        end
                        badTrials{iCell}(1,iEpoch)=inputData.analysis.badTrials(1,theInTrial);
                        badChans{iCell}(1,iEpoch,:)=inputData.analysis.badChans(1,theInTrial,:);
                    end
                    outputData.cellNames{end+1,1}=cellTable{iCell,1};
                    outputData.cellTypes{end+1,1}='SGL';
                    outputData.trialNames(end+1,1)=iEpoch+trialCountStart;
                    if isempty(TRSPspecNamesList)
                        outputData.trialSpecs{end+1,1}=[];
                    else
                        for iSpec=1:length(allTRSPspecValues{iEpoch})
                            if iSpec==1 %so also a new trial
                                outputData.trialSpecs{end+1,iSpec}=allTRSPspecValues{iEpoch}{iSpec};
                            else
                                outputData.trialSpecs{end,iSpec}=allTRSPspecValues{iEpoch}{iSpec};
                            end
                        end
                    end
                    if ~isempty(taskLabel)
                        theSpec=find(strcmp('EPtask',outTRSPnames));
                        if isempty(theSpec)
                            outTRSPnames{end+1}='EPtask';
                            theSpec=length(outTRSPnames);
                        end
                        outputData.trialSpecs{end,theSpec}=taskLabel;
                    end
                    outputData.recTime(end+1,1)=epochRecTime(iEpoch);
                end
            end
        end
    end%iCell
    if ~preview
        numTrials=length(outputData.cellNames);
        if numTrials==0
            disp(['Error: Segmented file for ' Name ' did not result in any segments and therefore will not be saved.']);
        else
            for iSpec=length(outTRSPnames):-1:1
                startChar=1;
                if strcmpi(outTRSPnames{iSpec}(1:3),'TS-') && ~any(strcmpi(outTRSPnames{iSpec}(4:end),outTRSPnames))
                    startChar=4;
                elseif any(strcmpi(outTRSPnames{iSpec}(1:4),{'EPM-','TSP-'}))
                    startChar=5;
                end
                outTRSPnames{iSpec}=outTRSPnames{iSpec}(startChar:end);
            end
            outputData.trialSpecNames=outTRSPnames;
            outputData.avgNum=ones(1,numTrials);
            outputData.subNum=ones(1,numTrials);
            outputData.covNum=ones(1,numTrials);
%             outputData.history{end+1}={'ep_segmentData',cellTable,inputFiles,importFormat,outputFormat,preview};
            outputData.analysis.blinkTrial=zeros(1,numTrials);
            outputData.analysis.saccadeTrial=zeros(1,numTrials);
            outputData.analysis.saccadeOnset=zeros(1,numTrials);
            outputData.analysis.moveTrial=zeros(1,numTrials);
            outputData.analysis.badTrials=zeros(1,numTrials);
            outputData.analysis.badChans=zeros(1,numTrials,numChans);
            count=0;
            for iCell=1:numOutCells
                outputData.analysis.badTrials(1,count+1:count+length(badTrials{iCell}))=badTrials{iCell};
                outputData.analysis.badChans(1,count+1:count+length(badTrials{iCell}),:)=badChans{iCell};
                count=count+length(badTrials{iCell});
            end
            if isempty(outputData.trialSpecNames)
                outputData.trialSpecs=cell(numTrials,0);
            end
            outputData.stims=struct('name',{},'image',{},'AOI',{});
            for iStim=1:length(inputData.stims)
                stimFlag=0;
                for iCell=1:length(outputData.cellNames)
                    for iEvent=1:length(outputData.events{iCell})
                        if any(strcmp(inputData.stims(iStim).name,{outputData.events{iCell}(iEvent).keys.data}))
                            stimFlag=1;
                        end
                    end
                end
                if stimFlag
                    outputData.stims(end+1)=inputData.stims(iStim); %keep only stim images whose events are still in the segmented data.
                end
            end
            
            [err]=ep_checkEPfile(outputData);
            if err
                disp(['Error: Segmented file for ' Name ' did not pass data integrity checks and therefore will not be saved.']);
            else
                [pathstr, name, ext] = fileparts(inputFiles{iFile});
                [fileSuffix,formatName]=ep_fileExtensions(outputFormat);
                
                sameName=1;
                theNumber=0;
                fileNameSuffix=[pathstr filesep name segSuffix fileSuffix{1}];
                while sameName
                    sameName=0;
                    if exist(fileNameSuffix,'file')
                        sameName=1;
                    end
                    if sameName
                        theNumber=theNumber+1;
                        fileNameSuffix=[pathstr filesep name segSuffix '-' num2str(theNumber) fileSuffix{1}];
                    end
                end
                if ~strcmp(importFormat,'ep_mat')
                    theDescription=['Imported the file ' inputFiles{iFile} '.'];
                    outputData.history=ep_addHistory(outputData.history, EPmain.preferences.records, theDescription, EPmain.handles.hMainWindow, EPmain.preferences, EPmain.preferences.EPver,inputFiles);
                end
                theDescription=['Segmented the file ' inputFiles{iFile} '.'];
                outputData.history=ep_addHistory(outputData.history, EPmain.preferences.records, theDescription, EPmain.handles.hMainWindow, EPmain.preferences, EPmain.preferences.EPver,inputFiles);
                [err]=ep_writeData(outputData,fileNameSuffix,EPmain.preferences.general.specSuffix,EPmain.preferences.general.subjectSpecSuffix,outputFormat);
            end
        end
    end
end
cellNums=cellCount;





