function [EPdataOut msgLog]=ep_readEyelink(EPdataIn,edfStruct,matchTable)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [EPdataOut msgLog]=ep_readEyelink(EPdataIn,edfStruct,matchTable);
% Reads in text output from EyeLink eye tracker and adds the sample-by-sample information to the data field and the events
% to the events field.  Assumes that the EEG is a continuous data file. Interpolates the EyeLink data to correct for uneven
% temporal sampling, using matching EEG-EyeLink events as the anchor points.
% Uses the mean of the two eye coordinates.  Uses the POR X [px] columns of the EyeLink files.
% Counts as fixation start run of samples where either eye registers as being fixations.
% Relies on the first such sample as fixation coordinate, ignoring subsequent drift.
%
%Inputs
%   EPdataIn:  Structured array with the data and accompanying information.  See readData.
%   edfStruct: Structured data from edfmex
%      .time   sample times in ms.
%      .gx:    eye gaze x-position (eye,time)
%      .gy:    eye gaze y-position (eye,time)
%      .pa:    pupil dilation (eye,ms)
%   matchTable: Table with first column of unique EEG event values and a second table of the EyeLink events they correspond
%   to, if any.  If they do not correspond, then the value will be 'none'.
%
%Outputs
%   EPdataOut: Structured array with the data and accompanying information.  See readData.
%   msgLog: messages.
%
% Whereas SMI recorded separate blink events for each eye, EyeLink is set
% for monocular (left eye) recording by default.
% In one recording, I see two STARTFIX events in a row.  The next ENDFIX
% event has the start time of the second of them.  The FIX
% and SACC periods appear to be exclusive of each other and start
% immediately (one ms later) after the ending of the prior period.
% The documentation indicates that blink events are always bracketed by
% artifactual saccade events that result from the distortion of eye
% position data around the blink period.  They should therefore be ignored
% as not really being blinks.  The blink period itself should be timed
% either according to the blinks themselves or by the artifactual saccade
% period.  The manual also suggests that fixation periods of 100ms or less
% before or after blinks should be ignored.  Will examine issue later.
% Unlike SMI data, the EyeLink data has a reliable 1000Hz sampling rate.
% Also, right now only the event times are being drift corrected.  EEG needs
% drift correcting as well?
%
% History:
%
% by Joseph Dien (12/31/24)
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

global EPtictoc

EPdataOut=[];
msgLog=cell(0);

if isempty(matchTable)
    disp('Error: no match table specified.');
    return
end

if size(matchTable,2) ~= 2
    disp('Error: match table defective.');
    return
end

if isempty(EPdataIn)
    disp('Error: no EEG data file.');
    return
end

if ~strcmp(EPdataIn.dataType,'continuous')
    disp('Error: EEG data file is not continuous.');
    return
end

tic

%match up the EEG and SMI events and add the stimulus info to the events
matchTable=matchTable(~strcmp('none',matchTable(:,2)),:);%drop EEG events for which there are no matching SMI events
EEGevents={EPdataIn.events{1}.value}';
ELeventsList=edfStruct.FEVENT;
ELdata.time=edfStruct.FSAMPLE.time;
ELdata.gx=edfStruct.FSAMPLE.gx;
ELdata.gy=edfStruct.FSAMPLE.gy;
ELdata.pa=edfStruct.FSAMPLE.pa;
EEGsampLength=(1000/EPdataIn.Fs);
numpoints=length(EPdataIn.timeNames);

%prepare to correct for timing synch drift
MSGrows=find(strcmp({ELeventsList.codestring},'MESSAGEEVENT')); %EL events that are potentially paired with EEG events.
EEGtimingList=[]; %list of timing events (provided by match table) in EEG events array
ELtimingList=[]; %list of corresponding rows in EyeLink event array
for iMSG=1:length(MSGrows)
    theMSG=MSGrows(iMSG);
    if any(strcmp(ELeventsList(theMSG).message,matchTable(:,2)))
        matchEEG=find(strcmp(matchTable{find(strcmp(ELeventsList(theMSG).message,matchTable(:,2))),1},EEGevents));
        if isempty(EEGtimingList)
            theEEGevent=matchEEG(1); %first EEG event
        else
            theEEGevent=matchEEG(min(find(matchEEG>EEGtimingList(end)))); %first EEG event that comes after what is already in the EEG list
        end

        if isempty(theEEGevent)
            msg='Warning: EEG-EyeLink event mismatch.  Have run out of EyeLink events to match up with.  Will integrate EyeLink events up to this point.  Be sure to verify the event matching.';
            msgLog{end+1}=msg;
            disp(msg);
            EEGtimingList=[];
            break
        elseif ~isempty(EEGtimingList) && (abs(EPdataIn.events{1}(theEEGevent).sample*EEGsampLength-double(ELeventsList(theMSG).sttime))-abs(EPdataIn.events{1}(EEGtimingList(1)).sample*EEGsampLength-double(ELeventsList(ELtimingList(1)).sttime))) > 1000
            disp(['Timing mismatch: ' num2str(abs(EPdataIn.events{1}(theEEGevent).sample*EEGsampLength-double(ELeventsList(theMSG).sttime))-abs(EPdataIn.events{1}(EEGtimingList(1)).sample*EEGsampLength-double(ELeventsList(ELtimingList(1)).sttime)))])
            disp(['Dropping event: ' num2str(iMSG)])
        else
            ELtimingList(end+1)=theMSG;
            EEGtimingList(end+1)=theEEGevent;
        end
    end
end

if length(ELtimingList) ~= length(EEGtimingList)
    disp('Programming error: contact the developer Joseph Dien for assistance.');
    return
end

if isempty(EEGtimingList)
    disp('Error: no timing events from match table present.  Aborting import of eye-tracker data.')
    return
end

driftTimes=[EPdataIn.events{1}(EEGtimingList).sample]*EEGsampLength-double([ELeventsList(ELtimingList).sttime]);
driftTimes=driftTimes-driftTimes(1);

for iMSG=1:size(matchTable,1)
    disp('')
    disp(['Mean drift ' matchTable{iMSG,2} ':' num2str(mean(driftTimes(strcmp({ELeventsList(ELtimingList).message},matchTable{iMSG,2}))))]);
end
% disp(['Synch drift by end of recording was: ' num2str(driftTimes(end)-driftTimes(1))])

%convert EyeLink times to have same time zero as EEG times.
offsetTime=ELeventsList(ELtimingList(1)).sttime-(EPdataIn.events{1}(EEGtimingList(1)).sample-1)*EEGsampLength-EPdataIn.timeNames(1);
for iEvent=1:length(ELeventsList)
    ELeventsList(iEvent).sttime=ELeventsList(iEvent).sttime-offsetTime;
    ELeventsList(iEvent).entime=ELeventsList(iEvent).entime-offsetTime;
end
for iTime=1:length(ELdata.time)
    ELdata.time(iTime)=ELdata.time(iTime)-offsetTime;
end

%correct EyeLink times for synch drift (arbitrarily treating EEG times as correct).
driftBeta=double([ELeventsList(ELtimingList).sttime])/([EPdataIn.events{1}(EEGtimingList).sample]*EEGsampLength);
ELdata.time=double(ELdata.time)*driftBeta;
for iEvent=1:length(ELeventsList)
    ELeventsList(iEvent).sttime=double(ELeventsList(iEvent).sttime)*driftBeta;
    ELeventsList(iEvent).entime=double(ELeventsList(iEvent).entime)*driftBeta;
end

%combine the data from the two eyes
ELdata.gx(1,ELdata.gx(1,:)==-32768)=NaN;
ELdata.gx(2,ELdata.gx(2,:)==-32768)=NaN;
ELdata.gy(1,ELdata.gy(1,:)==-32768)=NaN;
ELdata.gy(2,ELdata.gy(2,:)==-32768)=NaN;
ELdata.pa(1,ELdata.pa(1,:)==-32768)=NaN;
ELdata.pa(2,ELdata.pa(2,:)==-32768)=NaN;

ELdata.gx=mean(ELdata.gx,1,'omitnan');
ELdata.gy=mean(ELdata.gy,1,'omitnan');
ELdata.pa=mean(ELdata.pa,1,'omitnan');

%resample EyeLink data so that it is consistent with the EEG data
ELdata.gx=interp1(ELdata.time,ELdata.gx,EPdataIn.timeNames);
ELdata.gy=interp1(ELdata.time,ELdata.gy,EPdataIn.timeNames);
ELdata.pa=interp1(ELdata.time,ELdata.pa,EPdataIn.timeNames);
ELdata.time=EPdataIn.timeNames;

%remove artifactual saccade events around blinks and expand blinks to cover
%these false saccades as well as being partial blinks.
blinkStarts=find(strcmp('STARTBLINK ',{ELeventsList.codestring}));
for iBlink=1:length(blinkStarts)
    blinkStart=blinkStarts(iBlink);
    endBlinks=blinkStart+find(strcmp('ENDBLINK',{ELeventsList(blinkStart:end).codestring}))-1;
    if isempty(endBlinks)
        disp('Missing ENDBLINK');
    else
        blinkEnd=min(find(strcmp('ENDBLINK',{ELeventsList(blinkStart:endBlinks(1)).codestring})))+blinkStart-1;
        if ELeventsList(blinkStart).sttime ~= ELeventsList(blinkEnd).sttime
            disp('Mismatched ENDBLINK time');
        else
            %drop any fixation or saccade events within a blink period
            if (blinkEnd-blinkStart) > 1
                for iEvent=blinkStart:blinkEnd
                    if any(strcmp(edfStruct(iEvent).FEVENT.codestring,{'STARTFIX';'ENDFIX';'STARTSACC';'ENDSACC'}))
                        ELeventsList(iEvent).codestring='delete';
                    end
                end
            end
            %drop artifactual saccades bracketing blinks and expand the
            %blink period out to cover these partial blink periods
            falseSaccStart=max(find(strcmp('STARTSACC',{ELeventsList(1:blinkStart).codestring})));
            falseSaccEnd=blinkEnd+min(find(strcmp('ENDSACC',{ELeventsList(blinkEnd:end).codestring})))-1;
            if ~isempty(falseSaccStart) && ~isempty(falseSaccEnd) && ELeventsList(blinkStart).sttime ~= ELeventsList(blinkEnd).sttime
                disp('Mismatched false saccade times');
            else
                if ~isempty(falseSaccStart)
                    ELeventsList(blinkStart).codestring='delete';
                    ELeventsList(falseSaccStart).codestring='STARTBLINK ';
                end
                if ~isempty(falseSaccEnd)
                    ELeventsList(blinkEnd).codestring='delete';
                    ELeventsList(falseSaccEnd).codestring='ENDBLINK';
                end
            end
        end
    end
end

ELeventsList(strcmp('delete',{ELeventsList.codestring}))=[];

%Add saccade events
saccadeStarts=find(strcmp('STARTSACC',{ELeventsList.codestring}));
for iSaccade=1:length(saccadeStarts)
    %calculate saccade direction based on 50 ms following each saccade
    sacTime=ELeventsList(saccadeStarts(iSaccade)).sttime;
    sacSamples=find(([ELdata.time] <= sacTime+50) & ([ELdata.time] >= sacTime));
    sacX=ELdata.gx;
    sacY=ELdata.gy;
    warning('off','MATLAB:rankDeficientMatrix')
    Bx=[ones(length(sacSamples),1),sacSamples]\sacX(sacSamples); %increases rightward
    By=[ones(length(sacSamples),1),sacSamples]\(-sacY(sacSamples)); %increases downward
    warning('on','MATLAB:rankDeficientMatrix')
    sacAng=atand(By(2)/Bx(2)); %will be 360 degrees with zero upwards

    if (Bx(2)==0) || (By(2)==0)
        if (Bx(2)==0) && (By(2)==0)
            sacAng=NaN;
        elseif (Bx(2)==0) && (By(2)~=0)
            if By(2)>0
                sacAng=0;
            else
                sacAng=180;
            end
        elseif (Bx(2)~=0) && (By(2)==0)
            if Bx(2)>0
                sacAng=90;
            else
                sacAng=270;
            end
        end
    else
        if (Bx(2)>0) && (By(2)>0)
            sacAng=90-sacAng;
        elseif (Bx(2)<0) && (By(2)>0)
            sacAng=270-sacAng;
        elseif (Bx(2)>0) && (By(2)<0)
            sacAng=90-sacAng;
        elseif (Bx(2)<0) && (By(2)<0)
            sacAng=270-sacAng;
        end
    end

    EPdataIn.events{1}(end+1).type='eye-tracker';
    EPdataIn.events{1}(end).sample=ceil(sacTime/EEGsampLength);
    EPdataIn.events{1}(end).value='saccadeET';
    EPdataIn.events{1}(end).duration=1;
    EPdataIn.events{1}(end).keys(1).code='ELevent';
    EPdataIn.events{1}(end).keys(1).datatype='short';
    EPdataIn.events{1}(end).keys(1).data=num2str(saccadeStarts(iSaccade));
    EPdataIn.events{1}(end).keys(1).description='';
    EPdataIn.events{1}(end).keys(2).code='angle';
    EPdataIn.events{1}(end).keys(2).datatype='char';
    EPdataIn.events{1}(end).keys(2).data=num2str(round(sacAng));
    EPdataIn.events{1}(end).keys(2).description='';
end

%Add fixation events
fixationStarts=find(strcmp('STARTFIX',{ELeventsList.codestring}));
for iFixation=1:length(fixationStarts)
    fixTime=ELeventsList(fixationStarts(iFixation)).sttime;
    EPdataIn.events{1}(end+1).type='eye-tracker';
    EPdataIn.events{1}(end).sample=ceil(fixTime/EEGsampLength);
    EPdataIn.events{1}(end).value='fixationET';
    EPdataIn.events{1}(end).duration=1;
    EPdataIn.events{1}(end).keys(1).code='ELevent';
    EPdataIn.events{1}(end).keys(1).datatype='short';
    EPdataIn.events{1}(end).keys(1).data=num2str(fixationStarts(iFixation));
    EPdataIn.events{1}(end).keys(1).description='';
end

%Add blink events
blinkStarts=find(strcmp('STARTBLINK ',{ELeventsList.codestring}));
blinkEnds=find(strcmp('ENDBLINK',{ELeventsList.codestring}));
if length(blinkStarts) ~= length(blinkEnds)
    disp('Mismatched blink events');
else
    for iBlink=1:length(blinkStarts)
        blinkStartSamp=ceil(ELeventsList(blinkStarts(iBlink)).sttime/EEGsampLength);
        EPdataIn.events{1}(end+1).type='eye-tracker';
        EPdataIn.events{1}(end).sample=blinkStartSamp;
        EPdataIn.events{1}(end).value='blinkStartET';
        EPdataIn.events{1}(end).duration=1;
        EPdataIn.events{1}(end).keys(1).code='ELevent';
        EPdataIn.events{1}(end).keys(1).datatype='short';
        EPdataIn.events{1}(end).keys(1).data=num2str(blinkStarts(iBlink));
        EPdataIn.events{1}(end).keys(1).description='';

        blinkEndSamp=ceil(ELeventsList(blinkEnds(iBlink)).entime/EEGsampLength);
        EPdataIn.events{1}(end+1).type='eye-tracker';
        EPdataIn.events{1}(end).sample=blinkEndSamp;
        EPdataIn.events{1}(end).value='blinkEndET';
        EPdataIn.events{1}(end).duration=1;
        EPdataIn.events{1}(end).keys(1).code='ELevent';
        EPdataIn.events{1}(end).keys(1).datatype='short';
        EPdataIn.events{1}(end).keys(1).data=num2str(blinkEnds(iBlink));
        EPdataIn.events{1}(end).keys(1).description='';

        ELdata.pa(blinkStartSamp:blinkEndSamp)=NaN;
        ELdata.gx(blinkStartSamp:blinkEndSamp)=NaN;
        ELdata.gy(blinkStartSamp:blinkEndSamp)=NaN;
    end
end

%add the eye-tracker event information to the EEG dataset
EPadd=[];
EPadd.chanNames={'pupil';'x-eye';'y-eye'};
EPadd.chanTypes={'PPL';'XEY';'YEY'};
[EPdataIn]=ep_addData(EPdataIn,EPadd,'channels');
if isempty(EPdataIn)
    return;
end
EPdataIn.data(end-2,:,:,:,:,:,:)=ELdata.pa;
EPdataIn.data(end-1,:,:,:,:,:,:)=ELdata.gx;
EPdataIn.data(end,:,:,:,:,:,:)=ELdata.gy;

EPdataOut=EPdataIn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from Matlab documentation (8/26/19) at https://www.mathworks.com/help/matlab/ref/xmlread.html
% and covered by their copyrights.

function theStruct = parseXML(filename)
% PARSEXML Convert XML file to a MATLAB structure.
try
   tree = xmlread(filename);
catch
   error('Failed to read XML file %s.',filename);
end

% Recurse over child nodes. This could run into problems 
% with very deeply nested trees.
try
   theStruct = parseChildNodes(tree);
catch
   error('Unable to parse XML file %s.',filename);
end


% ----- Local function PARSECHILDNODES -----
function children = parseChildNodes(theNode)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
   childNodes = theNode.getChildNodes;
   numChildNodes = childNodes.getLength;
   allocCell = cell(1, numChildNodes);

   children = struct(             ...
      'Name', allocCell, 'Attributes', allocCell,    ...
      'Data', allocCell, 'Children', allocCell);

    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
    end
end

% ----- Local function MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
% Create structure of node info.

nodeStruct = struct(                        ...
   'Name', char(theNode.getNodeName),       ...
   'Attributes', parseAttributes(theNode),  ...
   'Data', '',                              ...
   'Children', parseChildNodes(theNode));

if any(strcmp(methods(theNode), 'getData'))
   nodeStruct.Data = char(theNode.getData); 
else
   nodeStruct.Data = '';
end

% ----- Local function PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = [];
if theNode.hasAttributes
   theAttributes = theNode.getAttributes;
   numAttributes = theAttributes.getLength;
   allocCell = cell(1, numAttributes);
   attributes = struct('Name', allocCell, 'Value', ...
                       allocCell);

   for count = 1:numAttributes
      attrib = theAttributes.item(count-1);
      attributes(count).Name = char(attrib.getName);
      attributes(count).Value = char(attrib.getValue);
   end
end