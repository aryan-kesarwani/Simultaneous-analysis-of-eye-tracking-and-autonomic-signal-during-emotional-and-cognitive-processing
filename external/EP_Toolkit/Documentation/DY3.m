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

addpath('trigger');

% Clear the workspace
close all;
clearvars;
sca;

%experiment parameters
latinSquare={'aAAbBBdDDcCC';'bBBcCCaAAdDD';'cCCdDDbBBaAA';'dDDaAAcCCbBB'}; %pseudo-Latin square design for four tasks.
expName='DY3';
taskNames={'Naming';'Masked';'Priming';'Rhyming'};
numTasks=length(taskNames);
numBlocks=length(latinSquare{1});
trialCounter=zeros(numTasks,1);

% Setup PTB with some default values
PsychDefaultSetup(2);

% Seed the random number generator. Here we use the an older way to be
% compatible with older systems. Newer syntax would be rng('shuffle'). Look
% at the help function of rand "help rand" for more information
rand('seed', sum(100 * clock));

% Set the screen number to the external secondary monitor if there is one
% connected
info.screenNumber = max(Screen('Screens'));

% Define colors
info.white = WhiteIndex(info.screenNumber);
info.grey = info.white / 2;
info.black = BlackIndex(info.screenNumber);
info.red=[1 0 0];
info.green= [0 1 0];
info.blue= [0 0 1];

% defining trigger port
config_io;
% optional step: verify that the inpoutx64 driver was successfully initialized
global cogent;
if( cogent.io.status ~= 0 )
    error('inp/outp installation failed');
end

% write a value to the default LPT1 printer output port (at 0x378)
info.address = hex2dec('2040');
info.sleepDur = .015;
info.CMU=0;

%Get session information and calculate counterbalances

sessionInfo=[];
doneFlag=0;
while ~doneFlag
    sessionInfo.subject=input('Subject?','s');
    sessionInfo.subject=floor(str2double(sessionInfo.subject));
    if isnumeric(sessionInfo.subject) && (sessionInfo.subject>0)
        doneFlag=1;
    end
end

doneFlag=0;
while ~doneFlag
    sessionInfo.gender=input('Gender (m or f or o)?','s');
    if any(strcmp(sessionInfo.gender,{'f','m','o'}))
        doneFlag=1;
    end
end

doneFlag=0;
while ~doneFlag
    sessionInfo.age=input('Age?','s');
    if isnumeric(str2double(sessionInfo.age))
        sessionInfo.age=str2double(sessionInfo.age);
        doneFlag=1;
    end
end

sessionInfo.session='';
sessionInfo.handedness='';
theTime=clock;
sessionInfo.time=[num2str(theTime(4)) ':' num2str(theTime(5)) ':' num2str(round(theTime(6)))];
sessionInfo.date=date;
sessionInfo.researcher='';
sessionInfo.group=[];
sessionInfo.expName=expName;

fingerCB=rem(sessionInfo.subject-1,2)+1; %1 or 2
stimCB=rem(floor((sessionInfo.subject-1)/2),2)+1; %1 or 2
orderCB=rem(floor((sessionInfo.subject-1)/4),4)+1; %1-4

disp(['Finger counterbalance:' num2str(fingerCB)])
disp(['Stimulus counterbalance:' num2str(stimCB)])
disp(['Order counterbalance:' num2str(orderCB)])

info.fingerCB=fingerCB;

for iBlock=1:length(latinSquare{orderCB})
    taskCode=latinSquare{orderCB}(iBlock);
    if double(taskCode) > 96
        taskNum=double(taskCode)-96;
        blockName=['practice ' taskNames{taskNum}];
    else
        taskNum=double(taskCode)-64;
        blockName=taskNames{taskNum};
    end
    disp([sprintf('%02d',iBlock) ': ' blockName]);
end

doneFlag=0;
while ~doneFlag
    startBlock=input(['Start block number (1-' num2str(numBlocks) '), type 1 unless restarting experiment?'],'s');
    startBlock=floor(str2double(startBlock));
    if isnumeric(startBlock) && (startBlock>0) && (startBlock<=numBlocks)
        doneFlag=1;
    end
end

%load in the stimulus lists
stimLists=cell(numTasks,1);
stimHeaders=cell(numTasks,1);
for iTask=1:numTasks
    disp([taskNames{iTask} '-' num2str(stimCB) '.txt'])
    fid=fopen([taskNames{iTask} '-' num2str(stimCB) '.txt'],'r');
    tempVar=fgetl(fid);
    delim='\t';
    if ~isempty(strfind(tempVar,',')) && isempty(strfind(tempVar,'\t'))
        delim=','; %if there are commas and no tabs, assume it is a comma-delimited file.
    elseif ~isempty(strfind(tempVar,' ')) && isempty(strfind(tempVar,'\t'))
        delim=' '; %if there are spaces and no tabs, assume it is a space-delimited file.
    end
    
    numcols=length(regexp(tempVar,delim))+1; %determine number of columns based on tab markers
    if regexp(tempVar,[delim '$'])
        numcols=numcols-1; %if there is an extra tab at the end of the line, drop it.
    end
    
    frewind(fid);
    theData=textscan(fid, [repmat('%s',1,numcols)],'Delimiter',delim);
    fclose(fid);
    numSpecs=length(theData);
    theHeaders=cell(numSpecs,1);
    trialData=cell(length(theData{1})-1,numSpecs);
    for iSpec=1:numSpecs
        theHeaders{iSpec}=theData{iSpec}{1};
        trialData(:,iSpec)=theData{iSpec}(2:end);
    end
    stimLists{iTask}=trialData;
    stimHeaders{iTask}=theHeaders;
    numTrials=0;
    for iRow=1:size(trialData,1)
        numTrials=numTrials+str2double(trialData{iRow,find(strcmp('numReps',stimHeaders{iTask}))});
    end
    numReps=length(find(findstr(char(65+iTask),latinSquare{orderCB})));
    if rem(numTrials,numReps) >0
        disp(['number of trials is not an even multiple of number of blocks for task ' taskNames{iTask}]);
        keyboard
    end
end

%Generate the trial sequences
trialSequences=cell(numTasks,1);
trialHeaders=cell(numTasks,1);
for iTask=1:numTasks
    repsCol=find(strcmp('numReps',stimHeaders{iTask}));
    if isempty(repsCol)
        trialSequences{iTask}=stimLists{iTask};
        trialHeaders{iTask}=stimHeaders{iTask};
    else
        numRows=size(stimLists{iTask},1);
        for iRow=1:numRows
            numReps=str2double(stimLists{iTask}{iRow,repsCol});
            for iRep=1:numReps
                trialSequences{iTask}(end+1,:)=stimLists{iTask}(iRow,[1:repsCol-1 repsCol+1:end]);
            end
        end
        trialHeaders{iTask}=stimHeaders{iTask}([1:repsCol-1 repsCol+1:end]);
        trialSequences{iTask}=Shuffle(trialSequences{iTask},2);
    end
end

%load in the practice lists
pracLists=cell(numTasks,1);
pracHeaders=cell(numTasks,1);
for iTask=1:numTasks
    fid=fopen([taskNames{iTask} 'Prac.txt'],'r');
    tempVar=fgetl(fid);
    delim='\t';
    if ~isempty(strfind(tempVar,',')) && isempty(strfind(tempVar,'\t'))
        delim=','; %if there are commas and no tabs, assume it is a comma-delimited file.
    elseif ~isempty(strfind(tempVar,' ')) && isempty(strfind(tempVar,'\t'))
        delim=' '; %if there are spaces and no tabs, assume it is a space-delimited file.
    end
    
    numcols=length(regexp(tempVar,delim))+1; %determine number of columns based on tab markers
    if regexp(tempVar,[delim '$'])
        numcols=numcols-1; %if there is an extra tab at the end of the line, drop it.
    end
    
    frewind(fid);
    theData=textscan(fid, [repmat('%s',1,numcols)],'Delimiter',delim);
    fclose(fid);
    numSpecs=length(theData);
    theHeaders=cell(numSpecs,1);
    trialData=cell(length(theData{1})-1,numSpecs);
    for iSpec=1:numSpecs
        theHeaders{iSpec}=theData{iSpec}{1};
        trialData(:,iSpec)=theData{iSpec}(2:end);
    end
    pracLists{iTask}=trialData;
    pracHeaders{iTask}=theHeaders;
end

try
    
    % Open the screen
    [info.window, info.windowRect] = PsychImaging('OpenWindow', info.screenNumber, info.grey, [], 32, 2);
    
    % Flip to clear
    Screen('Flip', info.window);
    
    % Hide cursor
    HideCursor;
    
    % Query the frame duration
    info.ifi = Screen('GetFlipInterval', info.window);
    
    % Set the text size
    Screen('TextSize', info.window, 24);
    
    % Query the maximum priority level
    topPriorityLevel = MaxPriority(info.window);
    
    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(info.windowRect);
    
    % Set the blend funciton for the screen
    Screen('BlendFunction', info.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    %----------------------------------------------------------------------
    %                       Keyboard information
    %----------------------------------------------------------------------
    
    % Define the keyboard keys that are listened for. We will be using the left
    % and right arrow keys as response keys for the task and the escape key as
    % a exit/reset key
    info.keys.leftShift = KbName('LeftShift');
    info.keys.leftAlt = KbName('LeftAlt');
    info.keys.leftCtl = KbName('LeftControl');
    info.keys.yes = KbName('y');
    info.keys.no = KbName('n');
    
    info.keys.voiceKey = KbName('0');
    if fingerCB==1
        info.keys.oneKey = [KbName('1') KbName('a') 49];
        info.keys.twoKey = [KbName('4') KbName('f') 52];
    else
        info.keys.oneKey = [KbName('4') KbName('f') 52];
        info.keys.twoKey = [KbName('1') KbName('a') 49];
    end
    
    status='';
    for iBlock=startBlock:numBlocks
        if any(strcmp(status,{'abort','error'}))
            break
        end
        taskCode=latinSquare{orderCB}(iBlock);
        if double(taskCode) > 96
            taskNum=double(taskCode)-96;
        else
            taskNum=double(taskCode)-64;
        end
        
        ch='';
        if strcmp(taskNames{taskNum},'Naming') && ~info.CMU
            DrawFormattedText(info.window, 'Experimenters need to set up mike','center', 'center', info.black);
            Screen('Flip',info.window);
            while ~strcmp(ch,'C')
                [ch, when]=GetChar;
            end
            Screen('Flip',info.window);
            info.pstHandle=CMUBox('Open','pst','COM5','ftdi');
            info.CMU=1;
        end
        
        if ~strcmp(taskNames{taskNum},'Naming') && info.CMU
            CMUBox('Close',info.pstHandle);
            info.CMU=0;
            DrawFormattedText(info.window, 'Experimenters need to put away mike','center', 'center', info.black);
            Screen('Flip',info.window);
            while ~strcmp(ch,'C')
                [ch, when]=GetChar;
            end
            Screen('Flip',info.window);
        end
        
        if double(taskCode) > 96
            taskNum=double(taskCode)-96;
            pracDone=0;
            while ~pracDone
                pracLists{taskNum}=Shuffle(pracLists{taskNum},2);
                eval(['[status, info]=' taskNames{taskNum} 'Task(pracLists{taskNum},pracHeaders{taskNum},info,''practice'',[]);']);
                if any(strcmp(status,{'abort','error'}))
                    break
                elseif strcmp(status,'repeat')
                    pracDone=0;
                else
                    pracDone=1;
                end
            end
        else
            numTaskTrials=size(trialSequences{taskNum},1);
            numReps=length(find(findstr(taskCode,latinSquare{orderCB})));
            numTrials=numTaskTrials/numReps;
            blockNum=(trialCounter(taskNum)+numTrials)/numTrials;
            specsFile=[expName '-' sprintf('%03d',sessionInfo.subject) '-' taskNames{taskNum}  'Task' num2str(blockNum)];
            
            eval(['status=' taskNames{taskNum} 'Task(trialSequences{taskNum}(trialCounter(taskNum)+1:trialCounter(taskNum)+numTrials,:),trialHeaders{taskNum},info,specsFile,sessionInfo);']);
            if ~isempty(status)
                break
            end
            trialCounter(taskNum)=trialCounter(taskNum)+numTrials;
        end
    end
    
    if isempty(status)
        % end of experiment screen. We clear the screen once they have made their
        % response
        DrawFormattedText(info.window, 'Experiment Finished \n\n Press Any Key To Exit',...
            'center', 'center', info.black);
        Screen('Flip', info.window);
        KbStrokeWait;
        sca;
    end
    
    if exist('info','var') && isfield(info,'CMU') && (info.CMU==1)
        CMUBox('Close',info.pstHandle);
        info.CMU=0;
    end
    
catch ME
    Screen('CloseAll');
    ShowCursor;
    if exist('info','var') && isfield(info,'CMU') && (info.CMU==1)
        CMUBox('Close',info.pstHandle);
        info.CMU=0;
    end
    rethrow(ME)
    return;
end


%----------------------------------------------------------------------------------------------
function [status, info]=NamingTask(trialSequences,trialHeaders,info,specsFile,sessionInfo)
% Task in which participant names a word stimulus and the onset time is measured using a voice trigger.

status='';
stimCol=find(strcmp('stim',trialHeaders));
condCol=find(strcmp('cond',trialHeaders));
corrCol=find(strcmp('correct',trialHeaders));
% trialSpecNames={'stim';'cond';'RT';'ACC';'wt';'rt';'ft'};
trialSpecNames={'stim';'cond';'RT';'ACC'};

numTrials=size(trialSequences,1);
respMat = cell(numTrials,length(trialSpecNames));
eventList=[];
eventList(1)=500; %word for 500 ms
eventList(2)=1000; %response fixation for 1000 ms

if ~info.CMU
    DrawFormattedText(info.window, 'Experimenters need to set up mike','center', 'center', info.black);
    Screen('Flip',info.window);
    ch='';
    while ~strcmp(ch,'C')
        [ch, when]=GetChar;
    end
    Screen('Flip',info.window);
    info.pstHandle=CMUBox('Open','pst','COM5','ftdi');
    info.CMU=1;
end

%directions

if strcmp(specsFile,'practice')
    instructionsText='This section is practice to learn the task.  You will repeat it until you reach 90% accuracy.\n\n';
else
    instructionsText='';
end

instructionsText=[instructionsText 'You will be presented with both words and non-words.\n\n'...
    'Speak the words out loud as quickly as you can while still being accurate.\n\n'...
    'Do not speak the non-words.\n\n'...
    'Wait until the plus sign turns into a question mark before speaking.\n\n'...
    'At the end of every trial, the plus sign will turn green for correct and red for incorrect.\n\n'...
    'Correct is speaking to a word and not speaking to a non-word.\n\n'...
    'The computer will not judge the correctness of the pronunciation but please try to be accurate.\n\n'...
    'You may need to speak loudly but please try to minimize head movement.\n\n'...
    'Keep your eyes on the center of the screen.\n\n'];

if strcmp(specsFile,'practice')
    instructionsText=[instructionsText 'Press a Button To Begin\n'];
else
    instructionsText=[instructionsText 'Wait for Experimenter To Begin\n'];
end

DrawFormattedText(info.window, instructionsText,'center', 'center', info.black);
Screen('Flip',info.window);
ch='';
FlushEvents
if strcmp(specsFile,'practice')
    KbStrokeWait;
else
    while ~strcmp(ch,'C')
        [ch, when]=GetChar;
    end
end

if ~strcmp(specsFile,'practice')
    DrawFormattedText(info.window, 'Get Ready',...
        'center', 'center', info.black);
    Screen('Flip',info.window);
    writeHeader(info,sessionInfo,trialSpecNames)
end

DrawFormattedText(info.window, '+', 'center', 'center', info.black);
vbl = Screen('Flip', info.window);
% Draw fixation for 1000 ms
theDuration=1;
DrawFormattedText(info.window, '+', 'center', 'center', info.black);
vbl = Screen('Flip', info.window, vbl + (round(theDuration / info.ifi) - 0.5 -1) * info.ifi);

for iTrial = 1:numTrials
    theStim=trialSequences{iTrial,stimCol};
    theCond=trialSequences{iTrial,condCol};
    theCorr=str2double(trialSequences{iTrial,corrCol});
    ACC=2;
    RT=0;
    eventsDone=0;
    eventCounter=1;
    
    %     % Draw the word
    %     DrawFormattedText(info.window, theStim, 'center', 'center', info.black);
    %     tarTime = Screen('Flip', info.window);
    %     outp(info.address,1);
    %     pause(.1);outp(info.address,0);
    
    evt=CMUBox('GetEvent',info.pstHandle);
    if ~isempty(evt)
        while ~isempty(evt) %flush the queue
            evt=CMUBox('GetEvent',info.pstHandle);
        end
    end
    
    % Draw the word
    DrawFormattedText(info.window, theStim, 'center', 'center', info.black);
    tarTime = Screen('Flip', info.window);
    outp(info.address,1);
    pause(.1);outp(info.address,0);
    
    % Draw fixation after 500 ms
    theDuration=.5;
    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
    waitTime = Screen('Flip', info.window, tarTime + (round(theDuration / info.ifi) -1) * info.ifi);
    
    % Draw ? after 1000 ms
    theDuration=1;
    DrawFormattedText(info.window, '?', 'center', 'center', info.black);
    respTime = Screen('Flip', info.window, waitTime + (round(theDuration/ info.ifi) -1) * info.ifi);
    nextEvent=respTime+((round((eventList(eventCounter)/1000) / info.ifi) -3) * info.ifi);
    
    %     nextEvent=tarTime+((round((eventList(eventCounter)/1000) / info.ifi) -1.5) * info.ifi);
    while ~eventsDone
        if ACC==2
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyIsDown
                if keyCode(info.keys.leftShift) && keyCode(info.keys.leftAlt) && keyCode(info.keys.leftCtl)
                    DrawFormattedText(info.window, 'Abort experiment?', 'center', 'center', info.black);
                    vbl=Screen('Flip', info.window);
                    keyCode=zeros(1,256);
                    while (~any(keyCode(info.keys.yes)) && ~any(keyCode(info.keys.no)))
                        [secs, keyCode,deltaSecs]=KbWait;
                    end
                    if any(keyCode(info.keys.yes))
                        if ~strcmp(specsFile,'practice')
                            err=saveTrialSpecs(trialSpecNames,respMat,specsFile);
                        end
                        ShowCursor;
                        sca;
                        status='abort';
                        return
                    end
                end
            end
            evt=CMUBox('GetEvent',info.pstHandle);
            if ~isempty(evt) && bitand(evt.state,32)
                %         elseif any(keyCode(info.keys.voiceKey))
                if theCorr==1
                    ACC=1;
                else
                    ACC=0;
                end
                outp(info.address,4);
                pause(.1);outp(info.address,0);
                RT = round((evt.time - tarTime)*1000);
                while ~isempty(evt) %flush the queue
                    evt=CMUBox('GetEvent',info.pstHandle);
                end
            end
        end
        
        if nextEvent < GetSecs-info.ifi
            eventCounter=eventCounter+1;
            switch eventCounter
                case 2 %response fixation for 1000 ms
                    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
                    vbl=Screen('Flip', info.window,nextEvent);
                    nextEvent=GetSecs+((round((eventList(eventCounter)/1000) / info.ifi) -1.5) * info.ifi);
                    fixTime=vbl;
                case 3
                    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
                    vbl=Screen('Flip', info.window,nextEvent);
                    eventsDone=1;
            end
        end
    end
    
    % Draw feedback
    if (theCorr==2) && (ACC==2)
        ACC=1;
    elseif (theCorr==1) && (ACC==2)
        ACC=0;
    end
    
    switch ACC
        case 0
            theColor=info.red;
        case 1
            theColor= info.green;
        case 2
            theColor= info.blue;
    end
    DrawFormattedText(info.window, '+', 'center', 'center', theColor);
    vbl = Screen('Flip', info.window);
    
    % Draw fixation after 250 ms ISI
    theDuration=.25;
    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
    vbl = Screen('Flip', info.window, vbl + (round(theDuration/ info.ifi) - 0.5 -1) * info.ifi);
    
    % Record the trial data into out data matrix
    respMat{iTrial,1} = theStim;
    respMat{iTrial,2} = theCond;
    respMat{iTrial,3} = RT;
    respMat{iTrial,4} = ACC;
    %     respMat{iTrial,5} = round((waitTime-tarTime)*1000);
    %     respMat{iTrial,6} = round((respTime-tarTime)*1000);
    %     respMat{iTrial,7} = round((fixTime-tarTime)*1000);
    trialSpecs=respMat(iTrial,:);
    if length(trialSpecs) > length(trialSpecNames)
        status='error';
        return
    end
    
    % ~1000 ms ISI
    theDuration=1;
    if strcmp(specsFile,'practice')
        pause(theDuration)
    else
        writeTrialSpecs(trialSpecs,info,theDuration);
    end
    fprintf('Trial #%d, Accuracy %03d\n',iTrial,round(100*length(find(cell2mat(respMat(1:iTrial,4))==1))/iTrial));
end

if strcmp(specsFile,'practice')
    pracAcc=length(find(cell2mat(respMat(:,4))==1))/numTrials;
    
    DrawFormattedText(info.window, sprintf('Accuracy: %d',pracAcc*100),'center', 'center', info.black);
    vbl = Screen('Flip', info.window);
    % leave message up for one second
    theDuration=1;
    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
    vbl = Screen('Flip', info.window, vbl + (round(theDuration/ info.ifi) - 0.5 -1) * info.ifi);
    
    if pracAcc < .9
        status='repeat';
    end
else
    DrawFormattedText(info.window, 'Rest Break','center', 'center', info.black);
    vbl = Screen('Flip', info.window);
    
    % leave message up for one second
    theDuration=1;
    DrawFormattedText(info.window, 'Rest Break','center', 'center', info.black);
    vbl = Screen('Flip', info.window, vbl + (round(theDuration/ info.ifi) - 0.5 -1) * info.ifi);
    
    err=saveTrialSpecs(trialSpecNames,respMat,specsFile);
end

end

%----------------------------------------------------------------------------------------------
function [status, info]=MaskedTask(trialSequences,trialHeaders,info,specsFile,sessionInfo)
% Task in which participant performs lexical decision on stimulus degraded by supraliminal masks and visual noise.
status='';

stimCol=find(strcmp('stim',trialHeaders));
condCol=find(strcmp('cond',trialHeaders));
corrCol=find(strcmp('correct',trialHeaders));
mask1Col=find(strcmp('mask1',trialHeaders));
mask2Col=find(strcmp('mask2',trialHeaders));

%trialSpecNames={'stim';'cond';'RT';'ACC';'resp';'mask1';'mask2';'tt';'m2t';'ft'};
trialSpecNames={'stim';'cond';'RT';'ACC';'resp'};

if strcmp(specsFile,'practice')
    stimFolder='pracs';
else
    stimFolder='stims';
end


numTrials=size(trialSequences,1);
respMat = cell(numTrials,length(trialSpecNames));

eventList=[];
eventList(1)=250; %target for 250 ms
eventList(2)=250; %backward mask for 250 ms
eventList(3)=500; %response fixation for 500 ms
% eventSchedule=zeros(round(sum(eventList/1000)/info.ifi)-1,1); %already flipped for one frame
% eventCounter=1;
% for iFrame=1:length(eventSchedule)
%     if iFrame>round(sum(eventList(1:eventCounter)/1000)/info.ifi)-1 %already flipped for one frame
%         eventCounter=min(eventCounter+1,length(eventList));
%     end
%     eventSchedule(iFrame)=eventCounter;
% end

%directions
if info.fingerCB==1
    wordButton='left';
    nonwordButton='right';
else
    wordButton='right';
    nonwordButton='left';
end

if strcmp(specsFile,'practice')
    instructionsText='This section is practice to learn the task.  You will repeat it until you reach 90% accuracy.\n\n';
else
    instructionsText='';
end

instructionsText=[instructionsText 'You will be presented with both words and non-words.\n\n'...
    'Press the ' wordButton ' button if it is a word and the ' nonwordButton ' if it is not.\n\n'...
    'Please respond as quickly as you can while still being accurate.\n\n'...
    'At the end of every trial, the plus sign will turn green for correct and red for incorrect.\n\n'...
    'If you do not respond in time, it will turn blue.\n\n'...
    'Please try to minimize head movement.\n\n'...
    'Keep your eyes on the center of the screen.\n\n'];

if strcmp(specsFile,'practice')
    instructionsText=[instructionsText 'Press a Button To Begin\n'];
else
    instructionsText=[instructionsText 'Wait for Experimenter To Begin\n'];
end

DrawFormattedText(info.window, instructionsText,'center', 'center', info.black);
Screen('Flip', info.window);
ch='';
FlushEvents
if strcmp(specsFile,'practice')
    KbStrokeWait;
else
    while ~strcmp(ch,'C')
        [ch, when]=GetChar;
    end
end

if ~strcmp(specsFile,'practice')
    DrawFormattedText(info.window, 'Get Ready',...
        'center', 'center', info.black);
    Screen('Flip',info.window);
    writeHeader(info,sessionInfo,trialSpecNames)
end

[screenXpixels, screenYpixels] = Screen('WindowSize', info.window);

DrawFormattedText(info.window, '+', 'center', 'center', info.black);
vbl = Screen('Flip', info.window);
% Draw fixation for 1000 ms
theDuration=1;
DrawFormattedText(info.window, '+', 'center', 'center', info.black);
vbl = Screen('Flip', info.window, vbl + (round(theDuration / info.ifi) - 0.5 -1) * info.ifi);

for iTrial = 1:numTrials
    theStim=trialSequences{iTrial,stimCol};
    theCond=trialSequences{iTrial,condCol};
    theCorr=str2double(trialSequences{iTrial,corrCol});
    theMask1=trialSequences{iTrial,mask1Col};
    theMask2=trialSequences{iTrial,mask2Col};
    ACC=2;
    theResp=0;
    RT=0;
    eventsDone=0;
    eventCounter=1;
    
    theImageLocation = [pwd filesep stimFolder filesep theMask1];
    mask1Image = imread(theImageLocation);
    [s1, s2, s3] = size(mask1Image);
    if s1 > screenYpixels || s2 > screenYpixels
        disp('ERROR! Image is too big to fit on the screen');
        sca;
        status='error';
        return;
    end
    mask1Texture = Screen('MakeTexture', info.window, mask1Image);
    
    theImageLocation = [pwd filesep stimFolder filesep theStim];
    stimImage = imread(theImageLocation);
    [s1, s2, s3] = size(stimImage);
    if s1 > screenYpixels || s2 > screenYpixels
        disp('ERROR! Image is too big to fit on the screen');
        sca;
        status='error';
        return;
    end
    stimTexture = Screen('MakeTexture', info.window, stimImage);
    
    theImageLocation = [pwd filesep stimFolder filesep theMask2];
    mask2Image = imread(theImageLocation);
    [s1, s2, s3] = size(mask2Image);
    if s1 > screenYpixels || s2 > screenYpixels
        disp('ERROR! Image is too big to fit on the screen');
        sca;
        status='error';
        return;
    end
    mask2Texture = Screen('MakeTexture', info.window, mask2Image);
    
    % Draw the premask
    Screen('DrawTexture', info.window, mask1Texture, [], [], 0);
    vbl = Screen('Flip', info.window);
    m1Time=vbl;
    m2Time=0;
    fixTime=0;
    
    % Draw stimulus after 250 ms
    theDuration=.25;
    Screen('DrawTexture', info.window, stimTexture, [], [], 0);
    tarTime = Screen('Flip', info.window, vbl + (round(theDuration / info.ifi) -1) * info.ifi);
    outp(info.address,1);
    pause(.1);outp(info.address,0);
    nextEvent=tarTime+((round((eventList(eventCounter)/1000) / info.ifi) -1.5) * info.ifi);
    while ~eventsDone
        if ACC==2
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyIsDown
                if keyCode(info.keys.leftShift) && keyCode(info.keys.leftAlt) && keyCode(info.keys.leftCtl)
                    DrawFormattedText(info.window, 'Abort experiment?', 'center', 'center', info.black);
                    vbl=Screen('Flip', info.window);
                    keyCode=zeros(1,256);
                    while (~any(keyCode(info.keys.yes)) && ~any(keyCode(info.keys.no)))
                        [secs, keyCode,deltaSecs]=KbWait;
                    end
                    if any(keyCode(info.keys.yes))
                        if ~strcmp(specsFile,'practice')
                            err=saveTrialSpecs(trialSpecNames,respMat,specsFile);
                        end
                        ShowCursor;
                        sca;
                        status='abort';
                        return
                    end
                elseif any(keyCode(info.keys.oneKey))
                    outp(info.address,4);
                    pause(.1);outp(info.address,0);
                    if theCorr==1
                        ACC=1;
                    else
                        ACC=0;
                    end
                    RT = round((secs - tarTime)*1000);
                    theResp=1;
                elseif any(keyCode(info.keys.twoKey))
                    outp(info.address,4);
                    pause(.1);outp(info.address,0);
                    if theCorr==2
                        ACC=1;
                    else
                        ACC=0;
                    end
                    RT = round((secs - tarTime)*1000);
                    theResp=2;
                end
            end
        end
        if nextEvent < GetSecs-info.ifi
            eventCounter=eventCounter+1;
            switch eventCounter
                %             case 1 %target for 250 ms
                %                 Screen('DrawTexture', info.window, stimTexture, [], [], 0);
                case 2 %backward mask for 250 ms
                    Screen('DrawTexture', info.window, mask2Texture, [], [], 0);
                    vbl=Screen('Flip', info.window,nextEvent);
                    nextEvent=GetSecs+((round((eventList(eventCounter)/1000) / info.ifi) -1.5) * info.ifi);
                    m2Time=vbl;
                case 3 %response fixation for 500 ms
                    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
                    vbl=Screen('Flip', info.window,nextEvent);
                    nextEvent=GetSecs+((round((eventList(eventCounter)/1000) / info.ifi) -1.5) * info.ifi);
                    fixTime=vbl;
                case 4
                    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
                    vbl=Screen('Flip', info.window,nextEvent);
                    eventsDone=1;
            end
        end
    end
    
    % Record the trial data into out data matrix
    respMat{iTrial,1} = theStim;
    respMat{iTrial,2} = theCond;
    respMat{iTrial,3} = RT;
    respMat{iTrial,4} = ACC;
    respMat{iTrial,5} = theResp;
    %     respMat{iTrial,6} = theMask1;
    %     respMat{iTrial,7} = theMask2;
    %     respMat{iTrial,8} = round((tarTime-m1Time)*1000);
    %     respMat{iTrial,9} = round((m2Time-m1Time)*1000);
    %     respMat{iTrial,10} = round((fixTime-m1Time)*1000);
    trialSpecs=respMat(iTrial,:);
    if length(trialSpecs) > length(trialSpecNames)
        status='error';
        return
    end
    
    % Draw feedback
    switch ACC
        case 0
            theColor=info.red;
        case 1
            theColor= info.green;
        case 2
            theColor= info.blue;
    end
    DrawFormattedText(info.window, '+', 'center', 'center', theColor);
    vbl = Screen('Flip', info.window);
    
    % Draw fixation after 250 ms feedback
    theDuration=.25;
    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
    vbl = Screen('Flip', info.window, vbl + (round(theDuration/ info.ifi) - 0.5 -1) * info.ifi);
    
    % Draw fixation after 1000 ms ISI
    theDuration=1.6;
    if strcmp(specsFile,'practice')
        pause(theDuration)
    else
        writeTrialSpecs(trialSpecs,info,theDuration);
    end
    fprintf('Trial #%d, Accuracy %03d\n',iTrial,round(100*length(find(cell2mat(respMat(1:iTrial,4))==1))/iTrial));
    Screen('Close') %clear out the textures so they don't clog up the memory
end

if strcmp(specsFile,'practice')
    pracAcc=length(find(cell2mat(respMat(:,4))==1))/numTrials;
    
    DrawFormattedText(info.window, sprintf('Accuracy: %d',pracAcc*100),'center', 'center', info.black);
    vbl = Screen('Flip', info.window);
    % leave message up for one second
    theDuration=1;
    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
    vbl = Screen('Flip', info.window, vbl + (round(theDuration/ info.ifi) - 0.5 -1) * info.ifi);
    
    if pracAcc < .9
        status='repeat';
    end
else
    DrawFormattedText(info.window, 'Rest Break','center', 'center', info.black);
    vbl = Screen('Flip', info.window);
    
    % leave message up for one second
    theDuration=1;
    DrawFormattedText(info.window, 'Rest Break','center', 'center', info.black);
    vbl = Screen('Flip', info.window, vbl + (round(theDuration/ info.ifi) - 0.5 -1) * info.ifi);
    
    err=saveTrialSpecs(trialSpecNames,respMat,specsFile);
end

end

%----------------------------------------------------------------------------------------------
function [status, info]=RhymingTask(trialSequences,trialHeaders,info,specsFile,sessionInfo)
% Task in which participant performs rhyming decision.
status='';

targetCol=find(strcmp('target',trialHeaders));
condCol=find(strcmp('cond',trialHeaders));
corrCol=find(strcmp('correct',trialHeaders));
primeCol=find(strcmp('prime',trialHeaders));

%trialSpecNames={'target';'cond';'RT';'ACC';'resp';'prime';'tt';'ft'};
trialSpecNames={'target';'cond';'RT';'ACC';'resp';'prime'};

numTrials=size(trialSequences,1);
respMat = cell(numTrials,length(trialSpecNames));

eventList=[];
eventList(1)=200; %target for 200 ms
eventList(2)=800; %response fixation for 800 ms
% eventSchedule=zeros(round(sum(eventList/1000)/info.ifi)-1,1); %already flipped for one frame
% eventCounter=1;
% for iFrame=1:length(eventSchedule)
%     if iFrame>round(sum(eventList(1:eventCounter)/1000)/info.ifi)-1  %already flipped for one frame
%         eventCounter=min(eventCounter+1,length(eventList));
%     end
%     eventSchedule(iFrame)=eventCounter;
% end

%directions
if info.fingerCB==1
    rhymeButton='left';
    nonrhymeButton='right';
else
    rhymeButton='right';
    nonrhymeButton='left';
end
if strcmp(specsFile,'practice')
    instructionsText='This section is practice to learn the task.  You will repeat it until you reach 90% accuracy.\n\n';
else
    instructionsText='';
end

instructionsText=[instructionsText 'You will be presented with two words in a row.\n\n'...
    'Press the ' rhymeButton ' button if they rhyme and the ' nonrhymeButton ' if they do not.\n\n'...
    'Please respond as quickly as you can while still being accurate.\n\n'...
    'At the end of every trial, the plus sign will turn green for correct and red for incorrect.\n\n'...
    'If you do not respond in time, it will turn blue.\n\n'...
    'Please try to minimize head movement.\n\n'...
    'Keep your eyes on the center of the screen.\n\n'];

if strcmp(specsFile,'practice')
    instructionsText=[instructionsText 'Press a Button To Begin\n'];
else
    instructionsText=[instructionsText 'Wait for Experimenter To Begin\n'];
end

DrawFormattedText(info.window,instructionsText,'center', 'center', info.black);
Screen('Flip', info.window);
ch='';
FlushEvents
if strcmp(specsFile,'practice')
    KbStrokeWait;
else
    while ~strcmp(ch,'C')
        [ch, when]=GetChar;
    end
end

if ~strcmp(specsFile,'practice')
    DrawFormattedText(info.window, 'Get Ready',...
        'center', 'center', info.black);
    Screen('Flip',info.window);
    writeHeader(info,sessionInfo,trialSpecNames)
end

[screenXpixels, screenYpixels] = Screen('WindowSize', info.window);

DrawFormattedText(info.window, '+', 'center', 'center', info.black);
vbl = Screen('Flip', info.window);
% Draw fixation for 1000 ms
theDuration=1;
DrawFormattedText(info.window, '+', 'center', 'center', info.black);
vbl = Screen('Flip', info.window, vbl + (round(theDuration / info.ifi) - 0.5 -1) * info.ifi);

for iTrial = 1:numTrials
    theTarget=trialSequences{iTrial,targetCol};
    theCond=trialSequences{iTrial,condCol};
    theCorr=str2double(trialSequences{iTrial,corrCol});
    thePrime=trialSequences{iTrial,primeCol};
    ACC=2;
    theResp=0;
    RT=0;
    eventsDone=0;
    eventCounter=1;
    
    % Draw the prime
    DrawFormattedText(info.window, thePrime, 'center', 'center', info.black);
    vbl = Screen('Flip', info.window);
    primeTime=vbl;
    fixTime=0;
    
    % Draw the fixation after 200 ms
    theDuration=.2;
    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
    vbl = Screen('Flip', info.window, vbl + (round(theDuration / info.ifi) -1) * info.ifi);
    
    % Draw the target after 800 ms
    theDuration=.8;
    DrawFormattedText(info.window, theTarget, 'center', 'center', info.black);
    tarTime = Screen('Flip', info.window, vbl + (round(theDuration / info.ifi) -1) * info.ifi);
    outp(info.address,1);
    pause(.1);outp(info.address,0);
    nextEvent=tarTime+((round((eventList(eventCounter)/1000) / info.ifi) -1.5) * info.ifi);
    while ~eventsDone
        if ACC==2
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyIsDown
                if keyCode(info.keys.leftShift) && keyCode(info.keys.leftAlt) && keyCode(info.keys.leftCtl)
                    DrawFormattedText(info.window, 'Abort experiment?', 'center', 'center', info.black);
                    vbl=Screen('Flip', info.window);
                    keyCode=zeros(1,256);
                    while (~any(keyCode(info.keys.yes)) && ~any(keyCode(info.keys.no)))
                        [secs, keyCode,deltaSecs]=KbWait;
                    end
                    if any(keyCode(info.keys.yes))
                        if ~strcmp(specsFile,'practice')
                            err=saveTrialSpecs(trialSpecNames,respMat,specsFile);
                        end
                        ShowCursor;
                        sca;
                        status='abort';
                        return
                    end
                elseif any(keyCode(info.keys.oneKey))
                    outp(info.address,4);
                    pause(.1);outp(info.address,0);
                    if theCorr==1
                        ACC=1;
                    else
                        ACC=0;
                    end
                    RT = round((secs - tarTime)*1000);
                    theResp=1;
                elseif any(keyCode(info.keys.twoKey))
                    outp(info.address,4);
                    pause(.1);outp(info.address,0);
                    if theCorr==2
                        ACC=1;
                    else
                        ACC=0;
                    end
                    RT = round((secs - tarTime)*1000);
                    theResp=2;
                end
            end
        end
        
        if nextEvent < GetSecs-info.ifi
            eventCounter=eventCounter+1;
            switch eventCounter
                %             case 1 %target for 200 ms
                %                 DrawFormattedText(info.window, theTarget, 'center', 'center', info.black);
                case 2 %response fixation for 800 ms
                    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
                    vbl=Screen('Flip', info.window,nextEvent);
                    nextEvent=GetSecs+((round((eventList(eventCounter)/1000) / info.ifi) -1.5) * info.ifi);
                    fixTime=vbl;
                case 3
                    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
                    vbl=Screen('Flip', info.window,nextEvent);
                    eventsDone=1;
            end
        end
    end
    
    % Record the trial data into out data matrix
    respMat{iTrial,1} = theTarget;
    respMat{iTrial,2} = theCond;
    respMat{iTrial,3} = RT;
    respMat{iTrial,4} = ACC;
    respMat{iTrial,5} = theResp;
    respMat{iTrial,6} = thePrime;
    %     respMat{iTrial,7} = round((tarTime-primeTime)*1000);
    %     respMat{iTrial,8} = round((fixTime-primeTime)*1000);
    trialSpecs=respMat(iTrial,:);
    if length(trialSpecs) > length(trialSpecNames)
        status='error';
        return
    end
    
    % Draw feedback
    switch ACC
        case 0
            theColor=info.red;
        case 1
            theColor= info.green;
        case 2
            theColor= info.blue;
    end
    DrawFormattedText(info.window, '+', 'center', 'center', theColor);
    vbl = Screen('Flip', info.window);
    
    % Draw fixation after 250 ms ISI
    theDuration=.25;
    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
    vbl = Screen('Flip', info.window, vbl + (round(theDuration/ info.ifi) - 0.5 -1) * info.ifi);
    
    % Draw fixation after 1300 ms ISI
    theDuration=1.3;
    if strcmp(specsFile,'practice')
        pause(theDuration)
    else
        writeTrialSpecs(trialSpecs,info,theDuration);
    end
    fprintf('Trial #%d, Accuracy %03d\n',iTrial,round(100*length(find(cell2mat(respMat(1:iTrial,4))==1))/iTrial));
end

if strcmp(specsFile,'practice')
    pracAcc=length(find(cell2mat(respMat(:,4))==1))/numTrials;
    
    DrawFormattedText(info.window, sprintf('Accuracy: %d',pracAcc*100),'center', 'center', info.black);
    vbl = Screen('Flip', info.window);
    % leave message up for one second
    theDuration=1;
    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
    vbl = Screen('Flip', info.window, vbl + (round(theDuration/ info.ifi) - 0.5 -1) * info.ifi);
    
    if pracAcc < .9
        status='repeat';
    end
else
    DrawFormattedText(info.window, 'Rest Break','center', 'center', info.black);
    vbl = Screen('Flip', info.window);
    
    % leave message up for one second
    theDuration=1;
    DrawFormattedText(info.window, 'Rest Break','center', 'center', info.black);
    vbl = Screen('Flip', info.window, vbl + (round(theDuration/ info.ifi) - 0.5 -1) * info.ifi);
    
    err=saveTrialSpecs(trialSpecNames,respMat,specsFile);
end

end

%----------------------------------------------------------------------------------------------
function [status, info]=PrimingTask(trialSequences,trialHeaders,info,specsFile,sessionInfo)
% Task in which participant performs associative primed lexical decision.
status='';

targetCol=find(strcmp('target',trialHeaders));
condCol=find(strcmp('cond',trialHeaders));
corrCol=find(strcmp('correct',trialHeaders));
primeCol=find(strcmp('prime',trialHeaders));

%trialSpecNames={'target';'cond';'RT';'ACC';'resp';'prime';'tt';'ft'};
trialSpecNames={'target';'cond';'RT';'ACC';'resp';'prime'};

numTrials=size(trialSequences,1);
respMat = cell(numTrials,length(trialSpecNames));

eventList=[];
eventList(1)=150; %target for 150 ms
eventList(2)=850; %response fixation for 850 ms
% eventSchedule=zeros(round(sum(eventList/1000)/info.ifi)-1,1); %already flipped for one frame
% eventCounter=1;
% for iFrame=1:length(eventSchedule)
%     if iFrame>round(sum(eventList(1:eventCounter)/1000)/info.ifi)-1 %already flipped for one frame
%         eventCounter=min(eventCounter+1,length(eventList));
%     end
%     eventSchedule(iFrame)=eventCounter;
% end

%directions
if info.fingerCB==1
    wordButton='left';
    nonwordButton='right';
else
    wordButton='right';
    nonwordButton='left';
end
if strcmp(specsFile,'practice')
    instructionsText='This section is practice to learn the task.  You will repeat it until you reach 90% accuracy.\n\n';
else
    instructionsText='';
end

instructionsText=[instructionsText 'You will be presented with two stimuli in a row.\n\n'...
    'Ignore the first stimulus.\n\n'...
    'The second stimulus will be words half the time and non-words half the time.\n\n'...
    'Press the ' wordButton ' button if it is a word and the ' nonwordButton ' if it is not.\n\n'...
    'Please respond as quickly as you can while still being accurate.\n\n'...
    'At the end of every trial, the plus sign will turn green for correct and red for incorrect.\n\n'...
    'If you do not respond in time, it will turn blue.\n\n'...
    'Please try to minimize head movement.\n\n'...
    'Keep your eyes on the center of the screen.\n\n'];

if strcmp(specsFile,'practice')
    instructionsText=[instructionsText 'Press a Button To Begin\n'];
else
    instructionsText=[instructionsText 'Wait for Experimenter To Begin\n'];
end

DrawFormattedText(info.window,instructionsText,'center', 'center', info.black);
Screen('Flip', info.window);
ch='';
FlushEvents
if strcmp(specsFile,'practice')
    KbStrokeWait;
else
    while ~strcmp(ch,'C')
        [ch, when]=GetChar;
    end
end

if ~strcmp(specsFile,'practice')
    DrawFormattedText(info.window, 'Get Ready',...
        'center', 'center', info.black);
    Screen('Flip',info.window);
    writeHeader(info,sessionInfo,trialSpecNames)
end


[screenXpixels, screenYpixels] = Screen('WindowSize', info.window);

DrawFormattedText(info.window, '+', 'center', 'center', info.black);
vbl = Screen('Flip', info.window);
% Draw fixation for 1000 ms
theDuration=1;
DrawFormattedText(info.window, '+', 'center', 'center', info.black);
vbl = Screen('Flip', info.window, vbl + (round(theDuration / info.ifi) - 0.5 -1) * info.ifi);

for iTrial = 1:numTrials
    theTarget=trialSequences{iTrial,targetCol};
    theCond=trialSequences{iTrial,condCol};
    theCorr=str2double(trialSequences{iTrial,corrCol});
    thePrime=trialSequences{iTrial,primeCol};
    ACC=2;
    theResp=0;
    RT=0;
    eventsDone=0;
    eventCounter=1;
    
    % Draw the prime
    DrawFormattedText(info.window, thePrime, 'center', 'center', info.black);
    vbl = Screen('Flip', info.window);
    primeTime=vbl;
    fixTime=0;
    
    % Draw the target after 150 ms
    theDuration=.15;
    DrawFormattedText(info.window, theTarget, 'center', 'center', info.black);
    tarTime = Screen('Flip', info.window, vbl + (round(theDuration / info.ifi) -1) * info.ifi);
    outp(info.address,1);
    pause(.1);outp(info.address,0);
    nextEvent=tarTime+((round((eventList(eventCounter)/1000) / info.ifi) -1.5) * info.ifi);
    while ~eventsDone
        if ACC==2
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyIsDown
                if keyCode(info.keys.leftShift) && keyCode(info.keys.leftAlt) && keyCode(info.keys.leftCtl)
                    DrawFormattedText(info.window, 'Abort experiment?', 'center', 'center', info.black);
                    vbl=Screen('Flip', info.window);
                    keyCode=zeros(1,256);
                    while (~any(keyCode(info.keys.yes)) && ~any(keyCode(info.keys.no)))
                        [secs, keyCode,deltaSecs]=KbWait;
                    end
                    if any(keyCode(info.keys.yes))
                        if ~strcmp(specsFile,'practice')
                            err=saveTrialSpecs(trialSpecNames,respMat,specsFile);
                        end
                        ShowCursor;
                        sca;
                        status='abort';
                        return
                    end
                elseif any(keyCode(info.keys.oneKey))
                    outp(info.address,4);
                    pause(.1);outp(info.address,0);
                    if theCorr==1
                        ACC=1;
                    else
                        ACC=0;
                    end
                    RT = round((secs - tarTime)*1000);
                    theResp=1;
                elseif any(keyCode(info.keys.twoKey))
                    outp(info.address,4);
                    pause(.1);outp(info.address,0);
                    if theCorr==2
                        ACC=1;
                    else
                        ACC=0;
                    end
                    RT = round((secs - tarTime)*1000);
                    theResp=2;
                end
            end
            %             evt=CMUBox('GetEvent',info.pstHandle);
        end
        
        %         currentFrame=round((GetSecs-tarTime)/info.ifi);
        %         if currentFrame>length(eventSchedule)
        %             break
        %         else
        if nextEvent < GetSecs-info.ifi
            eventCounter=eventCounter+1;
            switch eventCounter
                %                 case 1 %target for 150 ms
                %                     DrawFormattedText(info.window, theTarget, 'center', 'center', info.black);
                case 2 %response fixation for 1350 ms
                    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
                    vbl=Screen('Flip', info.window,nextEvent);
                    nextEvent=GetSecs+((round((eventList(eventCounter)/1000) / info.ifi) -1.5) * info.ifi);
                    fixTime=vbl;
                case 3
                    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
                    vbl=Screen('Flip', info.window,nextEvent);
                    eventsDone=1;
            end
        end
    end
    
    % Record the trial data into out data matrix
    respMat{iTrial,1} = theTarget;
    respMat{iTrial,2} = theCond;
    respMat{iTrial,3} = RT;
    respMat{iTrial,4} = ACC;
    respMat{iTrial,5} = theResp;
    respMat{iTrial,6} = thePrime;
    %     respMat{iTrial,7} = round((tarTime-primeTime)*1000);
    %     respMat{iTrial,8} = round((fixTime-primeTime)*1000);
    trialSpecs=respMat(iTrial,:);
    if length(trialSpecs) > length(trialSpecNames)
        status='error';
        return
    end
    
    % Draw feedback
    switch ACC
        case 0
            theColor=info.red;
        case 1
            theColor= info.green;
        case 2
            theColor= info.blue;
    end
    DrawFormattedText(info.window, '+', 'center', 'center', theColor);
    vbl = Screen('Flip', info.window);
    
    % Draw fixation after 250 ms ISI
    theDuration=.25;
    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
    vbl = Screen('Flip', info.window, vbl + (round(theDuration/ info.ifi) - 0.5 -1) * info.ifi);
    
    % Draw fixation after 1000 ms ISI
    theDuration=1.6;
    if strcmp(specsFile,'practice')
        pause(theDuration)
    else
        writeTrialSpecs(trialSpecs,info,theDuration);
    end
    fprintf('Trial #%d, Accuracy %03d\n',iTrial,round(100*length(find(cell2mat(respMat(1:iTrial,4))==1))/iTrial));
end

if strcmp(specsFile,'practice')
    pracAcc=length(find(cell2mat(respMat(:,4))==1))/numTrials;
    
    DrawFormattedText(info.window, sprintf('Accuracy: %d',pracAcc*100),'center', 'center', info.black);
    vbl = Screen('Flip', info.window);
    % leave message up for one second
    theDuration=1;
    DrawFormattedText(info.window, '+', 'center', 'center', info.black);
    vbl = Screen('Flip', info.window, vbl + (round(theDuration/ info.ifi) - 0.5 -1) * info.ifi);
    
    if pracAcc < .9
        status='repeat';
    end
else
    DrawFormattedText(info.window, 'Rest Break','center', 'center', info.black);
    vbl = Screen('Flip', info.window);
    
    % leave message up for one second
    theDuration=1;
    DrawFormattedText(info.window, 'Rest Break','center', 'center', info.black);
    vbl = Screen('Flip', info.window, vbl + (round(theDuration/ info.ifi) - 0.5 -1) * info.ifi);
    
    err=saveTrialSpecs(trialSpecNames,respMat,specsFile);
end

end

%----------------------------------------------------------------------------------------------
function err=saveTrialSpecs(trialSpecNames,respMat,specsFile)
% Saves a text file with the trial specs at the end of the block.

err=0;
sameName=1;
theNumber=0;
fileNameStem=specsFile;
while sameName
    sameName=0;
    if exist([pwd filesep specsFile '.txt'],'file')
        sameName=1;
    end;
    if sameName
        theNumber=theNumber+1;
        specsFile=[fileNameStem '-' num2str(theNumber)];
    end;
end;

outFID=fopen([specsFile '.txt'],'w');
if (outFID == -1)
    err=1;
    return
end

for iSpec = 1:length(trialSpecNames)
    fprintf(outFID,'%s\t',trialSpecNames{iSpec});
end
fprintf(outFID,'\r');
for iTrial=1:size(respMat,1)
    for iSpec = 1:size(respMat,2)
        theSpec=respMat{iTrial,iSpec};
        if isnumeric(theSpec)
            theSpec=num2str(theSpec);
        end
        fprintf(outFID,'%s\t',theSpec);
    end
    fprintf(outFID,'\r');
end
fclose(outFID);
end

%----------------------------------------------------------------------------------------------
function writePortByte(theInteger,info)
if isempty(theInteger)
    theInteger=0;
end
if theInteger > 255
    theInteger = 0;
elseif theInteger < 0
    theInteger = 0;
end
if theInteger == 0
    outp(info.address, 4)
    pause(info.sleepDur)
    outp(info.address, 0)
    pause(info.sleepDur)
else
    outp(info.address, 1)
    pause(info.sleepDur)
    outp(info.address, 0)
    pause(info.sleepDur)
    outp(info.address, theInteger)
    pause(info.sleepDur)
    outp(info.address, 0)
    pause(info.sleepDur)
end
end

%----------------------------------------------------------------------------------------------
function writePortWord(theInteger,info)
if isempty(theInteger)
    theInteger=0;
end
loByte=bitand(theInteger,255,'int16');
hiByte=bitand(bitshift(theInteger,-8,'int16'),255,'int16');
if theInteger == 0
    outp(info.address, 7)
    pause(info.sleepDur)
    outp(info.address, 0)
    pause(info.sleepDur)
elseif loByte == 0
    outp(info.address, 5)
    pause(info.sleepDur)
    outp(info.address, 0)
    pause(info.sleepDur)
    outp(info.address, hiByte)
    pause(info.sleepDur)
    outp(info.address, 0)
    pause(info.sleepDur)
elseif hiByte == 0
    outp(info.address, 6)
    pause(info.sleepDur)
    outp(info.address, 0)
    pause(info.sleepDur)
    outp(info.address, loByte)
    pause(info.sleepDur)
    outp(info.address, 0)
    pause(info.sleepDur)
else
    outp(info.address, 2)
    pause(info.sleepDur)
    outp(info.address, 0)
    pause(info.sleepDur)
    outp(info.address, loByte)
    pause(info.sleepDur)
    outp(info.address, 0)
    pause(info.sleepDur)
    outp(info.address, hiByte)
    pause(info.sleepDur)
    outp(info.address, 0)
    pause(info.sleepDur)
end
end

%----------------------------------------------------------------------------------------------
function writePortString(theString,info)
stringNum = length(theString);
if stringNum == 0
    outp(info.address, 8)
    pause(info.sleepDur)
    outp(info.address, 0)
    pause(info.sleepDur)
else
    outp(info.address, 3)
    pause(info.sleepDur)
    outp(info.address, 0)
    pause(info.sleepDur)
    if stringNum > 255
        stringNum = 255;
    end
    outp(info.address, stringNum)
    pause(info.sleepDur)
    outp(info.address, 0)
    pause(info.sleepDur)
    for iChar = 1:stringNum
        theChar = double(theString(iChar));
        writePortWord(theChar, info)
    end
end
end

%----------------------------------------------------------------------------------------------
function writeHeader(info,sessionInfo,trialSpecNames)

%Sends out header information for EP Toolkit via the trigger line
%1=unsigned integer byte (1-255), 2=signed two byte integer (-32768+32767), 3=text
%sleep after each writeport, 4=zero integer byte, 5=low two byte zero, 6=high two byte zero
%7=both bytes zero, 8=empty text

outp(info.address, double('h'))
pause(info.sleepDur)
outp(info.address, double('d'))
pause(info.sleepDur)
outp(info.address, double('r'))
pause(info.sleepDur)
outp(info.address, 10) %number of subject specs
pause(info.sleepDur)

theStringName = 'ExperimentName';
writePortString(theStringName,info) %write the field name
theString=sessionInfo.expName;
writePortString(theString,info) %write the text

theStringName = 'SessionDate';
writePortString(theStringName,info) %write the field name
theString=sessionInfo.date;
writePortString(theString,info) %write the text

theStringName = 'SessionTime';
writePortString(theStringName,info) %write the field name
theString=sessionInfo.time;
writePortString(theString,info) %write the text

theStringName = 'Subject';
writePortString(theStringName,info) %write the field name
theString=sessionInfo.subject;
writePortWord(theString,info) %write the text

theStringName = 'Session';
writePortString(theStringName,info) %write the field name
theString=sessionInfo.session;
writePortWord(theString,info) %write the text

theStringName = 'Handedness';
writePortString(theStringName,info) %write the field name
theString=sessionInfo.handedness;
writePortString(theString,info) %write the text

theStringName = 'Gender';
writePortString(theStringName,info) %write the field name
theString=sessionInfo.gender;
writePortString(theString,info) %write the text

theStringName = 'Age';
writePortString(theStringName,info) %write the field name
theString=sessionInfo.age;
writePortByte(theString,info) %write the text

theStringName = 'ResearcherID';
writePortString(theStringName,info) %write the field name
theString=sessionInfo.researcher;
writePortWord(theString,info) %write the text

theStringName = 'Group';
writePortString(theStringName,info) %write the field name
theString=sessionInfo.group;
writePortByte(theString,info) %write the text

%trial spec names
outp(info.address, length(trialSpecNames)) %number of trial specs
pause(info.sleepDur)

for iSpec=1:length(trialSpecNames)
    theStringName = trialSpecNames{iSpec};
    writePortString(theStringName,info) %write the trial spec name
end

outp(info.address, double('r'))
pause(info.sleepDur)
outp(info.address, double('d'))
pause(info.sleepDur)
outp(info.address, double('h'))
pause(info.sleepDur)
end
%----------------------------------------------------------------------------------------------
function [timeElapsed thePause]=writeTrialSpecs(trialSpecs,info,theDuration)

tic

outp(info.address, 255);
pause(info.sleepDur)

outp(info.address, length(trialSpecs)) %number of trial specs
pause(info.sleepDur)

for iSpec=1:length(trialSpecs)
    theSpec=trialSpecs{iSpec};
    if ischar(theSpec)
        writePortString(theSpec,info)
    elseif (theSpec==floor(theSpec)) && (theSpec>=0) && (theSpec<256)
        writePortByte(theSpec,info)
    elseif (theSpec~=floor(theSpec)) || (theSpec<-32767) || (theSpec>32767)
        writePortString(num2str(theSpec),info)
    else
        writePortWord(theSpec,info)
    end
end

outp(info.address, 255);
pause(info.sleepDur)

timeElapsed=toc;

thePause=theDuration-timeElapsed-.1+(rand*.2);
pause(thePause)
end


%practice
%timing errors
%countdown?

