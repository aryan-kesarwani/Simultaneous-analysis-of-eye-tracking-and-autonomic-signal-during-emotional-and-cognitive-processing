function [EPdata]=ep_alleeg2ep(ALLEEG)
% ep_alleeg2ep - [EPdata]=ep_alleeg2ep(ALLEEG) -
% Convert EEGlab's ALLEEG format to EPdata format.  Doesn't convert factor files.
%
%Inputs
%  ALLEEG         : EEGlab's ALLEG format.
%
%Outputs:
%   EPdata        : Structured array with the data and accompanying information in EP file format.  See readData.
%
%History
%  by Joseph Dien (4/30/20)
%  jdien07@mac.com
%
%  EEGlab files don't really have set conventions but to the extent that there is one,
%  the events should have both a 'value' field to denote the generic type of event,
%  as in 'trigger', and a 'type' field to denote the nature of this generic event,
%  as in the condition of the experiment.
%  Note also that this is the reverse of the FieldTrip convention.
%
% bugfix 1/24/21 JD
% Fixed crashes.  Hadn't actually been finished yet.
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

%note this is a quick and dirty conversion routine I've put together for a specific purpose so it hasn't been fully fleshed out or tested.

EPdata=ep_newFile;

subList=unique({ALLEEG.subject});
numSubs=length(subList);

sessList=unique({ALLEEG.session});
numSess=length(sessList);
EPdata.sessNames=sessList;
EPdata.sessNames=EPdata.sessNames(:);

if numSess > 0
    subList=cell(0);
    for iEEG=1:length(ALLEEG)
        subList{end+1}=[ALLEEG(iEEG).subject '-' ALLEEG(iEEG).session];
    end
    subList=unique(subList);
    numSubs=length(subList);
end

EPdata.subNames=cell(numSubs,1);

cellList=unique({ALLEEG.condition});
cellList=cellList(:);
numCells=length(cellList);

if (numSubs==1) && any([ALLEEG.trials]>1)
    EPdata.dataType='single_trial';
elseif (numSubs==1) && (numCells==1)
    EPdata.dataType='continuous';
else
    EPdata.dataType='average';    
end

if any(diff([ALLEEG.srate]))
    disp('Error: EP data format cannot handle differing sampling rates.');
    EPdata=[];
    return    
else
    EPdata.Fs=ALLEEG(1).srate;
end

timesList=[ALLEEG.times];
if any(any(timesList-timesList(:,1)))
    disp('Error: EP data format cannot handle differing time samples.');
    EPdata=[];
    return    
else
    EPdata.timeNames=ALLEEG(1).times;
    EPdata.timeNames=EPdata.timeNames(:);
end

EPdata.eloc=ALLEEG(1).chanlocs;
for iEEG=2:length(ALLEEG)
    if ~isequal(ALLEEG(1).chanlocs,ALLEEG(iEEG).chanlocs)
        disp('Error: EP data format cannot handle differing chanlocs.');
        EPdata=[];
        return
    end
end

numTotalEpochs=sum([ALLEEG.trials]);
EPdata.data=nan(ALLEEG(1).nbchan,ALLEEG(1).pnts,numTotalEpochs);
EPdata.cellNames=cell(numTotalEpochs,1);
if strcmp(EPdata.dataType,'single_trial')
    EPdata.trialNames=zeros(numTotalEpochs,1);
end
condCounter=zeros(numSubs,numCells);
for iEEG=1:length(ALLEEG)
    if numSess > 0
        theSub=find(strcmp([ALLEEG(iEEG).subject '-' ALLEEG(iEEG).session],subList));
        EPdata.sessNums=find(strcmp(ALLEEG(iEEG).session,sessList));
    else
        theSub=find(strcmp(ALLEEG(iEEG).subject,subList));
    end
    EPdata.subNames{theSub}=ALLEEG(iEEG).subject;
    if strcmp(EPdata.dataType,'average')
        EPdata.subTypes{theSub}='AVG';
    else
        EPdata.subTypes{theSub}='RAW';
    end
    theCell=find(strcmp(ALLEEG(iEEG).condition,cellList));
    numEpochs=ALLEEG(iEEG).trials;
    lastEpoch=sum(condCounter(theSub,:));
    condCounter(theSub,theCell)=condCounter(theSub,theCell)+numEpochs;
    if (condCounter(theSub,theCell) > 1) && strcmp(EPdata.dataType,'average')
        disp('Error: In EP average files, a subject cannot have more than one epoch with the same condition name.');
        EPdata=[];
        return
    end
    for iEpoch=1:numEpochs
        EPdata.data(:,:,lastEpoch+iEpoch,theSub)=ALLEEG(iEEG).data(:,:,iEpoch);
        if strcmp(EPdata.dataType,'single_trial')
            EPdata.cellNames{lastEpoch+iEpoch}=cellList{theCell};
            EPdata.trialNames(lastEpoch+iEpoch)=iEpoch;
        else
            EPdata.cellNames{lastEpoch+iEpoch}=ALLEEG(iEEG).condition;
        end
        EPdata.cellTypes{lastEpoch+iEpoch}='SGL';
        if isempty(ALLEEG(iEEG).epoch)
            EPdata.events{theSub,lastEpoch+iEpoch}=ALLEEG(iEEG).event;
        else
            EPdata.events{theSub,lastEpoch+iEpoch}=ALLEEG(iEEG).event(ALLEEG(iEEG).epoch(iEpoch).event);
        end
    end
     
    if isempty(EPdata.reference.type) %just assume these are all the same.  Not worth throwing an error if not.
        if strcmp(ALLEEG(iEEG).ref,'common')
            EPdata.reference.type='REG';
        elseif strcmp(ALLEEG(iEEG).ref,'averef')
            EPdata.reference.type='AVG';
        end
    end
end

[err]=ep_checkEPfile(EPdata);
if err
    disp('EP data format defective after ALEEG to EP conversion performed.');
    EPdata=[];
    return
end
