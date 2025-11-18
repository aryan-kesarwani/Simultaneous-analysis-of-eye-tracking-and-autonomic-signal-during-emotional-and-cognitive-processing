function [ALLEEG]=ep_ep2alleeg(EPdata)
% ep_ep2alleeg - [ALLEEG]=ep_ep2alleeg(EPdata) -
% Convert EPdata format to EEGlab's ALLEEG format.  Doesn't convert factor files.
%
%Inputs
%   EPdata        : Structured array with the data and accompanying information in EP file format.  See readData.
%
%Outputs:
%  ALLEEG         : EEGlab's ALLEG format.
%
%History
%  by Joseph Dien, with help by Grega Repovs (3/28/09)
%  jdien07@mac.com
%
%  EEGlab files don't really have set conventions but to the extent that there is one,
%  the events should have both a 'value' field to denote the generic type of event,
%  as in 'trigger', and a 'type' field to denote the nature of this generic event,
%  as in the condition of the experiment.
%  Note also that this is the reverse of the FieldTrip convention.


%
% bugfix 5/24/10 JD
% Now separates each cell into a separate dataset, not just each subject.
%
% bugfix 10/20/13 JD
% Fixed crash when saving single-trial data.
% Fixed including only first trial of each cell type.
% Fixed not including correct condition name.
% Fixed not setting up events field correctly.
% Fixed not filling out ref field.
% Fixed epoch field to correspond to epoch from entire session, not from just the one condition.
% Fixed not filling out rejected field based on bad channel and bad trial fields.
%
% bugfix 10/29/13 JD
% Fixed crash when translating single-trial EP file to EEGlab format and there are no events.
%
% modified 10/29/13 JD
% Added .keys field to events.
%
% bugfix 6/30/14 JD
% Fixed crash when file has an empty .keys field.
%
% modified 7/16/14 JD
% Simplified keys code field structure.
%
% bugfix 7/21/14 JD
% Fixed labels for REG channels being blank.
% For single_trial data, the event latency values now conform to "for epoched datasets the event latencies are also
% encoded in sample points with respect to the beginning of the data (as if the data were continuous)"
% . http://sccn.ucsd.edu/wiki/Chapter_03:_Event_Processing rather than being in terms of the beginning of the original continuous data.
% Fixed adding 'trigger' events if they are already present.
% Fixed epoch event field referring to trial number from complete dataset rather than in terms of the EEG file (the one
% condition).
% Fixed urevent event field reflecting event numbering of full dataset rather than just the one condition in the EEG
% file.
% Fixed sometimes adding too many epoch entries when exporting single_trial .set files, resulting in aborted export process.
%
% bugfix 1/13/19 JD
% Events field simply empty when there are no events rather than having empty subfields of its own, to be compatible with mffio code.
%
% modified 1/10/20 JD
% Added sessNums sessNames fields.
%
% bugfix 3/21/20 JD
% Fixed crash when events structure is empty.
% Added handles factor data by expanding it into full back-projected form.
%
% modified 5/25/20 JD
% Saves channels using labels in chanNames rather than elocs.
%
% bugfix 11/9/20 JD
% Fixed eloc labels were being saved as cells containing the labels, resulting in errors down the line.
%
% modified 10/29/21 JD
% Cells are now listed in the same order as the original EP file rather than being reordered in alphabetical order.
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

ALLEEG=[];

uniqueCells=unique(EPdata.cellNames,'stable');
numPoints=length(EPdata.timeNames);

for iEloc=1:length(EPdata.eloc)
    if isempty(EPdata.eloc(iEloc).labels)
        EPdata.eloc(iEloc).labels=EPdata.chanNames{iEloc};
    end
end

count=0;
for iSub=1:length(EPdata.subNames)
    for iCell=1:length(uniqueCells)
        count=count+1;
        
        ALLEEG(count).setname=[EPdata.ename '-' EPdata.subNames{iSub} '-' EPdata.cellNames{iCell}];
        ALLEEG(count).filename=[];
        ALLEEG(count).filepath=[];
        ALLEEG(count).subject=EPdata.subNames{iSub};
        ALLEEG(count).group='';
        ALLEEG(count).condition=uniqueCells{iCell};
        if ~isempty(EPdata.sessNames)
            ALLEEG(count).session=EPdata.sessNames{EPdata.sessNums(iSub)};
        else
            ALLEEG(count).session='';
        end
        ALLEEG(count).comments=[];
        ALLEEG(count).nbchan=length(EPdata.chanNames);
        ALLEEG(count).trials=length(strmatch(uniqueCells{iCell},EPdata.cellNames,'exact'));
        ALLEEG(count).pnts=length(EPdata.timeNames);
        ALLEEG(count).srate=EPdata.Fs;
        ALLEEG(count).xmin=min(EPdata.timeNames)/1000;
        ALLEEG(count).xmax=max(EPdata.timeNames)/1000;
        ALLEEG(count).times=EPdata.timeNames;
        ALLEEG(count).data=squeeze(ep_expandFacs(EPdata,[],[],strmatch(uniqueCells{iCell},EPdata.cellNames,'exact'),iSub,[],[],[]));
        ep_tictoc;if EPtictoc.stop;return;end
        ALLEEG(count).icaact=[];
        ALLEEG(count).icawinv=[];
        ALLEEG(count).icasphere=[];
        ALLEEG(count).icaweights=[];
        ALLEEG(count).icachansind=[];
        ALLEEG(count).chanlocs=EPdata.eloc;
        for iChan=1:length(ALLEEG(count).chanlocs)
            ALLEEG(count).chanlocs(iChan).labels=EPdata.chanNames{iChan};
        end
        ALLEEG(count).urchanlocs=[];
        ALLEEG(count).chaninfo.icachansind=[];
        if strcmp(EPdata.reference.type,'REG')
            ALLEEG(count).ref='common';
        elseif strcmp(EPdata.reference.type,'AVG')
            ALLEEG(count).ref='averef';
        else
            ALLEEG(count).ref=[];
        end
        
        trialList=strmatch(uniqueCells{iCell},EPdata.cellNames,'exact');
        clear events;
        if ~isempty(EPdata.events)
            for trial=1:length(trialList)
                theTrial=trialList(trial);
                for theEvent=1:length(EPdata.events{iSub,theTrial})
                    EPdata.events{iSub,theTrial}(theEvent).epoch=trial;
                    EPdata.events{iSub,theTrial}(theEvent).sample=(trial-1)*ALLEEG(count).pnts+EPdata.events{iSub,theTrial}(theEvent).sample;
                    if exist('events', 'var')
                        events(end+1)=EPdata.events{iSub,theTrial}(theEvent);
                    else
                        events(1)=EPdata.events{iSub,theTrial}(theEvent);
                    end
                end
            end
            if exist('events', 'var')
                tempEvents=events;
                for i=1:length(events)
                    events(i).value=tempEvents(i).type;
                    events(i).type=tempEvents(i).value;
                    events(i).latency=events(i).sample;
                    events(i).urevent=i;
                    for iKey=1:length(events(i).keys)
                        keyName=events(i).keys(iKey).code;
                        if ~isempty(keyName)
                            if ~any(strcmp(keyName,{'epoch','urevent'}))
                                if ~isempty(strfind(keyName,'#'))
                                    keyName=strrep(keyName,'#','hash_');
                                end
                                if ~isempty(strfind(keyName,'+'))
                                    keyName=strrep(keyName,'+','plus_');
                                end
                                eval(['events(i).' keyName '=events(i).keys(iKey).data;']);
                            end
                        end
                    end
                end
                events = rmfield(events,'sample');
                events = rmfield(events,'keys');
            end
        end
        
        if ~exist('events', 'var')
            %events(1)=struct('type',[],'value',[],'latency',[],'duration',[],'epoch',[]);
            events=[];
        end
        
        %add trigger events to mark each trial if single-trial data.
        if strcmp(EPdata.dataType,'single_trial')
            for trial=1:length(trialList)
                theTrial=trialList(trial);
                if isempty(events) || (length(find(strcmp('trigger',{events.value}))) ~= length(trialList))
                    if ~isempty(events) && isempty(events(end).value) && length(events)==1
                        events(end).value='trigger';
                    else
                        events(end+1).value='trigger';
                    end
                    events(end).type=EPdata.cellNames{theTrial};
                    events(end).latency=(trial-1)*ALLEEG(count).pnts+EPdata.baseline+1;
                    events(end).duration=0;
                    events(end).epoch=trial;
                end
                
                if ~isempty(EPdata.trialSpecNames)
                    events(end+1).value='TRSP';
                    events(end).type='';
                    events(end).latency=(trial-1)*ALLEEG(count).pnts+EPdata.baseline+1;
                    events(end).duration=0;
                    events(end).epoch=trial;
                    events(end).keys.names=EPdata.trialSpecNames;
                    events(end).keys.specs=EPdata.trialSpecs(theTrial,:);
                end
            end
        end
        
        ALLEEG(count).event=events;
        ALLEEG(count).urevent=events;
        
        ALLEEG(count).eventdescription=[];
        
        clear epoch;
        epoch(1)=struct('event',[],'eventtype',[],'eventvalue',[],'eventlatency',[],'eventduration',[]);
        
        for trial=1:length(trialList)
            if ~isempty(events)
                eventList=find([events.epoch]==trial);
                epoch(trial).event=eventList;
                epoch(trial).eventduration={events(eventList).duration};
                epoch(trial).eventlatency={events(eventList).latency};
                epoch(trial).eventtype={events(eventList).type};
                epoch(trial).eventvalue={events(eventList).value};
            end
            if strcmp(EPdata.dataType,'single_trial')
                ALLEEG(count).reject.manualE=(squeeze(EPdata.analysis.badChans(iSub,trialList,:))' ~=0);
                ALLEEG(count).reject.manual=EPdata.analysis.badTrials(iSub,trialList)~=0;
            else
                ALLEEG(count).reject.manualE=false(size(squeeze(EPdata.analysis.badChans(iSub,trialList,:))'));
                ALLEEG(count).reject.manual=false(EPdata.analysis.badTrials(iSub,trialList));
            end
        end
        
        ALLEEG(count).epoch=epoch;
        ALLEEG(count).epochdescription=cell(0);
        ALLEEG(count).stats=[];
        ALLEEG(count).specdata=[];
        ALLEEG(count).specicaact=[];
        ALLEEG(count).splinefile='';
        ALLEEG(count).icasplinefile='';
        ALLEEG(count).dipfit=[];
        ALLEEG(count).history=[];
        ALLEEG(count).saved='yes';
        ALLEEG(count).etc=[];
        ALLEEG(count).datfile=[];
    end
end