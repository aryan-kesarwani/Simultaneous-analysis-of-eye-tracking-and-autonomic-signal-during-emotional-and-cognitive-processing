function ep_saveEPdataset(EPdata,theDataset,fileSaved)
%  ep_saveEPdataset(EPdata,theDataset,fileSaved);
%       Adds a new dataset to the work directory.  If existing entry specified, it is deleted and replaced.
%
%Inputs:
%  EPdata         : Structured array with the data and accompanying information.  See readData.
%  theDataset    : Where to enter the dataset in the dataset list.  Defaults to the end of the list.
%  fileSaved          : Has the latest version of the file been saved?  Defaults to yes.
%
%Outputs:
%  EPdataset        : Structured array with the list of files in the work directory

%History:
%  by Joseph Dien (7/29/09)
%  jdien07@mac.com
%
%  bugfix 8/3/09 JD
%  Was dropping suffixes even though readData had already dropped the file suffix, resulting in truncation of file names
%  with multiple periods in the name.  Also, file names with multiple periods would confuse Matlab about the file type,
%  resulting in it saving the file as a ASCII file rather than as a .mat file.
%
%  modified 10/31/09 JD
%  Added more information to EPdataset to speed up main pane and preprocess pane.
%
% modified 6/15/10 JD
% Marks when a file isn't saved yet.
%
% bugfix 3/24/14 JD
% Fixed when loading in a new file that had the same name as multiple existing files, appending dashed number to prior dashed
% number instead of replacing it (e.g., name-1-2)
%
% modified 5/1/18 JD
% If trying to save a file and preferences are set to v6 or v7 and it is over 2GB, then instead of not saving it will now save in v7.3 format.
%
% modified 6/15/19 JD
% Now maintains a copy of the most recent datasets from the working set in RAM to minimize disk access time.
% Now uses EPdataset via global variable or will load it in from disk if not available.
%
% bugfix 8/19/19 JD
% Fixed no longer marking unsaved changes to data in working set with an asterisk.
%
% bugfix 2/2/20 JD
% Fixed not saving very large EPdataset caches.
%
% bugfix & modified 4/13/20 JD
% Don't change EPeeg if the dataset was already there, resulting in error that required EPeeg to be reinitialized.
% If existing entry specified, it is deleted and replaced.
%
% bugfix 4/22/20 JD
% Fixed EPeeg RAM cache being reinitialized (slowing things down) after saving changes via Edit function.
%
% bugfix 6/3/20 JD
% No longer adds -1 suffix when updating an existing entry in the working set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global EPeeg EPdataset EPmain;

if ~exist('theDataset','var')
    theDataset=length(EPdataset.dataset)+1;
end

if ~exist('fileSaved','var')
    fileSaved='yes';
end

if length(theDataset) > 1
    msg{1}='Multiple datasets were specified for deletion.  Only one at a time is allowed.';
    [msg]=ep_errorMsg(msg);
    return
end

if isempty(EPdata.dataName)
    EPdata.dataName='PCA';
end
sameName=1;
suffix=0;
baseName=EPdata.dataName;
tempName=baseName;
while sameName
    sameName=0;
    for iDataset=1:length(EPdataset.dataset)
        if strcmp(EPdataset.dataset(iDataset).dataName,tempName) && (iDataset ~= theDataset)
            sameName=1;
        end
    end
    if sameName
        suffix=suffix+1;
        tempName=[baseName '-' num2str(suffix)];
    end
end
EPdata.dataName=tempName;

newDataset=ep_addToEPworkCache(EPdata);
if isempty(newDataset)
    msg{1}='Aborting effort to read file.';
    [msg]=ep_errorMsg(msg);
    return
end

if ~exist([EPdataset.EPwork filesep 'EPwork'],'dir')
    msg{1}='The work directory cannot be found.';
    [msg]=ep_errorMsg(msg);
    return
end

ep_tictoc('ioStart');
ep_checkEPworkCache;

if theDataset <= length(EPdataset.dataset) %replace existing entry
    delete([EPdataset.EPwork filesep 'EPwork' filesep EPdataset.dataset(theDataset).dataName '.mat']);
    inMemory=EPdataset.dataset(theDataset).inMemory;
    EPdataset.dataset=[EPdataset.dataset(1:theDataset-1) newDataset EPdataset.dataset(theDataset+1:end)];
    EPdataset.dataset(theDataset).inMemory=inMemory;
else
    EPdataset.dataset=[EPdataset.dataset newDataset];
end

EPdataset.dataset(theDataset).saved = fileSaved;

lastwarn('');
warning('off','MATLAB:save:sizeTooBigForMATFile')
eval(['save ''' EPdataset.EPwork filesep 'EPwork' filesep 'EPdataset'' EPdataset']);
warning('on')
if ~isempty(lastwarn)
    [msg,msgID] = lastwarn;
    if strcmp(msgID,'MATLAB:save:sizeTooBigForMATFile')
        save('-mat', [EPdataset.EPwork filesep 'EPwork' filesep 'EPdataset.mat'], 'EPdataset', '-v7.3');
    end
end

lastwarn('');
warning('off','MATLAB:save:sizeTooBigForMATFile')
eval(['save ''' EPdataset.EPwork filesep 'EPwork' filesep EPdata.dataName '.mat'' EPdata;']);
warning('on')
if ~isempty(lastwarn)
    [msg,msgID] = lastwarn;
    if strcmp(msgID,'MATLAB:save:sizeTooBigForMATFile')
        save('-mat', [EPdataset.EPwork filesep 'EPwork' filesep EPdata.dataName '.mat'], 'EPdata','-v7.3');
    end
end

%update what datasets are being maintained in RAM in the EPeeg structure.
if EPdataset.dataset(theDataset).inMemory==0
    EEGlist=find([EPdataset.dataset(:).inMemory]);
    if ~isempty(EEGlist)
        for iEEG=1:length(EEGlist)
            theEEG=EEGlist(iEEG);
            theSlot=EPdataset.dataset(theEEG).inMemory;
            if (theSlot>length(EPeeg)) || isempty(EPeeg(theSlot))
                EPdataset.dataset(theEEG).inMemory=0;
            else
                EPdataset.dataset(theEEG).inMemory=EPdataset.dataset(theEEG).inMemory+1;
            end
        end
        EPdataset.dataset(theDataset).inMemory=1;
        EEGlist=find([EPdataset.dataset(:).inMemory]);
        numEEG=length(EEGlist);
        if numEEG>EPmain.read.numEEG
            for iEEG=1:length(EEGlist)
                theEEG=EEGlist(iEEG);
                if EPdataset.dataset(theEEG).inMemory>EPmain.read.numEEG
                    EPdataset.dataset(theEEG).inMemory=0;
                    EPeeg{theEEG}=[];
                end
            end
        end
    else
        EPdataset.dataset(theDataset).inMemory=1;
    end
    
    EEGlist=find([EPdataset.dataset(:).inMemory]);
    numEEG=length(EEGlist);
    for iEEG=numEEG:-1:2
        EPeeg(iEEG)=EPeeg(iEEG-1);
    end
    EPeeg{1}=EPdata;
else
    EPeeg{EPdataset.dataset(theDataset).inMemory}=EPdata;
end

ep_tictoc('ioFinish');


