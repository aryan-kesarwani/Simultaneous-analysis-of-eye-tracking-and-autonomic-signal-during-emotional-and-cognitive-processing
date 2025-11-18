function EPdata=ep_loadEPdataset(theDataset,regenFlag)
%  EPdata=ep_loadEPdataset(theDataset,regenFlag);
%       Loads a dataset from the work directory.
%
%Inputs:
%  theDataset     : Which dataset to load in.
%  regenFlag      : flag to indicate that cache is already known to be bad and in need of regenerating.
%
%Outputs:
%  EPdata         : Structured array with the data and accompanying information.  See readData.

%History:
%  by Joseph Dien (7/11/09)
%  jdien07@mac.com
%
%  bugfix 8/3/09 JD
%  File names with multiple periods would confuse Matlab about the file type,
%  resulting in it trying to load the file as an ASCII file rather than as a .mat file.
%
%  modified 2/24/10 JD
%  Backward compatibility for changes in EPdata format.
%
%  modified 10/12/10 JD
%  Backward compatibility for shift to epochwise format for analysis fields for continuous data files.
%
%  modified 10/17/10 JD
%  Added support for saccadeTrials and saccadeOnset fields.
%
% modified 4/19/11 JD
% Added support for transforms field
%
% modified 2/5/12 JD
% Eliminated transforms field and added instead freqNames and facVecF fields
%
% modified 2/22/12 JD
% Backward compatibility for noise and std no longer optional.
%
% modified 1/11/13 JD
% Added option to do internal calculations of frequency data in either amplitude or power form.
%
% modified 10/13/13 JD
% Added recTime field.
%
% bugfix 10/21/13 JD
% Ensures that power comes after analysis field so order of fields is always the same.
%
% modified 10/29/13 JD
% Added .keys field to events.
%
% bugfix 1/7/14 JD
% Added check for .recTime field in .pca field.
%
% modified 2/26/14 JD
% Made .pca field obligatory.  Reorders fields into a standard order.
%
% bufix 3/12/14 JD
% Handles decimal sampling rates gracefully.
%
% modified 3/12/14 JD
% Changed MEG channel type to MGA (axial gradiometer) and MGP (planar gradiometer) and MGM (magnetometer)
%
% bufix 3/19/14 JD
% Fixed recTime field not including space for the 'all' cell in PCA output, resulting in crashes when edited.
%
% modified 3/24/14 JD
% Added cov field.
%
% bufix 4/8/14 JD
% Fixed not putting factor variance information in correct location when loading PCA .ept files from older versions of
% the EP Toolkit, resulting in "There are 0 factors that meet the minimum variance criterion" error messages when trying
% to autoPCA them.
%
% modified 4/9/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure.
%
% modified 4/17/14 JD
% Added .covNum field.
%
% bufix 4/25/14 JD
% Added conversion of REF channel type to EEG for older EP files.
%
% bufix 6/12/14 JD
% Fixed blank keys field of events being produced without .key (e.g., .keys.keyCod instead of .keys.key.keyCode)
%
% modified 7/16/14 JD
% Simplified keys code field structure.
%
% bufix 8/27/14 JD
% Fixed crash when accessing a dataset in the working set with a different number of fields
% than that of the current EP format.
%
% modified 3/8/16 JD
% Eliminated support for .power field.
%
% modified 10/16/16 JD
% Added .stims field.
%
% modified 11/13/16 JD
% Added .calibration field.
%
% modified 6/13/17 JD
% Added .timeUnits field.
%
% modified 11/22/17 JD
% Added support for impedances field.
%
% modified 2/11/18 JD
% Added support for stdCM field.
%
% bugfix 3/18/18 JD
% Fixed crash when .ept file has event with empty key.
%
% bugfix 6/12/18 JD
% Fixed not ensuring subject spec names and trial spec names are column vectors.
%
% modified 1/8/19 JD
% FacVar and FacVarQ now include CMB factors.
%
% modified 3/30/19 JD
% Added support for task level performance measures.
%
% modified 6/15/19 JD
% Now maintains a copy of the most recent datasets from the working set in RAM to minimize disk access time.
% Now uses EPdataset via global variable or will load it in from disk if not available.
%
% bugfix & modified 11/4/19 JD
% Added sessNums sessNames fields.
% Fixed endless loop when regenerating the cache.
%
% bugfix 1/21/20 JD
% Fixed error check unnecessarily resetting EPeeg cache.
%
% bugfix 3/12/20 JD
% Fixed crash in some situations when deleting member of EPeeg cache.
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

EPdata=[];

global EPeeg EPdataset EPmain;

if ~exist('regenFlag','var')
    regenFlag=0;
end

if isempty(theDataset)
    msg{1}='Error: which dataset to load in was not specified.';
    [msg]=ep_errorMsg(msg);
    return
end

if length(theDataset) > 1
    msg{1}='Error: more than one dataset specified.';
    [msg]=ep_errorMsg(msg);
    return
end

if ~exist('EPdataset','var')
    ep_tictoc('ioStart');
    try
        eval(['load ''' EPdataset.EPwork filesep 'EPwork' filesep 'EPdataset.mat'';']);
    catch ME
        if strcmp(ME.identifier,'MATLAB:load:unableToReadMatFile')
            EPdataset=regenerateEPload([EPdataset.EPwork filesep 'EPwork' filesep 'EPdataset.mat']);
        else
            ep_tictoc('error');
            msg{1}='Error: Working directory has been corrupted.  Delete the contents of the EPwork directory, except for the EPprefs file, and restart the EP Toolkit.';
            [msg]=ep_errorMsg(msg);
            beep
            return
        end
    end
    ep_tictoc('ioFinish');
end

if ~regenFlag
    ep_checkEPworkCache
end

if theDataset > length(EPdataset.dataset)
    msg{1}='Error: there is no dataset corresponding to the one requested.';
    [msg]=ep_errorMsg(msg);
    return
end

if ~exist('EPeeg','var')
    EPeeg=cell(0);
else
    for iDataset=1:length(EPdataset.dataset)
        if EPdataset.dataset(iDataset).inMemory ~=0
            eegDataset=EPdataset.dataset(iDataset).inMemory;
            if eegDataset > 0
                if (eegDataset> length(EPeeg)) || isempty(EPeeg{eegDataset}) || ~isfield(EPeeg{eegDataset},'dataName')|| ~strcmp(EPdataset.dataset(iDataset).dataName,EPeeg{eegDataset}.dataName)
                    disp('Reinitializing RAM copies, probably due to restarting Matlab.  Should not affect user other than to require a pause while fresh copy is loaded from working set.')
                    EPeeg=cell(0);
                    [EPdataset.dataset(:).inMemory]=deal(0);
                end
            end
        end
    end
end

eegSlot=EPdataset.dataset(theDataset).inMemory;
if isempty(eegSlot)
    EPdataset.dataset(theDataset).inMemory=0;
    eegSlot=0;
end
if (eegSlot==0) || (eegSlot>length(EPeeg)) || isempty(EPeeg(eegSlot))
    if ~regenFlag
        ep_tictoc('ioStart');
    end
    try
        eval(['load ''' EPdataset.EPwork filesep 'EPwork' filesep EPdataset.dataset(theDataset).dataName '.mat'';']);
    catch
        if ~regenFlag
            ep_tictoc('error');
        end
        msg{1}='Error: Working directory has been corrupted.  Delete the contents of the EPwork directory, except for the EPprefs file, and restart the EP Toolkit.';
        [msg]=ep_errorMsg(msg);
        beep
        return
    end
    if ~regenFlag
        ep_tictoc('ioFinish');
    end

    EPdata=ep_updateEPfile(EPdata,regenFlag);
    
    refChans=find(strcmp('REF',EPdata.chanTypes)); %assume that REF normally indicates original reference channels
    if ~isempty(refChans)
        EPdata.chanTypes{refChans}='EEG'; %assume all REF channels are EEG channels
    end
            
    %update what datasets are being maintained in RAM in the EPeeg structure.
    EEGlist=find([EPdataset.dataset(:).inMemory]);
    if ~isempty(EEGlist)
        for iEEG=1:length(EEGlist)
            theEEG=EEGlist(iEEG);
            eegSlot=EPdataset.dataset(theEEG).inMemory;
            if (length(EPeeg) < eegSlot) || isempty(EPeeg(eegSlot))
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
    EPdata = EPeeg{EPdataset.dataset(theDataset).inMemory};
    EPdata=ep_updateEPfile(EPdata,regenFlag);
end








