function EPdataOut=ep_addEloc(ced,eloc,fileFormat,EPdataIn, silentMode, readMode)
% EPdataOut=ep_addEloc(ced,eloc,fileFormat,EPdataIn, silentMode, readMode) -
% Adds electrode coordinate information to data file.  Mostly dealing with divergences in the electrode lists between the data file and the CED file and also sorting out reference sites.
%
%Input:
%    ced       : The name of the .ced file for electrode coordinates.  Can also have the path.
%    eloc      : The electrode location information, one for each channel (see readlocs header)
%    fileFormat: The current file format.
%    EPdataIn:  Structured array with the data and accompanying information.  See readData.
%    silentMode: 'on' means don't issue messages.
%    readMode  : Use case for this function - "read" for reading the data for the first time, "remap" for interpolating new electrode locations, and "replace" for replacing the CED information.
%Outputs
%	EPdataOut: Structured array with the data and accompanying information.  See readData.
%
%  Recognized types in the CED files include the ones listed in ep_chanTypes plus BAD (channel should be deleted), REF (EEG channel that is a reference), and FID (fiducial).
%
% This function performs the following operations:
% 1) add appropriate references.
% 2) handle BAD channels (dropping from data files).
% 3) handles mismatches between data and eloc, as in dropped channels in the data, or where the entire set of channel names is different.
% 4) handles EGI lack of standardization for reference channel name.
% 5) handles specification of REF channels (reference channels), including implicit ones.
% 6) handles specification of FID channels (fiducial channels).
%
% "replace" means change out all the coordinate information.  For the case where one used one CED and then later replaced them with a different CED (as in an improved one or fixing an error).
% If there is any kind of mismatch between the two, then drop ones that don't match.  Don't make any guesses about channel names or references.  Either they match or they don't.
% BAD is irrelevant and non-recognized Types are ignored.  Don't add implicit REF channels.  Replace FIDs.  Non-EEG channels are unaffected in the original data and ignored in the CED.  
% REG unaffected.
%
% "remap" means keep only the channels in the new CED (after interpolation) and drop all the old ones (both EEG and REG).  
% Add implicit references.  No need to guess about channel names as it is assumed they will all be different and even if they match, 
% such matches should be ignored (as in "E1" for different EGI montages).  BAD is irrelevant and non-recognized Types are ignored.  Replace FIDs.  
% Non-EEG channels are unaffected in the original data and ignored in the CED.  The new channels will be zeros, so the interpolating will have to happen elsewhere.

%History
%  by Joseph Dien (1/15/14)
%  jdien07@mac.com
%
% bugfix 2/17/14 JD
% Fixed crash when reading ced file with REF channel type indicated.
%
% bugfix 3/2/14 JD
% Fixed REF channel type not being changed to EEG for files with one reference channel.
%
% modified 3/12/14 JD
% Changed MEG channel type to MGA (axial gradiometer) and MGP (planar gradiometer) and MGM (magnetometer)
%
% bugfix 3/24/14 JD
% Fixed crash when loading file type that has internal channel names (like Neuroscan) and the only mismatch between it
% and the ced file is a single implicit REF channel or there is no mismatch and there is a single explicit REF channel.
%
% modified 4/9/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure.
%
% bugfix 5/28/14 JD
% Fixed crash when the ced file has no type field.
%
% modified 1/7/15 JD
% Handle alternate EGI names for the vertex channel when matching channel names with the ced channel names.
%
% bugfix 5/21/15 JD
% Fixed losing electrode coordinates of eeglab files when type field is empty.
%
% modified 10/9/15 JD
% When none of the CED channel names match but there are the same number of EEG
% channels, it will be assumed that they are the same channels and in the
% same order.  A warning message is provided.
%
% modified 11/22/17 JD
% If none of the channel names match but there are the same number of non-ref EEG channels, will
% not only use the CED channel names, will assume the ref channels are implicit and add them.
% Eliminated restrictions on location of CED files.
% Added support for impedances field.
%
% modified 2/11/18 JD
% Added support for stdCM field.
%
% bufix 3/26/18 JD
% Fixed crash when adding a channel and there are impedance values.
%
% bufix 11/4/18 JD
% Fixed crash if there are impedances listed but no reference.
%
% bufix & modified 12/24/19 JD
% Upgraded support of std information by adding .covAVE and .GAVsubs fields and eliminating .std and .stdCM fields.
% Fixed implicit mastoid channel a mirror of the wrong channel for case where it is mean mastoid reference with one explicit and one implicit channel and the order of the channels has bee changed by the CED file.
%
% modified 2/4/20 JD
% Converted variable structure to EPdata.
%
% modified 3/18/20 JD
% Now uses ep_chanTypes function to keep track of valid CED types.
%
% bugfix 6/7/20 JD
% Fixed crash if there is a channel in the data file that is not in the CED file.
%
% bugfix 7/25/20 JD
% Fixed fiducial channels not being included in the output file.
% Fixed crash when reading in .ced file with fiducials.
%
% bugfix 11/11/20 JD
% Fixed crash when reading in file format with fixed channels (e.g., text) and REF channels are specified in the CED file.
%
% bugfix 2/10/21 JD
% Fixed crash when reading in ced file with BAD channels.
%
% modified 4/8/21 JD
% Now handles case where just some channels names are being changed, if the number of channels missing from the CED is equal to the number that are extra.
%
% bugfix 7/25/21 JD
% Fixed not handling case where remapping all the channels to a new montage.
%
% bugfix 6/20/22 JD
% Fixed not handling case where eloc is missing and ced is specified and the file is a spatial PCA file.
% Fixed crash when there is a BAD channel specified by a non .ept file.
% Fixed not including regional types in list of acceptable eloc types for non .ept files.
% Fixed crash when two REF channels specified in the CED file.
% If there are already implicit channels, replaces them rather than adding to them.
%
% bugfix 8/1/22 JD
% Fixed not properly matching up CED with data channels when some of the channel names do not match and at least one of them was a BAD channel.
%
% bugfix 8/18/22 JD
% Fixed crash when there are more than one unique unrecognized channel types.
%
% bugfix 9/26/22 JD
% Fixed crash when there is no electrode coordinate information.
% Fixed not recognizing lowercase versions of recognized channel types.
%
% modified 2/8/23 JD
% Attempts to accomodate non-standard EGI channel name conventions like "chan001" rather than "E1"
%
% modified 5/17/24 JD
% Added readMode field to accommodate different use cases for this function.

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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global EPtictoc;

EPdataOut=[];

if ~exist('silentMode','var')
    silentMode='off';
end

if ~exist('readMode','var')
    readMode='read';
end

if isempty(eloc)
    if strcmp(ced,'internal')
        disp('Error: Electrode information is designated as being provided by the file but no electrode coordinates have been obtained.');
        return;
    end
    if exist(ced,'file')
        whichCED=ced;
        [pathstr, name, fileSuffix] = fileparts(whichCED);
        ced=[name fileSuffix];
    else
        whichCED=which(ced);
    end
    if isempty(whichCED)
        disp(['Could not find the file ' ced '.  Please put it either in the electrodes folder of the EP Toolkit or in the active directory.']);
        ced = [];
        eloc = [];
        implicit=[];
    else
        try
            if ~strcmp(silentMode,'on')
                disp(['Loading the ced file: ' whichCED]);
            end
            eloc = ep_readlocsWrapper([whichCED],'filetype','chanedit');
            if EPtictoc.stop;return;end
            if ~isempty(eloc)
                implicit=eloc(1);
                implicit(1)=[];
            end
        catch
            disp(['The ced file ' ced ' did not work for some reason.  The error message was:']);
            disp(lasterr)
            ced = [];
            eloc = [];
            implicit=[];
        end
        if isempty(eloc)
            ced = [];
            implicit=[];
        end
    end
else
    if length(eloc) ~= (length(EPdataIn.chanNames)+length(find(strcmp('FID',{eloc.type})))) && strcmp(readMode,'rad')
        disp('Error: Length of eloc array does not match number of channels in the data.');
        return;
    end
    implicit=eloc(1);
    implicit(1)=[];
end

inEGInumChans=[];
if ~isempty(EPdataIn.montage)
    switch EPdataIn.montage
        case 'Adult GSN 64-channel 2.0'
            inEGInumChans=65;
        case 'Adult GSN200 128-channel 2.1'
            inEGInumChans=129;
        case 'Adult GSN 256-channel 2.1'
            inEGInumChans=257;
        case 'Adult Hydrocel 32-channel 1.0'
            inEGInumChans=33;
        case 'Adult Hydrocel 64-channel 1.0'
            inEGInumChans=65;
        case 'Adult Hydrocel 128-channel 1.0'
            inEGInumChans=129;
        case {'Adult Hydrocel 256-channel 1.0'}
            inEGInumChans=257;
    end
end

if ~isempty(ced) && ~isempty(inEGInumChans)
    %if is an EGI ced and there is a Cz channel in the eloc but not in the
    %set of channel names, then before adding Cz to the channels, check
    %first to make sure there isn't a VREF or an Ennn (where nnn is the number of channels including the reference),
    %which are alternate EGI names for the vertex channel.  If there is,
    %then change the Cz eloc label to the alternate channel name.
    theCzChan=find(strcmp('Cz',{eloc.labels}));
    if ~isempty(theCzChan)
        if ~any(strcmp('Cz',EPdataIn.chanNames))
            vertexChan=find(strcmp('VREF',EPdataIn.chanNames));
            if ~isempty(vertexChan)
                eloc(theCzChan).labels=EPdataIn.chanNames{vertexChan};
            else
                altName=['E' num2str(inEGInumChans)];
                vertexChan=find(strcmp(altName,EPdataIn.chanNames));
                if ~isempty(vertexChan) && isempty(find(strcmp(altName,{eloc.labels})))
                    eloc(theCzChan).labels=EPdataIn.chanNames{vertexChan};
                end
            end
        end
    end
end

if isfield(eloc,'type')
    for i=1:length(eloc)
        if ~isempty(eloc(i).labels) && ~isempty(eloc(i).theta) && isempty(eloc(i).type) && ~isnumeric(eloc(i).theta)
            eloc(i).type=eloc(i).theta; %if a CED file has just the label and the type filled out, the type info migrates over to the theta column for some reason.
            eloc(i).theta=[];
        end
    end
end

if ~isempty(eloc)
    if length(eloc) ~= length(unique({eloc.labels}))
        if ~strcmp(silentMode,'on')
            disp(['The ced file ' ced ' had channel names that were not unique so ignoring it.']);
        end
        ced = [];
        eloc = [];
    end
end

if ~isempty(ced)
    if ~strcmp(silentMode,'on')
        if strcmp(ced,'internal')
            disp('Adding electrode coordinates information contained inside the data file.')
        else
            disp(['Adding electrode coordinates from CED file ' ced '.']);
        end
    end
    if ~isfield(eloc,'type')
        eloc(1).type=[];
        implicit(1).type=[];
        if ~strcmp(silentMode,'on')
            disp('Type field missing from ced file.  Will assume all channels are EEG channels.');
        end
    end
    typeFlag=0;
    for i=1:length(eloc)
        if isempty(eloc(i).type)
            eloc(i).type='EEG'; %if type fields are empty, just assume they are EEG channels.
        end
        if isnumeric(eloc(i).type)
            eloc(i).type='EEG';
            typeFlag=1;
        end
    end
    if typeFlag
        if ~strcmp(silentMode,'on')
            disp('Warning: the ''type'' field from the CED file contained numbers rather than labels like EEG or REF.');
        end
    end
    
    [chanTypes, chanModes, chanRegs]=ep_chanTypes;
    chanTypes = unique([chanTypes; chanRegs(~isempty(chanRegs))]); %add regional types to the list of acceptable eloc types.
    numREF=length(find(strcmp('REF',{eloc.type})));
    numFID=length(find(strcmp('FID',{eloc.type})));
    otherList=find(~ismember({eloc.type},[chanTypes;'BAD';'REF';'FID']));
    %change type to uppercase if needed
    if length(otherList) > 0
        fixList=[];
        for iChan=1:length(otherList)
            if any(strcmp(upper(eloc(otherList(iChan)).type),chanTypes))
                eloc(otherList(iChan)).type=upper(eloc(otherList(iChan)).type);
                fixList=[fixList;iChan];
            end
        end
        otherList(fixList)=[];
    end
    if ~isempty(otherList)
        if ~strcmp(silentMode,'on')
            disp('#########');
            disp('Warning: Unrecognized channel type in CED file; the corresponding channels will be deleted:');
            if any(strcmp('EOG',{eloc.type}))
                disp('Note - EOG channels should just be specified as being EEG channels unless they are different in kind, as in bipolar arrays, in which case they should be REG.');
            end
            badTypes=unique({eloc(otherList).type});
            for iOther=1:length(badTypes)
                disp(badTypes{iOther});
            end
            disp('The following channels have these unrecognized types and will therefore be dropped:')
            for iOther=1:length(otherList)
                disp(eloc(otherList(iOther)).labels);
            end
            if strcmp(ced,'internal')
                disp('If you wish to keep these channels, you will need to check the "no internal" box under General Preferences and then use a CED file to specify what to keep and what to drop.')
            end
            disp('#########');
        end
        for iOther=1:length(otherList)
            eloc(otherList(iOther)).type='BAD'; %drop channels with unrecognized channel types
        end
    end

    if any(strcmpi('BAD',{eloc.type})) && ~strcmp(readMode,'read')
        disp('Warning: Since this data is not the initial import, it will be assumed BAD channels have already been dropped and these channel Types will be ignored.')
    end

    if any(strcmpi('EMG',{eloc(strcmp('EEG',{eloc.type})).labels}))
        disp('Warning: EMG channels should not be typed as EEG channels.')
    end
    if any(strcmpi('EMG',{eloc(strcmp('REG',{eloc.type})).labels}))
        disp('Warning: EMG channels should not be typed as REG channels.')
    end
    if any(strcmpi('EKG',{eloc(strcmp('EEG',{eloc.type})).labels}))
        disp('Warning: EKG channels should not be typed as EEG channels.')
    end
    if any(strcmpi('EKG',{eloc(strcmp('REG',{eloc.type})).labels}))
        disp('Warning: EKG channels should not be typed as REG channels.')
    end
    if any(strcmpi('ECG',{eloc(strcmp('EEG',{eloc.type})).labels}))
        disp('Warning: ECG channels should not be typed as EEG channels.')
    end
    if any(strcmpi('ECG',{eloc(strcmp('REG',{eloc.type})).labels}))
        disp('Warning: ECG channels should not be typed as REG channels.')
    end
   
    if numFID > 0
        FIDs=find(strcmp('FID',{eloc.type}));
        temp=eloc(FIDs);
        [temp.type] = deal('FID');
        EPdataIn.implicit = temp;
        eloc=eloc(setdiff([1:length(eloc)],FIDs));
    end

    %Formats that use a fixed channel order and have no channel labels in the header
    if any(strcmp(fileFormat,{'egi_egia','egi_egis','egi_sbin','text','ns_mat'}))
        nonFID=find(ismember({eloc.type},[chanTypes;'BAD';'REF']));
        temp={eloc.labels};
        if (length(EPdataIn.chanNames) ~= length(nonFID)) && (length(EPdataIn.chanNames)+1 ~= length(nonFID))
            disp(['Error: This is a file format where the number of channels should be fixed and yet the number of channels (' num2str(length(EPdataIn.chanNames)) ') differs from the number of EEG channels in the CED file (' num2str(length(nonFID)) ').']);
            if strcmp(fileFormat,'text') && ((length(EPdataIn.timeNames)==length(nonFID)) || ((length(EPdataIn.timeNames)+1) == length(nonFID)))
                disp('Is it possible that the rows and columns were swapped?  Try changing the setting in the Files Preferences.');
            end
            disp('Giving up on electrode coordinates.');
            eloc = [];
            return
        end
        EPdataIn.chanNames=temp(nonFID(1:length(EPdataIn.chanNames)));
        temp={eloc.type};
        EPdataIn.chanTypes=temp(nonFID(1:length(EPdataIn.chanNames)));
        EPdataIn.chanNames=EPdataIn.chanNames(:);
        EPdataIn.chanTypes=EPdataIn.chanTypes(:);
    end
    
    %Formats that include channel labels in the header and where channels have flexible ordering or may even be absent
    
    %For such files channels that were left out due to being bad data or due to being an implicit reference
    %but are known to exist from the CED file will be added back to the data file.  Missing data channels
    %will be marked bad and be set to zero.  An implicit reference channel will be set to zero but not be marked
    %bad.  If there are two references (one implicit and one explicit) then it will be assumed to be a mean
    %mastoid reference or equivalent and the implicit one will be set to be the negative of the explicit one (as
    %is appropriate for such a reference scheme).  If there are two reference channels (according to the CED)
    %and both are missing (as is sometimes done) then both channels will be set to be bad as their contents
    %cannot be determined.  Fiducial channels will always be set to be implicit and there will be no
    %corresponding channel in the voltage data.
    
    chanCount=0; %channels that are in the data and also in the CED (both EEG and BAD and anything else)
    nonBadChanCount=0; %channels that are in the data and also in the CED and are not BAD
    nonCED=[]; %channels that are in the data but not in the CED
    elocLabels={eloc.labels};
    elocTypes={eloc.type};
    extraDataNames=cell(0); %channels that are in the data that are not in the CED
    elocIndex=[]; %list of data channels that are in both the data and the CED and are an actual channel rather than BAD or FID
    elocCEDindex=[]; %list of which CED channels correspond to each entry in elocIndex.
    for iChan=1:length(EPdataIn.chanNames)
        if find(strcmpi(EPdataIn.chanNames(iChan),elocLabels(find(ismember(elocTypes,[chanTypes;'REF'])))))
            chanCount=chanCount+1;
            nonBadChanCount=nonBadChanCount+1;
            elocIndex(end+1)=iChan;
            elocCEDindex(end+1)=find(strcmpi(EPdataIn.chanNames(iChan),elocLabels));
        elseif find(strcmpi(EPdataIn.chanNames(iChan),elocLabels(find(strcmp('BAD',elocTypes)))))
            chanCount=chanCount+1;
        else
            nonCED(end+1)=iChan;
            extraDataNames{end+1}=EPdataIn.chanNames{iChan};
        end
    end
    
    if strcmp(readMode,'read')
        if (chanCount ==0) && (length(EPdataIn.chanNames) == length(elocLabels(find(ismember(elocTypes,{'EEG','REF'})))))
            if ~strcmp(silentMode,'on')
                disp(['None of the channel names in the ced file ' ced ' match those in the file.']);
                disp('But the number of EEG channels are the same.  Will make the assumption that they are the same channels and in the same order.')
                disp('If this assumption is incorrect, you will need to change the channel names in either the CED file or in the data file so they match.')
                disp('Will keep the channel names in the CED file.');
            end
            chanCount=length(EPdataIn.chanNames);
            EPdataIn.chanNames=elocLabels(find(ismember(elocTypes,{'EEG','REF'})))';
            nonCED=[];
            extraDataNames=cell(0);
        elseif (chanCount ==0) && (length(EPdataIn.chanNames) == length(elocLabels(find(ismember(elocTypes,'EEG')))))
            if ~strcmp(silentMode,'on')
                disp(['None of the channel names in the ced file ' ced ' match those in the file.']);
                disp('But the number of EEG channels are the same as the non-reference channels.  Will make the assumption that they are the same channels and in the same order.')
                disp('If this assumption is incorrect, you will need to change the channel names in either the CED file or in the data file so they match.')
                disp('Will keep the channel names in the CED file and also add the missing reference channels.');
            end
            chanCount=length(EPdataIn.chanNames);
            EPdataIn.chanNames=elocLabels(find(ismember(elocTypes,'EEG')))';
            nonCED=[];
            extraDataNames=cell(0);
        elseif (chanCount ==0) && (length(EPdataIn.chanNames) == length(elocLabels(find(ismember(elocTypes,{'EEG','BAD'})))))
            if ~strcmp(silentMode,'on')
                disp(['None of the channel names in the ced file ' ced ' match those in the file.']);
                disp('But the number of EEG and BAD channels are the same as the non-reference channels.  Will make the assumption that they are the same channels and in the same order.')
                disp('If this assumption is incorrect, you will need to change the channel names in either the CED file or in the data file so they match.')
                disp('Will keep the channel names in the CED file and also add the missing reference channels.');
            end
            chanCount=length(EPdataIn.chanNames);
            EPdataIn.chanNames=elocLabels(find(ismember(elocTypes,'EEG')))';
            nonCED=[];
            extraDataNames=cell(0);
        elseif (nonBadChanCount ==0) && (length(EPdataIn.chanNames) == length(elocLabels(find(ismember(elocTypes,{'EEG','BAD'})))))
            if ~strcmp(silentMode,'on')
                disp(['None of the non-BAD channel names in the ced file ' ced ' match those in the file.']);
                disp('But the number of EEG channels are the same as the non-reference channels.  Will make the assumption that they are the same channels and in the same order.')
                disp('If this assumption is incorrect, you will need to change the channel names in either the CED file or in the data file so they match.')
                disp('Will keep the channel names in the CED file and also add the missing reference channels.');
            end
            chanCount=length(EPdataIn.chanNames);
            EPdataIn.chanNames=elocLabels(find(ismember(elocTypes,{'EEG','BAD'})))';
            nonCED=[];
            extraDataNames=cell(0);
        elseif length(extraDataNames) == (length(elocLabels(find(ismember(elocTypes,{'EEG';'BAD'}))))-chanCount)
            if ~strcmp(silentMode,'on')
                disp(['The number of extra channel names in the ced file ' ced ' match the number of missing channels names.']);
                disp('Will make the assumption that they are the same channels and in the same order.')
                disp('If this assumption is incorrect, you will need to add the missing channels to the CED file.')
                disp('Will keep the channel names in the CED file and also add the missing reference channels.');
            end
            chanCount=length(EPdataIn.chanNames);
            EPdataIn.chanNames=elocLabels(find(ismember(elocTypes,{'EEG','BAD'})))';
            nonCED=[];
            extraDataNames=cell(0);
        else
            if ~strcmp(silentMode,'on')
                if chanCount > 0
                    disp(['Some of the channel names in the ced file ' ced ' match those in the file but not all']);
                else
                    disp(['None of the channel names in the ced file ' ced ' match those in the file']);
                end
                disp('and the number of channels in the ced file and in the data are different.');
                disp('The new channels in the ced file will simply be added to the existing ones.');
            end
        end

        %drop BAD channels
        badChans=find(strcmp('BAD',elocTypes));
        if ~isempty(badChans)
            goodChans=find(~strcmp('BAD',elocTypes));
            [EPdataIn]=ep_selectData(EPdataIn,{goodChans,[],[],[],[],[]});
            elocTypes=elocTypes(goodChans);
            eloc=eloc(goodChans);
            elocLabels=elocLabels(goodChans);
        end

    elseif strcmp(readMode,'remap')
        if ~strcmp(silentMode,'on')
            disp('Remapping EEG data so all channels (even if same names) in new CED file are being added iva interpolation and all original EEG data are being dropped.')
        end
    elseif strcmp(readMode,'replace')
        if ~strcmp(silentMode,'on')
            disp('Replacing original electrode coordinate information.  The EEG data are not being modified.  EEG channels with names that do not match up are dropped.')
        end
    else
        disp('Programmer Error: readMode not recognized.');
        return;
    end

    EPadd=[];

    extraCED=0; %number of channels only in the CED file
    extraCEDNames=cell(0);
    chanIndex=find(ismember(elocTypes,[chanTypes;'REF'])); %list of CED channels that are relevant, as opposed to locations like FID that are not part of the data.
    extraElocIndex=[];
    for newElocChan=1:length(chanIndex)
        oldElocChan=chanIndex(newElocChan);
        oldFileChan=find(strcmpi(eloc(oldElocChan).labels,EPdataIn.chanNames));
        if isempty(oldFileChan)
            extraCED=extraCED+1;
            extraCEDNames{end+1}=eloc(oldElocChan).labels;
            extraElocIndex(end+1)=oldElocChan;
        end
    end

    if strcmp(readMode,'read')
        %Set up new data arrays with the full set of channels
        EPadd.chanNames=extraCEDNames;
        EPadd.eloc=eloc(extraElocIndex);
        if isempty(eloc)
            blankEloc=eloc;
        else
            blankEloc=eloc(1);
            blankEloc(1)=[];
        end
        EPdataIn.eloc=blankEloc;
        for iChan=1:length(EPdataIn.chanNames)
            elocSlot=find(strcmpi(EPdataIn.chanNames{iChan},elocLabels));
            if ~isempty(elocSlot)
                EPdataIn.eloc(iChan)=eloc(elocSlot); %add elocs for those present in both data and CED
            else
                blankEloc(1).labels=EPdataIn.chanNames{iChan};
                EPdataIn.eloc(iChan)=blankEloc;
            end
        end

        if extraCED > 0
            EPdataIn=ep_addData(EPdataIn,EPadd,'channels'); %add channels that are present only in the CED
            if isempty(EPdataIn) || isempty(EPdataIn.data)
                return
            end
        end
        %reorder data back to order in the eloc file
        sortOrder=zeros(length(EPdataIn.chanNames),1);
        dataIndex=0;
        for iChan=1:length(EPdataIn.chanNames)
            elocSlot=find(strcmpi(EPdataIn.chanNames{iChan},elocLabels));
            if isempty(elocSlot)
                dataIndex=dataIndex+1;
                elocSlot=dataIndex+length(eloc);
            end
            sortOrder(elocSlot)=iChan;
        end
        EPdataIn=ep_reorderData(EPdataIn,'channels',sortOrder);
        if isempty(EPdataIn)
            disp('Error: Data malformed in some manner.');
            return;
        end
    elseif strcmp(readMode,'remap')
        if ~strcmp(silentMode,'on')
            disp('Remapping EEG data so all channels (even if same names) in new CED file are being added iva interpolation and all original EEG data are being dropped.')
        end

        chanIndex=find(ismember(elocTypes,{'EEG','REF'})); %list of remapped EEG channels being added
        if isempty(chanIndex)
            disp('Error: no EEG channels in the new CED file.');
            return;
        end

        extraCEDNames=cell(length(chanIndex),1);
        for newChan=1:length(chanIndex)
            extraCEDNames{newChan}=eloc(chanIndex(newChan)).labels;
        end

        keepChans=find(~ismember(EPdataIn.chanTypes,['EEG';'REG']));
        dropChans=setdiff([1:length(EPdataIn.chanNames)],keepChans);

        EPadd.chanNames=extraCEDNames;
        EPadd.eloc=eloc(chanIndex);
        EPdataIn.eloc=ep_elocFormat('initialize');
        EPdataIn.eloc(length(EPdataIn.chanNames),1).type='';
        EPdataIn=ep_addData(EPdataIn,EPadd,'channels'); %add new remapped EEG channels
        if isempty(EPdataIn) || isempty(EPdataIn.data)
            return
        end
        [EPdataIn]=ep_selectData(EPdataIn,{setdiff([1:length(EPdataIn.chanNames)],dropChans),[],[],[],[],[]}); %drop existing EEG channels

        %reorder data back to order in the eloc file
        sortOrder=zeros(length(EPdataIn.chanNames),1);
        dataIndex=0;
        for iChan=1:length(EPdataIn.chanNames)
            elocSlot=find(strcmpi(EPdataIn.chanNames{iChan},elocLabels));
            if isempty(elocSlot)
                dataIndex=dataIndex+1;
                elocSlot=dataIndex+length(eloc);
            end
            sortOrder(elocSlot)=iChan;
        end
        EPdataIn=ep_reorderData(EPdataIn,'channels',sortOrder);
        if isempty(EPdataIn)
            disp('Error: Data malformed in some manner.');
            return;
        end
    elseif strcmp(readMode,'replace')
        if ~strcmp(silentMode,'on')
            disp('Replacing original electrode coordinate information.  The EEG data are not being modified.  EEG channels with names that do not match up are dropped.')
        end

        chanIndex=find(ismember(elocTypes,{'EEG','REF','REG'})); %list of EEG channels in the CED file
        if isempty(chanIndex)
            disp('Error: no EEG channels in the new CED file.');
            return;
        end

        keepChans=find(ismember(EPdataIn.chanNames(find(ismember(EPdataIn.chanTypes,{'EEG','REG'}))),{eloc(chanIndex).labels})); %EEG channel names in data that are also in the CED.
        keepChans=unique([keepChans;find(~ismember(EPdataIn.chanTypes,{'EEG','REG'}))]);
        [EPdataIn]=ep_selectData(EPdataIn,{keepChans,[],[],[],[],[]}); %drop existing EEG channels that are not represented in the CED file.
        if isempty(EPdataIn) || isempty(EPdataIn.data)
            return
        end

        if isempty(EPdataIn.eloc)
            EPdataIn.eloc=ep_elocFormat('initialize');
        end
        for iChan=1:length(EPdataIn.chanNames)
            if any(strcmp({'EEG';'REG'},EPdataIn.chanTypes{iChan}))
                theEloc=find(strcmp({eloc.labels},EPdataIn.chanNames{iChan}));
                if ~isempty(theEloc)
                    EPdataIn.eloc(iChan)=eloc(theEloc);
                else
                    EPdataIn.eloc(iChan).type='';
                end
            else
                EPdataIn.eloc(iChan).type='';
            end
        end

    else
        disp('Programmer Error: readMode not recognized.');
        return;
    end

    if any(strcmp({'read';'remap'},readMode))

        %handle reference channels
        if numREF > 2
            if ~strcmp(silentMode,'on')
                disp('Warning: File has more than two reference channels marked in CED file.  Reference channel information ignored.');
            end
        end

        if numREF < 3
            impRefs=setdiff({eloc(find(strcmp('REF',{eloc.type}))).labels},EPdataIn.chanNames); %implicit refs
            impRefNums=zeros(length(impRefs),1);
            for i=1:length(impRefs)
                impRefNums(i)=find(strcmpi(impRefs{i},{eloc.labels}));
            end
            if (numREF==2) && (length(impRefs) == 1) %if mean mastoids with one explicit and one implicit
                impMastoid=find(strcmpi(eloc(impRefNums).labels,EPdataIn.chanNames));
                expMastoid=find(strcmpi(eloc(setdiff(find(strcmp('REF',{eloc.type})),impRefNums)).labels,EPdataIn.chanNames));
                [EPdataIn]=ep_combineData(EPdataIn,'channels',[expMastoid],[-1],EPdataIn.chanNames{impMastoid});
                if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
                if isempty(EPdataIn)
                    return
                end
                EPdataIn.analysis.badChans(:,:,impMastoid)=0; %implicit reference channels are not bad
                EPdataIn.reference.original=[impMastoid expMastoid];
                EPdataIn.reference.current=EPdataIn.reference.original;
                EPdataIn.chanTypes{impMastoid}='EEG';
                EPdataIn.chanTypes{expMastoid}='EEG';
                if ~strcmp(silentMode,'on')
                    disp('CED file indicates was originally mean mastoids reference and will assume is still so.');
                end
            elseif (numREF==2) && (length(impRefs) == 2) %if mean mastoids and both were implicit
                % both ref channels will be left as being marked bad as their voltage values will be missing
                M1=find(strcmpi(EPdataIn.eloc(impRefNums).labels,EPdataIn.chanNames));
                M2=find(strcmpi(EPdataIn.eloc(setdiff(find(strcmp('REF',{EPdataIn.eloc.type})),M1)).labels,EPdataIn.chanNames));
                EPdataIn.reference.original=[M1 M2];
                EPdataIn.reference.current=EPdataIn.reference.original;
                EPdataIn.chanTypes{M1}='EEG';
                EPdataIn.chanTypes{M2}='EEG';
                if ~isempty(EPdataIn.impedances.channels)
                    EPdataIn.impedances2.channels(M1,:)=EPdataIn.impedances.reference(1);
                    EPdataIn.impedances2.channels(M2,:)=EPdataIn.impedances.reference(2);
                end
                if ~strcmp(silentMode,'on')
                    disp('CED file indicates was originally mean mastoids reference and will assume is still so.');
                end
            elseif (numREF==2) && (length(impRefs) == 0) %if mean mastoids and both were explicit
                references=find(strcmp('REF',{eloc.type}));
                M1=references(1);
                M2=references(2);
                EPdataIn.reference.original=[M1 M2];
                EPdataIn.reference.current=EPdataIn.reference.original;
                EPdataIn.chanTypes{M1}='EEG';
                EPdataIn.chanTypes{M2}='EEG';
                if ~strcmp(silentMode,'on')
                    disp('CED file indicates was originally physically linked mastoids reference and will assume is still so (but incorrect if was rereferenced to mean mastoids after data collection).');
                end
            elseif (numREF==1) && (length(impRefs) == 1) %if a single implicit reference
                M1=find(strcmpi(EPdataIn.eloc(impRefNums).labels,EPdataIn.chanNames));
                EPdataIn.reference.original=M1;
                EPdataIn.reference.current=EPdataIn.reference.original; %assume that if there are implicit references then they were the original reference and still are
                EPdataIn.chanTypes{M1}='EEG';
                EPdataIn.analysis.badChans(:,:,M1)=0; %implicit reference channels are not bad
                if ~isempty(EPdataIn.impedances.channels)
                    if isfield(EPdataIn.impedances,'reference')
                        EPdataIn.impedances.channels(M1,:)=EPdataIn.impedances.reference(1);
                    end
                end
                if ~isempty(EPdataIn.covAVE)
                    if size(EPdataIn.covAVE,7)==1
                        EPdataIn.covAVE(end+1,:,:,:,:,:,:)=zeros(size(EPdataIn.covAVE,2),size(EPdataIn.covAVE,3),size(EPdataIn.covAVE,4),size(EPdataIn.covAVE,5),size(EPdataIn.covAVE,6),1);
                    else
                        EPdataIn.covAVE(end+1,:,:,:,:,:,:)=zeros(1,size(EPdataIn.covAVE,2),size(EPdataIn.covAVE,3),size(EPdataIn.covAVE,4),size(EPdataIn.covAVE,5),size(EPdataIn.covAVE,6),size(covAVE,1)+1);
                    end
                end
            elseif (numREF==1) && isempty(impRefs) %if a single explicit reference
                theRef=find(strcmpi(EPdataIn.eloc(find(strcmp('REF',{EPdataIn.eloc.type}))).labels,EPdataIn.chanNames));
                EPdataIn.reference.original=theRef;

                if strcmpi(EPdataIn.chanNames{theRef},'CMS')
                    disp('Warning: The CMS channel from BioSemi systems is not a conventional reference channel and should be dropped entirely from the data.');
                end

                EEGchans=find(ismember(EPdataIn.chanTypes,{'EEG','REF'}));
                [expandedData]=ep_expandFacs(EPdataIn);
                ep_tictoc;if EPtictoc.stop;return;end
                if ~squeeze(any(any(any(any(any(expandedData(theRef,:,:,:,:,:)))))))
                    EPdataIn.reference.current=EPdataIn.reference.original;
                elseif ~squeeze(any(any(any(any(sum(expandedData(EEGchans,:,:,:,:,:,:),1))))))
                    EPdataIn.reference.current=[];
                    EPdataIn.reference.type='AVG';
                else
                    EPdataIn.reference.current=[];
                    EPdataIn.reference.type='REG';
                end
                EPdataIn.chanTypes{theRef}='EEG';
            elseif numREF==0
                if ~strcmp(silentMode,'on')
                    disp('Note: No reference channel was found.');
                end
                EEGchans=find(ismember(EPdataIn.chanTypes,{'EEG','REF'}));
                [expandedData]=ep_expandFacs(EPdataIn);
                ep_tictoc;if EPtictoc.stop;return;end
                if ~squeeze(any(any(any(any(sum(expandedData(EEGchans,:,:,:,:,:,:),1))))))
                    EPdataIn.reference.current=[];
                    EPdataIn.reference.type='AVG';
                end
            else
                if ~strcmp(silentMode,'on')
                    disp('Warning: EP Toolkit was unable to match up contents of the CED file to the data file.');
                end
            end
        else
            impRefs=[];
            impRefNums=[];
        end
    end

    % if strcmp(readMode,'remap')
    %     if ~isempty(nonCED)
    %         if ~strcmp(silentMode,'on')
    %             disp(['The following data file channels are not represented in the ced file:']);
    %             for i=1:length(extraDataNames)
    %                 disp(extraDataNames{i});
    %             end
    %         end
    %     end
    % end

    if any(strcmp({'read';'replace'},readMode))
        if extraCED > 0
            if ~strcmp(silentMode,'on')
                disp(['The following ced file channels are not represented in the data file:']);
                for i=1:length(extraCEDNames)
                    disp(extraCEDNames{i});
                end
            end
        end
    end
end
% end

EPdataOut=EPdataIn;
