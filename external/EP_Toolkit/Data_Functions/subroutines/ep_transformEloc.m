function [outEloc outFid] = ep_transformEloc(inEloc, inFID, targetEloc, targetFID, inChanNames, targetChanNames, inMontage, targetMontage)
%  [outEloc outFid] = ep_transformEloc(inEloc, inFID, targetEloc, targetFID, inChanNames, targetChanNames, inMontage, targetMontage)
%       Transforms one set of electrode coordinates to match canonical space or vice versa, based on shared electrode names.
%       For EGI files, identified from the CED names, uses official EGI Technical Note (Luu and Feree, 2000) to map some electrodes onto closest 10-10 equivalents.
%
%Inputs:
%  inEloc         : The electrode location information, one for each channel, for the data (see readlocs header).
%  inFID        : The fiducial location information, one for each location, for the data (see readlocs header).
%  targetEloc     : The electrode location information, one for each channel, for the target coordinate system (see readlocs header).  If empty, will load in canonical file.
%  targetFID    : The fiducial location information, one for each location, for the target coordinate system  (see readlocs header).
%  inChanNames    : The corresponding channel names of the data elocs (which may differ from the eloc labels) not including fiducials. optional.
%  targetChanNames: The corresponding channel names of the target elocs (which may differ from the eloc labels) not including fiducials. optional.
%  inMontage          : The montage code for the data electrode coordinates.  Helps identify special cases.
%  targetMontage      : The montage code for target electrode coordinates.  Helps identify special cases.
%
%Outputs:
%  outEloc        : The electrode location information, one for each channel, rescaled to match the target eloc (see readlocs header).
%  outFid         : The fiducial location information, one for each location, rescaled to match the target eloc (see readlocs header).
%
%for EEGlab CED files, the convention is Y is +left -right, Z is anterior+ posterior-, X is dorsal+ ventral -.

%History:
%  by Joseph Dien (2/14/20)
%  jdien07@mac.com
%
% bugfix 5/9/20 JD
% Modified auto-alignment of electrode coordinates to work with EGI montages even when the fiducial locations are missing.
%
% bugfix 6/25/20 JD
% Now removes sph_theta_besa and sph_phi_besa fields.
% Fixed crash when transforming electrode locations in a .set file with internal eloc information.
%
% modified 12/23/20 JD
% Can still recognize EGI montages even when some of the electrodes have been renamed and the ced field is not a GSN montage name if at least half still start with 'E'.
%
% bugfix 1/29/21 JD
% Can still recognize EGI montages even when some channels are missing as long as they are still labeled in the EGI manner (e.g., E11) or also an equivalent (e.g., EEG 11).
% Also, can still recognize even if the file was not .mff or .set.
%
% bugfix 2/9/21 JD
% Fixed sometimes crashing when trying to read an EGI data file.
%
% modified 8/8/21 JD
% Added support for canonical 10-05 coordinates in eloc structure.
%
% bugfix 7/26/22 JD
% Updated code for recognizing alternate names for EGI reference channel.
% Fixed code so that fiducials are used as well if transforming from canonical.
% Corrected function so that it now explicitly assumes that it is transforming from or to canonical coordinates.
% Now relies on fiducials plus Cz (plus any other channel name matches) for EGI montages and only uses approximate channel matches if they are not available.
%
% bugfix 11/11/22 JD
% Fixed failing to match EGI data channels with error message about needing to match at least three channels.
%
% modified 2/8/23 JD
% Attempts to accomodate non-standard EGI channel name conventions like "chan001" rather than "E1"
%
% bugfix 11/3/23 JD
% Fixed no longer recognizing standard EGI channel names and then crashing.
%
% bugfix 5/19/24 JD
% Fixed failing to match EGI montages to other montages.
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

global EPtictoc

outEloc=[];
outFid=[];

if ~isempty(inEloc) && ~isempty(targetEloc)
        msg{1}='Error: At present ep_transformEloc can only transform to or from canonical coordinate space.';
        [msg]=ep_errorMsg(msg);
        return
end

[EPdir, ~, ~] = fileparts(which('ep.m'));

EGIdir=[EPdir filesep 'templates' filesep 'EGI_1005' filesep];

if isempty(inEloc)
    inEloc = ep_readlocsWrapper([EPdir filesep 'templates' filesep 'Standard-10-5-Cap385-VEOG.ced'],'filetype','chanedit');
    if EPtictoc.stop;return;end
    if isempty(inEloc)
        return
    end
    canonicalSpace='in';
end

if isempty(targetEloc)
    targetEloc = ep_readlocsWrapper([EPdir filesep 'templates' filesep 'Standard-10-5-Cap385-VEOG.ced'],'filetype','chanedit');
    if EPtictoc.stop;return;end
    if isempty(targetEloc)
        return
    end
    canonicalSpace='target';
end

if exist('inChanNames','var') && ~isempty(inChanNames)
    if length(inEloc) ~= length(inChanNames)
        msg{1}='Error: Number of data channel names does not match number of elocs.';
        [msg]=ep_errorMsg(msg);
        return
    end
else
    inChanNames={inEloc.labels};
end

if exist('targetChanNames','var') && ~isempty(targetChanNames)
    if length(targetEloc) ~= length(targetChanNames)
        msg{1}='Error: Number of data channel names does not match number of elocs.';
        [msg]=ep_errorMsg(msg);
        return
    end
else
    targetChanNames={targetEloc.labels};
end

inEGInumChans=[];
EGImatchListIn=cell(0); %list of corresponding 10-05 electrodes for some of the EGI channels.
if ~isempty(inMontage)    
    switch inMontage
        case 'Adult GSN 64-channel 2.0'
            [theHeader, EGImatchListIn, theDelim] = ep_textScan([EGIdir 'GSN64.txt']);
            inEGInumChans=65;
        case 'Adult GSN200 128-channel 2.1'
            [theHeader, EGImatchListIn, theDelim] = ep_textScan([EGIdir 'GSN128.txt']);
            inEGInumChans=129;
        case 'Adult GSN 256-channel 2.1'
            [theHeader, EGImatchListIn, theDelim] = ep_textScan([EGIdir 'GSN256.txt']);
            inEGInumChans=257;
        case 'Adult Hydrocel 32-channel 1.0'
            [theHeader, EGImatchListIn, theDelim] = ep_textScan([EGIdir 'Hydro32.txt']);
            inEGInumChans=33;
        case 'Adult Hydrocel 64-channel 1.0'
            [theHeader, EGImatchListIn, theDelim] = ep_textScan([EGIdir 'Hydro64.txt']);
            inEGInumChans=65;
        case 'Adult Hydrocel 128-channel 1.0'
            [theHeader, EGImatchListIn, theDelim] = ep_textScan([EGIdir 'Hydro128.txt']);
            inEGInumChans=129;
        case {'Adult Hydrocel 256-channel 1.0'}
            [theHeader, EGImatchListIn, theDelim] = ep_textScan([EGIdir 'Hydro256.txt']);
            inEGInumChans=257;
    end
    ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
end

if ~isempty(inEGInumChans)
    %if all the channel names do not conform to EGI standard (e.g., E1) or just channel n umber (e.g., 1) but all but one do conform to some kind of regular naming convention (e.g., chan001) 
    %then match the channels based on the first channel.
    eCounter=0;
    for iChan=1:length(inChanNames)
        if strcmp(inChanNames{iChan}(1),'E') || ~isnan(str2double(inChanNames{iChan}(1)))
            eCounter=eCounter+1;
        end
    end
    if eCounter == 0
        doneFlag =0;
        matchNum=0;
        while ~doneFlag
            eCounter=0;
            matchNum=matchNum+1;
            if matchNum == length(inChanNames{1})
                doneFlag=1;
            else
                for iChan=2:length(inChanNames)
                    if strcmp(inChanNames{1}(1:matchNum),inChanNames{iChan}(1:matchNum))
                        eCounter=eCounter+1;
                    end
                end
                if eCounter < (length(inChanNames)-1)
                    %finish looking for prefix if they are not pretty much all the same prefix (not including possible ref)
                    matchNum=matchNum-1;
                    doneFlag=1;
                elseif ~isnan(str2double(inChanNames{1}(matchNum+1:end)))
                    %finish looking for prefix if the remainder of the first channel translates to a number
                    doneFlag=1;
                end
            end
        end
        if matchNum > 0
            altLabel=inChanNames{1}(1:matchNum);
            disp(['Attempting to accomodate alternative EGI channel naming convention detected using prefix: ' altLabel])
            for iChan=1:length(inChanNames)
                if (length(inChanNames{iChan}) > matchNum) && strcmp(altLabel,inChanNames{iChan}(1:matchNum))
                    matchCand=str2double(inChanNames{iChan}(matchNum+1:end));
                    if ~isnan(matchCand)
                        inChanNames{iChan}=[altLabel num2str(matchCand)];
                        %change the eloc labels from the E1 convention to this alternative convention so they match
                    end
                end
            end
        end
    end
    
    %if is an EGI montage and there is no channel named "Cz" then if there is a "VREF" channel, change its name to Cz.
    %If there is neither then check for a channel of the form Ennn where nnn is the number of channels including the reference
    %and if present change it to Cz.
    theCzChan=find(strcmp('Cz',inChanNames));
    if isempty(theCzChan)
        vertexChan=find(strcmp('VREF',inChanNames));
        if ~isempty(vertexChan)
            inChanNames{vertexChan}='Cz';
        else
            altName=['E' num2str(inEGInumChans)];
            vertexChan=find(strcmp(altName,inChanNames));
            if ~isempty(vertexChan)
                inChanNames{vertexChan}='Cz';
            end
        end
    end
else
    inEGInumChans=0; %number of EEG channels in the EGI montage
end

targetEGInumChans=[];
EGImatchListTarget=cell(0); %list of corresponding 10-05 electrodes for some of the EGI channels.
if ~isempty(targetMontage)    
    switch targetMontage
        case 'Adult GSN 64-channel 2.0'
            [theHeader, EGImatchListTarget, theDelim] = ep_textScan([EGIdir 'GSN64.txt']);
            targetEGInumChans=65;
        case 'Adult GSN200 128-channel 2.1'
            [theHeader, EGImatchListTarget, theDelim] = ep_textScan([EGIdir 'GSN128.txt']);
            targetEGInumChans=129;
        case 'Adult GSN 256-channel 2.1'
            [theHeader, EGImatchListTarget, theDelim] = ep_textScan([EGIdir 'GSN256.txt']);
            targetEGInumChans=257;
        case 'Adult Hydrocel 32-channel 1.0'
            [theHeader, EGImatchListTarget, theDelim] = ep_textScan([EGIdir 'Hydro32.txt']);
            targetEGInumChans=33;
        case 'Adult Hydrocel 64-channel 1.0'
            [theHeader, EGImatchListTarget, theDelim] = ep_textScan([EGIdir 'Hydro64.txt']);
            targetEGInumChans=65;
        case 'Adult Hydrocel 128-channel 1.0'
            [theHeader, EGImatchListTarget, theDelim] = ep_textScan([EGIdir 'Hydro128.txt']);
            targetEGInumChans=129;
        case {'Adult Hydrocel 256-channel 1.0'}
            [theHeader, EGImatchListTarget, theDelim] = ep_textScan([EGIdir 'Hydro256.txt']);
            targetEGInumChans=257;
    end
    ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
end

if ~isempty(targetEGInumChans)
    %if all the channel names do not conform to EGI standard (e.g., E1) or just channel n umber (e.g., 1) but all but one do conform to some kind of regular naming convention (e.g., chan001) 
    %then match the channels based on the first channel.
    eCounter=0;
    for iChan=1:length(targetChanNames)
        if strcmp(targetChanNames{iChan}(1),'E') || ~isnan(str2double(targetChanNames{iChan}(1)))
            eCounter=eCounter+1;
        end
    end
    if eCounter == 0
        doneFlag =0;
        matchNum=0;
        while ~doneFlag
            eCounter=0;
            matchNum=matchNum+1;
            if matchNum == length(targetChanNames{1})
                doneFlag=1;
            else
                for iChan=2:length(targetChanNames)
                    if strcmp(targetChanNames{1}(1:matchNum),targetChanNames{iChan}(1:matchNum))
                        eCounter=eCounter+1;
                    end
                end
                if eCounter < (length(targetChanNames)-1)
                    %finish looking for prefix if they are not pretty much all the same prefix (not including possible ref)
                    matchNum=matchNum-1;
                    doneFlag=1;
                elseif ~isnan(str2double(targetChanNames{1}(matchNum+1:end)))
                    %finish looking for prefix if the remainder of the first channel translates to a number
                    doneFlag=1;
                end
            end
        end
        if matchNum > 0
            altLabel=targetChanNames{1}(1:matchNum);
            disp(['Attempting to accomodate alternative EGI channel naming convention detected using prefix: ' altLabel])
            for iChan=1:length(targetChanNames)
                if (length(targetChanNames{iChan}) > matchNum) && strcmp(altLabel,targetChanNames{iChan}(1:matchNum))
                    matchCand=str2double(targetChanNames{iChan}(matchNum+1:end));
                    if ~isnan(matchCand)
                        targetChanNames{iChan}=[altLabel num2str(matchCand)];
                        %change the eloc labels from the E1 convention to this alternative convention so they match
                    end
                end
            end
        end
    end
    
    %if is an EGI montage and there is no channel named "Cz" then if there is a "VREF" channel, change its name to Cz.
    %If there is neither then check for a channel of the form Ennn where nnn is the number of channels including the reference
    %and if present change it to Cz.
    theCzChan=find(strcmp('Cz',targetChanNames));
    if isempty(theCzChan)
        vertexChan=find(strcmp('VREF',targetChanNames));
        if ~isempty(vertexChan)
            targetChanNames{vertexChan}='Cz';
        else
            altName=['E' num2str(targetEGInumChans)];
            vertexChan=find(strcmp(altName,targetChanNames));
            if ~isempty(vertexChan)
                targetChanNames{vertexChan}='Cz';
            end
        end
    end
end

%add fiducials to the list of electrodes
if ~isempty(inFID)
    numFID=length(inFID);
    for iFID=1:numFID
        inEloc(end+1).X=inFID(iFID).X;
        inEloc(end).Y=inFID(iFID).Y;
        inEloc(end).Z=inFID(iFID).Z;
        inEloc(end).labels=inFID(iFID).labels;
        switch inFID(iFID).labels
            case 'FidNz'
                inChanNames{end+1}='Nz';
            case 'FidT9'
                inChanNames{end+1}='LPA';
            case 'FidT10'
                inChanNames{end+1}='RPA';
            otherwise
                inChanNames{end+1}=inFID(iFID).labels;
        end
    end
else
    numFID=0;
end

if ~isempty(targetFID)
    for iFID=1:length(targetFID)
        targetEloc(end+1).X=targetFID(iFID).X;
        targetEloc(end).Y=targetFID(iFID).Y;
        targetEloc(end).Z=targetFID(iFID).Z;
        targetEloc(end).labels=targetFID(iFID).labels;
        switch targetFID(iFID).labels
            case 'FidNz'
                targetChanNames{end+1}='Nz';
            case 'FidT9'
                targetChanNames{end+1}='LPA';
            case 'FidT10'
                targetChanNames{end+1}='RPA';
            otherwise
                targetChanNames{end+1}=targetFID(iFID).labels;
        end
    end
end

%compile list of channels with coordinate information.
chanList=zeros(length(inEloc),1);
for iChan=1:length(inEloc)
    if ~isempty(inEloc(iChan).X)
        chanList(iChan)=1;
    end
end
chanList=find(chanList); %in channels with coordinates
inEloc2=inEloc(chanList); %the in elocs with coordinates
numChans=length(inEloc2); %the number of in channels with coordinates

targetList=zeros(length(targetEloc),1);
for iChan=1:length(targetEloc)
    if ~isempty(targetEloc(iChan).X)
        targetList(iChan)=1;
    end
end
targetList=find(targetList); %target channels with coordinates
targetEloc2=targetEloc(targetList); %target channels with coordinates
numTargetChans=length(targetEloc2); %the number of target channels with coordinates

matchList=zeros(numChans,1);
for iChan=1:numChans
    theChan=chanList(iChan);
    theMatch=find(strcmpi(inChanNames{theChan},targetChanNames(targetList)));
    if~isempty(theMatch)
        matchList(iChan)=theMatch; %for the inchannels with coordinates, which of the target chans they match up to
    end
end

ep_tictoc;if EPtictoc.stop;return;end

if ~isempty(EGImatchListIn) && (length(find(matchList))<3) %if EGI data, guestimate based on best matches but only if fiducials plus Cz are too few.
    %disp('Finding best match of EGI channels to 10-05 locations.')
    EGIordered=[];
    inEchans=zeros(numChans,1); %which E channel each channel corresponds to, if any
    for iChan=1:numChans
        if (length(inChanNames{iChan})>1) && strcmpi(inChanNames{iChan}(1),'E')
            if (length(inChanNames{iChan})==2) || ~strcmpi(inChanNames{iChan}(1:3),'ECG')
                if isempty(EGIordered)
                    EGIordered=1; %at least one followed the E naming convention
                end
                if (length(inChanNames{iChan})>3) && strcmpi(inChanNames{iChan}(1:3),'EEG')
                    theChan=str2double(inChanNames{iChan}(4:end));
                    if ~isnan(theChan)
                        inEchans(iChan)=theChan;
                        if theChan~=iChan
                            EGIordered=0; %name did not correspond to the slot so out of order
                        end
                    end
                else
                    theChan=str2double(inChanNames{iChan}(2:end));
                    if ~isnan(theChan)
                        inEchans(iChan)=theChan;
                        if theChan~=iChan
                            EGIordered=0; %name did not correspond to the slot so out of order
                        end
                    end
                end
            end
        end
    end
    
    EGImatch=[]; %which of the canonical EGI channels each of the inchannels corresponds to
    if ~any(matchList) && ((numChans==inEGInumChans) || (numChans==inEGInumChans-1)) && isempty(EGIordered)
        %disp('None of the channel names matched the 10-05 names but the overall number matched, so assuming their orders correspond.')
        EGImatch=[1:numChans];
    elseif (length(find(matchList))<numChans) && ((numChans==inEGInumChans) || (numChans==inEGInumChans-1)) && EGIordered
        %disp('Only some of the channel names matched but the overall number matched and the ones that followed the E naming scheme were in the right slot, so assuming their orders correspond.')
        EGImatch=[1:numChans];
    else
        %disp('The channels appear to not be in order so will try to match them up based on the E channel naming convention.')
        for iChan=1:numChans
            if inEchans(iChan)
                EGImatch(iChan)=inEchans(iChan);
            end
        end
    end
    
    if ~isempty(EGImatch) 
        for iMatch=1:length(EGImatchListIn)
            theEloc=find(strcmpi(EGImatchListIn{iMatch,1},targetChanNames)); %the corresponding channel in the target chans
            EGIchan=str2double(EGImatchListIn{iMatch,2}); %which of the complete EGI montage it is
            theInChan=find(ismember(EGImatch,EGIchan));
            if ~isempty(theEloc) && ~isempty(theInChan)
                for iEGImatch=1:length(theInChan)
                    matchList(theInChan(iEGImatch))=theEloc; %which of the target list each of these channels correspond to
                end
            end
        end
    end
    if ~any(matchList)
        disp('It was not possible to match up the EGI channels to the 10-05 locations, perhaps because channel names and/or channel order were changed.')
    end
end

ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end

if ~isempty(EGImatchListTarget) && (length(find(matchList))<3) %if EGI data, guestimate based on best matches but only if fiducials plus Cz are too few.
    %disp('Finding best match of EGI channels to 10-05 locations.')
    EGIordered=[];
    targetEchans=zeros(numTargetChans,1); %which E channel each channel corresponds to, if any
    for iChan=1:numTargetChans
        if (length(targetChanNames{iChan})>1) && strcmpi(targetChanNames{iChan}(1),'E')
            if (length(targetChanNames{iChan})==2) || ~strcmpi(targetChanNames{iChan}(1:3),'ECG')
                if isempty(EGIordered)
                    EGIordered=1; %at least one followed the E naming convention
                end
                if (length(targetChanNames{iChan})>3) && strcmpi(targetChanNames{iChan}(1:3),'EEG')
                    theChan=str2double(targetChanNames{iChan}(4:end));
                    if ~isnan(theChan)
                        targetEchans(iChan)=theChan;
                        if theChan~=iChan
                            EGIordered=0; %name did not correspond to the slot so out of order
                        end
                    end
                else
                    theChan=str2double(targetChanNames{iChan}(2:end));
                    if ~isnan(theChan)
                        targetEchans(iChan)=theChan;
                        if theChan~=iChan
                            EGIordered=0; %name did not correspond to the slot so out of order
                        end
                    end
                end
            end
        end
    end
    
    EGImatch=[]; %which of the canonical EGI channels each of the target channels corresponds to
    if ~any(matchList) && ((numTargetChans==targetEGInumChans) || (numTargetChans==targetEGInumChans-1)) && isempty(EGIordered)
        %disp('None of the channel names matched but the overall number matched, so assuming their orders correspond.')
        EGImatch=[1:numTargetChans];
    elseif (length(find(matchList))<numTargetChans) && ((numTargetChans==targetEGInumChans) || (numTargetChans==targetEGInumChans-1)) && EGIordered
        %disp('Only some of the channel names matched but the overall number matched and the ones that followed the E naming scheme were in the right slot, so assuming their orders correspond.')
        EGImatch=[1:numTargetChans];
    else
        %disp('The channels appear to not be in order so will try to match them up based on the E channel naming convention.')
        for iChan=1:numTargetChans
            if targetEchans(iChan)
                EGImatch(iChan)=targetEchans(iChan);
            end
        end
    end
    
    if ~isempty(EGImatch) 
        for iMatch=1:length(EGImatchListTarget)
            theEloc=find(strcmpi(EGImatchListTarget{iMatch,1},inChanNames)); %the corresponding channel in the target chans
            EGIchan=str2double(EGImatchListTarget{iMatch,2}); %which of the complete EGI montage it is
            theTargetChan=find(ismember(EGImatch,EGIchan));
            if ~isempty(theEloc) && ~isempty(theTargetChan)
                for iEGImatch=1:length(theTargetChan)
                    matchList(theEloc)=theTargetChan(iEGImatch); %which of the target list each of these channels correspond to
                end
            end
        end
    end
    if ~any(matchList)
        disp('It was not possible to match up the EGI channels to the 10-05 locations, perhaps because channel names and/or channel order were changed.')
    end
end

matchData=find(matchList); %which of the in channels with coordinates have a match in the target chans
matchSfp=matchList(matchData);
numMatch=length(matchData);
if numMatch<3
    msg{1}='Error: Need at least three matching channel names to co-register the electrode coordinates.  If this is an EGI dataset, did you specify the correct montage?';
    [msg]=ep_errorMsg(msg);
    return
end
dataLocs=zeros(numChans,3);
for iChan=1:numChans
    dataLocs(iChan,1)=inEloc2(iChan).X;
    dataLocs(iChan,2)=inEloc2(iChan).Y;
    dataLocs(iChan,3)=inEloc2(iChan).Z;
end

sfpData=zeros(numTargetChans,3);
for iChan=1:numTargetChans
    sfpData(iChan,1)=targetEloc2(iChan).X;
    sfpData(iChan,2)=targetEloc2(iChan).Y;
    sfpData(iChan,3)=targetEloc2(iChan).Z;
end

%following discussion of the following citation by Nghia Ho: http://nghiaho.com/?page_id=671
% Arun, K. S., Huang, T. S., & Blostein, S. D. (1987). Least-squares fitting of two 3-D point sets. IEEE Transactions on pattern analysis and machine intelligence, (5), 698-700.
Xdata=dataLocs(matchData,:);
Ysfp=sfpData(matchSfp,:);
XdataM=mean(Xdata);
YsfpaM=mean(Ysfp);
H=(Xdata-repmat(XdataM,numMatch,1))'*(Ysfp-repmat(YsfpaM,numMatch,1));
[U,S,V]=svd(H);
R=V*U';
if det(R) < 0
    [U,S,V] = svd(R);
    V=V*diag([1 1 -1]);
    R = V * U';
end
t=YsfpaM-(R*XdataM')';

dataLocs2=(R*dataLocs'+repmat(t,numChans,1)')';
betaX = [ones(numMatch,1) dataLocs2(matchData,1)]\sfpData(matchSfp,1);
betaY = [ones(numMatch,1) dataLocs2(matchData,2)]\sfpData(matchSfp,2);
betaZ = [ones(numMatch,1) dataLocs2(matchData,3)]\sfpData(matchSfp,3);
dataLocs2=dataLocs2*diag([betaX(2) betaY(2) betaZ(2)])+[betaX(1) betaY(1) betaZ(1)];

outEloc=inEloc;
for iChan=1:numChans
    theChan=chanList(iChan);
    outEloc(theChan).X=dataLocs2(iChan,1);
    outEloc(theChan).Y=dataLocs2(iChan,2);
    outEloc(theChan).Z=dataLocs2(iChan,3);
end

outEloc = convertlocs(outEloc, 'cart2all');
[outEloc.type]=inEloc.type;

outEloc=ep_elocFormat(outEloc);

if numFID >0
    outFid=outEloc(end-numFID+1:end);
    outEloc=outEloc(1:end-numFID);
else
    outFid=[];
end