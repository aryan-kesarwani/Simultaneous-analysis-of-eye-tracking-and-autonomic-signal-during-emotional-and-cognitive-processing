function EPdataOut=ep_updateEPfile(EPdataIn,regenFlag)
%  EPdataOut=ep_updateEPfile(EPdataIn,regenFlag)
%       Updates EP file to current format.
%
%Inputs:
%  EPdataIn      : Structured array with the data and accompanying information.  See readData.
%  regenFlag      : flag to indicate that cache is already known to be bad and in need of regenerating.
%
%Outputs:
%  EPdataOut     : Structured array with the data and accompanying information.  See readData.

%History:
%  by Joseph Dien (1/9/20)
%  jdien07@mac.com
%
% bugfix 2/14/20 JD
% Fixed crash when loading old ept file without reference field.
%
% modified 4/9/20 JD
% Now enforces standard format for eloc and implicit fields.
%
% modified 10/15/20 JD
% Now handles case where member of working set has extra fields that are not part of the EPdata standard.
%
% modified 8/8/21 JD
% Added support for canonical 10-05 coordinates in eloc structure.
%
% modified 5/8/22 JD
% Added backward compatibility updating to new history field format.
%
% bugfix 11/20/24 JD
% Fixed problem with .std to .covAVE backward compatibility conversion that
% was preventing old .ept files from loading.
%
% modified 4/18/25 JD
% Added video field.
% Added support for virtual grand averages.
%
% bugfix 5/23/25 JD
% Fixed crash when updating ept files with 2.997 form of virtual subjects.
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

global EPtictoc EPmain

EPdataOut=EPdataIn;

if ~exist('regenFlag','var')
    regenFlag=0;
end



if ~isfield(EPdataOut,'freqNames')
    EPdataOut.freqNames=[];
end
if ~isfield(EPdataOut,'facVecF')
    EPdataOut.facVecF=[];
end

if ~isfield(EPdataOut,'pca')
    EPdataOut.pca=[];
end

if ~isfield(EPdataOut.pca,'stims')
    EPdataOut.pca.stims=struct('name',{},'image',{},'AOI',{});
end

if ~isfield(EPdataOut.pca,'calibration')
    EPdataOut.pca.calibration=[];
end

if ~isfield(EPdataOut,'facVar')
    if isfield(EPdataOut.pca,'facVar')
        EPdataOut.facVar=EPdataOut.pca.facVar;
    else
        EPdataOut.facVar=[];
    end
end

if ~isfield(EPdataOut,'facVarQ')
    if isfield(EPdataOut.pca,'facVarQ')
        EPdataOut.facVarQ=EPdataOut.pca.facVarQ;
    else
        EPdataOut.facVarQ=[];
    end
end

if ~isfield(EPdataOut,'relNames')
    EPdataOut.relNames=[];
end

if ~isfield(EPdataOut,'stims') || isempty(EPdataOut.stims)
    EPdataOut.stims=struct('name',{},'image',{},'AOI',{});
else
    if ~isfield(EPdataOut.stims,'AOI')
        EPdataOut.stims.AOI=cell(0);
    end
end

if ~isfield(EPdataOut,'timeUnits')
    EPdataOut.timeUnits=[];
end

refChans=find(strcmp('REF',EPdataOut.chanTypes)); %assume that REF normally indicates original reference channels
if ~isfield(EPdataOut,'reference')
    EPdataOut.reference.original=[];
    EPdataOut.reference.current=[];
    EPdataOut.reference.type='';
    
    if isempty(refChans) && ~isempty(EPdataOut.eloc)
        refChans=find(strcmp('REF',{EPdataOut.eloc.type}));
    end
    
    if ~isempty(EPdataOut.freqNames) && ~isempty(EPdataOut.facNames) %if not spectral data or factor data
        if ~any(mean(EPdataOut.data(strcmp('EEG',EPdataOut.chanTypes),:,:,:,:,:),2))
            EPdataOut.reference.type='AVG'; %channels sum to zero so must be average reference
        end
    end
    
    if ~isempty(refChans)
        if length(refChans) > 2
            disp('More than two reference channels indicated.  Will ignore those past first two.');
            refChans=refChans(1:2);
        end
        EPdataOut.reference.original=refChans;
        EPdataOut.chanTypes{refChans}='EEG'; %assume all REF channels are EEG channels
        if ~any(sum(ep_expandFacs(EPdataOut,refChans,[],[],[],[],[]),2))
            ep_tictoc;if EPtictoc.stop;return;end
            EPdataOut.reference.current=refChans; %reference channels still sum to zero so still being reference channels
        end
    end
else
    if ~isempty(refChans)
        EPdataOut.chanTypes{refChans}='EEG'; %assume all REF channels are EEG channels
    end
end

if isfield(EPdataOut,'power')
    if strcmp(EPdataOut.power,'Power')
        EPdataOut=sqrt(EPdataOut.data);
        disp('Converting power data to amplitude data.');
    end
    EPdataOut=rmfield(EPdataOut,'power');
end

if ~isfield(EPdataOut,'recTime')
    EPdataOut.recTime=[1:EPdataOut.Fs:EPdataOut.Fs*(length(EPdataOut.cellNames)-1)+1]';
elseif ~isempty(EPdataOut.facNames)
    if isfield(EPdataOut,'GAVsubs')
        vCells=max(0,size(EPdataOut.GAVsubs,2)-1);
    else
        vCells=0;
    end
    if length(EPdataOut.recTime) == (length(EPdataOut.cellNames)-1-vCells)
        EPdataOut.recTime(end+1)=0; %repair effects of bug
    end
end

if ~regenFlag
    ep_tictoc;if EPtictoc.stop;return;end
end
if ~isempty(EPdataOut.events)
    for i=1:size(EPdataOut.events,1)
        for k=1:size(EPdataOut.events,2)
            if ~isempty(EPdataOut.events{i,k})
                if ~isfield(EPdataOut.events{i,k},'keys')
                    EPdataOut.events{i,k}(1).keys=struct('code','','data','','datatype','','description','');
                elseif isfield(EPdataOut.events{i,k}(1).keys,'key')
                    if isfield(EPdataOut.events{i,k}(1).keys(1).key,'keyCode')
                        for iEvent=1:length(EPdataOut.events{i,k})
                            for iKey=1:length(EPdataOut.events{i,k}(iEvent).keys)
                                if ~isempty(EPdataOut.events{i,k}(iEvent).keys(iKey).key)
                                    EPdataOut.events{i,k}(iEvent).keys(iKey).code=EPdataOut.events{i,k}(iEvent).keys(iKey).key.keyCode;
                                    EPdataOut.events{i,k}(iEvent).keys(iKey).data=EPdataOut.events{i,k}(iEvent).keys(iKey).key.data.data;
                                    EPdataOut.events{i,k}(iEvent).keys(iKey).datatype=EPdataOut.events{i,k}(iEvent).keys(iKey).key.data.dataType;
                                    EPdataOut.events{i,k}(iEvent).keys(iKey).description='';
                                end
                            end
                        end
                    end
                elseif (~isempty(EPdataOut.events{i,k}(1).keys)) && isfield(EPdataOut.events{i,k}(1).keys(1),'keyCode')
                    for iEvent=1:length(EPdataOut.events{i,k})
                        newKeys=[];
                        for iKey=1:length(EPdataOut.events{i,k}(iEvent).keys)
                            newKeys(iKey).code=EPdataOut.events{i,k}(iEvent).keys(iKey).keyCode;
                            newKeys(iKey).data=EPdataOut.events{i,k}(iEvent).keys(iKey).data.data;
                            if isfield(EPdataOut.events{i,k}(iEvent).keys(iKey).data,'dataType')
                                newKeys(iKey).datatype=EPdataOut.events{i,k}(iEvent).keys(iKey).data.dataType;
                            else
                                newKeys(iKey).datatype='';
                            end
                            if isfield(EPdataOut.events{i,k}(iEvent).keys(iKey).data,'description')
                                newKeys(iKey).description=EPdataOut.events{i,k}(iEvent).keys(iKey).data.description;
                            else
                                newKeys(iKey).description='';
                            end
                        end
                        EPdataOut.events{i,k}(iEvent).keys=newKeys;
                    end
                end
            end
        end
    end
end

if ~regenFlag
    ep_tictoc;if EPtictoc.stop;return;end
end
MEGchans=find(strcmp('MEG',EPdataOut.chanTypes));
if ~isempty(MEGchans)
    EPdataOut.chanTypes(MEGchans)='MGA';
    disp('Converting MEG channel types to MGA (axial gradiometer MEG).  If any of the MEG channels are actually planar or magnetometers, you will need to use the Edit function to correct them.');
end

if ~isfield(EPdataOut,'covNum') || isempty(EPdataOut.covNum)
    EPdataOut.covNum=EPdataOut.avgNum;
end

if xor(isempty(EPdataOut.trialSpecNames),isempty(EPdataOut.trialSpecs))
    EPdataOut.trialSpecNames=[];
    EPdataOut.trialSpecs=[];
end

if ~isfield(EPdataOut,'calibration')
    EPdataOut.calibration=[];
end

if ~isfield(EPdataOut,'impedances')
    EPdataOut.impedances=[];
end
if ~isfield(EPdataOut.impedances,'channels')
    EPdataOut.impedances.channels=[];
end
if ~isfield(EPdataOut.impedances,'ground')
    EPdataOut.impedances.ground=[];
end
if ~isfield(EPdataOut,'video')
    EPdataOut.video=struct('frames',[],'times',[]);
    EPdataOut.video(1)=[];
elseif isscalar(EPdataOut.video) && isfield(EPdataOut.video,'frames') && isempty(EPdataOut.video.frames) && isfield(EPdataOut.video(1).frames,'cdata')
    EPdataOut.video=struct('frames',[],'times',[]);
    EPdataOut.video(1)=[];
end

if length(EPdataOut.facVar) < length(EPdataOut.facNames)
    EPdataOut.facVar(length(EPdataOut.facNames))=0;
end

if length(EPdataOut.facVarQ) < length(EPdataOut.facNames)
    EPdataOut.facVarQ(length(EPdataOut.facNames))=0;
end

if ~isfield(EPdataOut,'taskSpecs')
    EPdataOut.taskSpecs=[];
end
if ~isfield(EPdataOut,'taskNames')
    EPdataOut.taskNames=cell(0);
end
if ~isfield(EPdataOut,'taskMeasNames')
    EPdataOut.taskMeasNames=cell(0);
end


if ~isfield(EPdataOut,'sessNames')
    EPdataOut.sessNames=cell(0);
end
if ~isfield(EPdataOut,'sessNums')
    EPdataOut.sessNums=[];
end

if isfield(EPdataOut,'std')
    if ~isempty(EPdataOut.std) && ~isfield(EPdataOut,'covAVE')
        disp('Note: Converting contents of .std field to .covAVE field.  covAVE will not have covariance information so rereferencing and Edit operations will not update standard deviation information correctly.');
        for iChan=1:length(EPdataOut.chanNames)
            if ~regenFlag
                ep_tictoc;if EPtictoc.stop;return;end
            end
            if size(EPdataOut.std,6)==1
                EPdataOut.covAVE(iChan,:,:,:,:)=EPdataOut.std(iChan,:,:,:,:).^2;
            else
                EPdataOut.covAVE(iChan,:,:,:,:,:)=EPdataOut.std(iChan,:,:,:,:,:).^2;
            end
        end
    end
    EPdataOut=rmfield(EPdataOut,'std');
end

if isfield(EPdataOut,'stdCM')
    EPdataOut=rmfield(EPdataOut,'stdCM');
end

if ~isfield(EPdataOut,'covAVE')
    EPdataOut.covAVE=[];
end

if isfield(EPdataOut,'covGAV')
    EPdataOut=rmfield(EPdataOut,'covGAV');
end

if ~isfield(EPdataOut,'GAVsubs') || isempty(EPdataOut.GAVsubs)
    %standardize format of empty GAVsubs
    EPdataOut.GAVsubs=[];
end
%update GAVsubs from original format of cell array (subs,1) with each being empty except for GAVs, where it listed the subjects.
if ~isempty(EPdataOut.GAVsubs) && (size(EPdataOut.GAVsubs,2)==1) && (size(EPdataOut.GAVsubs,1)==length(EPdataOut.subNames))
    vGAVlist=find(~cellfun(@isempty,EPdataOut.GAVsubs) & strcmp(EPdataOut.subTypes,'GAV'));
    convertSubNumList=zeros(length(EPdataOut.subNames),1);
    subCounter=1;
    subNames=cell(0);
    subTypes=cell(0);
    GAVsubs=[];
    preNumSubs=length(EPdataOut.subNames);
    for iSub=1:preNumSubs
        convertSubNumList(iSub)=iSub;
        if ismember(vGAVlist,iSub)
            if any(EPdataOut.GAVsubs{iSub}==iSub)
                disp(['Error: ' EPdataOut.subNames{iSub} '''s list of subjects corrupted.  Will not convert to virtual GAVE.'])
                EPdataOut.GAVsubs{iSub}=cell(0);
                vGAVlist=setdiff(vGAVlist,iSub);
                convertSubNumList(iSub)=subCounter;
                subCounter=subCounter+1;
            else
                subNames{end+1}=EPdataOut.subNames{iSub};
                subTypes{end+1}='GAV';
                GAVsubs{end+1}=EPdataOut.GAVsubs{iSub};
                convertSubNumList(iSub)=length(EPdataOut.subNames)-length(vGAVlist)+length(subNames);
            end
        else
            convertSubNumList(iSub)=subCounter;
            subCounter=subCounter+1;
        end
    end

    if isempty(vGAVlist)
        EPdataOut.GAVsubs=[];
    else
        %move virtual grand averages to the end of the subject list and then delete the original entries.
        EPdataOut.subNames(end+1:end+length(subNames))=subNames;
        EPdataOut.subTypes(end+1:end+length(subTypes))=subTypes;
        %update GAVsubs to new format with one for each cell to accommodate trimmed GAVs.
        GAVsubs2=cell(length(GAVsubs)*2+1,1,length(EPdataOut.facNames));
        GAVsubs2(2:length(GAVsubs)+1,1,1)=GAVsubs;
        GAVsubs2(2+length(GAVsubs):end,1,1)=GAVsubs;
        GAVsubs3=GAVsubs2;
        for iGAV=2:size(GAVsubs2,1)
            for iFac=1:max(1,length(EPdataOut.facNames))
                %add default weighting of 1 to each entry
                if ~isempty(GAVsubs2{iGAV,1,iFac})
                    GAVsubs3{iGAV,1,iFac}=[GAVsubs2{iGAV,1,iFac}(:,1), ones(size(GAVsubs2{iGAV,1,iFac}(:,1),1),1)];
                else
                    GAVsubs3{iGAV,1,iFac}=[GAVsubs2{iGAV,1,1}(:,1), ones(size(GAVsubs2{iGAV,1,1}(:,1),1),1)];
                end
            end
        end
        EPdataOut.GAVsubs=GAVsubs3; %GAVsubs now contains only the virtual GAVEs
        [EPdataOut]=ep_selectData(EPdataOut,{[],[],[],setdiff([1:length(EPdataOut.subNames)],vGAVlist),[],[]});
        if isempty(EPdataOut)
            return;
        end

        %update the GAVsubs entries to the new numbering.
        for iGAV=2:size(EPdataOut.GAVsubs,1)
            for iFac=1:length(EPdataOut.facNames)
                if ~isempty(EPdataOut.GAVsubs{iGAV,1,iFac})
                    for iRow=1:length(EPdataOut.GAVsubs{iGAV,1,iFac})
                        EPdataOut.GAVsubs{iGAV,1,iFac}(iRow)=convertSubNumList(EPdataOut.GAVsubs{iGAV,1,iFac}(iRow));
                    end
                end
            end
        end

        %handle GAV referencing other GAVs
        for iGAV=2:size(EPdataOut.GAVsubs,1)
            for iFac=1:length(EPdataOut.facNames)
                subGAVE=[];
                deleteList=[];
                if ~isempty(EPdataOut.GAVsubs{iGAV,1,iFac})
                    for iRow=1:length(EPdataOut.GAVsubs{iGAV,1,iFac})
                        if EPdataOut.GAVsubs{iGAV,1,iFac}(iRow) > (length(EPdataOut.subNames)-size(EPdataOut.GAVsubs,1)+1)
                            %if a subGAVE entry is itself a subGAVE, then substitute its entries.  Assume can only be one recursion.
                            subGAVE=[subGAVE; EPdataOut.GAVsubs{length(EPdataOut.subNames)-EPdataOut.GAVsubs{iGAV,1,iFac}(iRow)}];
                            deleteList=[deleteList;iRow];
                        end
                    end
                    if ~isempty(deleteList)
                        EPdataOut.GAVsubs{iGAV,1,iFac}(deleteList)=[];
                        EPdataOut.GAVsubs{iGAV,1,iFac}=unique([EPdataOut.GAVsubs{iGAV,1,iFac}; subGAVE]);
                    end
                end
            end
        end
    end
end

if isfield(EPdataOut,'analysis')
    if isfield(EPdataOut.analysis,'blinkTrial')
        if strcmp(EPdataOut.dataType,'continuous')
            if (size(EPdataOut.data,2) >= EPdataOut.Fs*2) && size(EPdataOut.analysis.blinkTrial,2) ==1
                EPdataOut.analysis.blinkTrial = []; %backward compatibility conversion
            end
        end
    end
    if isfield(EPdataOut.analysis,'saccadeTrial')
        if strcmp(EPdataOut.dataType,'continuous')
            if (size(EPdataOut.data,2) >= EPdataOut.Fs*2) && size(EPdataOut.analysis.saccadeTrial,2) ==1
                EPdataOut.analysis.saccadeTrial = []; %backward compatibility conversion
            end
        end
    end
    if isfield(EPdataOut.analysis,'saccadeOnset')
        if strcmp(EPdataOut.dataType,'continuous')
            if (size(EPdataOut.data,2) >= EPdataOut.Fs*2) && size(EPdataOut.analysis.saccadeOnset,2) ==1
                EPdataOut.analysis.saccadeOnset = []; %backward compatibility conversion
            end
        end
    end
    if isfield(EPdataOut.analysis,'moveTrial')
        if strcmp(EPdataOut.dataType,'continuous')
            if (size(EPdataOut.data,2) >= EPdataOut.Fs*2) && size(EPdataOut.analysis.moveTrial,2) ==1
                EPdataOut.analysis.moveTrial = []; %backward compatibility conversion
            end
        end
    end
    if isfield(EPdataOut.analysis,'badTrials')
        if strcmp(EPdataOut.dataType,'continuous')
            if (size(EPdataOut.data,2) >= EPdataOut.Fs*2) && size(EPdataOut.analysis.badTrials,2) ==1
                EPdataOut.analysis.badTrials = []; %backward compatibility conversion
            end
        end
    end
    if isfield(EPdataOut.analysis,'badChans')
        if strcmp(EPdataOut.dataType,'continuous')
            if (size(EPdataOut.data,2) >= EPdataOut.Fs*2) && size(EPdataOut.analysis.badChans,2) ==1
                EPdataOut.analysis.badChans = []; %backward compatibility conversion
            end
        end
    end
end

%update montage
if ~isempty(EPdataOut.montage)
    switch EPdataOut.montage
        case 'GSN200-128-21'
            EPdataOut.montage= 'Adult GSN200 128-channel 2.1';
        case 'GSN200-128-1'
            EPdataOut.montage=  'Adult GSN 128-channel 1.0';
        case 'GSN200-128-R1'
            EPdataOut.montage=  'Adult Recumbent GSN 128-channel 1.0';
        case 'GSN200-64-A1'
            EPdataOut.montage=  'Adult GSN 64-channel 1.0';
        case 'GSN200-64-I1'
            EPdataOut.montage=  'Infant GSN 64-channel 1.0';
        case 'GSN200-64-T1'
            EPdataOut.montage=  'Toddler GSN 64-channel 1.0';
        case 'GSN200-64-A2'
            EPdataOut.montage=  'Adult GSN 64-channel 2.0';
        case 'preGSN-64'
            EPdataOut.montage=  'Adult preGSN 64-channel 1.0';
        case 'GSN200-256-2'
            EPdataOut.montage=  'Adult GSN 256-channel 2.0';
        case 'GSN200-256-21'
            EPdataOut.montage=  'Adult GSN 256-channel 2.1';
        case 'Hydrocel-32-1'
            EPdataOut.montage=  'Adult Hydrocel 32-channel 1.0';
        case 'Hydrocel-64-1'
            EPdataOut.montage=  'Adult Hydrocel 64-channel 1.0';
        case 'Hydrocel-128-1'
            EPdataOut.montage=  'Adult Hydrocel 128-channel 1.0';
        case 'Hydrocel-256-1'
            EPdataOut.montage=  'Adult Hydrocel 256-channel 1.0';
        case 'Generic'
            EPdataOut.montage=  'Generic 10-05';
    end
end
if isempty(EPdataOut.montage)
    switch EPdataOut.ced
        case 'GSN129.ced'
            EPdataOut.montage= 'Adult GSN200 128-channel 2.1';
        case 'GSN65v1_0.ced'
            EPdataOut.montage=  'Adult GSN 64-channel 1.0';
        case 'GSN65v2_0.ced'
            EPdataOut.montage=  'Adult GSN 64-channel 2.0';
        case 'preGSN65.ced'
            EPdataOut.montage=  'Adult preGSN 64-channel 1.0';
        case 'GSN257.ced'
            EPdataOut.montage=  'Adult GSN 256-channel 2.1';
        case 'GSN-Hydrocel-65 1.0.ced'
            EPdataOut.montage=  'Adult Hydrocel 64-channel 1.0';
        case 'GSN-Hydrocel-129.ced'
            EPdataOut.montage=  'Adult Hydrocel 128-channel 1.0';
        case 'GSN-Hydrocel-257.ced'
            EPdataOut.montage=  'Adult Hydrocel 256-channel 1.0';
        otherwise
            %make best guess
            Elist=cell2mat(strfind(cellstr(cellfun(@(v)v(1),EPdataOut.chanNames)'),'E'));
            inEchans=zeros(length(Elist),1);
            labelLength=1;
            for iChan=1:length(Elist)
                if (length(EPdataOut.chanNames{iChan})>3) && strcmp(EPdataOut.chanNames{iChan}(1:3),'EEG')
                    labelLength=3;
                end
                inEchans(iChan)=str2double(EPdataOut.chanNames{Elist(iChan)}(labelLength+1:end));
            end
            if ~isempty(inEchans)
                highestChan=max(inEchans,[],'omitnan');
                disp(['There are ' num2str(length(inEchans)) ' channels with names starting with E.']);
                if highestChan > 257
                    disp('no list of comparable electrodes available for Hydrocel montages with this channel count.');
                elseif highestChan > 129
                    EPdataOut.montage='Adult Hydrocel 256-channel 1.0';
                elseif highestChan > 65
                    EPdataOut.montage='Adult Hydrocel 128-channel 1.0';
                elseif highestChan > 33
                    EPdataOut.montage='Adult Hydrocel 64-channel 1.0';
                elseif highestChan > 17
                    EPdataOut.montage='Adult Hydrocel 32-channel 1.0';
                else
                    EPdataOut.montage=EPmain.preferences.general.defaultMontage;
                end
            else
                EPdataOut.montage=EPmain.preferences.general.defaultMontage;
            end
    end
    disp(['Making assumption this is the montage: ' EPdataOut.montage]);
end

if strcmp(EPdataOut.montage,'Generic 10-05')
    changeFlag=1;
    switch EPdataOut.ced
        case 'GSN129.ced'
            EPdataOut.montage= 'Adult GSN200 128-channel 2.1';
        case 'GSN65v1_0.ced'
            EPdataOut.montage=  'Adult GSN 64-channel 1.0';
        case 'GSN65v2_0.ced'
            EPdataOut.montage=  'Adult GSN 64-channel 2.0';
        case 'preGSN65.ced'
            EPdataOut.montage=  'Adult preGSN 64-channel 1.0';
        case 'GSN257.ced'
            EPdataOut.montage=  'Adult GSN 256-channel 2.1';
        case 'GSN-Hydrocel-65 1.0.ced'
            EPdataOut.montage=  'Adult Hydrocel 64-channel 1.0';
        case 'GSN-Hydrocel-129.ced'
            EPdataOut.montage=  'Adult Hydrocel 128-channel 1.0';
        case 'GSN-Hydrocel-257.ced'
            EPdataOut.montage=  'Adult Hydrocel 256-channel 1.0';
        otherwise
            changeFlag=0;
    end
    if changeFlag
        disp(['Based on the ced file, making assumption this should actually be the montage: ' EPdataOut.montage]);
    end
end

%Ensure that one-dimensional string arrays are column vectors
EPdataOut.chanNames=EPdataOut.chanNames(:);
EPdataOut.timeNames=EPdataOut.timeNames(:);
EPdataOut.subNames=EPdataOut.subNames(:);
EPdataOut.cellNames=EPdataOut.cellNames(:);
EPdataOut.trialNames=EPdataOut.trialNames(:);
EPdataOut.facNames=EPdataOut.facNames(:);
EPdataOut.freqNames=EPdataOut.freqNames(:);
EPdataOut.relNames=EPdataOut.relNames(:);
EPdataOut.sessNames=EPdataOut.sessNames(:);
EPdataOut.chanTypes=EPdataOut.chanTypes(:);
EPdataOut.subTypes=EPdataOut.subTypes(:);
EPdataOut.sessNums=EPdataOut.sessNums(:);
EPdataOut.cellTypes=EPdataOut.cellTypes(:);
EPdataOut.facTypes=EPdataOut.facTypes(:);
EPdataOut.recTime=EPdataOut.recTime(:);
EPdataOut.trialSpecNames=EPdataOut.trialSpecNames(:);
EPdataOut.subjectSpecNames=EPdataOut.subjectSpecNames(:);
EPdataOut.taskNames=EPdataOut.taskNames(:);
EPdataOut.taskMeasNames=EPdataOut.taskMeasNames(:);

%ensure fields are in standard order.
[EPfieldNames]=fieldnames(ep_newFile);

modelEPdata=[];
for iField=1:length(EPfieldNames)
    modelEPdata.(EPfieldNames{iField})=[];
end

dataFieldNames=fieldnames(EPdataOut);
if ~isequal(EPfieldNames,dataFieldNames)
    if length(EPfieldNames) < length(dataFieldNames)
        for iField = 1:length(dataFieldNames)
            if ~any(strcmp(dataFieldNames{iField},EPfieldNames))
                EPdataOut=rmfield(EPdataOut,dataFieldNames{iField});
            end
        end
    end
    if length(EPfieldNames) > length(dataFieldNames)
        for iField = 1:length(EPfieldNames)
            if ~any(strcmp(EPfieldNames{iField},dataFieldNames))
                eval(['EPdataOut.' EPfieldNames{iField} '=[];']);
            end
        end
        
    end
    outFields=fieldnames(EPdataOut);
    modelFields=fieldnames(modelEPdata);
    extraFields=outFields(~ismember(outFields,modelFields));
    if ~isempty(extraFields)
        for iField=1:length(extraFields)
            EPdataOut=rmfield(EPdataOut,extraFields{iField});
            disp(['Warning: Dropping excess field ' extraFields{iField} ' from dataset in working set named ' EPdataOut.dataName])
        end
    end
    EPdataOut = orderfields(EPdataOut, modelEPdata);
end

EPdataOut.eloc=ep_elocFormat(EPdataOut.eloc);
EPdataOut.implicit=ep_elocFormat(EPdataOut.implicit);

if ~isempty(EPdataOut.history)
    newHist=cell(0,5);
    for iHist=1:size(EPdataOut.history,1)
        if isstruct(EPdataOut.history{iHist,1})
            newHist(end+1,:)=EPdataOut.history(iHist,1:min(size(EPdataOut.history,2),5));
        end
    end
    EPdataOut.history=newHist;
end
if size(EPdataOut.history,2) == 4
    EPdataOut.history{1,end+1} = cell(0);
end
