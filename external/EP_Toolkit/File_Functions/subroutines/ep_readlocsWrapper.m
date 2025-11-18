function eloc = ep_readlocsWrapper(fileName,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% eloc = ep_readlocsWrapper(fileName)
%	Provides front-end for EEGlab's readlocs to make sure the ced fields are in the correct order, otherwise it scrambles the data.
%   This function only keeps the ten canonical eloc fields and drops
%   others, including the two besa fields that sometimes show up.
%
%Inputs
%	fileName:     filename including path.
%
%Outputs
%    eloc      : The electrode location information, one for each channel (see readlocs header)
%
% History:
%
% by Joseph Dien (10/17/18)
% jdien07@mac.com
%
% bugfix 11/4/18 JD
% Fixed function not actually making any difference for compensating for EEGlab bug.
%
% bugfix 11/4/18 JD
% Fixed readlocs taking contents of next column when type column is left empty by defaulting to EEG.
%
% bugfix 6/13/19 JD
% Fixed case where a previous error leaves behind a copy of the temporary CED file, resulting in errors.
% Fixed case where user does not have write permission for Matlab active directory, resulting in failure to read the CED file.
%
% modified 12/30/19 JD
% Enabled reading of text files with a greater range of variations in character encoding, end-of-line markers, and field separation markers.
%
% bugfix 2/7/20 JD
% Fixed defaulting to "EEG" for all empty fields rather than just for the type field.
% Fixed sometimes reading CED fields as unicode rather than numbers.
%
% modified 4/14/20 JD
% Removes besa fields if added by eeglab.
%
% bugfix 6/28/20 JD
% Fixed crash when a field was left blank in a .ced file.
% Fixed error when two eloc labels are the same.
%
% modified 8/8/21 JD
% Added support for canonical 10-05 coordinates in eloc structure.
%
% bugfix & modified 6/20/22 JD
% Added support for reading sfp files.
% Fixed check for non-unique label names.
%
% modified 8/1/22 JD
% Names that readlocs reserves for fiducials ('nz' 'lpa' 'rpa' 'nasion' 'left' 'right' 'nazion' 'fidnz' 'fidt9' 'fidt10' 'cms' 'drl' 'nas' 'lht' 'rht' 'lhj' 'rhj') may now be used as channel names.
%
% bugfix 10/10/22 JD
% Fixed type information going into theta field when coordinates information is empty and then subsequently being ignored and the type field defaulting to "EEG".
%
% bugfix 1/1/25 JD
% Fixed crash when reading sfp file due to no type field.
% Added types for sfp info, assumes first three rows are FID, and removed
% "Fid" from start of their names if present.
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

eloc=[];

global EPdataset EPtictoc;

startDataRow=1;
fileType=find(strcmp('filetype',varargin));
if ~isempty(fileType)
    if strcmp('chanedit',varargin{fileType+1})
        startDataRow=2;
        fileType='ced';
    elseif strcmp('sfp',varargin{fileType+1})
        startDataRow=1;
        fileType='sfp';
    else
        disp('Error - at present this function is only meant to work with .ced or .sfp files.')
        return
    end
end
[theHeader, theData, theDelim] = ep_textScan(fileName,startDataRow);
ep_tictoc;if EPtictoc.stop;return;end
if isempty(theData)
    disp('Problem with opening the CED file.  Perhaps you need to put it where it can be found, as in the current directory or the electrodes folder in EP Toolkit.')
    return
end

tempCED=[EPdataset.EPwork filesep 'EPwork' filesep ['ep_tempCEDfile.' fileType]];
if exist(tempCED,'file')
    delete(tempCED);
end
fidOut=fopen(tempCED,'w');
if fidOut==-1
    disp('Problem with writing temporary coordinates file.')
    return
end

if strcmp(fileType,'ced')
    listcolformat=ep_elocFormat;
    outData=cell(1,length(listcolformat)+1);
    cedFields=theHeader{1}(2:end); %skip number field
    indexArray=zeros(length(cedFields),1);
    otherCounter=1;
    for iField=1:length(cedFields)
        if any(strcmp(cedFields{iField},listcolformat))
            %CED field is part of the EP standard fields
            indexArray(iField)=find(strcmp(cedFields{iField},listcolformat));
        else
            %CED field is not part of the EP standard fields
            indexArray(iField)=length(listcolformat)+otherCounter;
            otherCounter=otherCounter+1;
        end

    end

    outData(1:size(theData,1),1)=theData(:,1);
    for iField=1:length(indexArray)
        outData(:,indexArray(iField)+1)=theData(:,iField+1);
    end

    for iField=1:length(theHeader{1})
        fprintf(fidOut, '%s\t',  theHeader{1}{iField}); %field labels
    end
    fprintf(fidOut, '\n');
    for iRow=1:size(theData,1)
        fprintf(fidOut, '%d\t',  str2double(outData{iRow,1}));
        for iField = 2:length(listcolformat)+1
            theVal=outData{iRow,iField};
            if isempty(deblank(theVal))
                if strcmp(listcolformat{iField-1},'labels')
                    theVal=['e' num2str(iRow)];
                elseif strcmp(listcolformat{iField-1},'type')
                    theVal='EEG';
                else
                    theVal=' ';
                end
            else
                if isnan(str2double(theVal)) && any(strcmp(listcolformat{iField-1},{'theta';'radius';'X';'Y';'Z';'sph_theta';'sph_phi';'sph_radius';'cX';'cY';'cZ'}))
                    disp(['Error in the CED file ' fileName ' - the value ' theVal ' cannot be present in the field ' listcolformat{iField-1} '.'])
                    return
                end
            end
            if strcmp(listcolformat{iField-1},'labels') && (length(find(strcmp(theVal,outData(:,2))))>1) %readlocs requires labels to be unique
                labelName=theVal;
                sameName=1;
                theNumber=0;
                while sameName
                    sameName=0;
                    if ~isempty(find(strcmp(theVal,outData(:,2))))
                        sameName=1;
                    end
                    if sameName
                        theNumber=theNumber+1;
                        labelName=[theVal '-' num2str(theNumber)];
                    end
                end
                theVal=labelName;
            end
            if isnumeric(theVal)
                fprintf(fidOut, '%3.3g\t',  theVal);
            else
                fprintf(fidOut, '%s\t',  theVal);
            end
        end
        fprintf(fidOut, '\n');
    end
elseif strcmp(fileType,'sfp')
    outData=theData;
    for iRow=1:size(theData,1)
        for iField = 1:size(theData,2)
            theVal=outData{iRow,iField};
            if isempty(theVal)
                if iField==1
                    theVal=['e' num2str(iRow)];
                else
                    theVal=' ';
                end
            end
            if (iField==1) && (length(find(strcmp(theVal,outData(:,1))))>1) %readlocs requires labels to be unique
                labelName=theVal;
                sameName=1;
                theNumber=0;
                while sameName
                    sameName=0;
                    if ~isempty(find(strcmp(theVal,outData(:,1))))
                        sameName=1;
                    end
                    if sameName
                        theNumber=theNumber+1;
                        labelName=[theVal '-' num2str(theNumber)];
                    end
                end
                theVal=labelName;
            end
            if isnumeric(theVal)
                fprintf(fidOut, '%3.3g\t',  theVal);
            else
                fprintf(fidOut, '%s\t',  theVal);
            end
        end
        fprintf(fidOut, '\n');
    end
end

fclose(fidOut);

theArgs='tempCED';
for iVar=1:length(varargin)
    theArgs=[theArgs ',''' varargin{iVar} ''''];
end
evalc(['eloc = readlocs(' theArgs ')']);
delete(tempCED);

%if readlocs has changed any of the types to FID based on the channel name, change it back.
if isfield(eloc,'type')
    for iRow=1:length(eloc)
        if isnan(str2double(eloc(iRow).theta)) && isempty(eloc(iRow).type)
            eloc(iRow).type=eloc(iRow).theta; %if the coordinate fields are empty, readlocs will put the type info in the first empty field (i.e., theta).
        end
        if strcmp(eloc(iRow).type,'FID') && ~strcmp(theData(iRow,find(strcmp('type',theHeader{1}))),'FID')
            eloc(iRow).type=theData{iRow,find(strcmp('type',theHeader{1}))};
        end
    end
end

eloc=ep_elocFormat(eloc);

if strcmp(fileType,'sfp')
    %sfp files typically start with three fiducial locations and the rest
    %will be assumed to be EEG channels.
    for iRow=1:length(eloc)
        if iRow<4
            eloc(iRow).type='FID';
            if strcmp(eloc(iRow).labels(1:3),'Fid')
                eloc(iRow).labels=eloc(iRow).labels(4:end);
            end
        else
            eloc(iRow).type='EEG';
        end
    end
end