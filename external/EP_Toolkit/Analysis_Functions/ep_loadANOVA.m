function [ANOVAheader, ANOVAdata] = ep_loadANOVA(fileName,betweenVars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [ANOVAheader, ANOVAdata] = ep_loadANOVA(fileName,betweenVars) -
% Loads ANOVA text file and organizes the contents.
%
%Inputs
%	fileName:   filename and sourcepath.
%   betweenVars: information about between subject variables.
%
%Outputs
%   ANOVAheader: cell array with cell arrays of the header labels, if any.  Else empty.
%   ANOVAdata:   the data, organized as a 2D cell array.

%History
%  by Joseph Dien (12/1/19)
%  jdien07@mac.com
%
% modified 12/27/19 JD
% Added support for multi-session data.
%
% bugfix 3/23/20 JD
% Fixed crash when Loading an example ANOVA file in ANOVA function that has only one session.
%
% modified 8/31/20 JD
% Added ability to read in text files of data without EP ANOVA header, treating them as behavioral data.
%
% bugfix 12/5/21 JD
% Fixed crash when trying to handle no-header file with between-group variables.
%
% bugfix 5/8/22 JD
% Fixed crash when trying to handle header file with between-group variables.
%
% bugfix 7/5/24 JD
% Fixed crash when loading an ANOVA file with more than one spec.
%
% modified 4/5/25 JD
% Added version 3 ANOVA file format.
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

ANOVAdata=[];
ANOVAheader=cell(0);

if ~exist('betweenVars','var')
    betweenVars=[];
end

%load in the ANOVA data
ANOVAdata.ANOVAver=1;
ANOVAdata.sessNames=cell(0);

[~, theData, theDelim] = ep_textScan(fileName,1);
ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
if isempty(theData)
    msg{1}=['There was an error trying to open ' fileName '.'];
    [msg]=ep_errorMsg(msg);
    ANOVAdata=[];
    return
end

numRows=size(theData,1);
numCols=size(theData,2);
numHeaders=7;
dataOnlyFlag=0;
if (numRows  < (numHeaders+1)) || ~isempty(theData{6,end}) || ~isempty(theData{7,end}) %just data, no EP header
    dataOnlyFlag=1;
    disp('Assuming data is behavioral and adding the appropriate EP ANOVA header.');
    theHeader=cell(numHeaders,1);
    theHeader{1}=cell(0);
    theHeader{2}{1}='behavioral';
    theHeader{3}=cell(0);
    theHeader{4}=cell(0);
    theHeader{5}=cell(0);
    theHeader{6}=cell(0);
    theHeader{8}=cell(0);
    for iCol=1:numCols
        theHeader{9}{iCol}=char(double('a')-1+iCol);
        theHeader{10}{iCol}='data';
    end
    newHeader=theHeader;
else
    [theHeader, theData, theDelim] = ep_textScan(fileName,numHeaders+1);
    ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
end

if (~isempty(theHeader{1}) && startsWith(theHeader{1,1},'EP Toolkit Version:')) || startsWith(theHeader{2}{1},'behavioral')
    ANOVAdata.ANOVAver=3;
    numHeaders=10;
    [theHeader, theData, theDelim] = ep_textScan(fileName,numHeaders+1);
    ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
    if isempty(theData)
        msg{1}=['There was an error trying to open ' fileName '.'];
        [msg]=ep_errorMsg(msg);
        return
    end    
elseif isnan(str2double(theData{1,1}))
    ANOVAdata.ANOVAver=2;
    numHeaders=8;
    [theHeader, theData, theDelim] = ep_textScan(fileName,numHeaders+1);
    ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
    if isempty(theData)
        msg{1}=['There was an error trying to open ' fileName '.'];
        [msg]=ep_errorMsg(msg);
        return
    end
else
    disp('This appears to be a Version 1 ANOVA file.');
end
ANOVAheader=cell(numHeaders,1);
numRows=size(theData,1);
numCols=size(theData,2);
if dataOnlyFlag
    theHeader=newHeader;
end

for iHeader=1:numHeaders
    fusedHeader='';
    for iField=1:length(theHeader{iHeader})
        if iField==1
            fusedHeader=deblank(theHeader{iHeader}{iField});
        else
            fusedHeader=deblank([fusedHeader theDelim theHeader{iHeader}{iField}]);
        end
    end
    ANOVAheader{iHeader}=fusedHeader;
    switch ANOVAdata.ANOVAver
        case 1
            switch iHeader
                case 1
                    ANOVAdata.name=fusedHeader;
                case 2
                    ANOVAdata.window=fusedHeader;
                case 3
                    ANOVAdata.measure=fusedHeader;
                case 4
                    ANOVAdata.changrp=fusedHeader;
                case 5
                    ANOVAdata.factor=fusedHeader;
                case 6
                    ANOVAdata.cellNames=theHeader{iHeader};
                case 7
                    ANOVAdata.areaNames=theHeader{iHeader};
                case 8
                    ANOVAdata.areaNames=theHeader{iHeader};
            end
        case 2
            switch iHeader
                case 1
                    ANOVAdata.name=fusedHeader;
                case 2
                    ANOVAdata.window=fusedHeader;
                case 3
                    ANOVAdata.measure=fusedHeader;
                case 4
                    ANOVAdata.changrp=fusedHeader;
                case 5
                    ANOVAdata.factor=fusedHeader;
                case 6
                    ANOVAdata.sessNames=theHeader{iHeader};
                case 7
                    ANOVAdata.cellNames=theHeader{iHeader};
                case 8
                    ANOVAdata.areaNames=theHeader{iHeader};
            end
        case 3
            switch iHeader
                case 1
                    ANOVAdata.EPver=fusedHeader;
                case 2
                    ANOVAdata.name=fusedHeader;
                case 3
                    ANOVAdata.window=fusedHeader;
                case 4
                    ANOVAdata.measure=fusedHeader;
                case 5
                    ANOVAdata.changrp=fusedHeader;
                case 6
                    ANOVAdata.factor=fusedHeader;
                case 7
                    ANOVAdata.inputCells=theHeader{iHeader};
                case 8
                    ANOVAdata.sessNames=theHeader{iHeader};
                case 9
                    ANOVAdata.cellNames=theHeader{iHeader};
                case 10
                    ANOVAdata.areaNames=theHeader{iHeader};
            end
        otherwise
            disp('Error: ANOVA file version not recognized.')
            return
    end

end

if isfield(ANOVAdata,'inputCells')
    if ~isempty(ANOVAdata.inputCells)
        inputCells=cell(0);
        for iCell=1:length(ANOVAdata.inputCells)
            inputCells{iCell,1}=split(ANOVAdata.inputCells{iCell}(2:end),'@');
        end
        ANOVAdata.inputCells=inputCells;
    end
else
    ANOVAdata.inputCells=cell(0);
end

if (isscalar(ANOVAdata.sessNames)) && isempty(ANOVAdata.sessNames{1})
    ANOVAdata.sessNames=cell(0);
end

if (length(ANOVAdata.cellNames)+1) == numCols %drop subject names column
    theData=theData(:,1:end-1);
    numCols=numCols-1;
end

specList=find(strcmp('spec',ANOVAdata.areaNames))';
if ~isempty(betweenVars)
    specList=[specList; betweenVars.between(find(~cellfun(@isempty,betweenVars.betweenName)))'];
end
numSpecs=length(unique(specList));
numDataCols=numCols-numSpecs;

if length(ANOVAdata.cellNames) ~= numCols
    msg{1}=['The number of cell names (' num2str(length(ANOVAdata.cellNames)) ') is different from the number of data columns (' num2str(numCols) ').'];
    [msg]=ep_errorMsg(msg);
    ANOVAdata=[];
    return
end
if length(ANOVAdata.areaNames) ~= numCols
    msg{1}=['The number of area names (' num2str(length(ANOVAdata.areaNames)) ') is different from the number of data columns (' num2str(numCols) ').'];
    [msg]=ep_errorMsg(msg);
    ANOVAdata=[];
    return
end
if (length(ANOVAdata.sessNames) ~= numDataCols) && ~isempty(ANOVAdata.sessNames)
    msg{1}=['The number of session names (' num2str(length(ANOVAdata.sessNames)) ') is different from the number of data columns (' num2str(numDataCols) ').'];
    [msg]=ep_errorMsg(msg);
    ANOVAdata=[];
    return
end

if ~isempty(betweenVars)
    if (betweenVars.totalLevels+size(betweenVars.betweenLvl,2))>length(ANOVAdata.cellNames)
        msg{1}=['The number of data columns in the file (' num2str(length(ANOVAdata.cellNames)) ') is smaller than the number of columns specified by the ANOVA table (' num2str(betweenVars.totalLevels+size(betweenVars.betweenLvl,2)) ').'];
        [msg]=ep_errorMsg(msg);
        ANOVAdata=[];
        return
    end
end

ANOVAdata.betweenLvl=[];
for iRow=1:numRows
    for iCol = 1:numDataCols
        ANOVAdata.data(iRow,iCol)=str2double(theData(iRow,iCol));
    end
    for iSpec = 1:numSpecs
        ANOVAdata.betweenLvl{iRow,iSpec} = theData{iRow,iSpec+numDataCols};
    end
end

% if ~all(cellfun(@isempty,betweenVars.betweenName)) && isempty(ANOVAdata.betweenLvl)
%     msg{1}=['Between group factors have been specified but none detected in the data.  Did you remember to use the ''spec'' label for the between group variables?'];
%     [msg]=ep_errorMsg(msg);
%     ANOVAdata=[];
%     return    
% end
%check each file to make sure the between group levels are fully crossed
betweenFactors=[];
totBetweenLevels=1;
badSubs=zeros(size(ANOVAdata.data,1),1);
for iFactor=1:6
    if ~isempty(betweenVars)
        if ~isempty(betweenVars.betweenName{iFactor})
            if betweenVars.between(iFactor) < numDataCols
                msg{1}=['Error: there are more data columns in the file than expected in ' fileName '.'];
                [msg]=ep_errorMsg(msg);
                ANOVAdata=[];
                return
            end
            theLevels={ANOVAdata.betweenLvl{:,betweenVars.between(iFactor)-numDataCols}};
            for iSub=1:length(theLevels)
                if isempty(theLevels{iSub})
                    badSubs(iSub)=1;
                end
            end
        end
    end
end

badSubList=find(badSubs);
goodSubList=find(~badSubs);
if ~isempty(badSubList)
    disp(['Warning: Dropping ' num2str(length(badSubList)) ' subject(s) due to missing between group information.']);
end
if ~isempty(betweenVars) && ~isempty(ANOVAdata.betweenLvl)
    ANOVAdata.betweenLvl=ANOVAdata.betweenLvl(goodSubList,:);
end
ANOVAdata.data=ANOVAdata.data(goodSubList,:);
ANOVAdata.origData=theData(goodSubList,:);

for iFactor=1:6
    if ~isempty(betweenVars) && ~isempty(betweenVars.betweenName{iFactor}) && (betweenVars.covariate(iFactor)~=2)
        theLevels={ANOVAdata.betweenLvl{:,betweenVars.between(iFactor)-numDataCols}};
        for iSub=1:length(theLevels)
            theLevelsFirstLetter{iSub}=theLevels{iSub}(1);
        end
        
        if (length(unique(theLevels)) ~= length(unique(theLevelsFirstLetter)))
            msg{1}=['The between group variable ' betweenVars.betweenName{iFactor} ' in ' fileName ' has different level labels with the same first letter.'];
            [msg]=ep_errorMsg(msg);
            ANOVAdata=[];
            return;
        end
        
        if isscalar(unique({ANOVAdata.betweenLvl{:,betweenVars.between(iFactor)-numDataCols}}))
            msg{1}=['The between group variable ' betweenVars.betweenName{iFactor} ' in ' fileName ' has only a single level.'];
            [msg]=ep_errorMsg(msg);
            ANOVAdata=[];
            return;
        end
        
        betweenFactors(end+1)=betweenVars.between(iFactor);
        totBetweenLevels=totBetweenLevels*length(unique(theLevels));
    end
end

if totBetweenLevels> 1
    %only use first letter of between levels
    for i=1:size(ANOVAdata.betweenLvl,1)
        for j=1:size(ANOVAdata.betweenLvl,2)
            if betweenVars.covariate(j)==2
                ANOVAdata.betweenLvl{i,j}='';
            else
                ANOVAdata.betweenLvl{i,j}=ANOVAdata.betweenLvl{i,j}(1);
            end
        end
    end
    
    [sortedBetween,index] = sortrows(ANOVAdata.betweenLvl(:,betweenFactors-numDataCols));
    ANOVAdata.data=ANOVAdata.data(index,:);
    ANOVAdata.index=index;
    betweenCombos=cellstr(cell2mat(sortedBetween));
    if length(unique(betweenCombos)) ~= totBetweenLevels
        msg{1}=['The between group factors are incompletely crossed.  Nested designs are not currently supported.'];
        [msg]=ep_errorMsg(msg);
        ANOVAdata=[];
        return;
    end
    
    uniqueCombos=unique(betweenCombos);
    ANOVAdata.subjects=zeros(totBetweenLevels,1);
    for i=1:totBetweenLevels
        ANOVAdata.subjects(i)=length(find(strcmp(uniqueCombos{i},betweenCombos)));
    end
else
    ANOVAdata.subjects=size(ANOVAdata.data,1);
    ANOVAdata.index=[1:ANOVAdata.subjects];
    betweenCombos=cellstr(repmat('gave',ANOVAdata.subjects,1));
end

if any(ANOVAdata.subjects < 2)
    msg{1}=['Each between group cell needs at least two subjects to be able to estimate error variance.'];
    [msg]=ep_errorMsg(msg);
    ANOVAdata=[];
    return;
end



