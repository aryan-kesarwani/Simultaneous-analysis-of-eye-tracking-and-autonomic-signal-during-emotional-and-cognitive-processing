function [EPdataOutFinal]=ep_combineData(EPdataIn,dataDimension,combineData,combineWeightsIn,combineName,trialWeights,negCell,posCell)
%  [EPdataOutFinal]=ep_combineData(EPdataIn,dataDimension,combineData,combineWeightsIn,combineName,trialWeights,negCell)
%       Combines data levels (like a subset of subjects) and adds the result to the data.  Additionally, one can specify that this operation will only be applied to certain levels of the other dimensions.
%       There is also a negCell option for subtracting a single cell of the negative subs from the positive subs.
%
%Inputs:
%  EPdataIn       : Structured array with the input data and accompanying information in EP file format.  See readData.
%  dataDimension  : Which data dimension to combine (i.e., 'channels', 'cells', 'subjects', 'factors')
%                   'convert' will change all the virtual cells and subjects to real ones.
%  combineData    : Which entries to use from each dimension (cell array with six cells, one per dimension).  Empty means keep everything.
%                   For the dataDimension, Which levels of the data dimension to combine (e.g., [1 3 4]).
%                   For the other dimensions, which levels are to be updated with this new combination.
%                   There's no separate entry for the 7th relationships dimension because it just mirrors the 1st channels dimension.
%                   If converting virtual grand average to real grand average, ignore this parameter if info in GAVsubs.
%  combineWeightsIn : Weights to use when combining the levels of the data (e.g., [.5 .25 .25]).
%                   If totals more than one, weights will be rescaled so that they do total to one.
%                   If weights total zero (as in a difference wave) then no rescaling.
%                   For factors, no rescaling is applied so simple addition rather than averages are computed.
%                   If empty, then equal weighting, unless converting virtual weights.
%                   If converting virtual grand average to real grand average, ignore this parameter if info in GAVsubs.
%  combineName     : Name of the new combined level (e.g., 'grand average').
%                    If already present, update the data (and convert to real if virtual cell/subject) rather than add a new level.
%  trialWeights    : Weight by number of trials in average if available, for cells (0=no, 1=yes)
%  negCell         : For subject difference waves, use this cell of the negative subjects to subtract from all of the positive subject cells.
%  posCell         : For subject difference waves, use this cell of the positive subjects to add to all of the negative subject cells.
%
%Outputs:
%  EPdataOutFinal      : Structured array with the output data and accompanying information in EP file format.  See readData.
%
% With respect to trimmed means: Things get a little complicated with trimmed means of factors since the set of trimmed subjects can differ depending on the time window for spatial factors and the channel for temporal factors.
% In order to keep things logically consistent, when adding a trimmed subject mean waveform, it will only be done for one factor.  The rest will be left as zero and filled in when that factor is later done, as in autoPCA.
% Ditto for the accompanying noise, STD, etc waveforms.
% When combining values, NaN values are ignored.  If there ends up being no data for a datapoint, it will end up being an NaN (for contrasts, it will end up being NaN if either the negative or positive terms are missing).
% Ditto for empty values in trial specs.
% For contrasts, the positive and the negative sides are given equal weighting.  Within each side, the values are weighted according to the specified weights.
% Trial specs (such as RT) are handled the same as the voltage data.  Meta-data like badTrials are simply added without respect to the signs but are weighted.
% AveNum is simply the addition of the involved cells without regard to weighting or sign.  SubNum as well for 'subjects' but just the max value for 'cells.'
%
% Averaging of coherence matrices (which are basically correlations) is performed via Fisher-Z transform.
% Corey, D. M., Dunlap ,William P., & and Burke, M. J. (1998). Averaging Correlations: Expected Values and Bias in Combined Pearson rs and Fisher’s z Transformations. The Journal of General Psychology, 125(3), 245–261. https://doi.org/10.1080/00221309809595548
% https://medium.com/@jan.seifert/averaging-correlations-part-i-3adab6995042

%History:
%  by Joseph Dien (5/24/12)
%  jdien07@mac.com
%
%  modified 7/17/12 JD
%  Added option to weight cell combinations by number of trials in averages.
%
%  bugfix 8/6/12 JD
%  Fixed edit's add cells trial weighting option not working correctly when the cells are not a consecutive series starting with the first.
%
%  bugfix 11/4/12 JD
%  Fixed crash when using combining cells or chans for spatial PCA data.
%
%  bugfix 1/10/13 JD
%  Fixed crash when combining channels under certain circumstances.
%
%  bugfix 1/30/13 JD
%  Fixed crash when combining channels for factor data under certain circumstances (presence of facData due to adds).
%
%  bugfix 3/25/13 JD
%  Fixed combining of subjects and chans not correct when weights not the same (as in difference wave).
%  Fixed weighting of difference waves for cells and chans and subjects incorrect (waves too small).
%
%  bugfix 7/14/13 JD
%  Fixed combining channels results in flat waveform.
%
% modified 10/9/13 JD
% Added recTime field.
%
%  bugfix 11/25/13 JD
%  Fixed noise and std fields set equal to the data when combining cells or subjects.
%
%  bugfix 11/28/13 JD
%  If new add fails data check then output empty matrix.
%
% modified 3/19/14 JD
% Added combining factors
%
% bugfix 3/20/14 JD
% Makes sure that additions to the name fields are added as a column vector.
% Fixed bad subtype when adding a single subject.
%
%  bugfix 3/22/14 JD
%  Fixed all but one channel is flat for grand average combined factors, as in the "all" factor from PCAs.
%  Fixed all but one channel is flat for combined cells if one already has a combined factor, as when one uses the Edit
%  function to combined cells on a factor cell containing an "all" cell.
%
% modified 3/24/14 JD
% Added .cov field.
%
% bugfix 4/2/14 JD
% Fixed crash when performing combination of subjects with file containing .cov information.
%
% modified 4/14/14 JD
% Added .covNum field.
%
% bugfix 4/2/14 JD
% Fixed weighting not correct when computing difference waves that do not sum to zero.
% Fixed calculation of the SubNum field (number of subjects going into averages) when combining cells.
%
% modified 4/24/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure.
%
% bugfix 6/1/14 JD
% Fixed crash when combining trials.
%
% modified 6/1/14 JD
% Added simple averaging of trials in cell subpane for single-trial data.
%
% bugfix 6/29/14 JD
% Fixed cov.Nq field not being formed correctly when combining subjects, resulting in crashes later on.
%
% bugfix 4/27/15 JD
% Fixed error when adding cells to single-trial data.
%
% modified 9/4/15 JD
% Added trial specs for average files.
%
% bugfix 10/24/15 JD
% Fixed crash when combining cells or subjects and there are trial specs with numbers.
%
% bugfix 6/18/17 JD
% Fixed not combining numeric trial specs correctly over subjects, resulting in crashes down the line in ANOVA function.
% Switch to amplitude scaling when adding freq data together other than channels.
%
% bugfix 7/2/17 JD
% Fixed crash when combining subjects for spatial PCA data.
%
% modified 11/22/17 JD
% Added support for impedances field.
%
% bugfix 11/25/17 JD
% Fixed crash when conducting PCA on data where there are trial spec fields that are identical across all the trials/cells and are characters.
%
% bugfix 12/9/17 JD
% Fixed crash when combining cells and the data is frequency-domain.
%
% modified 2/11/18 JD
% Changed std field of combined subjects (as via the Edit function) to be std of the newly generated grand average data rather than a combination of their std values.
% Changed noise field of combined subjects (as via the Edit function) to be noise of the newly generated grand average data rather than a combination of their noise values.
%
% modified 2/11/18 JD
% Added support for stdCM field.
%
% bugfix 4/8/18 JD
% Fixed bug when combining subjects with trial specs where they are all the same string value across the subjects, resulting in failure to combine.
%
% bugfix 6/24/18 JD
% Fixed bug in calculation of .covNum field for FIFF files.
% Corrected calculation of covNum field when using trial weighting option with combination of cell averages.
%
% bugfix 8/5/18 JD
% Fixed crash when using Edit function to generate a grand average of frequency-domain data.
%
% modified 8/7/18 JD
% Optimized stdCM field calculation so it wouldn't take so darn long.
% Changed Cousineau-Morey standard error computation to be specific to two-cell contrasts.
%
% modified & bugfix 1/14/19 JD
% Added factor specification.
% Now allows CMB and SGL factors to be intermixed.
% FacVar and FacVarQ now include CMB factors.
% Fixed not adding cell in case where only one cell is to be added.
%
% modified 4/15/19 JD
% Added support for task level performance measures.
%
% bugfix 4/15/19 JD
% Fixed crash for frequency data.
%
% modified 6/14/19 JD
% No longer maintains noise and std information when combining subjects in factor data as not really needed and requires substantial memory space.
%
% bugfix 8/13/19 JD
% Fixed not combining factors correctly.
%
% bugfix & modified 11/4/19 JD
% Added sessNums sessNames fields.
% Fixed covData being computed based only on subjects that were good in the last cell when combining subjects.
%
% modified 1/12/20 JD
% Upgraded support of std information by adding .covAVE and .GAVsubs fields and eliminating .std and .stdCM fields.
% Revised function call so that levels of other dimensions can be selectively updated other than the dataDimension.
% Also, if the combineName is already present, it will update the data rather than add a new level.
% Better support for adding regional channels to factor data.
% Fixed summed factors too large (not divided by number of factors going into waveform).
%
% bugfix 1/20/20 JD
% Fixed crash when adding a combined factor to a dataset with no combined factors yet.
%
% bugfix 3/22/20 JD
% Now correctly handles NaN values.
% Now handles .badChans field more correctly.
% Fixed crash bugs when combining channels.
% Fixed not weighting levels correctly when no weights were provided to the function.
%
% bugfix 4/12/20 JD
% Fixed negative weights being treated as positive.
%
% bugfix 5/6/20 JD
% Fixed crash when combining subjects into a grand average and there are task specs.
%
% bugfix 10/21/20 JD
% if new combined cell was meant to be a difference wave and enough cells are bad that the weights no longer balance to zero, then treat as invalid rather than just crash.
%
% bugfix 1/23/21 JD
% Fixed crash if the trial-wise variance information was retained during averaging and there is an attempt to generate a cell difference wave where cells are bad for a subject.
%
% modified 8/11/21 JD
% Added negCell and posCell options.
% Revised how trial-specs and meta-data are combined to align better with how the voltage data are being combined.
%
% modified 11/21/21 JD
% ACC and RT trial specs now produce average of the specs rather than the difference so that the RT option in the View function works properly.
%
% bugfix 12/12/21 JD
% Fixed crash when trial specs include empty cells.
% When trial specs are all empty values, then result will also be empty value.
%
% bugfix 5/19/22 JD
% Fixed crash when combining subjects and the optional standard deviation information is present.
% Fixed crash when combining cells and none are good and there are trial specs present.
%
% bugfix 6/27/22 JD
% Fixed crash when combining subjects and there are empty trial spec values.
%
% bugfix 8/19/22 JD
% Fixed crash when combining subjects and the subject list is a row vector by orienting all lists as column vectors.
%
% bugfix 11/18/23 JD
% Fixed crash when combining cells containing covariance data and all but one of the cells are empty.
%
% bugfix 5/30/24 JD
% Fixed unable to combine cells or subjects, as in Edit function, when all the values of a trial spec were missing data.
%
% bugfix 11/21/24 JD
% Fixed RT and ACC trial specs being set to NaN for cells and subjects combined via the Edit function.
%
% bugfix & modified 4/5/25 JD
% Added video field.
% Now handles subject specs in same way as trial specs, plus also generating summary string for gender/sex.
% Added support for virtual grand averages.
% Fixed regional channels computation not handling NaN data points correctly, sometimes resulting in blank waveforms.
% Fixed grand averaging sometimes not correctly computing averaged subject specs correctly.
% Fixed grand averaging sometimes not computing averaged task specs.
% Fixed crash in combining cells or subjects when non-consecutive.
% Fixed crash when combining cells or subjects and some are bad.
% Fixed crash when performing windowAdds with spatial PCA data.
% Fixed WindowAdds addition of channels to PCA files.
% Fixed crash when combining subjects in certain cases.
% Allows for combining of coherence data, using Fisher's Z transform.
% Fixed coherence data not converted to absolute form when combining, unlike normal FFT data.
%
% bugfix 4/30/25 JD
% Fixed crash when subtracting two cells and they have trial specs.
% Fixed output of aborted procedure not empty set, resulting in unpredictable behavior.
%
% bugfix 5/9/25 JD
% Fixed output being empty for 'convert' option when GAVsubs field is empty.
% Fixed crash when combining channels.
%
% bugfix 5/23/25 JD
% Fixed crash when a spatial PCA file has impedance values.
%
% bugfix 6/11/25 JD
% Fixed not converting virtual cell/subject to real cell/subject when in cells/subjects mode and a virtual cell/subject name is specified.
% Fixed crash when combining cells, as in Edit or SampleTest.
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

EPdataOutFinal=[];
EPdataOut=EPdataIn;

if ~isempty(EPdataIn.relNames)
    if strcmp('channels',dataDimension)
        disp('Error: EP Toolkit cannot combine channels of coherence data.  Have not worked out the math yet.')
        EPdataOut=[];
        return
    end
    %Fisher-Z transform
    EPdataOut.data=atanh(EPdataOut.data);
end

EPdataOut.video=struct('frames',[],'times',[]); %regardless of nature of combination, video would be set to empty.
EPdataOut.video(1)=[];

if ~exist('trialWeights','var') || isempty(trialWeights)
    trialWeights=0;
end

if ~exist('negCell','var')
    negCell=[];
end

if ~exist('posCell','var')
    posCell=[];
end

numChans=length(EPdataIn.chanNames);
facNumChans=numChans; %number of channels dimension for .facData if any
if ~isempty(EPdataIn.facVecS)
    numChans=1; %number of channels dimension for .data
end
numPoints=length(EPdataIn.timeNames);
facNumPoints=numPoints;
if ~isempty(EPdataIn.facVecT)
    numPoints=1;
end
if (numPoints==0 && ~isempty(EPdataIn.freqNames))
    numPoints=1;
    facNumPoints=1;
end
numSubs=length(EPdataIn.subNames);
numVsubs=max(0,size(EPdataIn.GAVsubs,1)-1);
numRsubs=numSubs-numVsubs;
numCells=length(EPdataIn.cellNames);
numVcells=max(0,size(EPdataIn.GAVsubs,2)-1);
numRcells=numCells-numVcells;
numFacs=length(EPdataIn.facNames);
numFreqs=length(EPdataIn.freqNames);
facNumFreqs=numFreqs;
if ~isempty(EPdataIn.facVecF)
    numFreqs=1;
end
if numFreqs==0
    numFreqs=1;
end
if facNumFreqs==0
    facNumFreqs=1;
end
numRels=length(EPdataIn.relNames);
if numRels==0
    numRels=1;
end
if ~isempty(EPdataIn.facData)
    numCMBfacs=size(EPdataIn.facData,5);
else
    numCMBfacs=0;
end
numSGLfacs=numFacs-numCMBfacs;
SGLfacs=find(strcmp(EPdataIn.facTypes,'SGL'));
CMBfacs=find(strcmp(EPdataIn.facTypes,'CMB'));

chanList=combineData{1};
if isempty(chanList)
    chanList=[1:numChans]'; %number of channels dimension to keep for .data
    facChanList=[1:facNumChans]'; %number of channels dimension to keep for .facData
elseif ~isempty(EPdataIn.facVecS)
    facChanList=chanList;
    chanList=1;
else
    facChanList=chanList;
end
chanList=chanList(:);
facChanList=facChanList(:);

pointList=combineData{2};
if isempty(pointList)
    pointList=[1:numPoints];
    facPointList=[1:facNumPoints];
elseif ~isempty(EPdataIn.facVecT)
    facPointList=pointList;
    pointList=size(EPdataIn.data,2);
else
    facPointList=pointList;
end
pointList=pointList(:);
facPointList=facPointList(:);

cellList=combineData{3};
if isempty(cellList)
    cellList=[1:numRcells]';
end
cellList=cellList(:);

subList=combineData{4};
if isempty(subList)
    subList=[1:numRsubs]';
end
subList=subList(:);

if ~isempty(EPdataOut.GAVsubs) && any(subList>numRsubs)
    %virtual grand averages not allowed to reference other virtual grand averages
    %replace with subs it referenced, using full list rather than trimmed lists.
    subList2=subList;
    for iSub=1:length(subList)
        theSub=subList(iSub);
        if subList(iSub) > numRsubs
            subList2=setdiff(subList,theSub);
            subList2=[subList2; EPdataOut.GAVsubs{theSub-numRsubs+1,1,1}(:,1)]; %assume fac entries all the same
        end
    end
    subList=unique(subList2);
end

facList=combineData{5};
if isempty(facList)
    facList=[1:numFacs]';
end
facList=facList(:);

freqList=combineData{6};
if isempty(freqList)
    freqList=[1:numFreqs];
    facFreqList=[1:facNumFreqs];
elseif ~isempty(EPdataIn.facVecF)
    facFreqList=freqList;
    freqList=size(EPdataIn.data,6);
else
    facFreqList=freqList;
end
freqList=freqList(:);
facFreqList=facFreqList(:);

if numFacs==0
    numFacs=1;
    facSGLlist=1;
    facCMBlist=[];
else
    facSGLlist=find(ismember(SGLfacs,facList));
    facCMBlist=find(ismember(CMBfacs,facList));
end
facList=facList(:);

relList=[1:numRels]';
relList=relList(:);

numChanList=length(chanList);
numPointList=length(pointList);
numCellList=length(cellList);
numSubList=length(subList);
numFacSGLlist=length(facSGLlist);
numFacCMBList=length(facCMBlist);
numFreqList=length(freqList);
numRelList=length(relList);
numFacChanList=length(facChanList);
numFacPointList=length(facPointList);
numFacFreqList=length(facFreqList);

if isempty(combineWeightsIn)
    switch dataDimension
        case 'channels'
            numCombs=length(chanList);
            if ~isempty(EPdataIn.facVecS)
                numCombs=length(facChanList);
            end
        case 'cells'
            numCombs=length(cellList);
        case 'subjects'
            numCombs=length(subList);
        case 'factors'
            numCombs=length(facList);
        case 'convert'
            combineName='';
            numCombs=[];
            if isempty(EPdataIn.GAVsubs)
                %nothing to convert
                EPdataOutFinal=EPdataIn;
                return
            end
        otherwise
            disp('Oops - programmer error in ep_combineData')
            EPdataOut=[];
            return
    end
    combineWeights=ones(numCombs,1);
else
    combineWeights=combineWeightsIn;
end
combineWeights=combineWeights(:);

if (sum(combineWeights)~=0) && any(diff(sign(combineWeights)))
    disp('Error: combination weights need to either sum to zero or they need to all be the same sign.')
    EPdataOut=[];
    return
end

if sum(combineWeights)==0
    contrastMode=1;
else
    contrastMode=0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp('channels',dataDimension)
    combineLevel=find(strcmp(combineName,EPdataOut.chanNames));
    if isempty(combineLevel)
        if ~isempty(EPdataOut.facVecS)
            facNumChans=facNumChans+1;
            combineLevel=facNumChans;
        else
            numChans=numChans+1;
            combineLevel=numChans;
        end
        EPdataOut.chanNames{end+1,1}=combineName;
        theChanType=unique(EPdataOut.chanTypes(chanList));
        if length(theChanType) > 1
            disp('Oops - cannot combine multiple types of channels')
            EPdataOut=[];
            return
        end
        switch theChanType{1}
            case 'EEG'
                theRegType='REG';
            case 'BSC'
                theRegType='RBS';
        end
        EPdataOut.chanTypes{end+1,1}=theRegType;
        if ~isempty(EPdataOut.eloc)
            EPdataOut.eloc(end+1)=cell2struct(cell(1,length(fieldnames(EPdataOut.eloc))),fieldnames(EPdataOut.eloc),2);
        end
    end

    if numRels>1
        newCellDataPos=zeros(1,numPointList,numCellList,numSubList,numFacSGLlist,numFreqList,numRelList+1);
        newCellDataNeg=zeros(1,numPointList,numCellList,numSubList,numFacSGLlist,numFreqList,numRelList+1);
    else
        newCellDataPos=zeros(1,numPointList,numCellList,numSubList,numFacSGLlist,numFreqList,1);
        newCellDataNeg=zeros(1,numPointList,numCellList,numSubList,numFacSGLlist,numFreqList,1);
    end
    if ~isempty(EPdataOut.facVecT)
        newCellDataPos=newCellDataPos(:,1,:,:,:,:,:);
        newCellDataNeg=newCellDataNeg(:,1,:,:,:,:,:);
    end
    if ~isempty(EPdataOut.facVecF)
        newCellDataPos=newCellDataPos(:,:,:,:,:,1,:);
        newCellDataNeg=newCellDataNeg(:,:,:,:,:,1,:);
    end

    if ~isempty(EPdataOut.noise)
        newCellNoisePos=zeros(1,numPointList,numCellList,numSubList,numFacSGLlist,numFreqList);
        newCellNoiseNeg=zeros(1,numPointList,numCellList,numSubList,numFacSGLlist,numFreqList);
    end
    if ~isempty(EPdataOut.covAVE)
        if size(EPdataOut.covAVE,7)==1
            newCovAVE=zeros(numChans,numPointList,numCellList,numSubList,numFacCMBList,numFreqList,1);
        else
            newCovAVE=zeros(numChans,numPointList,numCellList,numSubList,numFacCMBList,numFreqList,numChans);
        end
    else
        newCovAVE=[];
    end

    newFacVecSData=zeros(1,numFacSGLlist);
    newBadChanData=zeros(numSubList,numCellList,1);

    newFacDataPos=zeros(1,numFacPointList,numCellList,numSubList,numFacCMBList,numFacFreqList,numRelList);
    newFacDataNeg=zeros(1,numFacPointList,numCellList,numSubList,numFacCMBList,numFacFreqList,numRelList);

    if numRels>1
        totalWeightPos=zeros(1,numPointList,numCellList,numSubList,numFacSGLlist,numFreqList,numRelList+1);
        totalWeightNeg=zeros(1,numPointList,numCellList,numSubList,numFacSGLlist,numFreqList,numRelList+1);
    else
        totalWeightPos=zeros(1,numPointList,numCellList,numSubList,numFacSGLlist,numFreqList,1);
        totalWeightNeg=zeros(1,numPointList,numCellList,numSubList,numFacSGLlist,numFreqList,1);
    end

    totalNoiseWeightPos=zeros(1,numPointList,numCellList,numSubList,numFacSGLlist,numFreqList,1);
    totalNoiseWeightNeg=zeros(1,numPointList,numCellList,numSubList,numFacSGLlist,numFreqList,1);

    totalFacWeightPos=zeros(1,numFacPointList,numCellList,numSubList,numFacCMBList,numFacFreqList,1);
    totalFacWeightNeg=zeros(1,numFacPointList,numCellList,numSubList,numFacCMBList,numFacFreqList,1);

    %no need to skip bad cells and subjects since the EP Toolkit will ignore them anyway and maybe they'll still be of interest
    for iSub=1:length(subList)
        theSub=subList(iSub);
        ep_tictoc;if EPtictoc.stop;return;end
        for iCell=1:length(cellList)
            theCell=cellList(iCell);
            if strcmp(EPdataOut.dataType,'average')
                goodChans=find(~isnan(EPdataOut.analysis.badChans(theSub,theCell,:)));
            else
                goodChans=find(EPdataOut.analysis.badChans(theSub,theCell,:) >= 0);
            end
            if ~isempty(EPdataOut.facVecS)
                if ~isempty(goodChans)
                    %for the single spatial factor channel, as long as at least one channel was good, it's all good.  Normal behavior for REG spatial factor channels.
                    if ~isempty(intersect(find(strcmp('EEG',EPdataOut.chanTypes(facChanList))),goodChans))
                        chanGoodList=find(strcmp('EEG',EPdataOut.chanTypes(facChanList)),1,'first');
                    else
                        chanGoodList=[];
                    end
                    % goodRegChans=intersect(find(strcmp('REG',EPdataOut.chanTypes(chanList))),goodChans);
                    % if ~isempty(goodRegChans)
                    %     chanGoodList=[chanGoodList; goodRegChans];
                    % end
                    chanGoodWeights=1;
                    facChanGoodList=intersect(facChanList,goodChans);
                    facChanGoodWeights=combineWeights(ismember(facChanList,facChanGoodList));
                else
                    chanGoodList=[];
                    facChanGoodList=[];
                end
            else
                chanGoodList=intersect(chanList,goodChans);
                chanGoodWeights=combineWeights(ismember(chanList,chanGoodList));
            end
            if ~isempty(chanGoodList)
                if ~isempty(EPdataOut.facVecS)
                    facVecTotalWeight=sum(abs(facChanGoodWeights));
                end
                if numRels>1
                    for iFac=1:numFacSGLlist
                        theFac=facSGLlist(iFac);
                        for iPoint=1:numPointList
                            thePoint=freqList(iPoint);
                            for iFreq=1:numFreqList
                                theFreq=freqList(iFreq);
                                outCovMat=ep_covMat(squeeze(EPdataOut.data(:,thePoint,theCell,theSub,theFac,theFreq,:)), chanList, 'combineChan', combineWeights);
                                newCellDataPos(1,iPoint,iCell,iSub,theFac,iFreq,:)=outCovMat(end,:);
                            end
                        end
                    end
                else
                    for iChan=1:length(chanGoodList)
                        theChan=chanGoodList(iChan);
                        theWeight=chanGoodWeights(iChan);
                        if ~isempty(EPdataOut.facVecS)
                            %this loop will just run once for spatial PCA data
                            newFacVecSData=sum((1/facVecTotalWeight)*diag(facChanGoodWeights)*EPdataOut.facVecS(facChanList,facSGLlist),1);
                        else
                            theData=EPdataOut.data(theChan,pointList,theCell,theSub,facSGLlist,freqList,:);
                            if ~isempty(EPdataOut.freqNames) && any(strcmp(EPdataOut.chanTypes{theChan},{'EEG','REG'}))
                                theData=abs(theData); %convert to amplitude scaling when adding freq data together.
                            end
                            goodData=~isnan(theData);
                            theData(~goodData)=0;
                            if theWeight>0
                                newCellDataPos(1,:,iCell,iSub,:,:,:)=newCellDataPos(1,:,iCell,iSub,:,:,:)+(theWeight*theData);
                                totalWeightPos(1,:,iCell,iSub,:,:,:)=totalWeightPos(1,:,iCell,iSub,:,:,:)+(goodData*abs(theWeight));
                            elseif theWeight<0
                                newCellDataNeg(1,:,iCell,iSub,:,:,:)=newCellDataNeg(1,:,iCell,iSub,:,:,:)+(theWeight*theData);
                                totalWeightNeg(1,:,iCell,iSub,:,:,:)=totalWeightNeg(1,:,iCell,iSub,:,:,:)+(goodData*abs(theWeight));
                            end
                        end

                        if ~isempty(EPdataOut.noise)
                            theNoise=EPdataOut.noise(theChan,pointList,theCell,theSub,facSGLlist,freqList);
                            goodData=~isnan(theData);
                            theData(~goodData)=0;
                            if theWeight>0
                                newCellNoisePos(1,:,iCell,iSub,:,:)=newCellNoisePos(1,:,iCell,iSub,:,:)+(theWeight*theNoise);
                                totalNoiseWeightPos(1,:,iCell,iSub,:,:,:)=totalNoiseWeightPos(1,:,iCell,iSub,:,:,:)+(goodData*abs(theWeight));
                            elseif theWeight<0
                                newCellNoiseNeg(1,:,iCell,iSub,:,:)=newCellNoiseNeg(1,:,iCell,iSub,:,:)+(theWeight*theNoise);
                                totalNoiseWeightNeg(1,:,iCell,iSub,:,:,:)=totalNoiseWeightNeg(1,:,iCell,iSub,:,:,:)+(goodData*abs(theWeight));
                            end
                        end

                        theBadData=EPdataOut.analysis.badChans(theSub,theCell,theChan);
                        if strcmp(EPdataOut.dataType,'average')
                            if ~xor(sign(theBadData),sign(newBadChanData(iSub,iCell,1)))
                                newBadChanData(iSub,iCell,1)= newBadChanData(iSub,iCell,1)+theBadData;
                            elseif isnan(newBadChanData(iSub,iCell,1)) || isnan(theBadData)
                                newBadChanData(iSub,iCell,1)=NaN;
                            elseif ~sign(theBadData)
                                newBadChanData(iSub,iCell,1)=theBadData;
                            else
                                newBadChanData(iSub,iCell,1)=newBadChanData(iSub,iCell,1);
                            end
                        else
                            if ~sign(theBadData) || ~sign(newBadChanData(iSub,iCell,1))
                                newBadChanData(iSub,iCell,1)=-1;
                            else
                                newBadChanData(iSub,iCell,1)= 0;
                            end
                        end
                    end
                end
            else %no good data available for this regional channel
                if strcmp(EPdataOut.dataType,'average')
                    newBadChanData(iSub,iCell,1)=NaN;
                else
                    newBadChanData(iSub,iCell,1)=-1;
                end
            end

            %handle facData, which has a full set of channels even for spatial factors
            if ~isempty(EPdataOut.facData) && ~isempty(facChanGoodList)
                for iChan=1:length(facChanGoodList)
                    theChan=facChanGoodWeights(iChan);
                    theData=EPdataOut.facData(theChan,facPointList,theCell,theSub,facCMBlist,facFreqList,:);
                    theWeight=facChanGoodWeights(iChan);
                    if ~isempty(EPdataOut.freqNames) && any(strcmp(EPdataOut.chanTypes{theChan},{'EEG','REG'})) && isempty(EPdataOut.relNames)
                        theData=abs(theData); %convert to amplitude scaling when adding freq data together.
                    end
                    goodData=~isnan(theData);
                    theData(~goodData)=0;
                    if theWeight>0
                        newFacDataPos(1,:,iCell,iSub,:,:,:)=newFacDataPos(1,:,iCell,iSub,:,:,:)+(theWeight*theData);
                        totalFacWeightPos(1,:,iCell,iSub,:,:,:)=totalFacWeightPos(1,:,iCell,iSub,:,:,:)+(goodData*abs(theWeight));
                    elseif theWeight<0
                        newFacDataNeg(1,:,iCell,iSub,:,:,:)=newFacDataNeg(1,:,iCell,iSub,:,:,:)+(theWeight*theData);
                        totalFacWeightNeg(1,:,iCell,iSub,:,:,:)=totalFacWeightNeg(1,:,iCell,iSub,:,:,:)+(goodData*abs(theWeight));
                    end
                end
            end

            if ~isempty(newCovAVE)
                if size(newCovAVE,7)==1
                    newCovAVE=[]; %without covariance information, can't combine channel variances
                    EPdataOut.covAVE=[];
                else
                    for iPoint=1:numPointList
                        thePoint=pointList(iPoint);
                        for iFac=1:numFacSGLlist
                            theFac=facSGLlist(iFac);
                            for iFreq=1:numFreqList
                                theFreq=freqList(iFreq);
                                newCovAVE(:,iPoint,iCell,iSub,iFac,iFreq,:)=ep_covMat(EPdataOut.covAVE(:,thePoint,theCell,theSub,theFac,theFreq,:), chanGoodList, 'combineChan');
                            end
                        end
                    end
                end
            end
        end
    end
    if ~isempty(EPdataOut.facVecS)
        EPdataOut.facVecS(combineLevel,facSGLlist)=newFacVecSData;
    else
        if contrastMode
            %zero weights for any datapoint in either positive or negative matrix will result in an NaN.
            EPdataOut.data(combineLevel,pointList,cellList,subList,facSGLlist,freqList,relList)=newCellDataPos./totalWeightPos+newCellDataNeg./totalWeightNeg;
        else
            %zero weights for any datapoint will result in an NaN.
            if ~isempty(EPdataOut.relNames)
                numRels=numRels+1; %the new channel
                EPdataOut.relNames{end+1,1}=combineName;
                EPdataOut.data(end+1,pointList,cellList,subList,facSGLlist,freqList,[relList; numRels])=newCellDataPos;
                %permute(newCellDataPos,[7 2 3 4 5 6 1]);
                EPdataOut.data(end,:,:,:,:,:,end)=1; %coherence of channel with itself is real number one.
            else
                if sum(combineWeights)>0
                    newChanData=newCellDataPos./totalWeightPos;
                    EPdataOut.data(combineLevel,pointList,cellList,subList,facSGLlist,freqList,relList)=newChanData;
                elseif sum(combineWeights)<0
                    newChanData=newCellDataNeg./totalWeightNeg;
                    EPdataOut.data(combineLevel,pointList,cellList,subList,facSGLlist,freqList,relList)=newChanData;
                end
            end
        end
    end

    if ~isempty(EPdataOut.facData)
        if contrastMode
            %zero weights for any datapoint in either positive or negative matrix will result in an NaN.
            EPdataOut.facData(combineLevel,facPointList,cellList,subList,facCMBlist,freqList,relList)=newFacDataPos./totalFacWeightPos+newFacDataNeg./totalFacWeightNeg;
        else
            %zero weights for any datapoint will result in an NaN.
            if sum(combineWeights)>0
                EPdataOut.facData(combineLevel,facPointList,cellList,subList,facCMBlist,freqList,relList)=newFacDataPos./totalFacWeightPos;
            elseif sum(combineWeights)<0
                EPdataOut.facData(combineLevel,facPointList,cellList,subList,facCMBlist,freqList,relList)=newFacDataNeg./totalFacWeightNeg;
            end
        end
    end

    EPdataOut.analysis.badChans(subList,cellList,combineLevel)=newBadChanData;

    if ~isempty(EPdataOut.noise)
        if contrastMode
            %zero weights for any datapoint in either positive or negative matrix will result in an NaN.
            EPdataOut.noise(combineLevel,pointList,cellList,subList,facSGLlist,freqList,relList)=newCellNoisePos./totalNoiseWeightPos+newCellNoiseNeg./totalNoiseWeightNeg;
        else
            %zero weights for any datapoint will result in an NaN.
            if sum(combineWeights)>0
                EPdataOut.noise(combineLevel,pointList,cellList,subList,facSGLlist,freqList,relList)=newCellNoisePos./totalNoiseWeightPos;
            elseif sum(combineWeights)<0
                EPdataOut.noise(combineLevel,pointList,cellList,subList,facSGLlist,freqList,relList)=newCellNoiseNeg./totalNoiseWeightNeg;
            end
        end
    end
    if ~isempty(EPdataOut.covAVE)
        if ~isempty(newCovAVE)
            EPdataOut.covAVE=newCovAVE;
        else
            EPdataOut.covAVE(combineLevel,pointList,cellList,subList,facSGLlist,freqList,:)=NaN;
            if size(EPdataOut.covAVE,7)~=1
                EPdataOut.covAVE(:,pointList,cellList,subList,facSGLlist,freqList,combineLevel)=NaN;
            end
        end
    end
    if ~isempty(EPdataOut.cov)
        EPdataOut.cov.covMatrix(:,combineLevel,combineLevel)=NaN;
    end

    if ~isempty(EPdataOut.impedances.channels)
        EPdataOut.impedances.channels(end+1,:)=NaN;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(dataDimension,'cells') || (strcmp(dataDimension,'convert') && (size(EPdataOut.GAVsubs,2) > 1))
    if ~isempty(combineName)
        combineList=find(strcmp(combineName,EPdataOut.cellNames), 1);
        if isempty(combineList)
            %add a new virtual cell if the name is a new one,
            %otherwise an existing one will just be updated.
            %if a new cell, then just need to update GAVsubs structure.
            numCells=numCells+1;
            EPdataOut.cellNames{end+1,1}=combineName;
            EPdataOut.cellTypes{end+1,1}='CMB'; %even if technically SGL, doing so would often result in unexpected behavior for the user.
            if ~isempty(EPdataOut.trialNames)
                EPdataOut.trialNames(end+1,1)=1;
            end
            if isempty(EPdataOut.GAVsubs)
                %initialize a new GAVsubs structure
                EPdataOut.GAVsubs=cell(1,2,numFacs);
            else
                EPdataOut.GAVsubs{end,end+1,end}=[];
            end
            EPdataOut.GAVsubs{1,end,1}=[cellList combineWeights];
            for iGav=2:size(EPdataOut.GAVsubs,1)
                for iFac=1:size(EPdataOut.GAVsubs,3)
                    EPdataOut.GAVsubs{iGav,end,iFac}=[subList ones(size(subList))];
                end
            end
        end
    else
        %convert all the virtual cells
        combineList=[1:numVcells]+numRcells;
        cellList=[1:numCells]';
        numCellList=length(cellList);
    end
    for iCombine=1:length(combineList)
        combineLevel=combineList(iCombine);
        %either update an existing normal grand average or turn an existing virtual one into a real one.
        if strcmp(dataDimension,'convert') || (strcmp(dataDimension,'cells') && (combineLevel>numRcells) && (combineLevel<=numCells))
            cellList2=EPdataOut.GAVsubs{1,combineLevel-numRcells+1,1}(:,1);
            combineWeights=EPdataOut.GAVsubs{1,combineLevel-numRcells+1,1}(:,2);
            if sum(combineWeights)==0
                contrastMode=1;
            else
                contrastMode=0;
            end
            if strcmp(dataDimension,'cells')
                %convert the virtual cell into a real cell
                theGAVcell=combineLevel-numRcells;
                %real subjects need to all be listed prior to the virtual ones
                EPdataOut=ep_reorderData(EPdataOut,'cells',[1:numRcells combineLevel setdiff([1:numVcells],theGAVcell)+numRcells]);
                if isempty(EPdataOut)
                    EPdataOut=[];
                    return
                end
                EPdataOut.GAVsubs(:,theGAVcell+1,:)=[];
                numRcells=numRcells+1;
                combineLevel=numRcells;
                numVcells=numVcells-1;
            end
        else
            cellList2=cellList;
        end

        newCellDataPos=zeros(numChanList,numPointList,1,numSubList,numFacSGLlist,numFreqList,numRelList);
        newCellDataNeg=zeros(numChanList,numPointList,1,numSubList,numFacSGLlist,numFreqList,numRelList);

        if ~isempty(EPdataOut.facVecT)
            newCellDataPos=newCellDataPos(:,1,:,:,:,:,:);
            newCellDataNeg=newCellDataNeg(:,1,:,:,:,:,:);
        end
        if ~isempty(EPdataOut.facVecF)
            newCellDataPos=newCellDataPos(:,:,:,:,:,1,:);
            newCellDataNeg=newCellDataNeg(:,:,:,:,:,1,:);
        end
        if ~isempty(EPdataOut.facVecS)
            newCellDataPos=newCellDataPos(1,:,:,:,:,:,:);
            newCellDataNeg=newCellDataNeg(1,:,:,:,:,:,:);
        end

        if ~isempty(EPdataOut.noise)
            newCellNoisePos=zeros(numChanList,numPointList,1,numSubList,numFacSGLlist,numFreqList);
            newCellNoiseNeg=zeros(numChanList,numPointList,1,numSubList,numFacSGLlist,numFreqList);
        end
        if ~isempty(EPdataOut.covAVE)
            if size(EPdataOut.covAVE,7)==1
                newCovAVE=zeros(numChans,numPointList,1,numSubList,numFacCMBList,numFreqList,1);
            else
                newCovAVE=zeros(numChans,numPointList,1,numSubList,numFacCMBList,numFreqList,numChans);
            end
        else
            newCovAVE=[];
        end

        newBadChanData=zeros(numSubList,1,numChanList);

        newFacDataPos=zeros(numFacChanList,numFacPointList,1,numSubList,numFacCMBList,numFacFreqList,numRelList);
        newFacDataNeg=zeros(numFacChanList,numFacPointList,1,numSubList,numFacCMBList,numFacFreqList,numRelList);

        totalWeightPos=zeros(numChans,numPointList,1,numSubList,numFacSGLlist,numFreqList,numRelList);
        totalWeightNeg=zeros(numChans,numPointList,1,numSubList,numFacSGLlist,numFreqList,numRelList);

        totalNoiseWeightPos=zeros(numChans,numPointList,1,numSubList,numFacSGLlist,numFreqList,1);
        totalNoiseWeightNeg=zeros(numChans,numPointList,1,numSubList,numFacSGLlist,numFreqList,1);

        totalFacWeightPos=zeros(numFacChanList,numFacPointList,1,numSubList,numFacCMBList,numFacFreqList,1);
        totalFacWeightNeg=zeros(numFacChanList,numFacPointList,1,numSubList,numFacCMBList,numFacFreqList,1);

        for iSub=1:length(subList)
            theSub=subList(iSub);
            ep_tictoc;if EPtictoc.stop;return;end

            cellGoodList=cellList2(EPdataOut.avgNum(theSub,cellList2)>=0);
            cellGoodWeights=combineWeights(ismember(cellList2,cellGoodList));

            if trialWeights && ~any(EPdataOut.avgNum(theSub,cellGoodList) == 0) && strcmp(EPdataOut.dataType,'average')
                useTrialWeights=1;
            else
                useTrialWeights=0;
            end
            if ~sum(combineWeights) && sum(cellGoodWeights)
                cellGoodList=[]; %if this was meant to be a difference wave and enough cells are bad to no longer be balanced, then treat as invalid.
            end
            if ~isempty(cellGoodList)

                posCellList=cellList2(combineWeights>0);
                negCellList=cellList2(combineWeights<0);

                posCellGoodList=posCellList(EPdataOut.avgNum(theSub,posCellList)>=0);
                negCellGoodList=negCellList(EPdataOut.avgNum(theSub,negCellList)>=0);

                posGoodWeights=combineWeights(ismember(cellList2,posCellGoodList));
                negGoodWeights=combineWeights(ismember(cellList2,negCellGoodList));
                posGoodWeights=posGoodWeights./sum(posGoodWeights);
                negGoodWeights=negGoodWeights./sum(negGoodWeights);


                %just sum up size of sample involved in this combined waveform, except for subNum where it is the largest number involved.
                EPdataOut.avgNum(theSub,combineLevel)=sum(EPdataOut.avgNum(theSub,cellGoodList),'omitnan');
                EPdataOut.subNum(theSub,combineLevel)=max(EPdataOut.subNum(theSub,cellGoodList));

                if ~isempty(EPdataOut.covNum)
                    if useTrialWeights
                        %if using trialweights, then the new effective sample size is simply the total trials.
                        EPdataOut.covNum(theSub,combineLevel)=sum(EPdataOut.covNum(theSub,cellGoodList));
                    else
                        %assume cov matrix is the same so just need to figure out the new effective sample size for the new combination
                        %per p.128 of the 3.7.2 MNE manual, 1/Leff=Sigma weight-squared/L
                        covWeights=combineWeights(ismember(cellList2,cellGoodList));
                        covWeights=covWeights/sum(abs(covWeights));
                        EPdataOut.covNum(theSub,combineLevel)=sum([EPdataOut.covNum(theSub,cellGoodList).^-1]'.*covWeights.^2).^-1;
                    end
                end

                %for these, weighted mean of the involved data and saccadeOnset is the difference in the onsets for contrasts.
                EPdataOut.analysis.blinkTrial(theSub,combineLevel)= sum(EPdataOut.analysis.blinkTrial(theSub,posCellGoodList)'.*posGoodWeights,'omitnan')+sum(EPdataOut.analysis.blinkTrial(theSub,negCellGoodList)'.*negGoodWeights,'omitnan');
                EPdataOut.analysis.saccadeTrial(theSub,combineLevel)= sum(EPdataOut.analysis.saccadeTrial(theSub,posCellGoodList)'.*posGoodWeights,'omitnan')+sum(EPdataOut.analysis.saccadeTrial(theSub,negCellGoodList)'.*negGoodWeights,'omitnan');
                EPdataOut.analysis.saccadeOnset(theSub,combineLevel)= sum(EPdataOut.analysis.saccadeOnset(theSub,posCellGoodList)'.*posGoodWeights,'omitnan')+sum(EPdataOut.analysis.saccadeOnset(theSub,negCellGoodList)'.*negGoodWeights,'omitnan');
                EPdataOut.analysis.moveTrial(theSub,combineLevel)= sum(EPdataOut.analysis.moveTrial(theSub,posCellGoodList)'.*posGoodWeights,'omitnan')+sum(EPdataOut.analysis.moveTrial(theSub,negCellGoodList)'.*negGoodWeights,'omitnan');
                EPdataOut.analysis.badTrials(theSub,combineLevel)= sum(EPdataOut.analysis.badTrials(theSub,posCellGoodList)'.*posGoodWeights,'omitnan')+sum(EPdataOut.analysis.badTrials(theSub,negCellGoodList)'.*negGoodWeights,'omitnan');

                EPdataOut.recTime(combineLevel)=min(EPdataOut.recTime(cellGoodList));

                %for trial specs, for numerical data weighted differences for contrasts.
                if ~isempty(EPdataOut.trialSpecs)
                    %generate the averaged trial specs
                    for iSpec=1:length(EPdataOut.trialSpecNames)
                        posSpecs=squeeze(EPdataOut.trialSpecs(posCellGoodList,iSpec,theSub));
                        negSpecs=squeeze(EPdataOut.trialSpecs(negCellGoodList,iSpec,theSub));
                        allSpecs=[posSpecs; negSpecs];
                        if ischar(allSpecs{1}) && all(strcmp(allSpecs{1},allSpecs)) %if the spec is all just the same character string then just set it to that string.
                            EPdataOut.trialSpecs{combineLevel,iSpec,theSub}=allSpecs{1};
                        else
                            %convert char to numeric.  char strings that are not numbers will convert to NaN and then be ignored during the computations.
                            theStrSpecs=find(cellfun(@ischar,posSpecs));
                            for iStrSpec=1:length(theStrSpecs)
                                theStrSpec=theStrSpecs(iStrSpec);
                                posSpecs{theStrSpec}=str2double(posSpecs{theStrSpec});
                                if ~isnumeric(posSpecs{theStrSpec})
                                    posSpecs{theStrSpec}=NaN; %dealing with weird Matlab bug where 'cab' was being converted into a driver handle
                                end
                            end
                            theStrSpecs=find(cellfun(@ischar,negSpecs));
                            for iStrSpec=1:length(theStrSpecs)
                                theStrSpec=theStrSpecs(iStrSpec);
                                negSpecs{theStrSpec}=str2double(negSpecs{theStrSpec});
                                if ~isnumeric(negSpecs{theStrSpec})
                                    negSpecs{theStrSpec}=NaN; %dealing with weird Matlab bug where 'cab' was being converted into a driver handle
                                end
                            end

                            %reweight for cases where there are missing values
                            posGoodWeightsSpec=combineWeights(find(ismember(cellGoodList,posCellGoodList)));
                            posGoodWeightsSpec=posGoodWeightsSpec(~cellfun(@isempty,posSpecs));
                            negGoodWeightsSpec=combineWeights(find(ismember(cellGoodList,negCellGoodList)));
                            negGoodWeightsSpec=negGoodWeightsSpec(~cellfun(@isempty,negSpecs));
                            posGoodWeightsSpec=posGoodWeightsSpec./sum(posGoodWeightsSpec);
                            negGoodWeightsSpec=negGoodWeightsSpec./sum(negGoodWeightsSpec);

                            if contrastMode
                                %zero weights for any datapoint in either positive or negative specs will result in an NaN.  Positive and negative sides given equal weight to each other.
                                if all(cellfun(@isempty,posSpecs)) || all(cellfun(@isempty,negSpecs))
                                    EPdataOut.trialSpecs{combineLevel,iSpec,theSub}=NaN;
                                else
                                    if any(strcmp(EPdataOut.trialSpecNames(iSpec),{'RT','ACC'}))
                                        EPdataOut.trialSpecs{combineLevel,iSpec,theSub}=(sum(cell2mat(posSpecs).*posGoodWeightsSpec,'omitnan')+sum(cell2mat(negSpecs).*negGoodWeightsSpec,'omitnan'))/2;
                                    else
                                        EPdataOut.trialSpecs{combineLevel,iSpec,theSub}=sum(cell2mat(posSpecs).*posGoodWeightsSpec,'omitnan')-sum(cell2mat(negSpecs).*negGoodWeightsSpec,'omitnan');
                                    end
                                end
                            else
                                %zero weights for any datapoint will result in an NaN.
                                if sum(combineWeights)>0
                                    if all(cellfun(@isempty,posSpecs))
                                        EPdataOut.trialSpecs{combineLevel,iSpec,theSub}=NaN;
                                    else
                                        EPdataOut.trialSpecs{combineLevel,iSpec,theSub}=sum(cell2mat(posSpecs).*posGoodWeightsSpec,'omitnan');
                                    end
                                elseif sum(combineWeights)<0
                                    if all(cellfun(@isempty,negSpecs))
                                        EPdataOut.trialSpecs{combineLevel,iSpec,theSub}=NaN;
                                    else
                                        EPdataOut.trialSpecs{combineLevel,iSpec,theSub}=sum(cell2mat(negSpecs).*negGoodWeightsSpec,'omitnan');
                                    end
                                end
                            end
                        end
                    end
                else
                    EPdataOut.trialSpecs=cell(numCells,0,numSubs);
                end

            else
                EPdataOut.avgNum(theSub,combineLevel)=-1;
                EPdataOut.covNum(theSub,combineLevel)=-1;
                EPdataOut.subNum(theSub,combineLevel)=-1;
                EPdataOut.analysis.blinkTrial(theSub,combineLevel)= 0;
                EPdataOut.analysis.saccadeTrial(theSub,combineLevel)= 0;
                EPdataOut.analysis.saccadeOnset(theSub,combineLevel)= 0;
                EPdataOut.analysis.moveTrial(theSub,combineLevel)= 0;
                EPdataOut.analysis.badTrials(theSub,combineLevel)= 0;
                newBadChanData(iSub,1,:)= 0;
                EPdataOut.recTime(combineLevel)=1;
            end

            for iChan=1:length(chanList)
                theChan=chanList(iChan);
                if strcmp(EPdataOut.dataType,'average')
                    cellChanGoodList=cellGoodList(~isnan(EPdataOut.analysis.badChans(theSub,cellGoodList,theChan)));
                else
                    cellChanGoodList=cellGoodList(EPdataOut.analysis.badChans(theSub,cellGoodList,theChan) >= 0);
                end
                cellChanGoodWeights=combineWeights(ismember(cellList2,cellChanGoodList));
                if ~isempty(cellChanGoodList)
                    for iCell=1:length(cellChanGoodList)
                        theCell=cellChanGoodList(iCell);
                        if useTrialWeights
                            theWeight=cellChanGoodWeights(iCell)*EPdataOut.avgNum(theSub,theCell);
                        else
                            theWeight=cellChanGoodWeights(iCell);
                        end
                        if (iChan==1) || isempty(EPdataOut.facVecS)
                            theData=EPdataOut.data(theChan,pointList,theCell,theSub,facSGLlist,freqList,:);
                            if ~isempty(EPdataOut.freqNames) && any(strcmp(EPdataOut.chanTypes{theChan},{'EEG','REG'})) && isempty(EPdataOut.relNames)
                                theData=abs(theData); %convert to amplitude scaling when adding freq data together.
                            end
                            goodData=~isnan(theData);
                            theData(~goodData)=0;
                            if theWeight>0
                                newCellDataPos(iChan,:,1,iSub,:,:,:)=newCellDataPos(iChan,:,1,iSub,:,:,:)+(theWeight*theData);
                                totalWeightPos(iChan,:,1,iSub,:,:,:)=totalWeightPos(iChan,:,1,iSub,:,:,:)+(goodData*abs(theWeight));
                            elseif theWeight<0
                                newCellDataNeg(iChan,:,1,iSub,:,:,:)=newCellDataNeg(iChan,:,1,iSub,:,:,:)+(theWeight*theData);
                                totalWeightNeg(iChan,:,1,iSub,:,:,:)=totalWeightNeg(iChan,:,1,iSub,:,:,:)+(goodData*abs(theWeight));
                            end
                        end
                        if ~isempty(EPdataOut.facData)
                            if ~isempty(EPdataOut.facVecS)
                                theChans=facChanList;
                            else
                                theChans=theChan;
                            end
                            theData=EPdataOut.facData(theChans,facPointList,theCell,theSub,facCMBlist,facFreqList,:);
                            goodData=~isnan(theData);
                            theData(~goodData)=0;
                            if ~isempty(EPdataOut.freqNames) && any(strcmp(EPdataOut.chanTypes{theChan},{'EEG','REG'})) && isempty(EPdataOut.relNames)
                                theData=abs(theData); %convert to amplitude scaling when adding freq data together.
                            end
                            if theWeight>0
                                newFacDataPos(theChans,:,1,iSub,:,:,:)=newFacDataPos(theChans,:,1,iSub,:,:,:)+(theWeight*theData);
                                totalFacWeightPos(theChans,:,1,iSub,:,:,:)=totalFacWeightPos(theChans,:,1,iSub,:,:,:)+(goodData*abs(theWeight));
                            elseif theWeight<0
                                newFacDataNeg(theChans,:,1,iSub,:,:,:)=newFacDataNeg(theChans,:,1,iSub,:,:,:)+(theWeight*theData);
                                totalFacWeightNeg(theChans,:,1,iSub,:,:,:)=totalFacWeightNeg(theChans,:,1,iSub,:,:,:)+(goodData*abs(theWeight));
                            end
                        end
                        if ~isempty(EPdataOut.noise)
                            theNoise=EPdataOut.noise(theChan,pointList,theCell,theSub,facSGLlist,freqList);
                            goodData=~isnan(theData);
                            theData(~goodData)=0;
                            if theWeight>0
                                newCellNoisePos(iChan,:,1,iSub,:,:)=newCellNoisePos(iChan,:,1,iSub,:,:)+(theWeight*theNoise);
                                totalNoiseWeightPos(iChan,:,1,iSub,:,:)=totalNoiseWeightPos(iChan,:,1,iSub,:,:)+(goodData*abs(theWeight));
                            elseif theWeight<0
                                newCellNoiseNeg(iChan,:,1,iSub,:,:)=newCellNoiseNeg(iChan,:,1,iSub,:,:)+(theWeight*theNoise);
                                totalNoiseWeightNeg(iChan,:,1,iSub,:,:)=totalNoiseWeightNeg(iChan,:,1,iSub,:,:)+(goodData*abs(theWeight));
                            end
                        end
                        theBadData=EPdataOut.analysis.badChans(theSub,theCell,theChan);
                        if strcmp(EPdataOut.dataType,'average')
                            if ~xor(sign(theBadData),sign(newBadChanData(iSub,1,iChan)))
                                newBadChanData(iSub,1,iChan)= newBadChanData(iSub,1,iChan)+theBadData;
                            elseif isnan(newBadChanData(iSub,1,iChan)) || isnan(theBadData)
                                newBadChanData(iSub,1,iChan)=NaN;
                            elseif ~sign(theBadData)
                                newBadChanData(iSub,1,iChan)=theBadData;
                            else
                                newBadChanData(iSub,1,iChan)=newBadChanData(iSub,1,iChan);
                            end
                        else
                            if ~sign(theBadData) || ~sign(newBadChanData(iSub,1,iChan))
                                newBadChanData(iSub,1,iChan)=-1;
                            else
                                newBadChanData(iSub,1,iChan)= 0;
                            end
                        end
                    end
                else %none of the cells have good data for this channel
                    if strcmp(EPdataOut.dataType,'average')
                        newBadChanData(iSub,1,iChan)=NaN;
                    else
                        newBadChanData(iSub,1,iChan)=-1;
                    end
                end
            end

            if ~isempty(EPdataOut.covAVE)
                if size(EPdataOut.covAVE,7)==1
                    covRel=1;
                else
                    covRel=chanList;
                end
                for iPoint=1:numPointList
                    thePoint=pointList(iPoint);
                    for iFac=1:numFacSGLlist
                        theFac=facSGLlist(iFac);
                        for iFreq=1:numFreqList
                            theFreq=freqList(iFreq);
                            if isempty(cellGoodList) || (isscalar(cellGoodList))
                                newCovAVE(:,iPoint,1,iSub,iFac,iFreq,:)=NaN;
                            else
                                newCovAVE(:,iPoint,1,iSub,iFac,iFreq,:)=ep_covMat(EPdataOut.covAVE(chanList,thePoint,cellGoodList,theSub,theFac,theFreq,covRel), [], 'poolObs', cellGoodWeights, [], []);
                                ep_tictoc;if EPtictoc.stop;return;end
                            end
                        end
                    end
                end
            end
        end

        %events
        for iSub=1:length(subList)
            theSub=subList(iSub);
            for iCell=1:length(cellList2)
                theCell=cellList2(iCell);
                EPdataOut.events{theSub,combineLevel}=struct('type',{},'sample',{},'value',{},'duration',{},'keys',struct('code','','data','','datatype','','description',''));
                for iEvent=1:length(EPdataOut.events{theSub,theCell})
                    EPdataOut.events{theSub,combineLevel}(end+1)=EPdataOut.events{theSub,theCell}(iEvent);
                end
            end
        end

        if contrastMode
            %zero weights for any datapoint in either positive or negative matrix will result in an NaN.
            EPdataOut.data(chanList,pointList,combineLevel,subList,facSGLlist,freqList,relList)=newCellDataPos./totalWeightPos+newCellDataNeg./totalWeightNeg;
        else
            %zero weights for any datapoint will result in an NaN.
            if sum(combineWeights)>0
                EPdataOut.data(chanList,pointList,combineLevel,subList,facSGLlist,freqList,relList)=newCellDataPos./totalWeightPos;
            elseif sum(combineWeights)<0
                EPdataOut.data(chanList,pointList,combineLevel,subList,facSGLlist,freqList,relList)=newCellDataNeg./totalWeightNeg;
            end
        end

        EPdataOut.analysis.badChans(subList,combineLevel,chanList)=newBadChanData;

        if contrastMode
            %zero weights for any datapoint in either positive or negative matrix will result in an NaN.
            EPdataOut.facData(facChanList,facPointList,combineLevel,subList,facCMBlist,facFreqList,relList)=newFacDataPos./totalFacWeightPos+newFacDataNeg./totalFacWeightNeg;
        else
            %zero weights for any datapoint will result in an NaN.
            if sum(combineWeights)>0
                EPdataOut.facData(facChanList,facPointList,combineLevel,subList,facCMBlist,facFreqList,relList)=newFacDataPos./totalFacWeightPos;
            elseif sum(combineWeights)<0
                EPdataOut.facData(facChanList,facPointList,combineLevel,subList,facCMBlist,facFreqList,relList)=newFacDataNeg./totalFacWeightNeg;
            end
        end

        if ~isempty(EPdataOut.noise)
            if contrastMode
                %zero weights for any datapoint in either positive or negative matrix will result in an NaN.
                EPdataOut.noise(chanList,pointList,combineLevel,subList,facSGLlist,freqList)=newCellNoisePos./totalNoiseWeightPos+newCellNoiseNeg./totalNoiseWeightNeg;
            else
                %zero weights for any datapoint will result in an NaN.
                if sum(combineWeights)>0
                    EPdataOut.noise(chanList,pointList,combineLevel,subList,facSGLlist,freqList)=newCellNoisePos./totalNoiseWeightPos;
                elseif sum(combineWeights)<0
                    EPdataOut.noise(chanList,pointList,combineLevel,subList,facSGLlist,freqList)=newCellNoiseNeg./totalNoiseWeightNeg;
                end
            end
        end
        if ~isempty(EPdataOut.covAVE)
            if size(EPdataOut.covAVE,7)==1
                covRel=1;
            else
                covRel=chanList;
            end
            if ~isempty(newCovAVE)
                EPdataOut.covAVE(chanList,pointList,combineLevel,subList,facSGLlist,freqList,covRel)=newCovAVE;
            else
                EPdataOut.covAVE(chanList,pointList,combineLevel,subList,facSGLlist,freqList,covRel)=NaN;
            end
        end
    end
    if strcmp(dataDimension,'convert')
        EPdataOut.GAVsubs(:,2:end,:)=[];
        numRcells=numRcells+numVcells;
        numVcells=0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(dataDimension,'subjects') || (strcmp(dataDimension,'convert') && (size(EPdataIn.GAVsubs,1) > 1))
    if ~isempty(combineName)
        combineList=find(strcmp(combineName,EPdataOut.subNames), 1);
        if isempty(combineList)
            %add a new virtual grand average if the name is a new one,
            %otherwise an existing one will just be updated.
            %if a new virtual grand average, then just need to update GAVsubs structure.
            numSubs=numSubs+1;
            numVsubs=numVsubs+1;
            EPdataOut.subNames{end+1,1}=combineName;
            EPdataOut.subTypes{end+1,1}='GAV';
            if isempty(EPdataOut.GAVsubs)
                %initialize a new GAVsubs structure
                EPdataOut.GAVsubs=cell(2,1,numFacs);
            else
                EPdataOut.GAVsubs{end+1,end,end}=[];
            end
            for iCell=1:size(EPdataOut.GAVsubs,2)
                for iFac=1:size(EPdataOut.GAVsubs,3)
                    EPdataOut.GAVsubs{end,iCell,iFac}=[subList combineWeights];
                end
            end
        end
    else
        %convert all the virtual subs if no name is specified
        combineList=[1:numVsubs]+numRsubs;
        cellList=[1:numCells]'; %and all the cells
        numCellList=length(cellList);
    end
    for iCombine=1:length(combineList)
        combineLevel=combineList(iCombine);
        %either update an existing normal grand average or turn an existing virtual one into a real one.
        if strcmp(dataDimension,'convert') || (strcmp(dataDimension,'subjects') && (combineLevel>numRsubs) && (combineLevel<=numSubs))
            subList=EPdataOut.GAVsubs{combineLevel-numRsubs+1,1}(:,1);
            combineWeights=EPdataOut.GAVsubs{combineLevel-numRsubs+1,1}(:,2);
            if sum(combineWeights)==0
                contrastMode=1;
            else
                contrastMode=0;
            end
            if strcmp(dataDimension,'subjects')
                %convert the virtual subject into a real subject
                theGAV=combineLevel-numRsubs;
                %real subjects need to all be listed prior to the virtual ones
                EPdataOut=ep_reorderData(EPdataOut,'subjects',[1:numRsubs combineLevel setdiff([1:numVsubs],theGAV)+numRsubs]);
                if isempty(EPdataOut)
                    EPdataOut=[];
                    return
                end
                EPdataOut.GAVsubs(theGAV+1,:,:)=[];
                numRsubs=numRsubs+1;
                combineLevel=numRsubs;
                numVsubs=numVsubs-1;
            end
        end
        posSubList=subList(combineWeights>0);
        negSubList=subList(combineWeights<0);

        newCellDataPos=zeros(numChanList,numPointList,numCellList,1,numFacSGLlist,numFreqList,numRelList);
        newCellDataNeg=zeros(numChanList,numPointList,numCellList,1,numFacSGLlist,numFreqList,numRelList);

        if ~isempty(EPdataOut.facVecT)
            newCellDataPos=newCellDataPos(:,1,:,:,:,:,:);
            newCellDataNeg=newCellDataNeg(:,1,:,:,:,:,:);
        end
        if ~isempty(EPdataOut.facVecF)
            newCellDataPos=newCellDataPos(:,:,:,:,:,1,:);
            newCellDataNeg=newCellDataNeg(:,:,:,:,:,1,:);
        end
        if ~isempty(EPdataOut.facVecS)
            newCellDataPos=newCellDataPos(1,:,:,:,:,:,:);
            newCellDataNeg=newCellDataNeg(1,:,:,:,:,:,:);
        end

        if ~isempty(EPdataOut.noise)
            newCellNoisePos=zeros(numChanList,numPointList,numCellList,1,numFacSGLlist,numFreqList);
            newCellNoiseNeg=zeros(numChanList,numPointList,numCellList,1,numFacSGLlist,numFreqList);
        end
        if ~isempty(EPdataOut.covAVE)
            if size(EPdataOut.covAVE,7)==1
                newCovAVE=zeros(numChanList,numPointList,numCellList,1,numFacCMBList,numFreqList,1);
            else
                newCovAVE=zeros(numChanList,numPointList,numCellList,1,numFacCMBList,numFreqList,numChanList);
            end
        else
            newCovAVE=[];
        end

        newBadChanData=zeros(1,numCellList,numChanList);

        newFacDataPos=zeros(numFacChanList,numFacPointList,numCellList,1,numFacCMBList,numFacFreqList,numRelList);
        newFacDataNeg=zeros(numFacChanList,numFacPointList,numCellList,1,numFacCMBList,numFacFreqList,numRelList);

        newCellCov=zeros(1,numChans,numChans);

        totalWeightPos=zeros(numChans,numPointList,numCellList,1,numFacSGLlist,numFreqList,numRelList);
        totalWeightNeg=zeros(numChans,numPointList,numCellList,1,numFacSGLlist,numFreqList,numRelList);

        totalNoiseWeightPos=zeros(numChans,numPointList,numCellList,1,numFacSGLlist,numFreqList,1);
        totalNoiseWeightNeg=zeros(numChans,numPointList,numCellList,1,numFacSGLlist,numFreqList,1);

        totalFacWeightPos=zeros(numFacChanList,numFacPointList,numCellList,1,numFacCMBList,numFacFreqList,1);
        totalFacWeightNeg=zeros(numFacChanList,numFacPointList,numCellList,1,numFacCMBList,numFacFreqList,1);

        for iFac=1:length(facSGLlist)
            theFac=facSGLlist(iFac);
            for iCell=1:length(cellList)
                theCell=cellList(iCell);
                ep_tictoc;if EPtictoc.stop;return;end
                if strcmp(dataDimension,'convert')
                    if isempty(EPdataIn.GAVsubs{combineLevel-numRsubs+1,max(theCell-numRcells+1,1),theFac})
                        subCellList=combineData{4};
                        combineSubCellWeights=combineWeightsIn;
                        if ~isempty(subCellList)
                            subGoodList=subCellList(EPdataOut.avgNum(subCellList,theCell)>=0);
                            posSubCellList=subGoodList(combineSubCellWeights>0);
                            negSubCellList=subGoodList(combineSubCellWeights<0);
                        else
                            subGoodList=[];
                        end
                    else
                        subGoodList=EPdataIn.GAVsubs{combineLevel-numRsubs+1,max(theCell-numRcells+1,1),theFac}(:,1);
                        combineSubCellWeights=EPdataIn.GAVsubs{combineLevel-numRsubs+1,max(theCell-numRcells+1,1),theFac}(:,2);
                        posSubCellList=subGoodList(combineSubCellWeights>0);
                        negSubCellList=subGoodList(combineSubCellWeights<0);
                    end
                    if sum(combineSubCellWeights)==0
                        SubCellContrastMode=1;
                    else
                        SubCellContrastMode=0;
                    end
                    if iFac==1
                        EPdataOut.avgNum(combineLevel,theCell)=0;
                        EPdataOut.subNum(combineLevel,theCell)=0;
                        if ~isempty(EPdataOut.covNum)
                            EPdataOut.covNum(combineLevel,theCell)=NaN;
                        end
                        EPdataOut.analysis.blinkTrial(combineLevel,theCell)=0;
                        EPdataOut.analysis.saccadeTrial(combineLevel,theCell)=0;
                        EPdataOut.analysis.saccadeOnset(combineLevel,theCell)=0;
                        EPdataOut.analysis.moveTrial(combineLevel,theCell)=0;
                        EPdataOut.analysis.badTrials(combineLevel,theCell)=0;
                        if ~isempty(EPdataOut.trialSpecs)
                            for iSpec=1:length(EPdataOut.trialSpecNames)
                                EPdataOut.trialSpecs{theCell,iSpec,combineLevel}='';
                            end
                        end
                    end
                else
                    subGoodList=subList(EPdataOut.avgNum(subList,theCell)>=0);
                    SubCellContrastMode=contrastMode;
                    negSubCellList=negSubList;
                    posSubCellList=posSubList;
                    combineSubCellWeights=combineWeights;
                end

                if ~isempty(subGoodList)
                    if isempty(negCell)
                        theNegCell=theCell;
                    else
                        theNegCell=negCell;
                    end
                    if isempty(posCell)
                        thePosCell=theCell;
                    else
                        thePosCell=posCell;
                    end

                    negSubGoodList=negSubCellList(EPdataOut.avgNum(negSubCellList,theCell)>=0);
                    posSubGoodList=posSubCellList(EPdataOut.avgNum(posSubCellList,theCell)>=0);
                    posGoodWeights=combineSubCellWeights(ismember(subGoodList,posSubGoodList));
                    negGoodWeights=combineSubCellWeights(ismember(subGoodList,negSubGoodList));
                    posGoodWeights=posGoodWeights./sum(posGoodWeights);
                    negGoodWeights=negGoodWeights./sum(negGoodWeights);

                    % if ~strcmp(dataDimension,'convert')
                        %for virtual grand averages, if from trimmed ANOVAs, lists of subjects are different for each factor so these specs are meaningless.

                        %just sum up size of sample involved in this combined waveform.
                        EPdataOut.avgNum(combineLevel,theCell)=sum(EPdataOut.avgNum(posSubGoodList,thePosCell),'omitnan')+sum(EPdataOut.avgNum(negSubGoodList,theNegCell),'omitnan');
                        EPdataOut.subNum(combineLevel,theCell)=sum(EPdataOut.subNum(posSubGoodList,thePosCell),'omitnan')+sum(EPdataOut.subNum(negSubGoodList,theNegCell),'omitnan');

                        if ~isempty(EPdataOut.covNum)
                            %assume cov matrix is different so need to figure out new covariance matrix as well as the new effective sample size for the new combination
                            %per p.128 of the 3.7.2 MNE manual, 1/Leff=Sigma weight-squared/L
                            %for noise, which is what this is for correcting, subtraction and addition have the same effect.
                            EPdataOut.covNum(combineLevel,theCell)=(sum([EPdataOut.covNum(posSubGoodList,thePosCell).^-1].*posGoodWeights.^2)+sum([EPdataOut.covNum(negSubGoodList,theNegCell).^-1].*negGoodWeights.^2)).^-1;
                        end

                        %for these, weighted mean of the involved data and saccadeOnset is the difference in the onsets for contrasts.
                        EPdataOut.analysis.blinkTrial(combineLevel,theCell)= sum(EPdataOut.analysis.blinkTrial(posSubGoodList,thePosCell).*posGoodWeights,'omitnan')+sum(EPdataOut.analysis.blinkTrial(negSubGoodList,theNegCell).*negGoodWeights,'omitnan');
                        EPdataOut.analysis.saccadeTrial(combineLevel,theCell)=sum(EPdataOut.analysis.saccadeTrial(posSubGoodList,thePosCell).*posGoodWeights,'omitnan')+sum(EPdataOut.analysis.saccadeTrial(negSubGoodList,theNegCell).*negGoodWeights,'omitnan');
                        EPdataOut.analysis.saccadeOnset(combineLevel,theCell)=sum(EPdataOut.analysis.saccadeOnset(posSubGoodList,thePosCell).*posGoodWeights,'omitnan')-sum(EPdataOut.analysis.saccadeOnset(negSubGoodList,theNegCell).*negGoodWeights,'omitnan');
                        EPdataOut.analysis.moveTrial(combineLevel,theCell)=sum(EPdataOut.analysis.moveTrial(posSubGoodList,thePosCell).*posGoodWeights,'omitnan')+sum(EPdataOut.analysis.moveTrial(negSubGoodList,theNegCell).*negGoodWeights,'omitnan');
                        EPdataOut.analysis.badTrials(combineLevel,theCell)=sum(EPdataOut.analysis.badTrials(posSubGoodList,thePosCell).*posGoodWeights,'omitnan')+sum(EPdataOut.analysis.badTrials(negSubGoodList,theNegCell).*negGoodWeights,'omitnan');

                        %for trial specs, for numerical data weighted differences for contrasts.
                        if ~isempty(EPdataOut.trialSpecs)
                            for iSpec=1:length(EPdataOut.trialSpecNames)
                                posSpecs=squeeze(EPdataOut.trialSpecs(thePosCell,iSpec,posSubGoodList));
                                negSpecs=squeeze(EPdataOut.trialSpecs(theNegCell,iSpec,negSubGoodList));
                                allSpecs=[posSpecs; negSpecs];
                                if ischar(allSpecs{1}) && all(strcmp(allSpecs{1},allSpecs)) %if the spec is all just the same character string then just set it to that string.
                                    EPdataOut.trialSpecs{theCell,iSpec,combineLevel}=allSpecs{1};
                                else
                                    %convert char to numeric.  char strings that are not numbers will convert to NaN and then be ignored during the computations.
                                    theStrSpecs=find(cellfun(@ischar,posSpecs));
                                    for iStrSpec=1:length(theStrSpecs)
                                        theStrSpec=theStrSpecs(iStrSpec);
                                        posSpecs{theStrSpec}=str2double(posSpecs{theStrSpec});
                                        if ~isnumeric(posSpecs{theStrSpec})
                                            posSpecs{theStrSpec}=NaN; %dealing with weird Matlab bug where 'cab' was being converted into a driver handle
                                        end
                                    end
                                    theStrSpecs=find(cellfun(@ischar,negSpecs));
                                    for iStrSpec=1:length(theStrSpecs)
                                        theStrSpec=theStrSpecs(iStrSpec);
                                        negSpecs{theStrSpec}=str2double(negSpecs{theStrSpec});
                                        if ~isnumeric(negSpecs{theStrSpec})
                                            negSpecs{theStrSpec}=NaN; %dealing with weird Matlab bug where 'cab' was being converted into a driver handle
                                        end
                                    end

                                    %reweight for cases where there are missing values
                                    posGoodWeightsSpec=combineSubCellWeights(ismember(subGoodList,posSubGoodList));
                                    posGoodWeightsSpec=posGoodWeightsSpec(~cellfun(@isempty,posSpecs));
                                    negGoodWeightsSpec=combineSubCellWeights(ismember(subGoodList,negSubGoodList));
                                    negGoodWeightsSpec=negGoodWeightsSpec(~cellfun(@isempty,negSpecs));
                                    posGoodWeightsSpec=posGoodWeightsSpec./sum(posGoodWeightsSpec);
                                    negGoodWeightsSpec=negGoodWeightsSpec./sum(negGoodWeightsSpec);

                                    if SubCellContrastMode
                                        %zero weights for any datapoint in either positive or negative specs will result in an NaN.  Positive and negative sides given equal weight to each other.
                                        if all(cellfun(@isempty,posSpecs)) || all(cellfun(@isempty,negSpecs))
                                            EPdataOut.trialSpecs{theCell,iSpec,combineLevel}=NaN;
                                        else
                                            if any(strcmp(EPdataOut.trialSpecNames(iSpec),{'RT','ACC'}))
                                                EPdataOut.trialSpecs{theCell,iSpec,combineLevel}=(sum(cell2mat(posSpecs(~cellfun(@isempty,posSpecs))).*posGoodWeightsSpec,'omitnan')+sum(cell2mat(negSpecs(~cellfun(@isempty,negSpecs))).*negGoodWeightsSpec,'omitnan'))/2;
                                            else
                                                EPdataOut.trialSpecs{theCell,iSpec,combineLevel}=sum(cell2mat(posSpecs(~cellfun(@isempty,posSpecs))).*posGoodWeightsSpec,'omitnan')-sum(cell2mat(negSpecs(~cellfun(@isempty,negSpecs))).*negGoodWeightsSpec,'omitnan');
                                            end
                                        end
                                    else
                                        %zero weights for any datapoint will result in an NaN.
                                        if sum(combineSubCellWeights)>0
                                            if all(cellfun(@isempty,posSpecs))
                                                EPdataOut.trialSpecs{theCell,iSpec,combineLevel}=NaN;
                                            else
                                                EPdataOut.trialSpecs{theCell,iSpec,combineLevel}=sum(cell2mat(posSpecs(~cellfun(@isempty,posSpecs))).*posGoodWeightsSpec,'omitnan');
                                            end
                                        elseif sum(combineSubCellWeights)<0
                                            if all(cellfun(@isempty,negSpecs))
                                                EPdataOut.trialSpecs{theCell,iSpec,combineLevel}=NaN;
                                            else
                                                EPdataOut.trialSpecs{theCell,iSpec,combineLevel}=sum(cell2mat(negSpecs(~cellfun(@isempty,negSpecs)))*negGoodWeightsSpec,'omitnan');
                                            end
                                        end
                                    end
                                end
                            end
                        else
                            EPdataOut.trialSpecs=cell(numCells,0,numSubs);
                        end
                    % end

                    %the voltage data
                    for iChan=1:length(chanList)
                        theChan=chanList(iChan);
                        if strcmp(EPdataOut.dataType,'average')
                            subChanGoodList=union(posSubGoodList(~isnan(EPdataOut.analysis.badChans(posSubGoodList,thePosCell,theChan))),negSubGoodList(~isnan(EPdataOut.analysis.badChans(negSubGoodList,theNegCell,theChan))));
                        else
                            subChanGoodList=union(posSubGoodList((EPdataOut.analysis.badChans(posSubGoodList,thePosCell,theChan) >= 0)),negSubGoodList((EPdataOut.analysis.badChans(negSubGoodList,theNegCell,theChan) >= 0)));
                        end
                        subChanGoodList=intersect(subChanGoodList,subGoodList);
                        subChanGoodWeights=combineSubCellWeights(ismember(subGoodList,subChanGoodList));
                        if ~isempty(subChanGoodList)
                            %compute weighted mean average and so forth
                            for iSub=1:length(subChanGoodList)
                                theSub=subChanGoodList(iSub);
                                theWeight=subChanGoodWeights(iSub);
                                if (theWeight<0) && ~isempty(negCell)
                                    theSubCell=negCell;
                                else
                                    theSubCell=theCell;
                                end
                                if (iChan==1) || isempty(EPdataOut.facVecS)
                                    if theSubCell>numRcells
                                        %for now, just assuming all cells averaged together.
                                        theData=mean(EPdataOut.data(theChan,pointList,EPdataIn.GAVsubs{1,theSubCell-numRcells+1,1}(:,1),theSub,theFac,freqList,:),3);
                                    else
                                        theData=EPdataOut.data(theChan,pointList,theSubCell,theSub,theFac,freqList,:);
                                    end
                                    if ~isempty(EPdataOut.freqNames) && any(strcmp(EPdataOut.chanTypes{theChan},{'EEG','REG'})) && isempty(EPdataOut.relNames)
                                        theData=abs(theData); %convert to amplitude scaling when adding freq data together.
                                    end
                                    goodData=~isnan(theData);
                                    theData(~goodData)=0;
                                    if theWeight>0
                                        newCellDataPos(iChan,:,iCell,1,theFac,:,:)=newCellDataPos(iChan,:,iCell,1,theFac,:,:)+(theWeight*theData);
                                        totalWeightPos(iChan,:,iCell,1,theFac,:,:)=totalWeightPos(iChan,:,iCell,1,theFac,:,:)+(goodData*abs(theWeight));
                                    elseif theWeight<0
                                        newCellDataNeg(iChan,:,iCell,1,theFac,:,:)=newCellDataNeg(iChan,:,iCell,1,theFac,:,:)+(theWeight*theData);
                                        totalWeightNeg(iChan,:,iCell,1,theFac,:,:)=totalWeightNeg(iChan,:,iCell,1,theFac,:,:)+(goodData*abs(theWeight));
                                    end
                                end
                                if ~isempty(EPdataOut.facData) && (iFac==1)
                                    if ~isempty(EPdataOut.facVecS)
                                        theChans=facChanList;
                                    else
                                        theChans=theChan;
                                    end
                                    theData=EPdataOut.facData(theChans,facPointList,theSubCell,theSub,facCMBlist,facFreqList,:);
                                    goodData=~isnan(theData);
                                    theData(~goodData)=0;
                                    if ~isempty(EPdataOut.freqNames) && any(strcmp(EPdataOut.chanTypes{theChan},{'EEG','REG'})) && isempty(EPdataOut.relNames)
                                        theData=abs(theData); %convert to amplitude scaling when adding freq data together.
                                    end
                                    if theWeight>0
                                        newFacDataPos(theChans,:,iCell,1,:,:,:)=newFacDataPos(theChans,:,iCell,1,:,:,:)+(theWeight*theData);
                                        totalFacWeightPos(theChans,:,iCell,1,:,:,:)=totalFacWeightPos(theChans,:,iCell,1,:,:,:)+(goodData*abs(theWeight));
                                    elseif theWeight<0
                                        newFacDataNeg(theChans,:,iCell,1,:,:,:)=newFacDataNeg(theChans,:,iCell,1,:,:,:)+(theWeight*theData);
                                        totalFacWeightNeg(theChans,:,iCell,1,:,:,:)=totalFacWeightNeg(theChans,:,iCell,1,:,:,:)+(goodData*abs(theWeight));
                                    end
                                end
                                if ~isempty(EPdataOut.noise)
                                    theNoise=EPdataOut.noise(theChan,pointList,theSubCell,theSub,theFac,freqList);
                                    goodData=~isnan(theData);
                                    theData(~goodData)=0;
                                    if theWeight>0
                                        newCellNoisePos(iChan,:,iCell,1,:,:)=newCellNoisePos(iChan,:,iCell,1,theFac,:)+(theWeight*theNoise);
                                        totalNoiseWeightPos(iChan,:,iCell,1,:,:)=totalNoiseWeightPos(iChan,:,iCell,1,theFac,:)+(goodData*abs(theWeight));
                                    elseif theWeight<0
                                        newCellNoiseNeg(iChan,:,iCell,1,:,:)=newCellNoiseNeg(iChan,:,iCell,1,theFac,:)+(theWeight*theNoise);
                                        totalNoiseWeightNeg(iChan,:,iCell,1,:,:)=totalNoiseWeightNeg(iChan,:,iCell,1,theFac,:)+(goodData*abs(theWeight));
                                    end
                                end
                                if ~isempty(EPdataOut.cov)
                                    newCellCov(1,iChan,:)=newCellCov(1,iChan,:)+(EPdataOut.cov.covMatrix(theSub,theChan,:)*EPdataOut.cov.Nq(theSub)); %inefficient but whatever.  Not weighted as assumed equally good estimates of noise, except for sample size differences.
                                end

                                theBadData=EPdataOut.analysis.badChans(theSub,theSubCell,theChan);
                                if strcmp(EPdataOut.dataType,'average')
                                    if ~xor(sign(theBadData),sign(newBadChanData(1,iCell,iChan)))
                                        newBadChanData(1,iCell,iChan)= newBadChanData(1,iCell,iChan)+theBadData;
                                    elseif isnan(newBadChanData(1,iCell,iChan)) || isnan(theBadData)
                                        newBadChanData(1,iCell,iChan)=NaN;
                                    elseif ~sign(theBadData)
                                        newBadChanData(1,iCell,iChan)=theBadData;
                                    else
                                        newBadChanData(1,iCell,iChan)=newBadChanData(1,iCell,iChan);
                                    end
                                else
                                    if ~sign(theBadData) || ~sign(newBadChanData(1,iCell,iChan))
                                        newBadChanData(1,iCell,iChan)=-1;
                                    else
                                        newBadChanData(1,iCell,iChan)= 0;
                                    end
                                end
                            end
                        else %none of the subjects have good data for this channel
                            if strcmp(EPdataOut.dataType,'average')
                                newBadChanData(1,iCell,iChan)=NaN;
                            else
                                newBadChanData(1,iCell,iChan)=-1;
                            end
                        end
                    end

                    if ~isempty(EPdataOut.covAVE)
                        %computing covAVE based on only subjects that have good channels for each pairs of channels.  To avoid degeneracy, first the variances are computed.
                        %then the covariances are computed and then rescaled to these variances.  This is equivalent to imputing the missing data based on the intact groups.
                        for iPoint=1:numPointList
                            ep_tictoc;if EPtictoc.stop;return;end
                            thePoint=pointList(iPoint);
                            for iFreq=1:numFreqList
                                theFreq=freqList(iFreq);
                                %compute channel variances
                                for iChan=1:length(chanList)
                                    theChan=chanList(iChan);
                                    if strcmp(EPdataOut.dataType,'average')
                                        negSubChanGoodList=negSubGoodList(~isnan(EPdataOut.analysis.badChans(negSubGoodList,theNegCell,theChan)));
                                        posSubChanGoodList=posSubGoodList(~isnan(EPdataOut.analysis.badChans(posSubGoodList,thePosCell,theChan)));
                                    else
                                        negSubChanGoodList=subGoodList((EPdataOut.analysis.badChans(negSubGoodList,theNegCell,theChan) >= 0));
                                        posSubChanGoodList=subGoodList((EPdataOut.analysis.badChans(posSubGoodList,thePosCell,theChan) >= 0));
                                    end
                                    if size(EPdataOut.covAVE,7)==1
                                        theCovData=squeeze(EPdataOut.covAVE(theChan,thePoint,theNegCell,negSubChanGoodList,theFac,theFreq,1));
                                        theCovData=[theCovData; squeeze(EPdataOut.covAVE(theChan,thePoint,thePosCell,posSubChanGoodList,theFac,theFreq,1))];
                                    else
                                        theCovData=squeeze(EPdataOut.covAVE(theChan,thePoint,theNegCell,negSubChanGoodList,theFac,theFreq,theChan));
                                        theCovData=[theCovData; squeeze(EPdataOut.covAVE(theChan,thePoint,thePosCell,posSubChanGoodList,theFac,theFreq,theChan))];
                                    end
                                    theCovData=shiftdim(theCovData,-2);
                                    subChanGoodWeights=[combineSubCellWeights(ismember(subGoodList,negSubChanGoodList)); combineSubCellWeights(ismember(subGoodList,posSubChanGoodList))];
                                    if size(EPdataOut.covAVE,7)==1
                                        newCovAVE(theChan,thePoint,theCell,1,theFac,theFreq,1)=ep_covMat(theCovData, [], 'poolObs', subChanGoodWeights, [], []);
                                    else
                                        newCovAVE(theChan,thePoint,theCell,1,theFac,theFreq,theChan)=ep_covMat(theCovData, [], 'poolObs', subChanGoodWeights, [], []);
                                    end
                                end
                                if size(newCovAVE,7)>1
                                    %compute channel covariances
                                    for iChan=1:length(chanList)
                                        theChan=chanList(iChan);
                                        if strcmp(EPdataOut.dataType,'average')
                                            negSubChanGoodList=negSubGoodList(~isnan(EPdataOut.analysis.badChans(negSubGoodList,theNegCell,theChan)));
                                            posSubChanGoodList=posSubGoodList(~isnan(EPdataOut.analysis.badChans(posSubGoodList,thePosCell,theChan)));
                                        else
                                            negSubChanGoodList=subGoodList((EPdataOut.analysis.badChans(negSubGoodList,theNegCell,theChan) >= 0));
                                            posSubChanGoodList=subGoodList((EPdataOut.analysis.badChans(posSubGoodList,thePosCell,theChan) >= 0));
                                        end
                                        subChanGoodList=union(negSubChanGoodList,posSubChanGoodList);
                                        for iChan2=iChan+1:length(chanList)
                                            theChan2=chanList(iChan2);
                                            if strcmp(EPdataOut.dataType,'average')
                                                negSubChanGoodList2=negSubGoodList(~isnan(EPdataOut.analysis.badChans(negSubGoodList,theNegCell,theChan)));
                                                posSubChanGoodList2=posSubGoodList(~isnan(EPdataOut.analysis.badChans(posSubGoodList,thePosCell,theChan)));
                                            else
                                                negSubChanGoodList2=subGoodList((EPdataOut.analysis.badChans(negSubGoodList,theNegCell,theChan) >= 0));
                                                posSubChanGoodList2=subGoodList((EPdataOut.analysis.badChans(posSubGoodList,thePosCell,theChan) >= 0));
                                            end
                                            bothChanGoodWeights=[combineSubCellWeights(ismember(subGoodList,negSubChanGoodList2)); combineSubCellWeights(ismember(subCellList,posSubChanGoodList2))];
                                            theCovData=EPdataOut.covAVE(theChan,thePoint,theNegCell,negSubChanGoodList2,theFac,theFreq,theChan);
                                            theCovData=[theCovData; EPdataOut.covAVE(theChan,thePoint,thePosCell,posSubChanGoodList2,theFac,theFreq,theChan)];
                                            covData=ep_covMat(theCovData, [], 'poolObs', bothChanGoodWeights, [], []);
                                            theCor=covData(1,2)/(sqrt(covData(1,1))*sqrt(covData(2,2)));
                                            theCov=theCor*sqrt(newCovAVE(theChan,thePoint,theCell,1,theFac,theFreq,theChan))*sqrt(newCovAVE(theChan2,thePoint,theCell,1,theFac,theFreq,theChan2));
                                            newCovAVE(theChan,thePoint,theCell,1,theFac,theFreq,theChan2)=theCov;
                                            newCovAVE(theChan2,thePoint,theCell,1,theFac,theFreq,theChan)=theCov;
                                        end
                                    end
                                end
                            end
                        end
                    end
                else
                    EPdataOut.avgNum(combineLevel,theCell)=-1;
                    EPdataOut.covNum(combineLevel,theCell)=-1;
                    EPdataOut.subNum(combineLevel,theCell)=-1;
                    EPdataOut.analysis.blinkTrial(combineLevel,theCell)= 0;
                    EPdataOut.analysis.saccadeTrial(combineLevel,theCell)= 0;
                    EPdataOut.analysis.saccadeOnset(combineLevel,theCell)= 0;
                    EPdataOut.analysis.moveTrial(combineLevel,theCell)= 0;
                    EPdataOut.analysis.badTrials(combineLevel,theCell)= 0;
                    newBadChanData(1,iCell,:)= 0;
                end
            end
        end

        %final voltage computation
        if contrastMode
            %zero weights for any datapoint in either positive or negative matrix will result in an NaN.
            EPdataOut.data(chanList,pointList,cellList,combineLevel,facSGLlist,freqList,relList)=newCellDataPos./totalWeightPos+newCellDataNeg./totalWeightNeg;
        else
            %zero weights for any datapoint will result in an NaN.
            if sum(combineWeights)>0
                EPdataOut.data(chanList,pointList,cellList,combineLevel,facSGLlist,freqList,relList)=newCellDataPos./totalWeightPos;
            elseif sum(combineWeights)<0
                EPdataOut.data(chanList,pointList,cellList,combineLevel,facSGLlist,freqList,relList)=newCellDataNeg./totalWeightNeg;
            end
        end

        %subject specs (outside of cell loop as good subject list does not differ per cell)
        if ~isempty(EPdataOut.subjectSpecs)
            for iSpec=1:length(EPdataOut.subjectSpecNames)
                if any(strcmp(EPdataOut.subjectSpecNames{iSpec},{'Gender','Sex'}))
                    tempSpecs=EPdataOut.subjectSpecs(subList,iSpec);
                    specList=unique(tempSpecs(~cell2mat(cellfun(@isempty,tempSpecs,'UniformOutput',false))));
                    summaryString='';
                    for iSpecValue=1:length(specList)
                        summaryString=[summaryString num2str(length(find(strcmp(specList{iSpecValue},EPdataOut.subjectSpecs(subList,iSpec))))) specList{iSpecValue} ' '];
                    end
                    if ~isempty(summaryString)
                        EPdataOut.subjectSpecs{combineLevel,iSpec}=summaryString(1:end-1);
                    end
                else
                    posSpecs=squeeze(EPdataOut.subjectSpecs(posSubList,iSpec));
                    negSpecs=squeeze(EPdataOut.subjectSpecs(negSubList,iSpec));
                    allSpecs=[posSpecs; negSpecs];
                    if ischar(allSpecs{1}) && (((size(allSpecs{1},1)==1) && all(strcmp(allSpecs{1},allSpecs))) || ((size(allSpecs{1},1)==1) && all(strcmp(allSpecs(1),allSpecs))))
                        
                        
                        %if the spec is all just the same character string then just set it to that string.
                        EPdataOut.subjectSpecs{combineLevel,iSpec}=allSpecs{1};
                    else
                        %convert char to numeric.  char strings that are not numbers will convert to NaN and then be ignored during the computations.
                        theStrSpecs=find(cellfun(@ischar,posSpecs));
                        for iStrSpec=1:length(theStrSpecs)
                            theStrSpec=theStrSpecs(iStrSpec);
                            posSpecs{theStrSpec}=str2double(posSpecs{theStrSpec});
                            if ~isnumeric(posSpecs{theStrSpec})
                                posSpecs{theStrSpec}=NaN; %dealing with weird Matlab bug where 'cab' was being converted into a driver handle
                            end
                        end
                        theStrSpecs=find(cellfun(@ischar,negSpecs));
                        for iStrSpec=1:length(theStrSpecs)
                            theStrSpec=theStrSpecs(iStrSpec);
                            negSpecs{theStrSpec}=str2double(negSpecs{theStrSpec});
                            if ~isnumeric(negSpecs{theStrSpec})
                                negSpecs{theStrSpec}=NaN; %dealing with weird Matlab bug where 'cab' was being converted into a driver handle
                            end
                        end


                        if (isempty(posSpecs) || all(cell2mat(cellfun(@isnan,posSpecs,'UniformOutput',false)))) && (isempty(negSpecs) || all(cell2mat(cellfun(@isnan,negSpecs,'UniformOutput',false))))
                            EPdataOut.subjectSpecs{combineLevel,iSpec}='';
                        else
                            %reweight for cases where there are missing values
                            posGoodWeightsSpec=combineWeights(ismember(subList,posSubList));
                            posGoodWeightsSpec=posGoodWeightsSpec(~cellfun(@isempty,posSpecs));
                            negGoodWeightsSpec=combineWeights(ismember(subList,negSubList));
                            negGoodWeightsSpec=negGoodWeightsSpec(~cellfun(@isempty,negSpecs));
                            posGoodWeightsSpec=posGoodWeightsSpec./sum(posGoodWeightsSpec);
                            negGoodWeightsSpec=negGoodWeightsSpec./sum(negGoodWeightsSpec);

                            if contrastMode
                                %zero weights for any datapoint in either positive or negative specs will result in an NaN.  Positive and negative sides given equal weight to each other.
                                if all(cellfun(@isempty,posSpecs)) || all(cellfun(@isempty,negSpecs))
                                    EPdataOut.subjectSpecs{combineLevel,iSpec}=NaN;
                                else
                                    EPdataOut.subjectSpecs{combineLevel,iSpec}=sum(cell2mat(posSpecs(~cellfun(@isempty,posSpecs))).*posGoodWeightsSpec,'omitnan')-sum(cell2mat(negSpecs(~cellfun(@isempty,negSpecs))).*negGoodWeightsSpec,'omitnan');
                                end
                            else
                                %zero weights for any datapoint will result in an NaN.
                                if sum(combineWeights)>0
                                    if all(cellfun(@isempty,posSpecs))
                                        EPdataOut.subjectSpecs{combineLevel,iSpec}=NaN;
                                    else
                                        EPdataOut.subjectSpecs{combineLevel,iSpec}=sum(cell2mat(posSpecs(~cellfun(@isempty,posSpecs))).*posGoodWeightsSpec,'omitnan');
                                    end
                                elseif sum(combineWeights)<0
                                    if all(cellfun(@isempty,negSpecs))
                                        EPdataOut.subjectSpecs{combineLevel,iSpec}=NaN;
                                    else
                                        EPdataOut.subjectSpecs{combineLevel,iSpec}=sum(cell2mat(negSpecs(~cellfun(@isempty,negSpecs)))*negGoodWeightsSpec,'omitnan');
                                    end
                                end
                            end
                        end
                    end
                end
            end
        else
            EPdataOut.subjectSpecs=cell(numSubs,0);
        end

        %task specs
        if ~isempty(EPdataOut.taskSpecs)
            for iSpec=1:length(EPdataOut.taskNames)
                for iMeas=1:length(EPdataOut.taskMeasNames)
                    posSpecs=squeeze(EPdataOut.taskSpecs(posSubList,iSpec,iMeas));
                    negSpecs=squeeze(EPdataOut.taskSpecs(negSubList,iSpec,iMeas));

                    %reweight for cases where there are missing values
                    posGoodWeightsSpec=combineWeights(ismember(subList,posSubList));
                    posGoodWeightsSpec=posGoodWeightsSpec(~isnan(posGoodWeightsSpec));
                    negGoodWeightsSpec=combineWeights(ismember(subList,negSubList));
                    negGoodWeightsSpec=negGoodWeightsSpec(~isnan(negGoodWeightsSpec));
                    posGoodWeightsSpec=posGoodWeightsSpec./sum(posGoodWeightsSpec);
                    negGoodWeightsSpec=negGoodWeightsSpec./sum(negGoodWeightsSpec);

                    if contrastMode
                        %zero weights for any datapoint in either positive or negative specs will result in an NaN.  Positive and negative sides given equal weight to each other.
                        if all(isnan(posGoodWeightsSpec)) || all(isnan(negGoodWeightsSpec))
                            EPdataOut.taskSpecs(combineLevel,iSpec,iMeas)=NaN;
                        else
                            EPdataOut.taskSpecs(combineLevel,iSpec,iMeas)=sum(posSpecs.*posGoodWeightsSpec,'omitnan')-sum(negSpecs.*negGoodWeightsSpec,'omitnan');
                        end
                    else
                        %zero weights for any datapoint will result in an NaN.
                        if sum(combineWeights)>0
                            if all(isnan(posSpecs))
                                EPdataOut.taskSpecs(combineLevel,iSpec,iMeas)=NaN;
                            else
                                EPdataOut.taskSpecs(combineLevel,iSpec,iMeas)=sum(posSpecs.*posGoodWeightsSpec,'omitnan');
                            end
                        elseif sum(combineWeights)<0
                            if all(isnan(negSpecs))
                                EPdataOut.taskSpecs(combineLevel,iSpec,iMeas)=NaN;
                            else
                                EPdataOut.taskSpecs(combineLevel,iSpec,iMeas)=sum(negSpecs.*negGoodWeightsSpec,'omitnan');
                            end
                        end
                    end
                end
            end
        else
            EPdataOut.taskSpecs=cell(0);
        end

        %impedances
        if ~isempty(EPdataOut.impedances) && ~isempty(EPdataOut.impedances.channels)
            if ~isempty(EPdataIn.facVecS)
                impChanList=[1:length(EPdataIn.chanNames)];
            else
                impChanList=chanList;
            end
            for iChan=1:length(impChanList)+1
                if iChan==(length(impChanList)+1)
                    inImpedances=EPdataOut.impedances.ground';
                    theChan=1;
                else
                    inImpedances=EPdataOut.impedances.channels;
                    theChan=impChanList(iChan);
                end
                if isempty(inImpedances)
                    outImpedances=[];
                else
                    posSpecs=squeeze(inImpedances(theChan,posSubList))';
                    negSpecs=squeeze(inImpedances(theChan,negSubList))';

                    %reweight for cases where there are missing values
                    posGoodWeightsSpec=combineWeights(ismember(subList,posSubList));
                    posGoodWeightsSpec=posGoodWeightsSpec(~isnan(posGoodWeightsSpec));
                    negGoodWeightsSpec=combineWeights(ismember(subList,negSubList));
                    negGoodWeightsSpec=negGoodWeightsSpec(~isnan(negGoodWeightsSpec));
                    posGoodWeightsSpec=posGoodWeightsSpec./sum(posGoodWeightsSpec);
                    negGoodWeightsSpec=negGoodWeightsSpec./sum(negGoodWeightsSpec);

                    if contrastMode
                        %zero weights for any datapoint in either positive or negative specs will result in an NaN.  Positive and negative sides given equal weight to each other.
                        if all(isnan(posGoodWeightsSpec)) || all(isnan(negGoodWeightsSpec))
                            outImpedances=NaN;
                        else
                            outImpedances=sum(posSpecs.*posGoodWeightsSpec,'omitnan')-sum(negSpecs.*negGoodWeightsSpec,'omitnan');
                        end
                    else
                        %zero weights for any datapoint will result in an NaN.
                        if sum(combineWeights)>0
                            if all(isnan(posSpecs))
                                outImpedances=NaN;
                            else
                                outImpedances=sum(posSpecs.*posGoodWeightsSpec,'omitnan');
                            end
                        elseif sum(combineWeights)<0
                            if all(isnan(negSpecs))
                                outImpedances=NaN;
                            else
                                outImpedances=sum(negSpecs.*negGoodWeightsSpec,'omitnan');
                            end
                        end
                    end
                end
                if iChan==(length(impChanList)+1)
                    if ~isempty(EPdataOut.impedances.ground)
                        EPdataOut.impedances.ground(combineLevel)=outImpedances;
                    end
                else
                    EPdataOut.impedances.channels(iChan,combineLevel)=outImpedances;
                end
            end
        else
            EPdataOut.impedances(1)=struct('channels',[],'ground',[]);
        end

        %events
        EPdataOut.events{combineLevel,theCell}=struct('type',{},'sample',{},'value',{},'duration',{},'keys',struct('code','','data','','datatype','','description',''));
        % EPdataOut.events{combineLevel,theCell}(1)=[];
        for iSub=1:length(subList)
            theSub=subList(iSub);
            for iCell=1:length(cellList)
                theCell=cellList(iCell);
                for iEvent=1:length(EPdataOut.events{theSub,theCell})
                    EPdataOut.events{combineLevel,theCell}(end+1)=EPdataOut.events{theSub,theCell}(iEvent);
                end
            end
        end

        EPdataOut.analysis.badChans(combineLevel,cellList,chanList)= newBadChanData;

        if contrastMode
            %zero weights for any datapoint in either positive or negative matrix will result in an NaN.
            EPdataOut.facData(facChanList,facPointList,cellList,combineLevel,facCMBlist,facFreqList,relList)=newFacDataPos./totalFacWeightPos+newFacDataNeg./totalFacWeightNeg;
        else
            %zero weights for any datapoint will result in an NaN.
            if sum(combineWeights)>0
                EPdataOut.facData(facChanList,facPointList,cellList,combineLevel,facCMBlist,facFreqList,relList)=newFacDataPos./totalFacWeightPos;
            elseif sum(combineWeights)<0
                EPdataOut.facData(facChanList,facPointList,cellList,combineLevel,facCMBlist,facFreqList,relList)=newFacDataNeg./totalFacWeightNeg;
            end
        end

        if ~isempty(EPdataOut.noise)
            if contrastMode
                %zero weights for any datapoint in either positive or negative matrix will result in an NaN.
                EPdataOut.noise(chanList,pointList,cellList,combineLevel,facSGLlist,freqList)=newCellNoisePos./totalNoiseWeightPos+newCellNoiseNeg./totalNoiseWeightNeg;
            else
                %zero weights for any datapoint will result in an NaN.
                if sum(combineWeights)>0
                    EPdataOut.noise(chanList,pointList,cellList,combineLevel,facSGLlist,freqList)=newCellNoisePos./totalNoiseWeightPos;
                elseif sum(combineWeights)<0
                    EPdataOut.noise(chanList,pointList,cellList,combineLevel,facSGLlist,freqList)=newCellNoiseNeg./totalNoiseWeightNeg;
                end
            end
        end
        if ~isempty(EPdataOut.cov)
            %add together cov matrices weighted by their sample size, per MNE manual 2.7.3 p.90.
            Nq=sum(EPdataOut.cov.Nq(subList)); %computed based on all involved subjects without regard to whether they might be bad in some or all of the cells.
            EPdataOut.cov.covMatrix(combineLevel,chanList,chanList)=newCellCov/Nq;
            EPdataOut.cov.Nq(combineLevel)=Nq;
        end
        if ~isempty(EPdataOut.covAVE)
            if size(EPdataOut.covAVE,7)==1
                if ~isempty(newCovAVE)
                    EPdataOut.covAVE(chanList,pointList,cellList,combineLevel,facSGLlist,freqList,1)=newCovAVE;
                else
                    EPdataOut.covAVE(chanList,pointList,cellList,combineLevel,facSGLlist,freqList,1)=NaN;
                end
            else
                if ~isempty(newCovAVE)
                    EPdataOut.covAVE(chanList,pointList,cellList,combineLevel,facSGLlist,freqList,chanList)=newCovAVE;
                else
                    EPdataOut.covAVE(chanList,pointList,cellList,combineLevel,facSGLlist,freqList,chanList)=NaN;
                end
            end
        end

        if ~isempty(EPdataOut.sessNums)
            if ~any(diff(EPdataOut.sessNums(subList)))
                EPdataOut.sessNums(combineLevel)=EPdataOut.sessNums(subList(1));
            else
                EPdataOut.sessNums(combineLevel)=0;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(dataDimension,'factors')
    %combination of factors is straight addition rather than averages, unlike for the other dimensions.
    combineLevel=find(strcmp(combineName,EPdataOut.facNames));
    newFacVar=0;
    newFacVarQ=0;
    if isempty(combineLevel)
        numFacs=numFacs+1;
        numCMBfacs=numCMBfacs+1;
        combineLevel=numCMBfacs;
        EPdataOut.facNames{end+1,1}=combineName;
        EPdataOut.facTypes{end+1,1}='CMB';
    end

    for iFac=1:length(facList)
        theFac=facList(iFac);
        theWeight=combineWeights(iFac);
        newFacVar=newFacVar+(theWeight)*EPdataOut.facVar(theFac);
        newFacVarQ=newFacVarQ+(theWeight)*EPdataOut.facVarQ(theFac);
    end
    EPdataOut.facVar(end+1)=newFacVar;
    EPdataOut.facVarQ(end+1)=newFacVarQ;

    newFacData=zeros(numFacChanList,numFacPointList,numRcells,numRsubs,1,numFacFreqList,numRelList);

    for iFac=1:length(facList)
        theFac=facList(iFac);
        ep_tictoc;if EPtictoc.stop;return;end
        theWeight=combineWeights(iFac);
        theData=ep_expandFacs(EPdataOut,facChanList,facPointList,cellList(cellList<=numRcells),subList(subList<=numRsubs),theFac,facFreqList,relList);
        if isempty(theData)
            EPdataOut=[];
            return
        end
        ep_tictoc;if EPtictoc.stop;return;end
        newFacData(:,:,:,:,1,:,:)=newFacData(:,:,:,:,1,:,:)+(theWeight*theData);
    end

    EPdataOut.facData(facChanList,facPointList,cellList(cellList<=numRcells),subList(subList<=numRsubs),combineLevel,facFreqList,relList)=newFacData;

    if ~isempty(EPdataOut.GAVsubs)
        EPdataOut.GAVsubs{end,end,end+1}=[];
    end
end

if strcmp(dataDimension,'convert')
    EPdataOut.GAVsubs=[];
end

if ~isempty(EPdataIn.relNames)
    %undo Fisher-Z transform
    EPdataOut.data=tanh(EPdataOut.data);
end

[err]=ep_checkEPfile(EPdataOut);
if err
    EPdataOut=[];
end
EPdataOutFinal=EPdataOut;