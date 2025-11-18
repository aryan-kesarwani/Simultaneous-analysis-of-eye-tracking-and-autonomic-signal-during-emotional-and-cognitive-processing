function [legendNames, legendNumbers]=ep_waveNames
% ep_waveNames - legendNames=ep_waveNames
% Generates labels for a set of waveforms.
%
%Outputs
%   legendNames: Cell array with the labels of the displayed waveforms.
%   legendNumbers: Array with the corresponding color numbers.

%History
%  by Joseph Dien (12/24/24)
%  jdien07@mac.com
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global EPmain EPdataset

legendNames=cell(0);
legendNames{1}='';
legendNumbers=[];

%set up names for colors
numCols=length(find(EPmain.view.dataset(1:EPmain.numColors) <= length(EPdataset.dataset)));
colData=cell(0);
colCell=cell(0);
colSub=cell(0);
colTrial=cell(0);
colFactor=cell(0);

for iColor=1:EPmain.numColors
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)

        colData{end+1}=EPdataset.dataset(EPmain.view.dataset(iColor)).dataName;

        if EPmain.view.allTrials(iColor)==5
            colCell{end+1}='all';
        elseif EPmain.view.allTrials(iColor)==6
            colCell{end+1}='erpimage';
        else
            colCell(end+1)=EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor));
        end

        if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'average')
            if EPmain.view.allTrials(iColor)==1
                subName='all';
            elseif EPmain.view.allTrials(iColor)==2
                subName='erpimage';
            else
                subName=EPdataset.dataset(EPmain.view.dataset(iColor)).subNames{EPmain.view.subject(iColor)};
            end
        else
            subName=EPdataset.dataset(EPmain.view.dataset(iColor)).subNames{EPmain.view.subject(iColor)};
        end
        colSub{end+1}=subName;

        if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames)
            if ismember(EPmain.view.allTrials(iColor),[1 2])
                colTrial{end+1}='all trials';
            else
                colTrial{end+1}=num2str(EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames(EPmain.view.trial(iColor)));
            end
        elseif strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'average')
            switch EPmain.view.trial(iColor)
                case 1 %nothing
                    colTrial{end+1}='data';
                case 2 %GFP
                    colTrial{end+1}='GFP';
                case 3 %noise
                    colTrial{end+1}='noise';
                case 4 %StDev
                    colTrial{end+1}='StDev';
                case 5 %95CI
                    colTrial{end+1}='95%CI';
            end
        end

        facName='';
        if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)
            if EPmain.view.allTrials(iColor)==3
                facName='all';
            elseif EPmain.view.allTrials(iColor)==4
                facName='erpimage';
            else
                facName=EPdataset.dataset(EPmain.view.dataset(iColor)).facNames{EPmain.view.factor(iColor)};
            end
        end
        colFactor{end+1}=facName;
        legendNumbers(end+1)=iColor;
    end
end
if length(unique(colData))==numCols
    legendNames=colData;
elseif length(unique(colCell))==numCols
    legendNames=colCell;
elseif length(unique(colSub))==numCols
    legendNames=colSub;
elseif (length(unique(colTrial))==numCols) && ~isempty(colTrial)
    legendNames=colTrial;
elseif (length(unique(colFactor))==numCols) && ~isempty(colFactor)
    legendNames=colFactor;
else
    nameCounter=1;
    legendNames=cell(1,numCols);
    for iColor=1:EPmain.numColors
        if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
            if length(unique(colData))~=1
                legendNames{nameCounter}=[colData{nameCounter}];
            end
            if length(unique(colCell))~=1
                if ~isempty(legendNames{nameCounter})
                    legendNames{nameCounter}=[legendNames{nameCounter} '-'];
                end
                legendNames{nameCounter}=[legendNames{nameCounter} colCell{nameCounter}];
            end
            if length(unique(colSub))~=1
                if ~isempty(legendNames{nameCounter})
                    legendNames{nameCounter}=[legendNames{nameCounter} '-'];
                end
                legendNames{nameCounter}=[legendNames{nameCounter} colSub{nameCounter}];
            end
            if ~isempty(colTrial) && length(unique(colTrial))~=1
                if ~isempty(legendNames{nameCounter})
                    legendNames{nameCounter}=[legendNames{nameCounter} '-'];
                end
                legendNames{nameCounter}=[legendNames{nameCounter} colTrial{nameCounter}];
            end
            if ~isempty(colFactor) && length(unique(colFactor))~=1
                if ~isempty(legendNames{nameCounter})
                    legendNames{nameCounter}=[legendNames{nameCounter} '-'];
                end
                legendNames{nameCounter}=[legendNames{nameCounter} colFactor{nameCounter}];
            end
            if isempty(legendNames{nameCounter})
                legendNames{nameCounter}=sprintf('Dataset %d',iColor);
            end
            nameCounter=nameCounter+1;
        end
    end
end