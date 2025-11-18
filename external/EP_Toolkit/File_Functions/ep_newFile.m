function EPdata=ep_newFile(templateFile,numCells,numSubs)
% EPdata=ep_newFile(templateFile,numCells,numSubs)
% Constructs a new EP file with empty fields.
%
%Inputs:
%   templateFile         : Optional template file.
%   numCells             : Optional number of cells.
%   numSubs              : Optional number of subjects.
%
%Outputs:
%    EPdata      : blank EPdata file with empty fields except for date.  See ep_readData header for more information.
%                  When the optional template file input is used, the blank file will have the channel and time information of the template 
%                  plus sampling rate, baseline, and time units.  If cell and subject number are specified, then .data will be initialized to Nan
%                  and dummy values provided for the cell and subject fields.
%
%History
%  by Joseph Dien (12/30/19)
%  jdien07@mac.com
%
% bugfix 4/30/20 JD
% Fixed taskSpecs initialized as a cell array, resulting in crashes elsewhere.
%
% modified 11/4/21 JD
% Added optional template file input.
% Added numCells and numSubs optional inputs.
%
% modified 1/27/25 JD
% Added video field.
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

if (nargin>0) && isstruct(templateFile)
    templateFlag=1;
else
    templateFlag=0;
end

if ~exist('numCells','var')
    numCells=[];
end
if ~exist('numSubs','var')
    numSubs=[];
end

if templateFlag
    numChans=length(templateFile.chanNames);
    numPoints=length(templateFile.timeNames);
else
    numChans=[];
    numPoints=[];
end

if ~isempty(numCells) && ~isempty(numSubs)
    EPdata.data=nan(numChans,numPoints,numCells,numSubs);
else
    EPdata.data=[];
end
EPdata.noise=[];
EPdata.covAVE=[];
EPdata.GAVsubs=[];
EPdata.cov=[];
if templateFlag
    EPdata.montage=templateFile.montage;
else
    EPdata.montage='';
end
EPdata.fileFormat='';
if templateFlag
    EPdata.chanNames=templateFile.chanNames;
    EPdata.timeNames=templateFile.timeNames;
else
    EPdata.chanNames=cell(0);
    EPdata.timeNames=[];
end

if isempty(numCells) || (numCells==0)
    EPdata.cellNames=cell(0);
    EPdata.cellTypes=cell(0);
else
    for iCell=1:numCells
        EPdata.cellNames{iCell,1}=num2str(iCell);
        EPdata.cellTypes{iCell,1}='SGL';
    end
end

if isempty(numSubs) || (numSubs==0)
    EPdata.subNames=cell(0);
    EPdata.subTypes=cell(0);
else
    for iSub=1:numSubs
        EPdata.subNames{iSub,1}=num2str(iSub);
        EPdata.subTypes{iSub,1}='AVG';
    end
end

EPdata.trialNames=[];
EPdata.facNames=cell(0);
EPdata.freqNames=[];
EPdata.relNames=cell(0);
EPdata.sessNames=cell(0);
if templateFlag
    EPdata.chanTypes=templateFile.chanTypes;
else
    EPdata.chanTypes=cell(0);
end
if templateFlag
    EPdata.timeUnits=templateFile.timeUnits;
else
    EPdata.timeUnits='';
end

EPdata.facTypes=cell(0);
EPdata.sessNums=[];
try
    EPdata.EPver=ver('EP_Toolkit');
catch
    EPdata.EPver='unavailable'; %workaround for bug in earlier version of Matlab
end
EPdata.ver=ver;
EPdata.date=date;
if templateFlag
    EPdata.Fs=templateFile.Fs;
else
    EPdata.Fs=[];
end
if templateFlag
    EPdata.baseline=templateFile.baseline;
else
    EPdata.baseline=[];
end
EPdata.dataName='';
EPdata.ename='';
EPdata.dataType='';
EPdata.trialSpecs=cell(0);
EPdata.trialSpecNames=cell(0);
EPdata.subjectSpecs=cell(0);
EPdata.subjectSpecNames=cell(0);
EPdata.taskSpecs=[];
EPdata.taskNames=cell(0);
EPdata.taskMeasNames=cell(0);
EPdata.events=cell(0);
EPdata.avgNum=[];
EPdata.subNum=[];
EPdata.covNum=[];
EPdata.fileName='';
EPdata.history=cell(0,5);
if templateFlag
    EPdata.ced=templateFile.ced;
else
    EPdata.ced='';
end
if templateFlag
    EPdata.eloc=templateFile.eloc;
else
    EPdata.eloc=cell(0);
end
if templateFlag
    EPdata.implicit=templateFile.implicit;
else
    EPdata.implicit=cell(0);
end
EPdata.facVecT=[];
EPdata.facVecS=[];
EPdata.facVecF=[];
EPdata.facData=[];
EPdata.facVar=[];
EPdata.facVarQ=[];
EPdata.reference(1)=struct('original',[],'current',[],'type','');
EPdata.analysis(1)=struct('blinkTrial',[],'saccadeTrial',[],'saccadeOnset',[],'moveTrial',[],'badTrials',[],'badChans',[]);
EPdata.pca=[];
EPdata.recTime=[];
EPdata.stims=struct('name',{},'image',{},'AOI',{});
EPdata.calibration=[];
EPdata.impedances(1)=struct('channels',[],'ground',[]);
EPdata.video(1)=struct('times',[],'frames',[]);
EPdata.video(1)=[];

if ~isempty(numCells) && ~isempty(numSubs) && ~isempty(numChans)
    EPdata.events=cell(numSubs,numCells);
    EPdata.avgNum=zeros(numSubs,numCells);
    EPdata.covNum=zeros(numSubs,numCells);
    EPdata.subNum=zeros(numSubs,numCells);
    EPdata.recTime=zeros(numCells);
    EPdata.analysis.badChans=zeros(numSubs,numCells);
    EPdata.analysis.moveTrial=zeros(numSubs,numCells);
    EPdata.analysis.blinkTrial=zeros(numSubs,numCells);
    EPdata.analysis.saccadeTrial=zeros(numSubs,numCells);
    EPdata.analysis.saccadeOnset=zeros(numSubs,numCells);
    EPdata.analysis.badTrials=zeros(numSubs,numCells);
    EPdata.analysis.badChans=zeros(numSubs,numCells,numChans);
end
