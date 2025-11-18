function outEloc=ep_elocFormat(inEloc)
% ep_cedFormat - outEloc=ep_elocFormat(inEloc)
% Converts eloc structures to standard EP format.  Also provides standard eloc field name list when there are no input arguments.
% Will initialize an empty standard EP format eloc with the keyword 'initialize'.
%
%Input:
%   inEloc : The input eloc structure.
%
%Output:
%  outEloc : The output eloc structure.

%History
%  by Joseph Dien (8/9/21)
%  jdien07@mac.com
%
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

outEloc=[];

if nargin == 0
    %output fieldnames of EP standard eloc format.
    
    %canonical .ced field order from pop_chanedit.m file plus cX cY cZ canonical 10-05 coordinates
    outEloc = { 'labels' 'theta' 'radius' 'X' 'Y' 'Z' 'sph_theta' 'sph_phi' 'sph_radius' 'type' 'cX' 'cY' 'cZ'}';
    
elseif strcmp(inEloc,'initialize')
    %initialize EP standard eloc format.
    
    outEloc=struct('labels',{},'theta',{},'radius',{},'X',{},'Y',{},'Z',{},'sph_theta',{},'sph_phi',{},'sph_radius',{},'type',{},'cX',{},'cY',{},'cZ',{});
    
elseif isstruct(inEloc)
    %reformat the eloc into EP standard eloc format.
    
    outEloc=inEloc;
    listcolformat=ep_elocFormat;
    EPfieldNames=listcolformat;
    dataFieldNames=fieldnames(inEloc);
    if ~isequal(EPfieldNames,dataFieldNames)
        modelEPdata=[];
        for i=1:length(EPfieldNames)
            modelEPdata.(EPfieldNames{i})=[]; %form the model eloc
        end
        if ~isempty(setdiff(dataFieldNames,EPfieldNames))
            for iField = 1:length(dataFieldNames)
                if ~any(strcmp(dataFieldNames{iField},EPfieldNames))
                    outEloc=rmfield(outEloc,dataFieldNames{iField});
                end
            end
        end
        if ~isempty(setdiff(EPfieldNames,dataFieldNames))
            numEloc=length(outEloc);
            for iField = 1:length(EPfieldNames)
                if ~any(strcmp(EPfieldNames{iField},dataFieldNames))
                    
                    eval(['outEloc(1).' EPfieldNames{iField} '=[];']);
                end
            end
            if numEloc==0
                outEloc(1)=[];
            end
        end
        outEloc = orderfields(outEloc, modelEPdata);
    end
    numericFields={'theta';'radius';'X';'Y';'Z';'sph_theta';'sph_phi';'sph_radius';'cX';'cY';'cZ'};
    numericFields(~ismember(numericFields,fieldnames(inEloc)))=[];
    for iChan=1:length(outEloc)
        for iField=1:length(numericFields)
            evalc(['theVal = inEloc(' num2str(iChan) ').' numericFields{iField} ';']);
            if ~isnumeric(theVal)
                evalc(['outEloc(' num2str(iChan) ').' numericFields{iField} '=[];']);
            end
        end
    end
elseif isempty(inEloc)
    %do nothing
else
    disp('Error - did not recognize the input to ep_elocFormat function.')
end

