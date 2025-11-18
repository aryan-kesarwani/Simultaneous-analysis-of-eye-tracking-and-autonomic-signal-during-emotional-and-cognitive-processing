function theOutput=ep_enumVar(stringName,inString)
% theOutput=ep_enumVar(inString)
% Defines enumerated strings and allows for checks of them.
%
%Inputs:
%  stringName       : The name of the enumerated string.
%  inString         : Optional cell string array, signalling that it should be checked for consistency with the cell string array.  
%
%Outputs:
%    theOutput      : either the cell string array corresponding to the stringName or boolean for whether inString is consistent with the allowed values.
%
%History
%  by Joseph Dien (10/13/21)
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch stringName
    case 'postHoc'
        theArray={'tukey-kramer','dunn-sidak','scheffe'};
        
    case 'epsilon'
        theArray={'Greenhouse-Geisser','Huynd-Feldt','Lower Bound'};
        
    otherwise
        theArray=cell(0);
end

if exist('inString','var')
    theOutput=any(strcmp(inString,theArray));
else
    theOutput=theArray;
end