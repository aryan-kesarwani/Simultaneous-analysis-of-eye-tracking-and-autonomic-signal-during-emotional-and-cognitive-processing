function handleListOut = ep_disableGUI(figureHandle,handleListIn)
%  ep_disableGUI;
%       Disables or enables GUI controls.  If handleListIn provided, will enable GUI controls to their prior state.
%
%Input:
%    figureHandle   : handle of the GUI window
%    handleListIn
%           .Tag: tags of the GUI controls
%           .Enable : prior state of the GUI controls (enabled on or off)
%
%Outputs:
%    handleListIn
%           .Tag: tags of the GUI controls
%           .Enable : prior state of the GUI controls (enabled on or off)

%History:
%  by Joseph Dien (5/18/24)
%  jdien07@mac.com
%
% modified 12/27/24 JD
% No longer requires figure elements to be unchanged.
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

handleListOut=[];

if exist('handleListIn','var')
    % test=get(EPmain.handles.hMainWindow.Children,'Type');
    % if length(handleListIn) ~= length(find(ismember(get(figureHandle.Children,'Type'),{'uicontrol';'uitable'})))
    %     disp('Error: handle list does not match objects on figure.')
    %     return;
    % end
    endableMode='on';
else
    endableMode='off';
end

switch endableMode
    case 'on'
        for iChild=1:length(figureHandle.Children)
            if ishandle(figureHandle.Children(iChild))
                if ~any(strcmp({'T1';'T2';'T3';'T4';'Reset'},get(figureHandle.Children(iChild),'Tag')))
                    if isprop(figureHandle.Children(iChild),'Enable')
                        set(figureHandle.Children(iChild),'Enable',handleListIn(iChild).Enable);
                    end
                end
            end
        end

    case 'off'
        for iChild=1:length(figureHandle.Children)
            if ishandle(figureHandle.Children(iChild))
                if isprop(figureHandle.Children(iChild),'Enable')
                    handleListOut(iChild).Enable=get(figureHandle.Children(iChild),'Enable');
                else
                    handleListOut(iChild).Enable='on';
                end
                if ~any(strcmp({'T1';'T2';'T3';'T4';'Reset'},get(figureHandle.Children(iChild),'Tag')))
                    if isprop(figureHandle.Children(iChild),'Enable')
                        set(figureHandle.Children(iChild),'Enable','off');
                    end
                end
            end
        end
    otherwise
        disp('Programmer Error: Enable mode not recognized.')
end