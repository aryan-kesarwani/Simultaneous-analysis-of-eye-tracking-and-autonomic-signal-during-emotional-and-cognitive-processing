function historyList=ep_addHistory(historyList, theRecords, theDescription, figureHandle, thePreferences, EPversion, theFiles)
% historyList=ep_addHistory(historyList, theRecords, theDescription, figureHandle, thePreferences, EPversion, theFiles)
% Adds a new entry to the History array.
%
%Inputs:
%    historyList   : The cell array of the history of the dataset.
%    theRecords
%       .user             : The user making the change.
%       .lab              : The lab of the user making the change.
%       .institution      : The institution of the user making the change.
%       .project          : The project of the dataset.
%       .experiment       : The experiment of the dataset.
%    theDescription       : the description of the change in plain English.
%    figureHandle         : the handle of the figure from which the settings are to be harvested.
%    thePreferences       : the preferences structure.
%    EPversion            : the version of EP Toolkit.
%    theFiles             : cell array of the input files for the function.
%
%Outputs:
%    historyList      : The cell array of the history of the dataset.
%       1) who and when
%           .user             : The user making the change.
%           .lab              : The lab of the user making the change.
%           .institution      : The institution of the user making the change.
%           .project          : The project of the dataset.
%           .experiment       : The experiment of the dataset.
%           .time             : The time and date of the change.
%           .EPversion        : The version of the EP Toolkit.
%       2) description of change.
%       3)  .position         : position of the figure window.
%           .name             : name of the figure window.
%           .children         : the control settings of the figure window.
%       4)  relevant preferences settings.
%       5)  input files if any, as a cell array.
%
%History
%  by Joseph Dien (10/19/21)
%  jdien07@mac.com
%
% modified 8/27/24 JD
% When type of change is same as the last entry and changes entry to indicate how many times it was done, also updates time of the entry.
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

if ~exist('theFiles','var')
    theFiles=cell(0);
end

if ~iscell(theFiles)
    temp{1}=theFiles;
    theFiles=temp;
end

theTime=char(datetime("now"));

if ~isempty(historyList) && (length(theDescription) <= length(historyList{end,2})) && strcmp(theDescription(1:end-1),historyList{end,2}(1:length(theDescription)-1)) && strcmp(historyList{end,1}.time(1:11),theTime(1:11))
    if  strcmp(historyList{end,2}(end-5:end),'times.')
        newDescript=[theDescription(1:end-1) ' ' num2str(str2double(historyList{end,2}(length(theDescription):length(historyList{end,2})-7))+1) ' times.'];
    else
        newDescript=[theDescription(1:end-1) ' 2 times.'];
    end
    historyList{end,2}=newDescript;
    historyList{end,1}.time=theTime;
else
    historyList{end+1,1}=theRecords;
    historyList{end,1}.time=theTime;
    historyList{end,1}.EPversion=EPversion;
    historyList{end,2}=theDescription;
    historyList{end,3}.position=get(figureHandle,'Position');
    historyList{end,3}.name=get(figureHandle,'Name');
    for iChild=1:length(figureHandle.Children)
        if isprop(figureHandle.Children(iChild),'Type')
            historyList{end,3}.children(iChild).Type=figureHandle.Children(iChild).Type;
        end
        if isprop(figureHandle.Children(iChild),'Style')
            historyList{end,3}.children(iChild).Style=figureHandle.Children(iChild).Style;
        end
        if isprop(figureHandle.Children(iChild),'Units')
            historyList{end,3}.children(iChild).Units=figureHandle.Children(iChild).Units;
        end
        if isprop(figureHandle.Children(iChild),'HorizontalAlignment')
            historyList{end,3}.children(iChild).HorizontalAlignment=figureHandle.Children(iChild).HorizontalAlignment;
        end
        if isprop(figureHandle.Children(iChild),'FontSize')
            historyList{end,3}.children(iChild).FontSize=figureHandle.Children(iChild).FontSize;
        end
        if isprop(figureHandle.Children(iChild),'Position')
            historyList{end,3}.children(iChild).Position=figureHandle.Children(iChild).Position;
        end
        if isprop(figureHandle.Children(iChild),'String')
            historyList{end,3}.children(iChild).String=figureHandle.Children(iChild).String;
        end
        if isprop(figureHandle.Children(iChild),'Value')
            historyList{end,3}.children(iChild).Value=figureHandle.Children(iChild).Value;
        end
        if isprop(figureHandle.Children(iChild),'BackgroundColor')
            historyList{end,3}.children(iChild).BackgroundColor=figureHandle.Children(iChild).BackgroundColor;
        end
        if isprop(figureHandle.Children(iChild),'Data')
            historyList{end,3}.children(iChild).Data=figureHandle.Children(iChild).Data;
        end
        if isprop(figureHandle.Children(iChild),'ColumnFormat')
            historyList{end,3}.children(iChild).ColumnFormat=figureHandle.Children(iChild).ColumnFormat;
        end
        if isprop(figureHandle.Children(iChild),'ColumnName')
            historyList{end,3}.children(iChild).ColumnName=figureHandle.Children(iChild).ColumnName;
        end
        if isprop(figureHandle.Children(iChild),'ColumnWidth')
            historyList{end,3}.children(iChild).ColumnWidth=figureHandle.Children(iChild).ColumnWidth;
        end
        if isprop(figureHandle.Children(iChild),'Tag')
            historyList{end,3}.children(iChild).Tag=figureHandle.Children(iChild).Tag;
        end
        if isprop(figureHandle.Children(iChild),'Enable')
            historyList{end,3}.children(iChild).Enable=figureHandle.Children(iChild).Enable;
        end
    end
    historyList{end,4}=thePreferences;
    historyList{end,5}=theFiles;
end