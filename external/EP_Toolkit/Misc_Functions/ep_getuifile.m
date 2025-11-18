function [file, path, indx]=ep_getuifile(filter,title,defname)
% [file, path, indx]=ep_getuifile(filter,title,defname)
% workaround for bug on OS X where title is not presented in the requestor window
%
%Input:
%    filter     : The file suffix mask
%    title      : The title to present in the requestor window
%    defname    : The default file name.
%Outputs
%	file: The chosen file
%   path: The path to the chosen file
%   indx: Which of the presented files was chosen.

%History
%  by Joseph Dien (2/2/21)
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

file='';
path='';
indx=[];

if ismac
    d=dialog('Name',title);
    [file, path, indx]=uigetfile(filter,title,defname);
    if ishandle(d)
        close(d);
    end
else
    [file, path, indx]=uigetfile(filter,title,defname);
end