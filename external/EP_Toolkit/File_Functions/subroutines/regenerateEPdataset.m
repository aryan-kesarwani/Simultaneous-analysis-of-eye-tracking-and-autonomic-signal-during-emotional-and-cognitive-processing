function EPdataset=regenerateEPdataset(theEPdataset) 
%try to regenerate a corrupted EPdataset file.
%
%Inputs:
%  theEPdataset     : The corrupted EPdataset file.
%
%Outputs:
%  EPdataset         : The EPdataset structured array.

%History:
%  by Joseph Dien (4/9/20)
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

disp('EPdataset file corrupted.  Will attempt to regenerate based on contents of the EPwork directory.');
delete(theEPdataset);
[path, name, ext]=fileparts(theEPdataset);
EPdataset=[];
EPdataset.EPwork=path;
EPdataset.dataset=[];

fileList=dir(path);
for iFile=1:length(fileList)
    if ~any(strcmp(fileList(iFile).name,{'.','..','EPprefs.mat'}))
        disp(fileList(iFile).name)
        load([fileList(iFile).folder filesep fileList(iFile).name]);
        EPdata=ep_updateEPfile(EPdata,1);
        if isempty(EPdataset.dataset)
            EPdataset.dataset=ep_addToEPworkCache(EPdata);
        else
            EPdataset.dataset=[EPdataset.dataset ep_addToEPworkCache(EPdata)];
        end
    end
end

save('-mat', theEPdataset, 'EPdataset', '-v7.3');

