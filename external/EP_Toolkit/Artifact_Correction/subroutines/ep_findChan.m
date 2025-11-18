function [theChan theOrder] = ep_findChan(eloc, implicit, chanNames, ced, badChans, targetChan, montage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [theChan theOrder] = ep_findChan(eloc, implicit, chanNames, ced, badChans, targetChan, montage)
%	Given the electrode locations, figures out which non-bad channel is closest to the desired channel coordinates.
%
%Inputs
%  eloc          : electrode coordinate structure containing the channel names and locations.
%  implicit      : electrode coordinate structure containing the fiducial names and locations.
%  chanNames     : The names of the data channels.
%  ced           : The ced name of the data channels.
%                 Assumes only electrodes present, no fiducials.
%   .X           : x-coordinate  
%   .Y           : y-coordinate 
%   .Z           : z-coordinate 
%  badChans:   list of bad channels.
%  targetChan    : The 10-05 label of the desired channel.
%  montage       : The montage of the eloc channels.
%
%Outputs
%	theChan      : The channel closest to Cz without being a bad channel.  Zero if they are all bad.
%   theOrder     : This channel is which closest to Cz (e.g., 3 means 3rd, and 1st and 2nd were therefore bad channels and so skipped).
%
% History:
%
% by Joseph Dien (9/7/16)
% jdien07@mac.com
%
% modified 12/18/17 JD
% Generalized function to find any given channel by coordinates.  Changed name of function as well.
%
% modified 11/3/18 JD
% Added support for EEGlab derived Cartesian coordinates from standard-10-5-cap385.ced and using EEGlab's more precise transform for cartesian coordinates from .sfp file.
%
% modified 2/14/20 JD
% Uses improved method for co-registering eloc coordinates, based on Oostenveld's Standard-10-5-Cap385.sfp file.
%
% modified 3/18/20 JD
% Added fiducials structure.
%
% modified 8/8/21 JD
% Added support for montage information.
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

theChan=[];
theOrder=[];

%transform the canonical sfp electrode coordinates into the space of the data eloc coordinates
[sfpEloc, ~] = ep_transformEloc([], [], eloc, implicit, [], chanNames, [], montage);
if isempty(sfpEloc)
    return;
end

theEloc=find(strcmpi(targetChan,{sfpEloc.labels}));
numChans=length(eloc);

elecDistances = zeros(numChans,1);
for iChan = 1:numChans
    if ~isempty(eloc(iChan).X) && ~isempty(eloc(iChan).Y) && ~isempty(eloc(iChan).Z) && ~isnan(eloc(iChan).X) && ~isnan(eloc(iChan).Y) && ~isnan(eloc(iChan).Z)
        elecDistances(iChan)=sqrt((eloc(iChan).X-sfpEloc(theEloc).X)^2+(eloc(iChan).Y-sfpEloc(theEloc).Y)^2+(eloc(iChan).Z-sfpEloc(theEloc).Z)^2);
    else
        elecDistances(iChan)=inf;
    end
end
[B IX]=sort(elecDistances);

theOrder=1;
while ismember(IX(theOrder),badChans)
    theOrder=theOrder+1;
end
if theOrder > numChans
    theChan=0;
else
    theChan=IX(theOrder);
end