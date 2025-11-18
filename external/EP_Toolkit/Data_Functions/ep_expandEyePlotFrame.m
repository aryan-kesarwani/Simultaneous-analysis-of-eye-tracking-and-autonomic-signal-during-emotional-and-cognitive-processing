function ep_expandEyePlotFrame(thePoint)
% ep_expandEyePlotFrame - ep_expandEyePlotFrame -
% Renders a frame for a movie file of the topo data.
%
%Inputs
%   thePoint: the index for displayPoints for the frames to be rendered.

%History
%  by Joseph Dien (1/31/16)
%  jdien07@mac.com
%
% modified 12/25/24 JD
% Added popupmenu controls for selecting the content of the topomaps, which can now be any of the displayed datasets.
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

global EPscanCont EPmain EPtictoc

scrsz = EPmain.scrsz;

EPscanCont.handles.waves.frameFigure = figure('Name', 'Scan Plot', 'NumberTitle', 'off', 'Position',[scrsz(3)/2 scrsz(4)/2 800 800], 'Units','normalized');

topoLoc=zeros(EPscanCont.topoNum,4); %coordinates of the topos
topoLoc(1,:)=[.1 .5 .4 .4];
topoLoc(2,:)=[.5 .5 .4 .4];
topoLoc(3,:)=[.1 .1 .4 .4];
topoLoc(4,:)=[.5 .1 .4 .4];

for iTopo=1:EPscanCont.topoNum
    EPscanCont.handles.movie.topo(iTopo) = axes('Units','pixel','position',topoLoc(iTopo,:));
    set(EPscanCont.handles.movie.topo(iTopo),'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);
end









cmap=colormap(jet);

colorData=zeros(length(EPscanCont.displayPoints),3);
if EPscanCont.chanShow~=length(EPscanCont.chanList)
    theData=squeeze(EPscanCont.EPdata(theColor).data(EPscanCont.chanShow,EPscanCont.displayPoints,:,:,:,EPscanCont.freqPos))';
    if strcmp(EPscanCont.dataType,'TFT')
        theData=abs(theData); %convert complex number to real number
        theData=theData/mean(diff(EPscanCont.EPdata(theColor).freqNames)); %convert to spectral density
        theData=theData.^2; %convert amplitude to power
        theData=log10(abs(theData))*10; %convert to dB log scaling
        tempVar=theData;
        tempVar(isinf(tempVar))=-flintmax;
        theData=tempVar; %log10 of zero is -inf.  Replace with maximum possible double-precision negative number.
    end
    theData=theData-median(theData);
    %     cData=ceil((theData/max(theData))*size(cmap,1));
    %     for iPoint=1:length(cData)
    %         if cData(iPoint) > 0
    %             colorData(iPoint,:)=cmap(max(cData(iPoint),1),:);
    %         end
    %     end
    colorData(:,1)=((theData-mean(theData))/max(theData-mean(theData)))>.75;
    %colorData(:,3)=(theData-mean(theData))/min(theData-mean(theData));
    %colorData(colorData<0)=0;
%     if strcmp(EPscanCont.dataType,'TFT')
%         colorData(:,1)=(theData-mean(theData))/max(theData-mean(theData));
%         colorData(:,3)=(theData-mean(theData))/min(theData-mean(theData));
%         colorData(colorData<0)=0;
%     end
end

stimPict=[];
if ~isempty(EPscanCont.stimList{theColor})
    whichPict=max(find([EPscanCont.stimList{theColor}.sample] <= EPscanCont.displayPoints(thePoint)));
    if ~isempty(whichPict)
        stimPict=EPscanCont.EPdata(theColor).stims(EPscanCont.stimList{theColor}(whichPict).stim).image;
    end
end

fixationEvents=find(strcmp('fixationET',{EPscanCont.EPdata(theColor).events{1}.value}));
fixationEvents(~ismember(round([EPscanCont.EPdata(theColor).events{1}(fixationEvents).sample]),EPscanCont.displayPoints))=[];

Xdata=squeeze(EPscanCont.EPdata(theColor).data(EPscanCont.XEYchan,EPscanCont.displayPoints,:,:,:,1));
Ydata=squeeze(EPscanCont.EPdata(theColor).data(EPscanCont.YEYchan,EPscanCont.displayPoints,:,:,:,1));
Xdata=((Xdata-EPscanCont.eyeCenterX)/EPscanCont.eyeScaleX)+.5;
Ydata=((Ydata-EPscanCont.eyeCenterY)/EPscanCont.eyeScaleY)+.5;
EPscanCont.handles.topos.expandEye=axes(EPscanCont.handles.waves.frameFigure);
axis([0 1 0 1])
set(EPscanCont.handles.topos.expandEye,'Units','normalized')
set(EPscanCont.handles.topos.expandEye,'XTickLabel','','YTickLabel','','XTick',[],'YTick',[]);
hold on
if ~isempty(stimPict)
    imageSize=size(stimPict);
    imageMax=max(imageSize(1),imageSize(2));
    x=max(1,round((imageSize(1)-imageSize(2))/2))/imageMax;
    y=max(1,round((imageSize(2)-imageSize(1))/2))/imageMax;
    image([x 1-x],[1-y y],stimPict);
end

%plot AOIs
for iAOI=1:length(EPscanCont.EPdata(theColor).stims(EPscanCont.stimList{theColor}(whichPict).stim).AOI)
    ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;return;end
    if isfield(EPscanCont.EPdata(theColor).stims(EPscanCont.stimList{theColor}(whichPict).stim).AOI(iAOI),'tags') &&...
        ~isempty(EPscanCont.EPdata(theColor).stims(EPscanCont.stimList{theColor}(whichPict).stim).AOI(iAOI).tags{1})
        AOIx1=EPscanCont.EPdata(theColor).stims(EPscanCont.stimList{theColor}(whichPict).stim).AOI(iAOI).Coords(1);
        AOIx2=EPscanCont.EPdata(theColor).stims(EPscanCont.stimList{theColor}(whichPict).stim).AOI(iAOI).Coords(3);
        AOIy1=EPscanCont.EPdata(theColor).stims(EPscanCont.stimList{theColor}(whichPict).stim).AOI(iAOI).Coords(2);
        AOIy2=EPscanCont.EPdata(theColor).stims(EPscanCont.stimList{theColor}(whichPict).stim).AOI(iAOI).Coords(4);
        AOIx1Plot=((AOIx1-EPscanCont.eyeCenterX)/EPscanCont.eyeScaleX)+.5;
        AOIx2Plot=((AOIx2-EPscanCont.eyeCenterX)/EPscanCont.eyeScaleX)+.5;
        AOIy1Plot=((AOIy1-EPscanCont.eyeCenterY)/EPscanCont.eyeScaleY)+.5;
        AOIy2Plot=((AOIy2-EPscanCont.eyeCenterY)/EPscanCont.eyeScaleY)+.5;
        switch EPscanCont.EPdata(theColor).stims(EPscanCont.stimList{theColor}(whichPict).stim).AOI(iAOI).tags{1}
            case 'T'
                theColor='blue';
            case 'S'
                theColor='yellow';
            case 'A'
                theColor='cyan';
            case 'C'
                theColor='magenta';
        end
        s=patch([AOIx1Plot AOIx2Plot AOIx2Plot AOIx1Plot], [AOIy1Plot AOIy1Plot AOIy2Plot AOIy2Plot],theColor,'EdgeColor','none');
        alpha(s,.1);
    end
end

% if EPscanCont.chanShow~=length(EPscanCont.chanList)
%     for iPoint=1:length(EPscanCont.displayPoints)-1
% %        if any(colorData(iPoint,:)>.5)
%             %plot(Xdata(iPoint:iPoint+1),Ydata(iPoint:iPoint+1),'color',colorData(iPoint+1,:));
%             plot(Xdata(iPoint:iPoint+1),Ydata(iPoint:iPoint+1),'color','black');
% %        end
%     end
% else
%     plot(Xdata,Ydata,'color','black');
% end
% for iFix=1:length(fixationEvents)
%     theEvent=fixationEvents(iFix);
%     thePoint=find(ismember(EPscanCont.displayPoints,round(EPscanCont.EPdata(theColor).events{1}(theEvent).sample)));
%     scatter(Xdata(thePoint),Ydata(thePoint),240,'MarkerEdgeColor','black');
% %     if thePoint < length(EPscanCont.displayPoints)
% %         scatter(Xdata(thePoint),Ydata(thePoint),240,'filled','MarkerFaceAlpha',2/8,'MarkerFaceColor',colorData(thePoint+1,:));
% %     end
% end
% if strcmp(EPscanCont.dataType,'VLT')
%     scatter(Xdata(thePoint),Ydata(thePoint),30,'filled','MarkerFaceAlpha',1/8,'MarkerFaceColor',colorData(thePoint,:));
% else
if (thePoint+EPscanCont.ETlatency) > length(EPscanCont.displayPoints)
    scatter(Xdata(thePoint),Ydata(thePoint),240,'filled','MarkerFaceAlpha',1/8,'MarkerFaceColor',[0 0 0]);
else
    scatter(Xdata(thePoint),Ydata(thePoint),240,'filled','MarkerFaceAlpha',3/8,'MarkerFaceColor',colorData(thePoint+EPscanCont.ETlatency,:));
end
% end

%plot(Xdata(thePoint),Ydata(thePoint),'color','green','Marker','*','MarkerSize',2);
hold off
