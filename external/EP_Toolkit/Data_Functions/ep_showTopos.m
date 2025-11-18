function ep_showTopos(startTime,endTime,startHz,endHz,marker1,marker2,FFTunits,minVolt,maxVolt,eventLines)
% ep_showTopos - ep_showTopos(startTime,endTime,startHz,endHz,FFTunits,minVolt,maxVolt) -
% Displays the scalp topographies of the data and performs source analyses.
% Will choose the maximum channel, time point, and herz for display.  If multiple datasets are being displayed,
% will use the one in the leftmost column not set to "none" as the index dataset and only judge based on columns
% corresponding to that one dataset.  That way there is no problem of apples and orange comparison for max values,
% as in correlations for coherence versus power for FFT.
% Regional channels are not included in the peak calculations.
%
%Input:
%  startTime         : Start of epoch
%  endTime           : End of epoch in ms (right side unless both startTime and endTime are the same number)
%  startHz           : Start of Hz band
%  endHz             : End of Hz band
%  marker1           : Optional marker for plots
%  marker2           : Second optional marker for plots
%  FFTunits          : spectral data units (1=complex, 3=asd, 3=psd, 4=dB)
%  minVolt           : Minimum data value in raw units.
%  maxVolt           : Maximum data value in raw units.
%  eventLines        : Cell array of samples of events to be displayed (4 colors)(cells/subs)
%
%  Time is counted for an epoch as, say, -200 to 800 ms.  In this case the samples have some length to them depending on
% the digitization rate, such as 4ms for a 250Hz sampling rate.  The points prior to the baseline are counted from the
% left side of the baseline sample and the points following the baseline are counted from the right side of the baseline
% sample.  Thus, when the actual samples are calculated, one first adds the baseline to the numbers, yielding 0 to
% 1000ms.  One then must add the length of the baseline to the first number to obtain 4 to 1000ms.  Then one would
% divide the numbers by the sample size, resulting in 1 to 250 samples, which are the correct samples to use.
% For continuous data or other such data where there is no baseline, one must by this convention start with 0 ms,
% as in 0-1000ms and 1000-2000ms.

%History
%  by Joseph Dien (4/20/10)
%  jdien07@mac.com

% bugfix 5/9/10 JD
% Fixed crash when using Topo button of View EEG pane and not all four colors are being used.
% Fixed can only change channel and latency settings for number of rows equal to number of pages of factors.
% Fixed sometimes crashes when first color is set to none.
%
% modified 5/15/10 JD
% Added white marker to topos for electrode corresponding to the waveform figures.
% Small black dots indicate electrode locations in topographical plots.
% May click on electrode dots to move the waveform plot channel.
% May right-click on topographical plot to obtain expanded 2D plot.
% May right-click on topographical plot to obtain expanded 3D plot.
% May right-click on topographical plot to obtain basic dipole analysis.
% Eliminated sorting of channel names.
%
% bugfix 6/6/10 JD
% Fixed crash when trying to display 3D plot or dipole source using .ced file generated from a .elp file.
%
% modified 6/15/10 JD
% Added dipole analysis of jack-knifed PCA results to the Topos function of the View Pane.
% Marks when a file isn't saved yet.
%
% bugfix 7/23/10 JD
% Fixed crash when doing dipole or 3D function and ced is either empty or "none".
% Fixed 3D, dipole, and jackknife results incorrect for topos not on the first page (when there are multiple pages of
% topos).
%
% bugfix 12/26/10 JD
% Fixed bottom of Topos window being cut off on laptop screens.
%
% bugfix 1/9/11 JD
% Fixed channels not being unstandardized correctly in jack-knife, resulting in some inaccuracy in dipole results.
% Fixed crash when trying to display data where some channels are missing electrode coordinates.
%
% bugfix 2/12/11 JD
% Fixed contents of topos window getting shifted upwards off window when OS X Dock is at bottom of screen.
%
% bugfix 2/14/11 JD
% Fixed crash when displaying data with only one timepoint.
%
% bugfix 3/24/11 JD
% Fixed crash when displaying data with only two timepoints.
%
% modified 3/15/12 JD
% Added support for plotting FFT data, including psd and dB scaling.
%
% bugfix 6/3/12 JD
% Fixed crash when conducting dipole analysis on data with regional channels.
%
% bugfix 6/4/12 JD
% Fixed crash after clicking on Done when the final colors are not PCA data.
%
% bugfix 9/1/12 JD
% Fixed minimum voltage scaling being set to -1000 whenever minimum voltages are below 1000.
%
% bugfix 11/4/12 JD
% Fixed crash when invoking 2D plots.
%
% modified 1/16/13 JD
% Handles situation where FFT data has negative values (e.g., due to imprecision in PCA results) and transforming to dB
% would result in complex numbers, by taking absolute value first.
%
% modified 1/24/13 JD
% Added option to contextual menu to rescale figures according to selected topo map.
% Changed peak Hz calculation for FFT data to operate on absolute values to accommodate sometimes negative numbers from
% PCA results.
%
% modified 1/25/13 JD
% Added markers and expanding window to waveform figures.
%
% bugfix 1/31/13 JD
% Fixed a variety of issues with the automatic scaling in the topomap figures, especially for spectral data.
%
% bugfix 2/6/13 JD
% Fixed peak channels not being identified correctly.
%
% bugfix 2/19/13 JD
% Fixed crash when displaying data with a regional channel.
%
% bugfix 2/25/13 JD
% Fixed crash when displaying frequency data in dB scaling and the maximum value is negative.
%
% bugfix 3/24/13 JD
% Fixed peak samples of ERP data not being identified by absolute amplitude.
% Fixed Topos not allowing two datasets to be shown in parallel when they have different regional channels.
% Fixed Topos crashing when trying to display factor and non-factor data side by side.
%
% modified 3/24/13 JD
% Added display of topos at every 50 ms for ERP and TFT data and every Hz for FFT data for non-factor data.
%
% modified 4/2/13 JD
% Added peak point/Hz line to expanded waveform windows.
%
% bugfix 4/2/13 JD
% Markers in waveform plots can be set at zero ms.
%
% bugfix 4/12/13 JD
% Scaling of topos now obeys the values on the View pane.  Also, fixed manual changes to plotting range no longer
% working.
%
% modified 9/26/13 JD
% Restricted peak channels to EEG channels.
%
% bugfix 10/10/13 JD
% Fixed topoplot not showing the correct peak channel for non-factor data.
%
% bugfix 11/1/13 JD
% Fixed font sizes on Windows.
%
% bugfix 11/22/13 JD
% Fixed crash when trying to change scaling of frequency-domain data.
%
% bugfix 1/12/14 JD
% Workaround for Matlab bug periodically causing screen size to register as having zero size.
%
% bugfix 1/15/14 JD
% Fixed crash when trying to plot frequency data where electrode coordinate information is present but all coordinates are missing.
%
% bugfix 2/25/14 JD
% Fixed crash when trying to plot in dB data where the power equals zero.
%
% modified 3/18/14 JD
% Changed uses of "temp" as a variable name to "tempVar" due to other Matlab programmers often using it as a function
% name, resulting in collisions.
%
% modified 4/23/14 JD
% Added coherence and phase-locking options, including support for complex numbers.
% For peak channels and points and Hz, uses only the first dataset as the index for this.
%
% bufix 4/28/14 JD
% Fixed crash when changing pages in Topos view with frequency data.
% Fixed TFT color scale for expanded plots.
% Fixed 2D expanded head plots for Topos view for frequency data.
% Fixed jack-knife test possibly conducted on wrong factor or just crashing in Topos view.
%
% bufix 5/29/14 JD
% Fixed crash when there are electrodes without coordinates.
%
% modified 6/19/14 JD
% Added support for sample-by-sample t-tests, including STS chanType.
%
% bufix 8/5/14 JD
% Fixed not displaying frequency-domain data.
%
% modified 8/12/14 JD
% Added explore rereferencing option to contextual menus.
%
% bufix 12/31/14 JD
% Fixed crash with 3D head when original CED file not available.
%
% bufix 1/9/15 JD
% Fixed axis labels being shown for waveform plots under Matlab 2014b due
% to Matlab bug.
%
% bufix 3/22/15 JD
% Fixed peak chans sometimes not being computed correctly for voltage data.
%
%  modified 5/25/14 JD
%  Set colormap to jet even for Matlab 2014b onwards.
%
%  modified 5/29/14 JD
%  Changed min and max scale to be set by plotted data unless overriden.
%
%  bugfix 7/4/15 JD
%  Fixed crash when using rereference function to change the displayed
%  referencing.
%  Fixed crash when plotting two datasets with different numbers of
%  electrodes and the larger one is the first one.
%
% modified 7/5/15 JD
% Performs average reference prior to performing PARE correction.
%
% bugfix 8/30/15 JD
% Fixed amplitude spectral density calculated as divided by Hz rather than
% by square root of Hz.
% Fixed dB of amplitude data not converted to power first.
% Fixed dB units should be labeled as dBV since it is a ratio and therefore
% has no units.
%
% bugfix 1/15/16 JD
% Now allows power scaled data to be displayed as amplitudes.
% Now handles complex FFT numbers.
%
% modified 1/21/16 JD
% Consolidated spectral unit controls so just four options (cm, am, pw, dB).
%
% modified 3/8/16 JD
% Eliminated support for .power field.
%
% bugfix 10/9/16 JD
% Fixed FFT data on amplitude scaling not converted to real numbers.
%
% bugfix & modified 6/20/17 JD
% Fixed displaying imaginary component of FFT data even if FFT units not specified as being complex.
% Improved auto scaling for dB unit FFT data.
% Fixed conversion to spectral density dividing by bin width rather than sqrt(bin width).
% Fixed sometimes crashing when frequency band is just one Hz.
% Now uses ep_expandChan for expanded view as well.
% Added Compare Channels option to the View Topos pop-up menu.
%
% modified & bugfix 5/23/18 JD
% Fixed sometimes one too many pages when displaying factors, resulting in a crash when one visits the erroneous one.
% Fixed color bar missing.
% Added support for listing all trials or cells or subjects.
%
% modified 6/5/18 JD
% Added option to add marks for RT and selected events to View figures.
% Added Synch checkbox to Topos view so that when the channel or the time point is changed for one row, it is changed for all of them.
%
% bugfix 7/14/18 JD
% Now taking absolute value for calculating peak channel and timepoint, which can make a big difference for non average reference data.
% Improved auto-scaling for dB scaled frequency-domain data.
% Fixed dipole and jack-knife ricght-click options not available for frequency-domain data.
%
% modified 11/3/18 JD
% Added support for EEGlab derived Cartesian coordinates from standard-10-5-cap385.ced and using EEGlab's more precise transform for cartesian coordinates from .sfp file.
%
% modified & bugfix 12/1/18 JD
% Fixed View Topos sometimes crashing with -all- cells or -all- subjects.
% Added Electrodes option to contextual menu to see electrode locations on a 2D map.
% Fixed crash when displaying data using -all- factors option and the leftmost dataset is not PCA data.
%
% modified & bugfix 1/12/19 JD
% Improved auto-scaling of minimum and maximum data.
% Fixes to handling of coherence data.
%
% bugfix 4/11/19 JD
% Fixed crash in View Topos for continuous files with events.
%
% bugfix 5/7/19 JD
% Fixed crash for 3D heads due to changes in Matlab 2019a.
% Fixed crash when displaying continuous data in Matlab 2016a.
%
% modified & bugfix 6/17/19 JD
% Fixed crash when using fixed min and max voltages.
% Added popupmenu to change range to symmetrical numbers.
% Changed time mode so that every sample is displayed rather than just every 50 ms.
% Updated dipole analysis to support current FieldTrip code.
% Fixed Greek letter mu garbled in some figure labels.
%
% bugfix 8/19/19 JD
% Fixed crash for 3D heads.
%
% bugfix 9/17/19 JD
% Fixed not correctly labeling cells at top of columns for single-trial data.
%
% bugfix 9/23/19 JD
% Fixed crash when mixing -all- and single choices, as with factors.
% Fixed showing real rather than imaginary data in the imaginary lines when not using -all- option for subjects, trials, or factors.
% Fixed crash for 3D heads due to changes in EEGlab2019.
% Fixed crash when viewing continuous data with View Topos.
%
% modified 11/5/19 JD
% Added sessNums sessNames fields.
%
% bugfix 1/30/20 JD
% Fixed crash when sampling rate of voltage data is not an even multiple of 1000 (e.g., 512).
% Fixed crash when using -all- option and one of the datasets has less items than the first one.
%
% modified 2/14/20 JD
% 3D head and source analysis now uses improved method for co-registering eloc coordinates, based on Oostenveld's Standard-10-5-Cap385.sfp file.
%
% modified 3/13/20 JD
% Added support for BOSC data.
%
% modified 4/13/20 JD
% ep_saveEPdataset now handles replacing existing EPdataset entries.
%
% modified 4/22/20 JD
% Added support for up to eight colors in waveform figures.
%
% bugfix & modified 7/21/20 JD
% Displaying blank data when there is only one timepoint.
% Added Topo Grid.
%
% bugfix 1/27/21 JD
% Fixed error when trying to produce 3D graph and there are regional EEG channels present.
%
% bugfix 3/3/21 JD
% Fixed problems after using all eight colors to display factor results.
%
% modified 8/9/21 JD
% Added support for canonical 10-05 coordinates in eloc structure.
% Fixed crash when viewing single-trial data.
%
% bugfix 1/15/22 JD
% Fixed crash when colors are non-sequential.
%
% bugfix 12/5/23 JD
% Fixed View>Topos with continuous data providing popupmenu with timepoints for the entire dataset, causing crash when a timepoint outside the specified time window is chosen.
% Fixed crash with View>Topos with Topogrids option when viewing as timepoints.
%
% modified 8/9/21 JD
% Now color codes the channel menus with blue for the minimum voltage and
% red for the maximum voltage so both peak channels can be readily found.
%
% modified 12/9/24 JD
% Now uses more attractive round topo maps for 2D expanded figures in View>Topos function.
%
% modified 1/17/25 JD
% Added lineStyle control for the waveform lines.
% In View>Topos, factor names are grayed out if they are smaller than the minimum factor variance criterion.
%
% bugfix 4/15/25 JD
% Fixed rereference context menu option not marking zero voltage line correctly.
% Added support for virtual grand averages.
% The topos in the first row can now be averaged over a time points and/or Hz window.
% Fixed various things not restricting themselves to just the displayed colors.
% Updated jack-knife function to handle frequency data.
% Improved determination of display range for dB unit FFT data in View Topos.
% Restored black electrode dots in main display.
% Fixed bug in display of TFT coherence data.
%
% modified 5/11/25 JD
% Added preference for electrode markings in topo maps.
% Fixed 2D contextual menu map not exactly corresponding to topomap.
%
% bugfix 5/26/25 JD
% Fixed crash when a regional channel is not at the end of the channels.
%
% modified 6/12/25 JD
% Sample-by-sample results can be rendered as highlighted region as long as the adjoining colors are normal waves,
% not just for case with one sample-by-sample color and two normal waves.
% Fixed crash when skipping colors.
% Fixed crash for PCA datasets that are only single-step.
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

global EPdataset EPmain EPtopos EPtictoc

ep_tictoc('begin');

scrsz = EPmain.scrsz;

windowHeight=scrsz(4)-100;

if EPtopos.page==0 %set up topoplot for the first time
    set(groot,'defaultAxesToolbarVisible','off');
    if EPmain.numColors<EPmain.maxColors
        paneWidth=800;
    else
        paneWidth=1300;
    end
    mainPane=get(EPmain.handles.hMainWindow,'Position');
    h=findobj('name','TopoWindow');
    for i=1:length(h)
        close(h(i))
    end
    EPtopos.handles.topos.topoWindow = figure('Name', 'TopoWindow', 'NumberTitle', 'off', 'Position',[mainPane(1)+mainPane(3)+1 101 paneWidth windowHeight], 'CloseRequestFcn',@closeTopo);
    colormap jet;
    EPtopos.plotWidth=140;
    EPtopos.plotHeight=100;
    EPtopos.topoSize=100;
    EPtopos.margin=.1;

    EPtopos.done=0;
    EPtopos.changed=0;
    
    EPtopos.perPage=min(floor(windowHeight/120)-1,10); %how many rows of figures per page
    EPtopos.page=1; %current page
    EPtopos.chans=[];
    EPtopos.minChans=[];
    EPtopos.maxChans=[];
    EPtopos.points=[];
    EPtopos.synch=0;
    
    EPtopos.FFTunits=FFTunits;
    EPtopos.CSD=NaN;
    
    EPtopos.marker1=marker1;
    EPtopos.marker2=marker2;
    
    EPtopos.twoChan=[];
    EPtopos.handles.twoChan.figure=gobjects(0);
    
    EPtopos.eventLines=eventLines;
    
    EPtopos.theFirstColor=0;
    
    EPtopos.nonRegRelChans=cell(EPmain.numColors,1);
    EPtopos.numNonRegRelChans=1;
    
    EPtopos.BOSC.BOSCcolor=[]; %no BOSC mode in View Topos
    EPtopos.BOSC.dataColor=[];
    
    EPtopos.rangeList{1}=''; %list of values for popupmenu for setting symmetrical voltage range
    if any(EPmain.view.rel)
        for iRange=.1:.1:1
            EPtopos.rangeList{end+1}=['+/-' num2str(iRange)];
        end
    else
        for iRange=1:20
            EPtopos.rangeList{end+1}=['+/-' num2str(iRange)];
        end
    end

    [EPtopos.legendNames, EPtopos.legendNumbers]=ep_waveNames; %labels for the waveform colors
    
    %make sure the plotted datasets are compatible
    numChans=0;
    EPtopos.numRows=1;
    for iColor=1:EPmain.numColors
        ep_tictoc;if EPtictoc.stop;return;end
        if EPmain.view.dataset(iColor) <= length(EPdataset.dataset) %if not set to "none"
            if numChans == 0
                EPtopos.flexMode=strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).timeUnits,'per');
                %windowName=[EPdataset.dataset(EPmain.view.dataset(iColor)).dataName '-' EPdataset.dataset(EPmain.view.dataset(iColor)).subNames{EPmain.view.subject(iColor)}];
                windowName=EPdataset.dataset(EPmain.view.dataset(iColor)).dataName;
                set(EPtopos.handles.topos.topoWindow,'Name',windowName);
                EPtopos.firstTime=min([EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames]);
                EPtopos.lastTime=max([EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames]);
                EPtopos.firstHz=min([EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames]);
                EPtopos.lastHz=max([EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames]);
                if isscalar(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames)
                    EPtopos.sampleSize=1000/EPdataset.dataset(EPmain.view.dataset(iColor)).Fs;
                else
                    EPtopos.sampleSize=mean(diff(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames));
                end
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames)
                    EPtopos.binSize=mean(diff(EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames));
                end
                implicit=EPdataset.dataset(EPmain.view.dataset(iColor)).implicit; %use implicits of the first dataset
                EPtopos.chanNames=EPdataset.dataset(EPmain.view.dataset(iColor)).chanNames; %use channel names of the first dataset
                EPtopos.nonRegChans=find(~ismember(EPdataset.dataset(EPmain.view.dataset(iColor)).chanTypes,{'ANS','ECG','REG'}));
                
                EPtopos.eloc=EPdataset.dataset(EPmain.view.dataset(iColor)).eloc;
                if isempty(EPtopos.eloc) || all(isempty([EPtopos.eloc.radius])) || all(isempty([EPtopos.eloc.theta]))
                    msg{1}=['Error: The dataset ' EPdataset.dataset(EPmain.view.dataset(iColor)).dataName ' has no electrode coordinates.'];
                    [msg]=ep_errorMsg(msg);
                    done
                    return
                end
                hasLoc=[];
                for iChan=1:length(EPtopos.nonRegChans)
                    if ~isempty(EPtopos.eloc(EPtopos.nonRegChans(iChan)).theta)
                        hasLoc(end+1)=EPtopos.nonRegChans(iChan);
                    end
                end  
                EPtopos.nonRegChans=hasLoc;
                
                EPtopos.ced=EPdataset.dataset(EPmain.view.dataset(iColor)).ced;
                EPtopos.dataName=EPdataset.dataset(EPmain.view.dataset(iColor)).dataName;
                if ~isempty(implicit)
                    implicit=implicit(~strcmp('FID',{implicit.type}));
                    EPtopos.chanNames=[EPtopos.chanNames {implicit.labels}]; %use channel names of the first dataset
                    EPtopos.eloc=[EPtopos.eloc implicit];
                end
                numChans=length(EPtopos.nonRegChans);
                if numChans==0
                    msg{1}='Error: No channels with known locations.';
                    [msg]=ep_errorMsg(msg);
                    done
                    return
                end
                
                if any(EPmain.view.allTrials(1:EPmain.numColors)==1)
                    if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'average')
                        EPtopos.numRows=length(EPdataset.dataset(EPmain.view.dataset(iColor)).subNames);
                        EPtopos.type='subject';
                    else %single-trial
                        EPtopos.numRows=length(EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames);
                        EPtopos.type='trial';
                    end
                elseif any(EPmain.view.allTrials(1:EPmain.numColors)==3)
                    EPtopos.numRows=length(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames);
                    EPtopos.type='factor';
                elseif any(EPmain.view.allTrials(1:EPmain.numColors)==5)
                    EPtopos.numRows=length(EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames);
                    EPtopos.type='cell';
                elseif any(strcmp(EPmain.view.dataTransform,{'VLT','TFT'}))
                    EPtopos.type='time';
                elseif strcmp(EPmain.view.dataTransform,'FFT')
                    EPtopos.type='freq';
                else
                    error('Programming Error - aborting.');
                end
                
                if isfield(EPtopos.eloc,'labels')
                    eLabels={EPtopos.eloc.labels};
                else
                    eLabels=cell(numChans,1);
                end
                for iChan=1:length(eLabels)
                    if isempty(eLabels{iChan})
                        eLabels{iChan}=EPtopos.chanNames{iChan};
                    end
                end
                %                 [B,chanIX{iColor}] = sort(eLabels);
                %                 if ~isempty(EPtopos.eloc)
                %                     sortedEloc=EPtopos.eloc(chanIX{iColor});
                %                 end
                %EPtopos.chanNames=EPtopos.chanNames(chanIX{iColor});
                if ~strcmp(EPtopos.type,'factor')
                    EPtopos.theFirstColor=iColor; %the leftmost color that is not set to "none" and hence used as the index color
                end
                EPtopos.rowList=[];
            else
                newImplicit=EPdataset.dataset(EPmain.view.dataset(iColor)).implicit;
                if ~isempty(newImplicit)
                    newImplicit=newImplicit(find(~strcmp('FID',{newImplicit.type})));
                end
                newNonRegChans=find(~ismember(EPdataset.dataset(EPmain.view.dataset(iColor)).chanTypes,{'ANS','ECG','REG'}));
                hasLoc=[];
                for iChan=1:length(newNonRegChans)
                    if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).eloc(newNonRegChans(iChan)).theta)
                        hasLoc(end+1)=newNonRegChans(iChan);
                    end
                end  
                newNonRegChans=hasLoc;

                if length(EPtopos.nonRegChans) ~= length(newNonRegChans)
                    msg{1}='Error: The number of channels are not compatible.';
                    [msg]=ep_errorMsg(msg);
                    done
                    return
                end
                EPtopos.firstTime =max(min([EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames]),EPtopos.firstTime);
                EPtopos.lastTime = min(max([EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames]),EPtopos.lastTime);
                EPtopos.firstHz =max(min([EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames]),EPtopos.firstHz);
                EPtopos.lastHz = min(max([EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames]),EPtopos.lastHz);
                if ~isnan(EPtopos.sampleSize) && EPtopos.sampleSize ~= mean(diff(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames))
                    msg{1}='Error: The sampling rates are not compatible.';
                    [msg]=ep_errorMsg(msg);
                    done
                    return
                end
                
                switch EPtopos.type
                    case 'subject'
                        EPtopos.numRows=max(EPtopos.numRows,length(EPdataset.dataset(EPmain.view.dataset(iColor)).subNames));
                    case 'trial'
                        EPtopos.numRows=max(EPtopos.numRows,length(EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames));
                    case 'factor'
                        EPtopos.numRows=max(EPtopos.numRows,length(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames));
                    case 'cell'
                        EPtopos.numRows=max(EPtopos.numRows,length(EPdataset.dataset(EPmain.view.dataset(iColor)).cellNames));
                end
                        
                newEloc=[EPdataset.dataset(EPmain.view.dataset(iColor)).eloc newImplicit];
                if ~isempty(EPtopos.eloc) && ~isempty(newEloc)
                    eLabels={newEloc.labels};
                    for chan=1:length(eLabels)
                        if isempty(eLabels{chan})
                            eLabels{chan}=EPdataset.dataset(EPmain.view.dataset(iColor)).chanNames{chan};
                        end
                    end
                    %[B,chanIX{iColor}] = sort(eLabels);
                    %newEloc=newEloc(chanIX{iColor});
                    if any([EPtopos.eloc(hasLoc).theta]-[newEloc(hasLoc).theta]) || any([EPtopos.eloc(hasLoc).radius]-[newEloc(hasLoc).radius])
                        msg{1}='Error: The electrode coordinates are not compatible.';
                        [msg]=ep_errorMsg(msg);
                        done
                        return
                    end
                    %                 elseif xor(isempty(EPtopos.eloc),isempty(newEloc)) %if one but not the other has no electrode coordinates
                    %                     msg{1}='Error: The electrode coordinates are not compatible.';
                    %                     [msg]=ep_errorMsg(msg);
                    %                     done
                    %                     return
                else
                    if isempty(newEloc)
                        msg{1}=['Error: The dataset ' EPdataset.dataset(EPmain.view.dataset(iColor)).dataName ' has no electrode coordinates.'];
                        [msg]=ep_errorMsg(msg);
                        done
                        return
                    end
                end
            end
            if strcmp(EPtopos.type,'factor') && EPtopos.theFirstColor==0
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)
                    EPtopos.theFirstColor=iColor; %set default color to the leftmost dataset that is PCA data.
                end
            end
        end
    end

    if startTime==endTime
        EPtopos.firstTime=max(EPtopos.firstTime,startTime);
        EPtopos.lastTime=min(EPtopos.lastTime,endTime); %assume already  left side of sample
    else
        EPtopos.firstTime=max(EPtopos.firstTime,startTime);
        EPtopos.lastTime=min(EPtopos.lastTime,endTime-EPtopos.sampleSize); %shift endTime to left side of sample
    end

    EPtopos.firstHz=max(EPtopos.firstHz,startHz);
    EPtopos.lastHz=min(EPtopos.lastHz,endHz);
    
    %compute epoch samples for each color
    EPtopos.numPoints=1;
    EPtopos.numHz=1;
    EPtopos.startBins=ones(1,EPmain.numColors);
    EPtopos.lastBins=ones(1,EPmain.numColors);
    EPtopos.startSamp=ones(1,EPmain.numColors);
    EPtopos.lastSamp=ones(1,EPmain.numColors);
    
    if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
        for iColor=1:EPmain.numColors
            ep_tictoc;if EPtictoc.stop;return;end
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).relNames)
                    EPtopos.nonRegRelChans{iColor}=EPtopos.nonRegChans;
                    EPtopos.numNonRegRelChans=length(EPtopos.nonRegChans);
                else
                    EPtopos.nonRegRelChans{iColor}=1;
                end
                
                [~, EPtopos.startBins(iColor)]=min(abs(EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames-EPtopos.firstHz));
                [~, EPtopos.lastBins(iColor)]=min(abs(EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames-EPtopos.lastHz));
                EPtopos.lastBins(iColor)=max(EPtopos.lastBins(iColor),EPtopos.startBins(iColor));
                EPtopos.numHz=EPtopos.lastBins(iColor)-EPtopos.startBins(iColor)+1;
            end
        end
        EPtopos.freqNameList=cell(0);
        plotDataList=find(EPmain.view.dataset(1:EPmain.numColors) <= length(EPdataset.dataset));
        for iBin=min(EPtopos.startBins(plotDataList)):max(EPtopos.lastBins(plotDataList))
            EPtopos.freqNameList{end+1}=sprintf('%5.1f',round(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).freqNames(iBin)*10)/10);
        end
    end
    if any(strcmp(EPmain.view.dataTransform,{'VLT','TFT'}))
        for iColor=1:EPmain.numColors
            ep_tictoc;if EPtictoc.stop;return;end
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                [~, EPtopos.startSamp(iColor)]=min(abs(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames-EPtopos.firstTime));
                [~, EPtopos.lastSamp(iColor)]=min(abs(EPdataset.dataset(EPmain.view.dataset(iColor)).timeNames-EPtopos.lastTime));
                EPtopos.lastSamp(iColor)=max(EPtopos.lastSamp(iColor),EPtopos.startSamp(iColor));
                EPtopos.numPoints=EPtopos.lastSamp(iColor)-EPtopos.startSamp(iColor)+1;
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).relNames)
                    EPtopos.nonRegRelChans{iColor}=EPtopos.nonRegChans;
                    EPtopos.numNonRegRelChans=length(EPtopos.nonRegChans);
                else
                    EPtopos.nonRegRelChans{iColor}=1;
                end
            end
        end
    end
    
    %organize the data for plotting

    EPtopos.colorIndex=[];
    EPtopos.plotColors=[];
    EPtopos.thePlotColors=[];
    EPtopos.plotLineIndex=cell(0);
    EPtopos.plotColorIndex=[];
    theDataset=0;
    EPtopos.totalData=zeros(numChans,EPtopos.numPoints,EPtopos.numRows,EPtopos.numHz,EPtopos.numNonRegRelChans,EPmain.numColors); %default is zero voltage if, for example, two factor datasets with different numbers of factors
    EPtopos.totalImagData=zeros(numChans,EPtopos.numPoints,EPtopos.numRows,EPtopos.numHz,EPtopos.numNonRegRelChans,EPmain.numColors); %default is zero voltage if, for example, two factor datasets with different numbers of factors
    EPtopos.colsForMax=[];
    EPtopos.complexData=0;
    EPtopos.plotForm=EPmain.view.dataTransform;
    if strcmp('TFT',EPmain.view.dataTransform) && (EPtopos.firstHz==EPtopos.lastHz)
        EPtopos.plotForm='VLT';
        EPtopos.type='time';
    end
    if strcmp('TFT',EPmain.view.dataTransform) && (EPtopos.firstTime==EPtopos.lastTime)
        EPtopos.plotForm='FFT';
        EPtopos.type='freq';
    end
    EPtopos.eventWave=cell(EPmain.numColors,1);
    EPtopos.boundary=cell(EPmain.numColors,1);
    theMax=0;
    for iColor=1:EPmain.numColors
        ep_tictoc;if EPtictoc.stop;return;end
        if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
            EPtopos.colorIndex(end+1)=iColor;
            if isempty(EPtopos.rowList)
                switch EPtopos.type
                    case 'subject'
                        EPtopos.rowList=[1:EPtopos.numRows];
                    case 'trial'
                        EPtopos.rowList=[1:EPtopos.numRows];
                    case 'factor'
                        EPtopos.rowList=[1:EPtopos.numRows];
                    case 'cell'
                        EPtopos.rowList=[1:EPtopos.numRows];
                    case 'time'
%                         if EPtopos.flexMode
%                             EPtopos.rowList=ceil([1:(5/EPtopos.sampleSize):EPtopos.numPoints]);
%                         else
%                             EPtopos.rowList=ceil([1:(50/EPtopos.sampleSize):EPtopos.numPoints]);
%                         end
                        EPtopos.rowList=[1:EPtopos.numPoints];
                        EPtopos.numRows=length(EPtopos.rowList);
                    case 'freq'
                        EPtopos.rowList=[1:EPtopos.numHz];
                        EPtopos.numRows=length(EPtopos.rowList);
                    otherwise
                        error('Programming Error - aborting.');
                end
            end
            
            if theDataset ~= EPmain.view.dataset(iColor) %load in the dataset if different from the one already loaded in.
                theDataset=EPmain.view.dataset(iColor);
                EPdata=ep_loadEPdataset(EPmain.view.dataset(iColor));
                if ~isempty(EPdata.GAVsubs)
                    numSubs=length(EPdata.subNames);
                    numVsubs=size(EPdata.GAVsubs,1)-1;
                    numRsubs=numSubs-numVsubs;
                    %convert virtual GAVEs to normal form so waveforms are available.
                    EPdata=ep_combineData(EPdata,'convert',{[],[],[],[],[],[]},[],[],[]);
                    if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
                    if isempty(EPdata)
                        msg{1}='Error: Could not convert virtual grand average into real one.';
                        [msg]=ep_errorMsg(msg);
                        return
                    end
                end
            end

            %While the global min and max values are determined by all four columns, the maximum for each row
            %(that determines which channel/latency/Hz to display) is determined by the columns that belong to the index dataset,
            %which is the leftmost one.  EPtopos.colsForMax are the columns belonging to the index dataset.
            
            if EPmain.view.dataset(EPtopos.theFirstColor) == EPmain.view.dataset(iColor)
                EPtopos.colsForMax=[EPtopos.colsForMax; iColor];
            end
            
            if any(strcmp(EPtopos.type,{'subject','time','freq','factor'})) || ~EPmain.view.allTrials(iColor)
                if isempty(EPdata.trialNames) %averaged data
                    theCell=EPmain.view.cell(iColor);
                else  %if single_trial data
                    theCell=intersect(find(strcmp(EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),EPdata.cellNames)),...
                        find(EPmain.view.trial(iColor)==EPdata.trialNames));
                end
            else
                theCell=1;
            end
            
            refChans=EPdata.reference.current;
            if isempty(refChans) && ~any(ismember({'AVG','CSD'},EPdata.reference.type))
                refChans=EPdata.reference.original;
            end
            EPtopos.reference(iColor)=EPdata.reference;
            
            if EPmain.view.rel(iColor)
                if strcmp(EPdata.dataType,'average')
                    EPtopos.goodRelChans{iColor}=find(squeeze(any(~isnan(EPdata.analysis.badChans(EPmain.view.subject(iColor),theCell,:)),2)));
                else
                    EPtopos.goodRelChans{iColor}=find(squeeze(any((EPdata.analysis.badChans(EPmain.view.subject(iColor),theCell,:)~=-1),2)));
                end
                
                if isscalar(refChans)
                    EPtopos.goodRelChans{iColor}=setdiff(EPtopos.goodRelChans{iColor},refChans); %coherence with a single reference channel is NaN.
                end
            else
                EPtopos.goodRelChans{iColor}=1;
            end
            
            if strcmp(EPdata.reference.type,'CSD')
                if EPtopos.CSD==0
                    disp('Warning: Some data are CSD and some are not.  Labeling will be according to the last dataset.');
                end
                EPtopos.CSD=1;
            else
                if EPtopos.CSD==1
                    disp('Warning: Some data are CSD and some are not.  Labeling will be according to the last dataset.');
                end
                EPtopos.CSD=0;
            end
            EPtopos.plotColors=[EPtopos.plotColors iColor];
            EPtopos.thePlotColors=[EPtopos.thePlotColors; EPmain.preferences.view.color(iColor).RGB];
            EPtopos.plotLineIndex{end+1}=EPmain.preferences.view.color(iColor).lineStyle;
            EPtopos.plotColorIndex=[EPtopos.plotColorIndex; EPmain.preferences.view.color(iColor).RGB];
            
            ep_tictoc;if EPtictoc.stop;return;end
            switch EPtopos.type
                case 'subject'
                    if ~isempty(EPdata.facNames)
                        theData=ep_expandFacs(EPdata,EPtopos.nonRegChans,EPtopos.startSamp(iColor):EPtopos.lastSamp(iColor),theCell,[],EPmain.view.factor(iColor),EPtopos.startBins(iColor):EPtopos.lastBins(iColor),EPtopos.nonRegRelChans{iColor});
                        ep_tictoc;if EPtictoc.stop;return;end
                    else
                        theData=EPdata.data(EPtopos.nonRegChans,EPtopos.startSamp(iColor):EPtopos.lastSamp(iColor),theCell,:,EPmain.view.factor(iColor),EPtopos.startBins(iColor):EPtopos.lastBins(iColor),EPtopos.nonRegRelChans{iColor});
                    end
                case {'trial','cell'}
                    if EPmain.view.allTrials(iColor)
                        if ~isempty(EPdata.facNames)
                            theData=ep_expandFacs(EPdata,EPtopos.nonRegChans,EPtopos.startSamp(iColor):EPtopos.lastSamp(iColor),[],EPmain.view.subject(iColor),EPmain.view.factor(iColor),EPtopos.startBins(iColor):EPtopos.lastBins(iColor),EPtopos.nonRegRelChans{iColor});
                            ep_tictoc;if EPtictoc.stop;return;end
                        else
                            theData=EPdata.data(EPtopos.nonRegChans,EPtopos.startSamp(iColor):EPtopos.lastSamp(iColor),:,EPmain.view.subject(iColor),EPmain.view.factor(iColor),EPtopos.startBins(iColor):EPtopos.lastBins(iColor),EPtopos.nonRegRelChans{iColor});
                        end
                    else
                        if ~isempty(EPdata.facNames)
                            theData=ep_expandFacs(EPdata,EPtopos.nonRegChans,EPtopos.startSamp(iColor):EPtopos.lastSamp(iColor),theCell,EPmain.view.subject(iColor),EPmain.view.factor(iColor),EPtopos.startBins(iColor):EPtopos.lastBins(iColor),EPtopos.nonRegRelChans{iColor});
                            ep_tictoc;if EPtictoc.stop;return;end
                        else
                            theData=EPdata.data(EPtopos.nonRegChans,EPtopos.startSamp(iColor):EPtopos.lastSamp(iColor),theCell,EPmain.view.subject(iColor),EPmain.view.factor(iColor),EPtopos.startBins(iColor):EPtopos.lastBins(iColor),EPtopos.nonRegRelChans{iColor});
                        end
                    end
                case 'factor'
                    theData=ep_expandFacs(EPdata,EPtopos.nonRegChans,EPtopos.startSamp(iColor):EPtopos.lastSamp(iColor),theCell,EPmain.view.subject(iColor),[],EPtopos.startBins(iColor):EPtopos.lastBins(iColor),EPtopos.nonRegRelChans{iColor});
                    ep_tictoc;if EPtictoc.stop;return;end
                    if isfield(EPdata.pca,'PCAmode')
                        if isfield(EPdata.pca,'jack')
                            if isfield(EPdata.pca.jack,'FacPat')
                                for theJK=1:size(EPdata.pca.jack.FacPat,3)
                                    for theFac=1:size(EPdata.pca.jack.FacPat,2)
                                        EPtopos.jack(iColor).FacPat(:,theFac,theJK)=EPdata.pca.jack.FacPat(:,theFac,theJK).*EPdata.pca.jack.varSD(:,theJK);
                                        EPtopos.jack(iColor).PCAmode=EPdata.pca.PCAmode;
                                    end
                                end
                            elseif isfield(EPdata.pca.jack,'FacPatST')
                                for theJK=1:size(EPdata.pca.jack.FacPatST,3)
                                    for theFac=1:size(EPdata.pca.jack.FacPatST,2)
                                        EPtopos.jack(iColor).FacPat(:,theFac,theJK)=EPdata.pca.jack.FacPatST(:,theFac,theJK).*EPdata.pca.jack.varSDST(floor((theFac-1)/EPdata.pca.numFacs2)+1,:,theJK)';
                                        EPtopos.jack(iColor).PCAmode=EPdata.pca.PCAmode2;
                                    end
                                end
                            elseif isfield(EPdata.pca.jack,'FacPat3')
                                for theJK=1:size(EPdata.pca.jack.FacPat3,3)
                                    for theFac=1:size(EPdata.pca.jack.FacPat3,2)
                                        EPtopos.jack(iColor).FacPat(:,theFac,theJK)=EPdata.pca.jack.FacPat3(:,theFac,theJK).*EPdata.pca.jack.varSD3(floor((theFac-1)/EPdata.pca.numFacs3)+1,:,theJK)';
                                        EPtopos.jack(iColor).PCAmode=EPdata.pca.PCAmode3;
                                    end
                                end
                            end
                        end
                    end
                case {'time','freq'}
                    if ~isempty(EPdata.facNames)
                        theData=ep_expandFacs(EPdata,EPtopos.nonRegChans,EPtopos.startSamp(iColor):EPtopos.lastSamp(iColor),theCell,EPmain.view.subject(iColor),EPmain.view.factor(iColor),EPtopos.startBins(iColor):EPtopos.lastBins(iColor),EPtopos.nonRegRelChans{iColor});
                        ep_tictoc;if EPtictoc.stop;return;end
                    else
                        theData=EPdata.data(EPtopos.nonRegChans,EPtopos.startSamp(iColor):EPtopos.lastSamp(iColor),theCell,EPmain.view.subject(iColor),EPmain.view.factor(iColor),EPtopos.startBins(iColor):EPtopos.lastBins(iColor),EPtopos.nonRegRelChans{iColor});
                    end
                otherwise
                    error('Programming Error - aborting.');
            end
            
            %if the data has an imaginary component, as in spectral data, and units are not complex then will be combined together.
            %but if coherence data, then will need to keep the real and imaginary components separate regardless of units
            if ~isreal(theData) && ((FFTunits ==1) || ((FFTunits > 1) && EPmain.view.correl(iColor)))
                theDataImag=imag(theData);
                theData=real(theData);
                EPtopos.complexData=1;
                EPtopos.plotLineIndex{end+1}=':';
                EPtopos.plotColorIndex=[EPtopos.plotColorIndex; EPmain.preferences.view.color(iColor).RGB];
            else
                theDataImag=[];
            end
            
            if any(strcmp(EPtopos.type,{'subject','trial','cell','factor','time'}))
                for iRow=1:EPtopos.numRows
                    ep_tictoc;if EPtictoc.stop;return;end
                    switch EPtopos.type
                        case 'subject'
                            if EPmain.view.allTrials(iColor)
                                if iRow <=size(theData,4)
                                    theRowData=squeeze(theData(:,:,:,iRow,:,:,:));
                                else
                                    theRowData=[];
                                end
                            elseif any(EPmain.view.allTrials(1:EPmain.numColors))
                                theRowData=squeeze(theData(:,:,:,EPmain.view.subject(iColor),:,:,:));
                            else
                                theRowData=theData;
                            end
                        case 'trial'
                            if EPmain.view.allTrials(iColor)
                                if iRow <=size(theData,3)
                                    theRowData=squeeze(theData(:,:,iRow,:,:,:,:));
                                else
                                    theRowData=[];
                                end
                            elseif any(EPmain.view.allTrials(1:EPmain.numColors))
                                theRowData=squeeze(theData(:,:,theCell,:,:,:,:));
                            else
                                theRowData=theData;
                            end
                        case 'factor'
                            if EPmain.view.allTrials(iColor)
                                if iRow <=size(theData,5)
                                    theRowData=squeeze(theData(:,:,:,:,iRow,:,:));
                                else
                                    theRowData=[];
                                end
                            elseif any(EPmain.view.allTrials(1:EPmain.numColors))
                                theRowData=squeeze(theData(:,:,:,:,EPmain.view.factor(iColor),:,:));
                            else
                                theRowData=theData;
                            end
                        case 'cell'
                            if EPmain.view.allTrials(iColor)
                                if iRow <=size(theData,3)
                                    theRowData=squeeze(theData(:,:,iRow,:,:,:,:));
                                else
                                    theRowData=[];
                                end
                            elseif any(EPmain.view.allTrials(1:EPmain.numColors))
                                theRowData=squeeze(theData(:,:,theCell,:,:,:,:));
                            else
                                theRowData=theData;
                            end
                        case 'time'
                            if EPmain.view.allTrials(iColor)
                                if iRow <=size(theData,4)
                                    theRowData=squeeze(theData(:,EPtopos.rowList(iRow),:,:,:,:,:));
                                else
                                    theRowData=[];
                                end
                            else
                                theRowData=theData;
                            end
                    end
                    if ~isempty(theDataImag)
                        switch EPtopos.type
                            case 'subject'
                                if EPmain.view.allTrials(iColor)
                                    if iRow <=size(theDataImag,4)
                                        theRowDataImag=squeeze(theDataImag(:,:,:,iRow,:,:,:));
                                    else
                                        theRowDataImag=[];
                                    end
                                elseif any(EPmain.view.allTrials(1:EPmain.numColors))
                                    theRowDataImag=squeeze(theDataImag(:,:,:,EPmain.view.subject(iColor),:,:,:));
                                else
                                    theRowDataImag=theDataImag;
                                end
                            case 'trial'
                                if EPmain.view.allTrials(iColor)
                                    if iRow <=size(theDataImag,3)
                                        theRowDataImag=squeeze(theDataImag(:,:,iRow,:,:,:,:));
                                    else
                                        theRowDataImag=[];
                                    end
                                elseif any(EPmain.view.allTrials(1:EPmain.numColors))
                                    theRowDataImag=squeeze(theDataImag(:,:,theCell,:,:,:,:));
                                else
                                    theRowDataImag=theDataImag;
                                end
                            case 'factor'
                                if EPmain.view.allTrials(iColor)
                                    if iRow <=size(theDataImag,5)
                                        theRowDataImag=squeeze(theDataImag(:,:,:,:,iRow,:,:));
                                    else
                                        theRowDataImag=[];
                                    end
                                elseif any(EPmain.view.allTrials(1:EPmain.numColors))
                                    theRowDataImag=squeeze(theData(:,:,:,:,EPmain.view.factor(iColor),:,:));
                                else
                                    theRowDataImag=theDataImag;
                                end
                            case 'cell'
                                if EPmain.view.allTrials(iColor)
                                    if iRow <=size(theDataImag,3)
                                        theRowDataImag=squeeze(theDataImag(:,:,iRow,:,:,:,:));
                                    else
                                        theRowDataImag=[];
                                    end
                                elseif any(EPmain.view.allTrials(1:EPmain.numColors))
                                    theRowDataImag=squeeze(theDataImag(:,:,theCell,:,:,:,:));
                                else
                                    theRowDataImag=[];
                                end
                            case 'time'
                                if EPmain.view.allTrials(iColor)
                                    if iRow <=size(theDataImag,4)
                                        theRowDataImag=squeeze(theDataImag(:,EPtopos.rowList(iRow),:,:,:,:,:));
                                    else
                                        theRowDataImag=[];
                                    end
                                else
                                    theRowDataImag=theDataImag;
                                end
                        end
                    end
                    if ~isempty(theRowData)
                        EPtopos.totalData(1:size(theData,1),1:size(theData,2),iRow,1:size(theData,6),1:size(theData,7),iColor)=theRowData; %allows for combining datasets with differing dimensions
                    end
                    if ~isempty(theDataImag) && ~isempty(theRowDataImag)
                        EPtopos.totalImagData(1:size(theDataImag,1),1:size(theDataImag,2),iRow,1:size(theDataImag,6),1:size(theDataImag,7),iColor)=theRowDataImag; %allows for combining datasets with differing dimensions
                    end
                end
            else
                EPtopos.totalData(1:size(theData,1),1:size(theData,2),1,1:size(theData,6),1:size(theData,7),iColor)=theData; %allows for combining datasets with differing dimensions
                if ~isempty(theDataImag)
                    EPtopos.totalImagData(1:size(theDataImag,1),1:size(theDataImag,2),1,1:size(theDataImag,6),1:size(theDataImag,7),iColor)=theDataImag; %allows for combining datasets with differing dimensions
                end
            end
            
            %event marks
            EPtopos.eventWave{iColor}=cell(EPtopos.numRows,1);
            EPtopos.boundary{iColor}=cell(EPtopos.numRows,1);
            for iRow=1:EPtopos.numRows
                theRow=EPtopos.rowList(iRow);
                tempEvents=[];
                if strcmp(EPdata.dataType,'continuous')
                    if strcmp('VLT',EPmain.view.dataTransform) && (size(EPtopos.totalData,2)>1)
                        for theChan=1:size(EPtopos.totalData,1)
                            EPtopos.totalData(theChan,:,iRow,:,iColor)=EPtopos.totalData(theChan,:,iRow,:,iColor)-mean(EPtopos.totalData(theChan,:,iRow,:,iColor),2); %center the waveforms if continuous
                        end
                    end
                    if ~isempty(EPdata.events{EPmain.view.subject(iColor),1})
                        tempEvents=EPdata.events{EPmain.view.subject(iColor),1}(([EPdata.events{1}.sample]>=EPtopos.startSamp(iColor)) & ([EPdata.events{1}.sample]<=EPtopos.lastSamp(iColor)));
                    end
                    whichEvent=1;
                else
                    switch EPtopos.type
                        case 'subject'
                            tempEvents=EPdata.events{min(theRow,size(EPdata.events,1)),EPmain.view.cell(iColor)};
                            whichEvent=iRow;
                        case 'cell'
                            tempEvents=EPdata.events{EPmain.view.subject(iColor),min(theRow,size(EPdata.events,2))};
                            whichEvent=iRow;
                        otherwise
                            tempEvents=EPdata.events{EPmain.view.subject(iColor),EPmain.view.cell(iColor)};
                            whichEvent=1;
                    end
                end
                if ~isempty(EPtopos.eventLines{iColor}) && iscell(EPtopos.eventLines{iColor}) && (length(EPtopos.eventLines{iColor})>=whichEvent) && ~isempty(EPtopos.eventLines{iColor}{whichEvent})
                    EPtopos.eventWave{iColor}{iRow}=histc(EPtopos.eventLines{iColor}{whichEvent},[1:size(EPtopos.totalData,2)]);
                    theMax=max([theMax,max(EPtopos.eventWave{iColor}{whichEvent})]);
                end
                if ~isempty(tempEvents)
                    boundaryEvents=find(strcmp('boundary',{tempEvents.value}));
                    if ~isempty(boundaryEvents)
                        EPtopos.boundary{iColor}{iRow}=tempEvents(boundaryEvents).sample;
                    end
                end
            end
        end
    end
    for iColor=1:EPmain.numColors
        ep_tictoc;if EPtictoc.stop;return;end
        for iRow=1:EPtopos.numRows
            if ~isempty(EPtopos.eventWave{iColor})
                if length(EPtopos.eventWave{iColor})>0
                    if ~isempty(EPtopos.eventWave{iColor}{iRow})
                        EPtopos.eventWave{iColor}{iRow}=EPtopos.eventWave{iColor}{iRow}/theMax; %rescale event information to between zero and one.
                    end
                end
            end
        end
    end
    
    for iColor=1:EPmain.numColors
        ep_tictoc;if EPtictoc.stop;return;end
        if EPmain.view.dataset(iColor) <= length(EPdataset.dataset) && ~EPmain.view.correl(iColor)
            if (EPmain.view.cell(iColor) > length(EPdataset.dataset(EPmain.view.dataset(iColor)).cellTypes)) || ~strcmp('STS',EPdataset.dataset(EPmain.view.dataset(iColor)).cellTypes(EPmain.view.cell(iColor)))
                %if the data are not correlations or STS output then rescale to the chosen units if FFT data
                if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
                    if (FFTunits > 1)
                        EPtopos.totalData(:,:,:,:,:,iColor)=abs(EPtopos.totalData(:,:,:,:,:,iColor)); %convert complex number to real number
                    end
                    EPtopos.totalData(:,:,:,:,:,iColor)=EPtopos.totalData(:,:,:,:,:,iColor)/sqrt(mean(diff(EPdataset.dataset(EPmain.view.dataset(iColor)).freqNames))); %convert to spectral density
                    if FFTunits > 2
                        EPtopos.totalData(:,:,:,:,:,iColor)=EPtopos.totalData(:,:,:,:,:,iColor).^2; %convert amplitude to power
                    end
                    if (FFTunits == 4)
                        if ~all(EPtopos.totalData(:,:,:,:,:,iColor) >=0)
                            disp('Some values negative so will take absolute value prior to log transform.  This can happen with PCA results due to minor imprecisions in PCA.');
                        end
                        EPtopos.totalData(:,:,:,:,:,iColor)=log10(abs(EPtopos.totalData(:,:,:,:,:,iColor)))*10; %convert to dB log scaling
                        tempVar=EPtopos.totalData(:,:,:,:,:,iColor);
                        tempVar(isinf(tempVar))=-flintmax;
                        EPtopos.totalData(:,:,:,:,:,iColor)=tempVar; %log10 of zero is -inf.  Replace with maximum possible double-precision negative number.
                    end
                end
            end
        end
    end
end

%check for sample test output
plotDataList=find(EPmain.view.dataset(1:EPmain.numColors) <= length(EPdataset.dataset));
EPtopos.STSdata=[];
EPtopos.normCells=[];
for iColor=1:length(plotDataList)
    theColor=plotDataList(iColor);
    if (EPmain.view.cell(theColor) <= length(EPdataset.dataset(EPmain.view.dataset(theColor)).cellNames))
        theCell=find(strcmp(EPmain.view.theCells{EPmain.view.dataset(theColor)}(EPmain.view.cell(theColor)),EPdataset.dataset(EPmain.view.dataset(theColor)).cellNames));
        theCell=theCell(1);
        if strcmp('STS',EPdataset.dataset(EPmain.view.dataset(theColor)).cellTypes(theCell))
            EPtopos.STSdata=[EPtopos.STSdata; theColor];
        end
        if ~any(ismember({'STS';'BSC'},EPdataset.dataset(EPmain.view.dataset(theColor)).cellTypes{theCell}))
            EPtopos.normCells=[EPtopos.normCells; 1];
        else
            EPtopos.normCells=[EPtopos.normCells; 0];
        end
    end
end
if EPtopos.complexData
    EPtopos.perPage=floor(EPtopos.perPage/2)*2; %if complex, ensure even number of rows per page and reduce by half so can present both imaginary and real
end

if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
    totalRows=length(EPtopos.rowList);
else
    totalRows=length(EPtopos.rowList)+1;
end

if EPtopos.complexData
    EPtopos.numPages=ceil((totalRows*2)/EPtopos.perPage);
else
    EPtopos.numPages=ceil((totalRows)/EPtopos.perPage);
end

if isempty(EPtopos.chans) %if first time through
    ep_tictoc;if EPtictoc.stop;return;end
    %find max channels and latencies and Hz based on the selected data for the index dataset.
    maxData=EPtopos.totalData;
    if ~any(EPmain.view.rel(EPtopos.colsForMax))
        maxData=maxData(:,:,:,:,1,:); %if the index data are not relationship data, drop this dimension even if some of the other plotted data is.
    end
    
    if any(EPmain.view.rel(EPtopos.colsForMax)) 
        %in order to determine the max chan/latency/Hz for coherence data, need to first collapse over the relations.
        maxData(:,:,:,:,:,find(EPmain.view.rel(EPtopos.colsForMax)))=abs(maxData(:,:,:,:,:,find(EPmain.view.rel(EPtopos.plotColors))));
        maxData2=zeros([size(maxData,1),size(maxData,2),size(maxData,3),size(maxData,4),1,size(maxData,6)]);
        for iColor=1:length(EPtopos.plotColors)
            maxData2(:,:,:,:,:,iColor)=mean(maxData(:,:,:,:,EPtopos.goodRelChans{EPtopos.plotColors(iColor)},iColor),5); %skip bad channels including the NaN for reference channel
        end
        maxData=maxData2;
    end
    if EPtopos.FFTunits == 4
        %if data is in dB units, then one does not want to use absolute values to determine maximum measurements.
        absData=maxData;
    else
        absData=abs(maxData);
    end

    theMin=inf;
    theMax=-inf;

    switch EPtopos.type
        case {'subject','trial','cell','factor'}
            EPtopos.chans=ones(length(EPtopos.rowList),1); %default chan is one
            EPtopos.minChans=ones(length(EPtopos.rowList),1); %default chan is one
            EPtopos.maxChans=ones(length(EPtopos.rowList),1); %default chan is one
            EPtopos.points=ones(length(EPtopos.rowList),1); %default point is one
            EPtopos.freqs=ones(length(EPtopos.rowList),1); %default freq is one
            for theFactor=1:length(EPtopos.rowList)
                theMax=max(theMax,max(reshape(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),1,[])));
                theMin=min(theMin,min(reshape(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),1,[])));
                [~, EPtopos.points(theFactor)]=max(max(reshape(shiftdim(absData(:,:,theFactor,:,:,EPtopos.colsForMax),1),EPtopos.numPoints,[])'));
                [A, EPtopos.minChans(theFactor)]=min(min(reshape(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])'));
                [B, EPtopos.maxChans(theFactor)]=max(max(reshape(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])'));
                if abs(A) > abs(B)
                    EPtopos.chans(theFactor)=EPtopos.minChans(theFactor);
                else
                    EPtopos.chans(theFactor)=EPtopos.maxChans(theFactor);
                end
                [~, EPtopos.freqs(theFactor)]=max(max(reshape(shiftdim(absData(:,:,theFactor,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])'));
            end
        case 'time'
            EPtopos.chans=ones(length(EPtopos.rowList)+1,1); %default chan is one
            EPtopos.minChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
            EPtopos.maxChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
            EPtopos.points=ones(length(EPtopos.rowList)+1,1); %default point is one
            EPtopos.freqs=ones(length(EPtopos.rowList)+1,1); %default freq is one
            [~, EPtopos.points(1)]=max(max(reshape(shiftdim(absData(:,:,1,:,:,EPtopos.colsForMax),1),EPtopos.numPoints,[])'));
            EPtopos.firstRow.latency1=EPtopos.points(1);
            EPtopos.firstRow.latency2=EPtopos.points(1);
            EPtopos.firstRow.Hz1=EPtopos.freqs(1);
            EPtopos.firstRow.Hz2=EPtopos.freqs(1);
            EPtopos.points(2:end)=EPtopos.rowList;
            EPtopos.rowList=[0 EPtopos.rowList];
            for thePoint=1:length(EPtopos.rowList)
                theMax=max(theMax,max(reshape(maxData(:,EPtopos.points(thePoint),:,:,:,EPtopos.colsForMax),1,[])));
                theMin=min(theMin,min(reshape(maxData(:,EPtopos.points(thePoint),:,:,:,EPtopos.colsForMax),1,[])));
                [A, EPtopos.minChans(thePoint)]=min(min(reshape(maxData(:,EPtopos.points(thePoint),:,:,:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])',[],1));
                [B, EPtopos.maxChans(thePoint)]=max(max(reshape(maxData(:,EPtopos.points(thePoint),:,:,:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])',[],1));
                if abs(A) > abs(B)
                    EPtopos.chans(thePoint)=EPtopos.minChans(thePoint);
                else
                    EPtopos.chans(thePoint)=EPtopos.maxChans(thePoint);
                end
                [~, EPtopos.freqs(thePoint)]=max(max(reshape(shiftdim(absData(:,EPtopos.points(thePoint),:,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])',[],1));
            end
        case 'freq'
            EPtopos.chans=ones(length(EPtopos.rowList)+1,1); %default chan is one
            EPtopos.minChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
            EPtopos.maxChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
            EPtopos.points=ones(length(EPtopos.rowList)+1,1); %default point is one
            EPtopos.freqs=ones(length(EPtopos.rowList)+1,1); %default freq is one
            theMax=max(theMax,max(reshape(absData(:,:,1,:,:,EPtopos.colsForMax),1,[])));
            [~, EPtopos.freqs(1)]=max(max(reshape(shiftdim(absData(:,:,1,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])'));
            EPtopos.firstRow.latency1=EPtopos.points(1);
            EPtopos.firstRow.latency2=EPtopos.points(1);
            EPtopos.firstRow.Hz1=EPtopos.freqs(1);
            EPtopos.firstRow.Hz2=EPtopos.freqs(1);
            EPtopos.freqs(2:end)=EPtopos.rowList;
            EPtopos.rowList=[0 EPtopos.rowList];
            for theFreq=1:length(EPtopos.rowList)
                theMax=max(theMax,max(reshape(maxData(:,:,:,EPtopos.freqs(theFreq),:,EPtopos.colsForMax),1,[])));
                theMin=min(theMin,min(reshape(maxData(:,:,:,EPtopos.freqs(theFreq),:,EPtopos.colsForMax),1,[])));
                [~, EPtopos.points(theFreq)]=max(max(reshape(shiftdim(absData(:,:,:,EPtopos.freqs(theFreq),:,EPtopos.colsForMax),1),EPtopos.numPoints,[])',[],1));
                [A, EPtopos.minChans(theFreq)]=min(min(reshape(absData(:,:,:,EPtopos.freqs(theFreq),:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])',[],1));
                [B, EPtopos.maxChans(theFreq)]=max(max(reshape(absData(:,:,:,EPtopos.freqs(theFreq),:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])',[],1));
                if abs(A) > abs(B)
                    EPtopos.chans(theFreq)=EPtopos.minChans(theFreq);
                else
                    EPtopos.chans(theFreq)=EPtopos.maxChans(theFreq);
                end
            end
        case 'chan'
            EPtopos.chans=ones(length(EPtopos.rowList)+1,1); %default chan is one
            EPtopos.minChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
            EPtopos.maxChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
            EPtopos.points=ones(length(EPtopos.rowList)+1,1); %default point is one
            EPtopos.freqs=ones(length(EPtopos.rowList)+1,1); %default freq is one
            EPtopos.firstRow.latency1=EPtopos.points(1);
            EPtopos.firstRow.latency2=EPtopos.points(1);
            EPtopos.firstRow.Hz1=EPtopos.freqs(1);
            EPtopos.firstRow.Hz2=EPtopos.freqs(1);
            [~, EPtopos.chans(1)]=max(max(reshape(absData(:,:,1,:,EPtopos.colsForMax),EPtopos.numNonRegRelChans,[])'));
            EPtopos.chans(2:end)=EPtopos.rowList;
            EPtopos.minChans(2:end)=EPtopos.rowList;
            EPtopos.maxChans(2:end)=EPtopos.rowList;
            EPtopos.rowList=[0 EPtopos.rowList];
            for theChan=1:length(EPtopos.rowList)
                theMax=max(theMax,max(reshape(maxData(EPtopos.chans(theChan),:,:,:,:,EPtopos.colsForMax),1,[])));
                theMin=min(theMin,min(reshape(maxData(EPtopos.chans(theChan),:,:,:,:,EPtopos.colsForMax),1,[])));
                [~, EPtopos.points(theChan)]=max(max(reshape(absData(EPtopos.chans(theChan),:,:,:,:,EPtopos.colsForMax),EPtopos.numPoints,[])',[],1));
                [~, EPtopos.freqs(theChan)]=max(max(reshape(shiftdim(absData(EPtopos.chans(theChan),:,:,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])',[],1));
            end
        otherwise
            error('Programming Error - aborting.');
            return
    end
    if EPtopos.complexData
        maxData=EPtopos.totalImagData;
        
        if ~any(EPmain.view.rel(EPtopos.colsForMax))
            maxData=maxData(:,:,:,:,1,:); %if the index data are not relationship data, drop this dimension even if some of the other plotted data is.
        end
        
        if any(EPmain.view.rel(EPtopos.colsForMax)) %in order to determine the max chan/latency/Hz for coherence data, need to first collapse over the relations.
            maxData(:,:,:,:,:,find(EPmain.view.rel(EPtopos.colsForMax)))=abs(maxData(:,:,:,:,:,find(EPmain.view.rel(EPtopos.plotColors))));
            maxData2=zeros([size(maxData,1),size(maxData,2),size(maxData,3),size(maxData,4),1,size(maxData,6)]);
            for iColor=1:length(EPtopos.plotColors)
                maxData2(:,:,:,:,:,iColor)=mean(maxData(:,:,:,:,EPtopos.goodRelChans{EPtopos.plotColors(iColor)},iColor),5); %skip bad channels including the NaN for reference channel
            end
            maxData=maxData2;
        end
        if EPtopos.FFTunits == 4
            %if data is in dB units, then one does not want to use absolute values to determine maximum measurements.
            absData=maxData;
        else
            absData=abs(maxData);
        end

        switch EPtopos.type
            case {'subject','trial','cell','factor'}
                EPtopos.iChans=ones(length(EPtopos.rowList),1); %default chan is one
                EPtopos.iMinChans=ones(length(EPtopos.rowList),1); %default chan is one
                EPtopos.iMaxChans=ones(length(EPtopos.rowList),1); %default chan is one
                EPtopos.iPoints=ones(length(EPtopos.rowList),1); %default point is one
                EPtopos.iFreqs=ones(length(EPtopos.rowList),1); %default freq is one
                theMax=max(theMax,max(reshape(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),1,[])));
                theMin=min(theMin,min(reshape(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),1,[])));
                for theFactor=1:length(EPtopos.rowList)
                    [~, EPtopos.iPoints(theFactor)]=max(max(reshape(shiftdim(absData(:,:,theFactor,:,:,EPtopos.colsForMax),1),EPtopos.numPoints,[])'));
                    [A, EPtopos.iMinChans(theFactor)]=min(min(reshape(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])'));
                    [B, EPtopos.iMaxChans(theFactor)]=max(max(reshape(maxData(:,:,theFactor,:,:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])'));
                    if abs(A) > abs(B)
                        EPtopos.iChans(theFactor)=EPtopos.minChans(theFactor);
                    else
                        EPtopos.iChans(theFactor)=EPtopos.maxChans(theFactor);
                    end
                    [~, EPtopos.iFreqs(theFactor)]=max(max(reshape(shiftdim(absData(:,:,theFactor,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])'));
                end
            case 'time'
                EPtopos.rowList=EPtopos.rowList(2:end); %undo the addition of the first row so it doesn't end up getting added twice.
                EPtopos.iChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
                EPtopos.iMinChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
                EPtopos.iMaxChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
                EPtopos.iPoints=ones(length(EPtopos.rowList)+1,1); %default point is one
                EPtopos.iFreqs=ones(length(EPtopos.rowList)+1,1); %default freq is one
                [~, EPtopos.points(1)]=max(max(reshape(shiftdim(absData(:,:,1,:,:,EPtopos.colsForMax),1),EPtopos.numPoints,[])'));
                
                EPtopos.iPoints(2:end)=EPtopos.rowList;
                EPtopos.rowList=[0 EPtopos.rowList];
                for thePoint=1:length(EPtopos.rowList)
                    theMax=max(theMax,max(reshape(maxData(:,EPtopos.points(thePoint),:,:,:,EPtopos.colsForMax),1,[])));
                    theMin=min(theMin,min(reshape(maxData(:,EPtopos.points(thePoint),:,:,:,EPtopos.colsForMax),1,[])));
                    [A, EPtopos.iMinChans(thePoint)]=min(min(reshape(maxData(:,EPtopos.points(thePoint),:,:,:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])',[],1));
                    [B, EPtopos.iMaxChans(thePoint)]=max(max(reshape(maxData(:,EPtopos.points(thePoint),:,:,:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])',[],1));
                    if abs(A) > abs(B)
                        EPtopos.iChans(thePoint)=EPtopos.iMinChans(thePoint);
                    else
                        EPtopos.iChans(thePoint)=EPtopos.iMaxChans(thePoint);
                    end
                    [~, EPtopos.iFreqs(thePoint)]=max(max(reshape(shiftdim(absData(:,EPtopos.iPoints(thePoint),:,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])',[],1));
                end
            case 'freq'
                EPtopos.rowList=EPtopos.rowList(2:end); %undo the addition of the first row so it doesn't end up getting added twice.
                EPtopos.iChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
                EPtopos.iMinChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
                EPtopos.iMaxChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
                EPtopos.iPoints=ones(length(EPtopos.rowList)+1,1); %default point is one
                EPtopos.iFreqs=ones(length(EPtopos.rowList)+1,1); %default freq is one
                [~, EPtopos.iFreqs(1)]=max(max(reshape(shiftdim(absData(:,:,1,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])'));
                
                EPtopos.iFreqs(2:end)=EPtopos.rowList;
                EPtopos.rowList=[0 EPtopos.rowList];
                for theFreq=1:length(EPtopos.rowList)
                    theMax=max(theMax,max(reshape(maxData(:,:,:,EPtopos.freqs(theFreq),:,EPtopos.colsForMax),1,[])));
                    theMin=min(theMin,min(reshape(maxData(:,:,:,EPtopos.freqs(theFreq),:,EPtopos.colsForMax),1,[])));
                    [~, EPtopos.iPoints(theFreq)]=max(max(reshape(shiftdim(absData(:,:,:,EPtopos.iFreqs(theFreq),:,EPtopos.colsForMax),1),EPtopos.numPoints,[])',[],1));
                    [A, EPtopos.iMinChans(theFreq)]=min(min(reshape(absData(:,:,:,EPtopos.freqs(theFreq),:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])',[],1));
                    [B, EPtopos.iMaxChans(theFreq)]=max(max(reshape(absData(:,:,:,EPtopos.freqs(theFreq),:,EPtopos.colsForMax),length(EPtopos.nonRegChans),[])',[],1));
                    if abs(A) > abs(B)
                        EPtopos.iChans(theFreq)=EPtopos.iMinChans(theFreq);
                    else
                        EPtopos.iChans(theFreq)=EPtopos.iMaxChans(theFreq);
                    end
                end
            case 'chan'
                EPtopos.rowList=EPtopos.rowList(2:end); %undo the addition of the first row so it doesn't end up getting added twice.
                EPtopos.iChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
                EPtopos.iMinChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
                EPtopos.iMaxChans=ones(length(EPtopos.rowList)+1,1); %default chan is one
                EPtopos.iPoints=ones(length(EPtopos.rowList)+1,1); %default point is one
                EPtopos.iFreqs=ones(length(EPtopos.rowList)+1,1); %default freq is one
                [~, EPtopos.iChans(1)]=max(max(reshape(absData(:,:,1,:,:,EPtopos.colsForMax),EPtopos.numNonRegRelChans,[])'));
                EPtopos.iChans(2:end)=EPtopos.rowList;
                EPtopos.iMinChans(2:end)=EPtopos.rowList;                
                EPtopos.iMaxChans(2:end)=EPtopos.rowList;
                EPtopos.rowList=[0 EPtopos.rowList];
                for theChan=1:length(EPtopos.rowList)
                    theMax=max(theMax,max(reshape(maxData(EPtopos.chans(theChan),:,:,:,:,EPtopos.colsForMax),1,[])));
                    theMin=min(theMin,min(reshape(maxData(EPtopos.chans(theChan),:,:,:,:,EPtopos.colsForMax),1,[])));
                    [~, EPtopos.iPoints(theChan)]=max(max(reshape(shiftdim(absData(EPtopos.iChans(theChan),:,:,:,:,EPtopos.colsForMax),1),EPtopos.numPoints,[])',[],1));
                    [~, EPtopos.iFreqs(theChan)]=max(max(reshape(shiftdim(absData(EPtopos.iChans(theChan),:,:,:,:,EPtopos.colsForMax),3),EPtopos.numHz,[])',[],1));
                end
            otherwise
                disp('Programming Error in ep_showTopos - aborting.');
                return
        end
        
    end
    
    switch EPtopos.type
        case 'subject'
            EPtopos.pageNames=EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).subNames(1:EPtopos.perPage:end);
        case 'trial'
            EPtopos.pageNames=EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).trialNames(1:EPtopos.perPage:end);
        case 'cell'
            EPtopos.pageNames=EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).cellNames(1:EPtopos.perPage:end);
        case 'factor'
            EPtopos.pageNames=EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).facNames(1:EPtopos.perPage:end);
        case 'time'
            EPtopos.pageNames=num2cell(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).timeNames(EPtopos.startSamp(1)-1+[EPtopos.points(2); EPtopos.points(EPtopos.perPage+1:EPtopos.perPage:end)]));
        case 'freq'
            if EPtopos.complexData
                numFreqs=length(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).freqNames);
                tempNames=cell(numFreqs*2,1);
                for iFreq=1:numFreqs
                    tempNames{(iFreq*2)-1}=['i' num2str(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).freqNames(iFreq))];
                    tempNames{iFreq*2}=num2str(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).freqNames(iFreq));
                end
            else
                tempNames=num2cell([EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).freqNames(1); EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).freqNames]);
            end
            EPtopos.pageNames=tempNames(1:EPtopos.perPage:end);
        case 'chan'
            tempNames={EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).chanNames{1}; EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).chanNames};
            EPtopos.pageNames=tempNames(1:EPtopos.perPage:end);
        otherwise
            disp('Programming Error - aborting.');
            return
    end

    if all(EPmain.view.correl(EPtopos.plotColors))
        theMin=min(min(EPtopos.totalData,[],'all'),min(EPtopos.totalImagData,[],'all'));
        theMax=max(max(EPtopos.totalData,[],'all'),max(EPtopos.totalImagData,[],'all'));
        theAbsMax=max([abs(theMin) abs(theMax)]);
        EPtopos.plotMVmin=-theAbsMax;
        EPtopos.plotMVmax=theAbsMax;
    else
        theAbsMax=max([abs(theMin) abs(theMax)]);
        if ~isempty(minVolt)
            EPtopos.plotMVmin=minVolt;
        else
            if ((theMax <= 0) || (theMin >= 0)) && ~any(EPmain.view.rel(EPtopos.colsForMax))
                EPtopos.plotMVmin=theMin;
            else
                EPtopos.plotMVmin=-theAbsMax;
            end
        end
        if ~isempty(maxVolt)
            EPtopos.plotMVmax=maxVolt;
        else
            if ((theMax <= 0) || (theMin >= 0)) && ~any(EPmain.view.rel(EPtopos.colsForMax))
                EPtopos.plotMVmax=theMax;
            else
                EPtopos.plotMVmax=theAbsMax;
            end
        end

        if EPtopos.plotMVmin == EPtopos.plotMVmax
            EPtopos.plotMVmin=-1;
            EPtopos.plotMVmax=1;
        end

        if (EPtopos.FFTunits == 4) %dB scaling
            if ~strcmp(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).cellTypes(EPmain.view.cell(EPtopos.theFirstColor)),'STS')
                EPtopos.plotMVmax=theMax;
            end
        end
    end

    if any(EPmain.view.correl(EPtopos.plotColors)) && ~all(EPmain.view.correl(EPtopos.plotColors))
        %if not all the waveforms are correlations but some are, then rescale the correlations to match the plots
        EPtopos.totalData(:,:,:,:,:,find(EPmain.view.correl))=(EPtopos.totalData(:,:,:,:,:,find(EPmain.view.correl))*(EPtopos.plotMVmax-EPtopos.plotMVmin))+EPtopos.plotMVmin;
        if EPtopos.complexData
            EPtopos.totalImagData(:,:,:,:,:,find(EPmain.view.correl))=(EPtopos.totalImagData(:,:,:,:,:,find(EPmain.view.correl))*(EPtopos.plotMVmax-EPtopos.plotMVmin))+EPtopos.plotMVmin;
        end
    end
end

if EPtopos.page < EPtopos.numPages
    numRows=EPtopos.perPage;
else
    if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
        if EPtopos.complexData
            numRows=mod((EPtopos.numRows*2)-1,EPtopos.perPage)+1;
        else
            numRows=mod(EPtopos.numRows-1,EPtopos.perPage)+1;
        end
    else
        if EPtopos.complexData
            numRows=mod((length(EPtopos.rowList)*2)-1,EPtopos.perPage)+1;
        else
            numRows=mod(length(EPtopos.rowList)-1,EPtopos.perPage)+1;
        end
    end
end

clf(EPtopos.handles.topos.topoWindow)
figure(EPtopos.handles.topos.topoWindow)

EPtopos.handles.synch = uicontrol('Style', 'checkbox', 'String', 'Synch','Value',EPtopos.synch,'FontSize',EPmain.fontsize,'Position', [20 windowHeight-80 60 30],...
    'TooltipString','When checked, changes in time or channel to one row changes all of them.',...
    'Callback', ['global EPtopos;','EPtopos.synch=get(EPtopos.handles.synch,''Value'');']);

EPtopos.handles.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
    'Position', [20 windowHeight-110 50 30], 'Callback', @done);

EPtopos.handles.plotMVmin = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
    'String',num2str(EPtopos.plotMVmin),'HorizontalAlignment','left',...
    'Position',[120 windowHeight-110 30 20],'Callback',@changeMV);

EPtopos.handles.plotMVmax = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
    'String',num2str(EPtopos.plotMVmax),'HorizontalAlignment','left',...
    'Position',[230 windowHeight-110 30 20],'Callback',@changeMV);

EPtopos.handles.plotRange = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
    'String',EPtopos.rangeList,'HorizontalAlignment','left','Value',1,...
    'Position',[160 windowHeight-110 60 20],'Callback',@changeMV);

theLabel='';

if all(EPmain.view.correl(EPmain.view.dataset(1:EPmain.numColors) <= length(EPdataset.dataset)))
    theLabel='Correlation';
elseif strcmp(EPmain.view.dataTransform,'VLT')
    if strcmp(EPtopos.CSD,'CSD')
        theLabel='V/m^2';
    else
        theLabel=[char(181) 'v'];
    end
elseif any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
    if (EPtopos.FFTunits > 2)
        if (EPtopos.FFTunits == 4)
        	theLabel='dBV';
        else
            if strcmp(EPtopos.CSD,'CSD')
                theLabel='(V^2)/(Hz*m^4)';
            else
                theLabel=['(' char(181) 'v^2)/Hz'];
            end
        end
    else
        if strcmp(EPtopos.CSD,'CSD')
            theLabel='V/(sqrt(Hz)*m^2)';
        else
            theLabel=[char(181) 'v/sqrt(Hz)'];
        end
    end
end

uicontrol('Style','text','HorizontalAlignment','left','String',theLabel,'FontSize',EPmain.fontsize,...
    'Position',[130 windowHeight-65 80 20]);

EPtopos.handles.topos.colorbar = axes('Units','pixel','position',[120 windowHeight-65 109 20],'Visible','off');
colorbar('location','southoutside','Visible','on', 'XTick', [])
% set(EPtopos.handles.topos.colorbar,'ButtonDownFcn',@expandColorBar)

for iColor=1:EPmain.numColors
    ep_tictoc;if EPtictoc.stop;return;end
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        uicontrol('Style','text','ForegroundColor',EPmain.preferences.view.color(iColor).RGB,'FontSize',EPmain.fontsize,...
            'String',EPdataset.dataset(EPmain.view.dataset(iColor)).dataName,'HorizontalAlignment','left',...
            'Position',[100+EPtopos.plotWidth+50+(iColor-1)*(EPtopos.topoSize+20) windowHeight-70 100 20]);
    end
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        if EPmain.view.allTrials(iColor)==5
            uicontrol('Style','text','ForegroundColor',EPmain.preferences.view.color(iColor).RGB,'FontSize',EPmain.fontsize,...
                'String','-all cells-','HorizontalAlignment','left',...
                'Position',[100+EPtopos.plotWidth+50+(iColor-1)*(EPtopos.topoSize+20) windowHeight-90 100 20]);
        else
            uicontrol('Style','text','ForegroundColor',EPmain.preferences.view.color(iColor).RGB,'FontSize',EPmain.fontsize,...
                'String',EPmain.view.theCells{EPmain.view.dataset(iColor)}(EPmain.view.cell(iColor)),'HorizontalAlignment','left',...
                'Position',[100+EPtopos.plotWidth+50+(iColor-1)*(EPtopos.topoSize+20) windowHeight-90 100 20]);
        end
    end
    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
        if strcmp(EPdataset.dataset(EPmain.view.dataset(iColor)).dataType,'single_trial')
            if EPmain.view.allTrials(iColor)==1
                uicontrol('Style','text','ForegroundColor',EPmain.preferences.view.color(iColor).RGB,'FontSize',EPmain.fontsize,...
                    'String','-all trials-','HorizontalAlignment','left',...
                    'Position',[100+EPtopos.plotWidth+50+(iColor-1)*(EPtopos.topoSize+20) windowHeight-110 100 20]);
            else
                uicontrol('Style','text','ForegroundColor',EPmain.preferences.view.color(iColor).RGB,'FontSize',EPmain.fontsize,...
                    'String',EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames(EPmain.view.trial(iColor)),'HorizontalAlignment','left',...
                    'Position',[100+EPtopos.plotWidth+50+(iColor-1)*(EPtopos.topoSize+20) windowHeight-110 100 20]);
            end
        else
            if EPmain.view.allTrials(iColor)==1
                uicontrol('Style','text','ForegroundColor',EPmain.preferences.view.color(iColor).RGB,'FontSize',EPmain.fontsize,...
                    'String','-all subs-','HorizontalAlignment','left',...
                    'Position',[100+EPtopos.plotWidth+50+(iColor-1)*(EPtopos.topoSize+20) windowHeight-110 100 20]);
            else
                uicontrol('Style','text','ForegroundColor',EPmain.preferences.view.color(iColor).RGB,'FontSize',EPmain.fontsize,...
                    'String',EPdataset.dataset(EPmain.view.dataset(iColor)).subNames{EPmain.view.subject(iColor)},'HorizontalAlignment','left',...
                    'Position',[100+EPtopos.plotWidth+50+(iColor-1)*(EPtopos.topoSize+20) windowHeight-110 100 20]);
            end
        end
    end
    if (EPmain.view.dataset(iColor) <= length(EPdataset.dataset)) && ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)
        if EPmain.view.allTrials(iColor)==3
            uicontrol('Style','text','ForegroundColor',EPmain.preferences.view.color(iColor).RGB,'FontSize',EPmain.fontsize,...
                'String','-all facs-','HorizontalAlignment','left',...
                'Position',[100+EPtopos.plotWidth+50+(iColor-1)*(EPtopos.topoSize+20) windowHeight-130 100 20]);
        else
            uicontrol('Style','text','ForegroundColor',EPmain.preferences.view.color(iColor).RGB,'FontSize',EPmain.fontsize,...
                'String',EPdataset.dataset(EPmain.view.dataset(iColor)).facNames{EPmain.view.factor(iColor)},'HorizontalAlignment','left',...
                'Position',[100+EPtopos.plotWidth+50+(iColor-1)*(EPtopos.topoSize+20) windowHeight-130 100 20]);
        end
    end
end

if EPtopos.numPages > 1
    %page buttons
    if EPtopos.numPages <= 10
        for iPage=1:EPtopos.numPages
            EPtopos.handles.pageButton(iPage) = uicontrol('Style', 'pushbutton', 'String', EPtopos.pageNames{min(iPage,length(EPtopos.pageNames))},'FontSize',EPmain.fontsize,...
                'Position', [(20 + (iPage-1)*70) windowHeight-40 60 20], 'Callback', ['global EPtopos;','EPtopos.page=' num2str(iPage) ';','ep_showTopos(EPtopos.firstTime,EPtopos.lastTime)']);
        end
        set(EPtopos.handles.pageButton(EPtopos.page),'ForegroundColor','blue');
    else
        EPtopos.handles.pageMenu = uicontrol('Style', 'popupmenu', 'String', EPtopos.pageNames,'Value',EPtopos.page,'FontSize',EPmain.fontsize,...
            'Position', [20 windowHeight-40 100 20],...
            'Callback', ['global EPtopos;','tempVar=get(EPtopos.handles.pageMenu,''Value'');','if tempVar ~=0,EPtopos.page=tempVar;end;','ep_showTopos(EPtopos.firstTime,EPtopos.lastTime)']);
    end
end

%start drawing the figures
set(EPtopos.handles.topos.topoWindow ,'DefaultAxesXTickLabel',[])
set(EPtopos.handles.topos.topoWindow ,'DefaultAxesYTickLabel',[])
set(EPtopos.handles.topos.topoWindow ,'DefaultAxesXTick',[])
set(EPtopos.handles.topos.topoWindow ,'DefaultAxesYTick',[])
set(EPtopos.handles.topos.topoWindow,'DefaultAxesColorOrder',EPtopos.thePlotColors)

numChans=length(EPtopos.eloc);
%start adapted code from EEGlab's topoplot

locChans=true(numChans,1);
Rd=[];
Th=[];
for chan=1:numChans
    theRd=EPtopos.eloc(chan).radius;
    if isempty(theRd)
        locChans(chan)=false;
    end
    theTheta=EPtopos.eloc(chan).theta;
    if isempty(theTheta)
        locChans(chan)=false;
    end
    if locChans(chan)
        Rd=[Rd theRd];
        Th=[Th theTheta];
    end
end
EPtopos.locChans=find(locChans);

Th = pi/180*Th;
[x,y] = pol2cart(Th,Rd);
plotrad = min(1.0,max(Rd)*1.02);            % default: just outside the outermost electrode location
plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
squeezefac = .5/plotrad;
x    = x*squeezefac;
y    = y*squeezefac;
%end code from EEGlab's topoplot

if EPtopos.page ==1 && ~any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
    if EPtopos.complexData
        uicontrol('Style','frame',...
            'Position',[0 windowHeight-345 800 1]);
    else
        uicontrol('Style','frame',...
            'Position',[0 windowHeight-225 800 1]);
    end
end

if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
    if EPtopos.complexData 
        EPtopos.firstRow.Hz1=EPtopos.iFreqs(1);
        EPtopos.firstRow.Hz2=EPtopos.iFreqs(1);
    else
        EPtopos.firstRow.Hz1=EPtopos.freqs(1);
        EPtopos.firstRow.Hz2=EPtopos.freqs(1);
    end
else
    EPtopos.freqNameList='none';
    EPtopos.firstRow.Hz1=1;
    EPtopos.firstRow.Hz2=1;
end


ep_tictoc;if EPtictoc.stop;return;end
for iRow=1:numRows
    if EPtopos.complexData
        rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(iRow/2);
    else
        rowCounter=(EPtopos.page-1)*EPtopos.perPage+iRow;
    end
    theRow=1;
    if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
        theRow=rowCounter;
    end
    
    switch EPtopos.type
        case 'subject'
            uicontrol('Style','text',...
                'String',EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).subNames{min(rowCounter,length(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).subNames))},'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[20 windowHeight-(140+(iRow-1)*120) 60 20]);
        case 'trial'
            uicontrol('Style','text',...
                'String',num2str(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).trialNames(min(rowCounter,length(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).trialNames)))),'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[20 windowHeight-(140+(iRow-1)*120) 60 20]);
        case 'factor'
            theFactor=min(rowCounter,length(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).facNames));
            facNameHandle=uicontrol('Style','text',...
                'String',EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).facNames{theFactor},'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[20 windowHeight-(140+(iRow-1)*120) 60 20]);
            if strcmp(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).facTypes{theFactor},'SGL')
                %gray out factors that are smaller than the threshold value
                if isfield(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).pca,'facVar3')
                    if EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).pca.facVar3(theFactor) < EPmain.preferences.window.minFacVar
                        set(facNameHandle,'Enable','off');
                    end
                elseif isfield(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).pca,'facVarST')
                    if EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).pca.facVarST(theFactor) < EPmain.preferences.window.minFacVar
                        set(facNameHandle,'Enable','off');
                    end
                elseif isfield(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).pca,'facVar')
                    if EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).pca.facVar(theFactor) < EPmain.preferences.window.minFacVar
                        set(facNameHandle,'Enable','off');
                    end
                end
            end
        case 'cell'
            uicontrol('Style','text',...
                'String',EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).cellNames{min(rowCounter,length(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).cellNames))},'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[20 windowHeight-(140+(iRow-1)*120) 60 20]);
    end

    peakPoint=NaN;
    peakHz=NaN;
    if any(strcmp(EPmain.view.dataTransform,{'VLT','TFT'}))
        EPtopos.spacing=EPtopos.sampleSize;
        if EPtopos.complexData && rem(iRow,2)
            peakPoint=EPtopos.firstTime+(EPtopos.iPoints(rowCounter)-1)*EPtopos.sampleSize;
        else
            peakPoint=EPtopos.firstTime+(EPtopos.points(rowCounter)-1)*EPtopos.sampleSize;
        end
    end
    if any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'}))
        if EPtopos.complexData && rem(iRow,2)
            peakHz=EPtopos.iFreqs(rowCounter);
        else
            peakHz=EPtopos.freqs(rowCounter);
        end
    else
        EPtopos.freqNameList='none';
        peakHz=1;
    end
    if EPtopos.complexData && rem(iRow,2)
        peakChan=EPtopos.iChans(rowCounter);
    else
        peakChan=EPtopos.chans(rowCounter);
    end
    if (rowCounter==1) && any(strcmp(EPtopos.type,{'time';'freq';'chan'}))
        EPtopos.handles.Hz1=uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPtopos.freqNameList,'HorizontalAlignment','left',...
            'Value',EPtopos.firstRow.Hz1,...
            'Position',[10 windowHeight-(160+(iRow-1)*120) 50 20],...
            'Callback',{@changeHz,-1});
        EPtopos.handles.Hz2=uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPtopos.freqNameList,'HorizontalAlignment','left',...
            'Value',EPtopos.firstRow.Hz2,...
            'Position',[60 windowHeight-(160+(iRow-1)*120) 50 20],...
            'Callback',{@changeHz,-2});
        if  strcmp('VLT',EPmain.view.dataTransform)
            set(EPtopos.handles.Hz1,'enable','off');
            set(EPtopos.handles.Hz2,'enable','off');
        end
    else
        EPtopos.handles.Hz(iRow)=uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPtopos.freqNameList,'HorizontalAlignment','left',...
            'Value',peakHz,...
            'Position',[15 windowHeight-(160+(iRow-1)*120) 80 20],...
            'Callback',{@changeHz,iRow});
        EPtopos.handles.HzLabel(iRow)=uicontrol('Style','text','HorizontalAlignment','left','String', 'Hz','FontSize',EPmain.fontsize,...
            'Position',[90 windowHeight-(165+(iRow-1)*120) 20 20]);
        if  strcmp('VLT',EPmain.view.dataTransform)
            set(EPtopos.handles.Hz(iRow),'enable','off');
            set(EPtopos.handles.HzLabel(iRow),'enable','off');
        end
    end
    if (rowCounter==1) && any(strcmp(EPtopos.type,{'time';'freq';'chan'}))
        if any(strcmp(EPmain.view.dataTransform,{'VLT','TFT'}))
            EPtopos.spacing=EPtopos.sampleSize;
            peakPoint1=EPtopos.firstTime+(EPtopos.firstRow.latency1-1)*EPtopos.sampleSize;
            peakPoint2=EPtopos.firstTime+(EPtopos.firstRow.latency2-1)*EPtopos.sampleSize;
        else
            peakPoint1=NaN;
            peakPoint2=NaN;
        end
        EPtopos.handles.latency1=uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',num2str(round(peakPoint1)),'HorizontalAlignment','left',...
            'Position',[10 windowHeight-(180+(iRow-1)*120) 50 20],...
            'Callback',{@changeLatency,-1});
        EPtopos.handles.latency2=uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',num2str(round(peakPoint2)),'HorizontalAlignment','left',...
            'Position',[60 windowHeight-(180+(iRow-1)*120) 50 20],...
            'Callback',{@changeLatency,-2});
        if  strcmp('FFT',EPmain.view.dataTransform)
            set(EPtopos.handles.latency1,'enable','off');
            set(EPtopos.handles.latency2,'enable','off');
        end
    else
        EPtopos.handles.latency(iRow)=uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',num2str(round(peakPoint)),'HorizontalAlignment','left',...
            'Position',[20 windowHeight-(180+(iRow-1)*120) 50 20],...
            'Callback',{@changeLatency,iRow});
        if EPtopos.flexMode
            theUnit='%';
        else
            theUnit='ms';
        end
        EPtopos.handles.latencyLabel(iRow)=uicontrol('Style','text','HorizontalAlignment','left','String',theUnit,'FontSize',EPmain.fontsize,...
            'Position',[70 windowHeight-(185+(iRow-1)*120) 20 20]);
        if  strcmp('FFT',EPmain.view.dataTransform)
            set(EPtopos.handles.latency(iRow),'enable','off');
            set(EPtopos.handles.latencyLabel(iRow),'enable','off');
        end
    end
    
    chanNames=EPtopos.chanNames(EPtopos.nonRegChans)';
    if EPtopos.complexData && rem(iRow,2)
        chanNames{EPtopos.iMaxChans(rowCounter)}=['<HTML><FONT COLOR="red">' chanNames{EPtopos.iMaxChans(rowCounter)}];
        chanNames{EPtopos.iMinChans(rowCounter)}=['<HTML><FONT COLOR="blue">' chanNames{EPtopos.iMinChans(rowCounter)}];
    else
        chanNames{EPtopos.maxChans(rowCounter)}=['<HTML><FONT COLOR="red">' chanNames{EPtopos.maxChans(rowCounter)}];
        chanNames{EPtopos.minChans(rowCounter)}=['<HTML><FONT COLOR="blue">' chanNames{EPtopos.minChans(rowCounter)}];
    end
    EPtopos.handles.channel(iRow)=uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',chanNames',...
        'Value',peakChan,...
        'Position',[15 windowHeight-(200+(iRow-1)*120) 80 20],...
        'Callback',{@changeChannel,iRow});
    
    %draw waveforms
    
    if strcmp('FFT',EPtopos.plotForm)
        EPtopos.sampleSize=0;
        EPtopos.spacing=(EPtopos.lastHz-EPtopos.firstHz)/(EPtopos.numHz-1);
    else
        EPtopos.spacing=EPtopos.sampleSize;
    end
    
    switch EPtopos.plotForm
        case 'VLT'
            if EPtopos.numPoints > 1
                EPtopos.handles.waves.hWave(iRow) = axes('Units','pixel','position',[120 windowHeight-(220+(iRow-1)*120) EPtopos.plotWidth EPtopos.plotHeight],'XTickMode','manual','YTickMode','manual');
                
                if EPtopos.complexData && rem(iRow,2)
                    theData=EPtopos.totalImagData(EPtopos.iChans(rowCounter),:,theRow,:,:,EPtopos.plotColors);
                else
                    theData=EPtopos.totalData(EPtopos.chans(rowCounter),:,theRow,:,:,EPtopos.plotColors);
                end
                if any(EPmain.view.rel(EPmain.view.dataset(1:EPmain.numColors) <= length(EPdataset.dataset)))
                    theData(:,:,:,:,:,EPmain.view.rel(EPtopos.plotColors))=abs(theData(:,:,:,:,:,EPmain.view.rel(EPtopos.plotColors)));
                    theData=mean(theData,5); %collapse over relations if any
                end
                theData=squeeze(theData);
                if isscalar(EPtopos.plotColors)
                    theData=theData(:);
                end
                
                hold on
                for iWave=1:size(theData,2)
                    theWave=EPtopos.plotColors(iWave);
                    if any(EPtopos.STSdata==EPtopos.colorIndex(iWave)) && (iWave>1) && (iWave<size(theData,2)) && EPtopos.normCells(iWave-1) && EPtopos.normCells(iWave+1)
                        breakList=sort([find(diff([0 (theData(:,iWave)'>0) 0])<0)-1 find(diff([0 (theData(:,iWave)'>0) 0])>0)]);
                        if ~isempty(breakList)
                            theData1=theData(:,EPtopos.plotColors(iWave-1));
                            theData2=theData(:,EPtopos.plotColors(iWave+1));
                            for iSigArea=1:length(breakList)/2
                                theTimePoints=[breakList((iSigArea-1)*2+1):breakList(iSigArea*2)];
                                patch(([theTimePoints flip(theTimePoints)]*EPtopos.spacing)+(EPtopos.firstTime-EPtopos.spacing),[theData1(theTimePoints)' theData2(flip(theTimePoints))'],EPtopos.thePlotColors(iWave,:),'FaceColor',EPtopos.thePlotColors(iWave,:),'EdgeColor','none','FaceAlpha',.25);
                            end
                        end
                    else
                        EPtopos.handles.waves.hLines{iRow} = plot(EPdataset.dataset(EPmain.view.dataset(theWave)).timeNames(EPtopos.startSamp(theWave):EPtopos.lastSamp(theWave)),theData(:,iWave),'color',EPtopos.thePlotColors(iWave,:),'LineWidth',EPmain.preferences.view.lineSize);
                    end
                end
                hold off
                
                axis([EPtopos.firstTime EPtopos.lastTime EPtopos.plotMVmin EPtopos.plotMVmax]); %left side of first sample to left side of last sample
                if EPtopos.direction ==2
                    set(EPtopos.handles.waves.hWave(iRow),'YDir','reverse')
                end
                line([EPtopos.firstTime EPtopos.lastTime],[0 0],'Color','black','LineWidth',EPmain.preferences.view.lineSize) % zero line
                line([0 0],[0 EPtopos.plotMVmax],'Color','black','LineWidth',EPmain.preferences.view.lineSize) %stimulus onset
                if (rowCounter==1) && any(strcmp(EPtopos.type,{'time';'freq';'chan'}))
                    line(repmat(peakPoint1,2),[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',[.5 .5 .5],'LineWidth',EPmain.preferences.view.lineSize);
                    line(repmat(peakPoint2,2),[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',[.5 .5 .5],'LineWidth',EPmain.preferences.view.lineSize);
                else
                    line(repmat(peakPoint,2),[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',[.5 .5 .5],'LineWidth',EPmain.preferences.view.lineSize);
                end

                if ~isempty(EPtopos.marker1)
                    line(repmat(EPtopos.marker1,2),[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',[.5 .5 .5],'LineWidth',EPmain.preferences.view.lineSize);
                end
                if ~isempty(EPtopos.marker2)
                    line(repmat(EPtopos.marker2,2),[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',[.5 .5 .5],'LineWidth',EPmain.preferences.view.lineSize);
                end
                set(EPtopos.handles.waves.hWave(iRow),'ButtonDownFcn',{@ep_expandChan,'EPtopos'});
                for iLine=1:length(EPtopos.handles.waves.hLines{iRow})
                    set(EPtopos.handles.waves.hLines{iRow}(iLine),'ButtonDownFcn',{@ep_expandChan,'EPtopos'});
                end

                %plot event lines
                for iColor=1:length(EPtopos.plotColors)
                    theColor=EPtopos.plotColors(iColor);
                    switch EPtopos.type
                        case 'subject'
                            theEventWave=rowCounter;
                        case 'cell'
                            theEventWave=rowCounter;
                        otherwise
                            theEventWave=1;
                    end
                    if ~isempty(EPtopos.eventWave{theColor}{theEventWave}) && (EPtopos.plotMVmin < 0)
                        plotPoints=find(EPtopos.eventWave{theColor}{theEventWave}>min(EPtopos.eventWave{theColor}{theEventWave}));
                        plotTimes=[EPtopos.firstTime:EPtopos.spacing:EPtopos.firstTime+(EPtopos.spacing*(EPtopos.numPoints-1))];
                        if (EPtopos.plotMVmin < 0) && (EPtopos.plotMVmax >= 0)
                            if isscalar(plotPoints)
                                line([plotTimes(plotPoints) plotTimes(plotPoints)],[EPtopos.plotMVmin EPtopos.eventWave{theColor}{theEventWave}(plotPoints)*(EPtopos.plotMVmin/2)],'Color',EPtopos.thePlotColors(theColor,:),'LineWidth',2) %event line
                            else
                                hold on
                                EPtopos.handles.waves.eventLines{iRow,theColor} = plot([EPtopos.firstTime:EPtopos.spacing:EPtopos.firstTime+(EPtopos.spacing*(EPtopos.numPoints-1))],(EPtopos.eventWave{theColor}{theEventWave}*(abs(EPtopos.plotMVmin/2)))+EPtopos.plotMVmin,'LineWidth',5,'Color',EPtopos.thePlotColors(theColor,:));
                                hold off
                                for iLine=1:length(EPtopos.handles.waves.eventLines{iRow,theColor})
                                    set(EPtopos.handles.waves.eventLines{iRow,theColor}(iLine),'YDataSource',['EPtopos.eventWave{' num2str(theColor) '}(' num2str(iLine) ',:)']);
                                end
                            end
                        end
                    end
                    if ~isempty(EPtopos.boundary{theColor})
                        hold on
                        for iBoundary=1:length(EPtopos.boundary{theColor}{theEventWave})
                            theSample=EPtopos.boundary{theColor}{theEventWave}(iBoundary);
                            EPtopos.handles.waves.boundary{iRow,theColor} = line([theSample theSample],[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',EPtopos.thePlotColors(theColor,:),'LineWidth',EPmain.preferences.view.lineSize);
                        end
                        hold off
                    end
                end
            end
        case 'FFT'
            if EPtopos.numHz > 1
                EPtopos.handles.waves.hWave(iRow) = axes('Units','pixel','position',[120 windowHeight-(220+(iRow-1)*120) EPtopos.plotWidth EPtopos.plotHeight],'XTickMode','manual','YTickMode','manual');
                
                if EPtopos.complexData && rem(iRow,2)
                    theData=EPtopos.totalImagData(EPtopos.iChans(rowCounter),:,theRow,:,:,EPtopos.plotColors);
                    theLineStyle=':';
                else
                    theData=EPtopos.totalData(EPtopos.chans(rowCounter),:,theRow,:,:,EPtopos.plotColors);
                    theLineStyle='-';
                end
                
                if any(EPmain.view.rel(EPmain.view.dataset(1:EPmain.numColors) <= length(EPdataset.dataset)))
                    if strcmp(EPtopos.type,'freq')
                        if EPtopos.complexData && rem(iRow,2)
                            theFreq=EPtopos.freqs(rowCounter);
                        else
                            theFreq=EPtopos.iFreqs(rowCounter);
                        end
                        theData=theData(:,:,:,theFreq,:,find(EPmain.view.rel(EPtopos.plotColors)));
                    else
                        theData(:,:,:,:,:,find(EPmain.view.rel(EPtopos.plotColors)))=abs(theData(:,:,:,:,:,find(EPmain.view.rel(EPtopos.plotColors))));
                        theData2=zeros([size(theData,1),size(theData,2),size(theData,3),size(theData,4),1,size(theData,6)]);
                        for iColor=1:length(EPtopos.plotColors)
                            theData2(:,:,:,:,:,iColor)=mean(theData(:,:,:,:,EPtopos.goodRelChans{EPtopos.plotColors(iColor)},iColor),5); %collapse over relations if any
                        end
                        theData=theData2;
                    end
                end
                theData=squeeze(theData);
                if isscalar(EPtopos.plotColors)
                    theData=theData(:);
                end

                if strcmp(EPtopos.type,'freq')
                    EPtopos.handles.waves.hLines{iRow} = plot(1:numChans,theData,'LineStyle',theLineStyle,'LineWidth',EPmain.preferences.view.lineSize);
                    axis([1 numChans EPtopos.plotMVmin EPtopos.plotMVmax]);
                else
                    EPtopos.handles.waves.hLines{iRow} = plot(EPtopos.startBins(iColor):EPtopos.lastBins(iColor),theData,'LineStyle',theLineStyle,'LineWidth',EPmain.preferences.view.lineSize);
                    axis([EPtopos.firstHz EPtopos.lastHz EPtopos.plotMVmin EPtopos.plotMVmax]);
                end

                line([EPtopos.firstHz EPtopos.lastHz],[0 0],'Color','black','LineWidth',EPmain.preferences.view.lineSize) % zero line
                line(repmat(peakHz,2),[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',[.5 .5 .5],'LineWidth',EPmain.preferences.view.lineSize);
                if ~isempty(EPtopos.marker1)
                    line(repmat(EPtopos.marker1,2),[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',[.5 .5 .5],'LineWidth',EPmain.preferences.view.lineSize);
                end
                if ~isempty(EPtopos.marker2)
                    line(repmat(EPtopos.marker2,2),[EPtopos.plotMVmin EPtopos.plotMVmax],'Color',[.5 .5 .5],'LineWidth',EPmain.preferences.view.lineSize);
                end
                set(EPtopos.handles.waves.hWave(iRow),'ButtonDownFcn',{@ep_expandChan,'EPtopos'});
                for iLine=1:length(EPtopos.handles.waves.hLines{iRow})
                    set(EPtopos.handles.waves.hLines{iRow}(iLine),'ButtonDownFcn',{@ep_expandChan,'EPtopos'});
                end
            end
        case 'TFT'
            if EPtopos.numPoints > 1
                imageSpace=length(EPtopos.plotColors);
                imageCount=0;
                for iColor=1:EPmain.numColors
                    if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                        imageCount=imageCount+1;
                        EPtopos.handles.waves.hLines{iRow}(iColor) = axes('Units','pixel','position',[120 windowHeight-(220+(iRow-1)*120)+(EPtopos.plotHeight/imageSpace)*(imageSpace-imageCount) EPtopos.plotWidth (EPtopos.plotHeight/imageSpace)]);
                        if EPtopos.complexData && rem(iRow,2)
                            theData=EPtopos.totalImagData(EPtopos.iChans(rowCounter),:,theRow,:,:,iColor);
                        else
                            theData=EPtopos.totalData(EPtopos.chans(rowCounter),:,theRow,:,:,iColor);
                        end
                        if EPmain.view.rel(iColor)
                            theData=abs(theData);
                            theData=mean(theData(:,:,:,:,EPtopos.goodRelChans{iColor},1),5); %collapse over relations if any
                        end
                        theData=squeeze(theData)';
                        
                        EPtopos.handles.waves.hLines{iRow}(EPmain.numColors+iColor) = imagesc(EPtopos.firstTime:EPtopos.lastTime+EPtopos.sampleSize,EPtopos.firstHz:EPtopos.lastHz,theData,[EPtopos.plotMVmin EPtopos.plotMVmax]);
                        axis([EPtopos.firstTime EPtopos.lastTime+EPtopos.sampleSize EPtopos.firstHz EPtopos.lastHz]);
                        EPtopos.handles.waves.hLines{iRow}(iColor).Toolbar.Visible = 'off';
                        line([0 0],[EPtopos.firstHz EPtopos.lastHz],'Color','black','LineWidth',EPmain.preferences.view.lineSize) %stimulus onset
                        line([EPtopos.firstTime EPtopos.lastTime],[peakHz peakHz],'Color','white','LineWidth',EPmain.preferences.view.lineSize) % Hz line
                        line([peakPoint peakPoint],[EPtopos.firstHz EPtopos.lastHz],'Color',[.5 .5 .5],'LineWidth',EPmain.preferences.view.lineSize); %ms line
                        if ~isempty(EPtopos.marker1)
                            line(repmat(EPtopos.marker1,2),[EPtopos.firstHz EPtopos.lastHz],'Color',[.5 .5 .5],'LineWidth',EPmain.preferences.view.lineSize);
                        end
                        if ~isempty(EPtopos.marker2)
                            line(repmat(EPtopos.marker2,2),[EPtopos.firstHz EPtopos.lastHz],'Color',[.5 .5 .5],'LineWidth',EPmain.preferences.view.lineSize);
                        end
                        set(EPtopos.handles.waves.hLines{iRow}(iColor),'ButtonDownFcn',{@ep_expandChan,'EPtopos'});
                        set(EPtopos.handles.waves.hLines{iRow}(EPmain.numColors+iColor),'ButtonDownFcn',{@ep_expandChan,'EPtopos'});
                    end
                end
            end
    end
    
    %draw topos
    for iCol=1:EPmain.numColors
        if EPmain.view.dataset(iCol) <= length(EPdataset.dataset)
            EPtopos.handles.topos.topo(iRow,iCol) = axes('Units','pixel','position',[100+EPtopos.plotWidth+50+(iCol-1)*(EPtopos.topoSize+20) windowHeight-(220+(iRow-1)*120) EPtopos.topoSize EPtopos.topoSize]);
            if (rowCounter==1) && any(strcmp(EPtopos.type,{'time';'freq';'chan'}))
                if EPmain.view.rel(iCol) %if relational data
                    if EPtopos.complexData && rem(iRow,2)
                        theData=mean(mean(EPtopos.totalImagData(:,EPtopos.firstRow.latency1:EPtopos.firstRow.latency2,theRow,EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2,EPtopos.iChans(rowCounter),iCol),2),4);
                    else
                        theData=mean(mean(EPtopos.totalData(:,EPtopos.firstRow.latency1:EPtopos.firstRow.latency2,theRow,EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2,EPtopos.chans(rowCounter),iCol),2),4);
                    end
                else
                    if EPtopos.complexData && rem(iRow,2)
                        theData=mean(mean(EPtopos.totalImagData(:,EPtopos.firstRow.latency1:EPtopos.firstRow.latency2,theRow,EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2,1,iCol),2),4);
                    else
                        theData=mean(mean(EPtopos.totalData(:,EPtopos.firstRow.latency1:EPtopos.firstRow.latency2,theRow,EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2,1,iCol),2),4);
                    end
                end
                theData=squeeze(theData);
            else
                if EPmain.view.rel(iCol) %if relational data
                    if EPtopos.complexData && rem(iRow,2)
                        theData=EPtopos.totalImagData(:,EPtopos.iPoints(rowCounter),theRow,EPtopos.iFreqs(rowCounter),EPtopos.iChans(rowCounter),iCol);
                    else
                        theData=EPtopos.totalData(:,EPtopos.points(rowCounter),theRow,EPtopos.freqs(rowCounter),EPtopos.chans(rowCounter),iCol);
                    end
                else
                    if EPtopos.complexData && rem(iRow,2)
                        theData=EPtopos.totalImagData(:,EPtopos.iPoints(rowCounter),theRow,EPtopos.iFreqs(rowCounter),1,iCol);
                    else
                        theData=EPtopos.totalData(:,EPtopos.points(rowCounter),theRow,EPtopos.freqs(rowCounter),1,iCol);
                    end
                end
                theData=squeeze(theData);
            end
            if strcmp(EPdataset.dataset(EPmain.view.dataset(iCol)).cellTypes(EPmain.view.cell(iCol)),'STS')
                theData(isnan(theData))=0;
                STScorrect='theData(isnan(theData))=0;';
            else
                STScorrect='';
            end
            el_topoplot(theData, EPtopos.eloc(EPtopos.nonRegChans),'maplimits',[EPtopos.plotMVmin EPtopos.plotMVmax],'electrodes',EPmain.preferences.view.topoElectrodes,'efontsize',4);
            if EPtopos.complexData && rem(iRow,2)
                line([y(EPtopos.iChans(rowCounter)) y(EPtopos.iChans(rowCounter))],[x(EPtopos.iChans(rowCounter)) x(EPtopos.iChans(rowCounter))],'Marker','o','MarkerFaceColor','white'); %add marker for the electrode of the waveform plot
            else
                line([y(EPtopos.chans(rowCounter)) y(EPtopos.chans(rowCounter))],[x(EPtopos.chans(rowCounter)) x(EPtopos.chans(rowCounter))],'Marker','o','MarkerFaceColor','white'); %add marker for the electrode of the waveform plot
            end
            for iLine=1:length(EPtopos.locChans)
                EPtopos.handles.topos.line(iRow,iCol,iLine)=line([y(iLine) y(iLine)],[x(iLine) x(iLine)],'color','black','ButtonDownFcn',['global EPtopos;','sel_typ = get(gcbf,''SelectionType'');','if strcmp(sel_typ,''normal'')','EPtopos.chans(' num2str(rowCounter) ')=' num2str(EPtopos.locChans(iLine)) ';','ep_showTopos(EPtopos.firstTime,EPtopos.lastTime);','end']);
            end
            
            hcmenu = uicontextmenu;
            if (rowCounter==1) && any(strcmp(EPtopos.type,{'time';'freq';'chan'}))
                pointsRange='EPtopos.firstRow.latency1:EPtopos.firstRow.latency2';
                binsRange='EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2';
                if EPtopos.complexData
                    iPointsRange='EPtopos.firstRow.latency1:EPtopos.firstRow.latency2';
                    iBinsRange='EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2';
                end
            else
                pointsRange=num2str(EPtopos.points(rowCounter));
                binsRange=num2str(EPtopos.freqs(rowCounter));
                if EPtopos.complexData
                    iBinsRange=num2str(EPtopos.iFreqs(rowCounter));
                    iPointsRange=num2str(EPtopos.iPoints(rowCounter));
                end
            end

            if EPmain.view.rel(iCol) %if relational data
                if EPtopos.complexData && rem(iRow,2)
                    D2map1 = ['global EPtopos EPmain;','figure;','theData=squeeze(mean(mean(EPtopos.totalImagData(EPtopos.iChans(' num2str(rowCounter) '),' iPointsRange ',' num2str(theRow) ',' iBinsRange ',:,' num2str(iCol) '),2),4));',STScorrect,'el_topoplot(theData,EPtopos.eloc(EPtopos.nonRegChans),''maplimits'',[EPtopos.plotMVmin EPtopos.plotMVmax]'];
                    D2map2 = [',''electrodes'',EPmain.preferences.view.topoElectrodes'];
                    D2map3 = ['), disp(''Using EEGlab function topoplot to perform 2D head display.'');','colormap jet'];
                else
                    D2map1 = ['global EPtopos EPmain;','figure;','theData=squeeze(mean(mean(EPtopos.totalData(EPtopos.chans(' num2str(rowCounter) '),' pointsRange ',' num2str(theRow) ',' binsRange ',:,' num2str(iCol) '),2),4));',STScorrect,'el_topoplot(theData,EPtopos.eloc(EPtopos.nonRegChans),''maplimits'',[EPtopos.plotMVmin EPtopos.plotMVmax]'];
                    D2map2 = [',''electrodes'',EPmain.preferences.view.topoElectrodes'];
                    D2map3 = ['), disp(''Using EEGlab function topoplot to perform 2D head display.'');','colormap jet'];
                end
            else
                if EPtopos.complexData && rem(iRow,2)
                    D2map1 = ['global EPtopos EPmain;','figure;','theData=squeeze(mean(mean(EPtopos.totalImagData(:,' iPointsRange ',' num2str(theRow) ',' iBinsRange ',1,' num2str(iCol) '),2),4));',STScorrect,'el_topoplot(theData,EPtopos.eloc(EPtopos.nonRegChans),''maplimits'',[EPtopos.plotMVmin EPtopos.plotMVmax]'];
                    D2map2 = [',''electrodes'',EPmain.preferences.view.topoElectrodes'];
                    D2map3 = ['), disp(''Using EEGlab function topoplot to perform 2D head display.'');','colormap jet'];
                else
                    D2map1 = ['global EPtopos EPmain;','figure;','theData=squeeze(mean(mean(EPtopos.totalData(:,' pointsRange ',' num2str(theRow) ',' binsRange ',1,' num2str(iCol) '),2),4));',STScorrect,'el_topoplot(theData,EPtopos.eloc(EPtopos.nonRegChans),''maplimits'',[EPtopos.plotMVmin EPtopos.plotMVmax]'];
                    D2map2 = [',''electrodes'',EPmain.preferences.view.topoElectrodes'];
                    D2map3 = ['), disp(''Using EEGlab function topoplot to perform 2D head display.'');','colormap jet'];
                end
            end
            D2map=[D2map1 D2map2 D2map3];
            ElectrodeMap=[D2map1 ',''electrodes'',''labels''' D2map3];
            % Define the context menu items and install their callbacks
            item1 = uimenu(hcmenu, 'Label', '2D', 'Callback', D2map);
            item2 = uimenu(hcmenu, 'Label', '3D', 'Callback', @D3head);
            item3 = uimenu(hcmenu, 'Label', 'Electrodes',  'Callback', ElectrodeMap);
            item4 = uimenu(hcmenu, 'Label', 'Rescale',  'Callback', @rescaleFigures);
            item5 = uimenu(hcmenu, 'Label', '2-Channels',  'Callback', @twoChan);
            item6 = uimenu(hcmenu, 'Label', 'Topo Grid',  'Callback', {@topoMatrix,iRow,iCol});
            if strcmp('VLT',EPmain.view.dataTransform) && ~strcmp(EPtopos.CSD,'CSD')
                item7 = uimenu(hcmenu, 'Label', 'Rereference',  'Callback', @referenceChan);
            end
            if strcmp('VLT',EPmain.view.dataTransform) || ((any(strcmp(EPmain.view.dataTransform,{'FFT','TFT'})) && ~EPmain.view.rel(iCol))) && ~strcmp(EPtopos.CSD,'CSD') && (EPtopos.FFTunits == 1)
                item8 = uimenu(hcmenu, 'Label', 'Dipole',  'Callback', @dipoles);
                if isfield(EPtopos,'jack')
                    if length(EPtopos.jack) >=iCol
                        if strcmp(EPtopos.jack(iCol).PCAmode,'spat')
                            item8 = uimenu(hcmenu, 'Label', 'Jack-Knife',  'Callback', @jackknife);
                        end
                    end
                end
            end
            % Attach the context
            set(EPtopos.handles.topos.topo(iRow,iCol),'UIContextMenu',hcmenu)
            theChildren=get(EPtopos.handles.topos.topo(iRow,iCol),'Children');
            for iChild=1:length(theChildren)
                set(theChildren(iChild),'UIContextMenu',hcmenu)
            end
            for iLine=1:length(EPtopos.locChans)
                set(EPtopos.handles.topos.line(iRow,iCol,iLine),'UIContextMenu',hcmenu)
            end
        end
    end
end

ep_tictoc('end');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function done(src,eventdata)
global EPtopos EPmain EPdataset

%save results of jack-knife PCA dipole analyses
if EPtopos.changed
    for col=1:length(EPtopos.jack)
        if EPmain.view.dataset(col) <= length(EPdataset.dataset) %if not set to "none"
            if ~isempty(EPtopos.jack(col))
                if isfield(EPtopos.jack(col),'sources')
                    if ~isempty(EPtopos.jack(col).sources)
                        EPdata=ep_loadEPdataset(EPmain.view.dataset(col));
                        EPdata.pca.jack.sources=EPtopos.jack(col).sources;
                        
                        try
                            EPver=ver('EP_Toolkit');
                        catch
                            EPver='unavailable'; %workaround for bug in earlier version of Matlab
                        end
                        EPdata.EPver=EPver;
                        EPdata.ver=ver;
                        EPdata.date=date;
                        
                        [err]=ep_checkEPfile(EPdata);
                        
                        if ~exist([EPdataset.EPwork filesep 'EPwork'],'dir')
                            warndlg('The work directory cannot be found.')
                            return
                        end
                        
                        ep_saveEPdataset(EPdata,EPmain.view.dataset(col),'no');
                    end
                end
            end
        end
    end
end

EPtopos.changed=0;
EPtopos.done=1;

close(EPtopos.handles.topos.topoWindow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeMV(src,eventdata)
global EPtopos EPmain

if src==EPtopos.handles.plotRange
    theRange=get(EPtopos.handles.plotRange,'value')-1;
    if any(EPmain.view.rel)
        set(EPtopos.handles.plotMVmin,'string',num2str(-theRange/10));
        set(EPtopos.handles.plotMVmax,'string',num2str(theRange/10));
    else
        set(EPtopos.handles.plotMVmin,'string',num2str(-theRange));
        set(EPtopos.handles.plotMVmax,'string',num2str(theRange));
    end
    set(EPtopos.handles.plotRange,'value',1);
end
plotMVmin=str2double(get(EPtopos.handles.plotMVmin,'string'));
plotMVmax=str2double(get(EPtopos.handles.plotMVmax,'string'));

if plotMVmin < plotMVmax
    EPtopos.plotMVmin=plotMVmin;
    EPtopos.plotMVmax=plotMVmax;
end

ep_showTopos(EPtopos.firstTime,EPtopos.lastTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeLatency(src,eventdata,theRow)
global EPtopos

if theRow<0
    latency1=str2double(get(EPtopos.handles.latency1,'string'));
    latency2=str2double(get(EPtopos.handles.latency2,'string'));
    if latency1>latency2
        disp(['Error: First latency cannot be after second latency number'])
        return
    end
    if theRow==-1
        theLatency=latency1; %the input number in ms
        oldTime=EPtopos.firstTime+(EPtopos.firstRow.latency1-1)*EPtopos.sampleSize; %existing number number in ms
    else
        theLatency=latency2; %the input number in ms
        oldTime=EPtopos.firstTime+(EPtopos.firstRow.latency2-1)*EPtopos.sampleSize; %existing number number in ms
    end
    rowCounter=1;
else
    if EPtopos.complexData
        rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(theRow/2);
    else
        rowCounter=(EPtopos.page-1)*EPtopos.perPage+theRow;
    end
    theLatency=str2double(get(EPtopos.handles.latency(theRow),'string')); %the input number in ms
    if EPtopos.complexData && rem(theRow,2)
        oldTime=EPtopos.firstTime+(EPtopos.iPoints(rowCounter)-1)*EPtopos.sampleSize; %existing number number in ms
    else
        oldTime=EPtopos.firstTime+(EPtopos.points(rowCounter)-1)*EPtopos.sampleSize; %existing number number in ms
    end
end

newSample=round((theLatency-EPtopos.firstTime)/EPtopos.sampleSize)+1;
if (theLatency ~= oldTime) || EPtopos.synch
    if (theLatency >= EPtopos.firstTime) && (theLatency <= EPtopos.lastTime)
        if (newSample > 0) && (newSample <= EPtopos.numPoints)
            if theRow<0
                if theRow==-1
                    EPtopos.firstRow.latency1=newSample;
                else
                    EPtopos.firstRow.latency2=newSample;
                end
            end
            if EPtopos.synch
                EPtopos.points(:)=newSample; %update all rows to new number (samples or bins)
                EPtopos.iPoints(:)=newSample; %update all rows to new number (samples or bins)
            else
                if EPtopos.complexData && rem(theRow,2)
                    EPtopos.iPoints(rowCounter)=newSample; %update to new number (samples or bins)
                else
                    EPtopos.points(rowCounter)=newSample; %update to new number (samples or bins)
                end
            end
        end
    end
end

ep_showTopos(EPtopos.firstTime,EPtopos.lastTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeHz(src,eventdata,theRow)
%change the Hz settings.  Variables are in bins while the actual Hz are denoted by the freqName list.
global EPtopos

if theRow<0
    Hz1=get(EPtopos.handles.Hz1,'value');
    Hz2=get(EPtopos.handles.Hz2,'value');
    if Hz1>Hz2
        disp(['Error: First Hz cannot be after second Hz number'])
        return
    end
    if theRow==-1
        newBin=Hz1;
        oldHz=EPtopos.firstHz+(EPtopos.firstRow.Hz1-1)*EPtopos.binSize;
    else
        newBin=Hz2;
        oldHz=EPtopos.firstHz+(EPtopos.firstRow.Hz2-1)*EPtopos.binSize; %existing number number in Hz
    end
else
    if EPtopos.complexData
        rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(theRow/2);
    else
        rowCounter=(EPtopos.page-1)*EPtopos.perPage+theRow;
    end
    newBin=get(EPtopos.handles.Hz(theRow),'value');
    oldHz=EPtopos.firstHz+(EPtopos.freqs(rowCounter)-1)*EPtopos.binSize;
end

newHz=round((newBin-EPtopos.firstHz)/EPtopos.binSize)+1;
if (newHz ~= oldHz) || EPtopos.synch
    if (newHz >= EPtopos.firstHz) && (newHz <= EPtopos.lastHz)
        if (newBin > 0) && (newBin <= EPtopos.numHz)
            if theRow<0
                if theRow==-1
                    EPtopos.firstRow.Hz1=newBin;
                else
                    EPtopos.firstRow.Hz2=newBin;
                end
            else
                if EPtopos.synch
                    EPtopos.freqs(:)=newBin; %update all rows to new bins
                    EPtopos.iFreqs(:)=newBin; %update all rows to new bins
                else
                    if EPtopos.complexData && rem(theRow,2)
                        EPtopos.iFreqs(rowCounter)=newBin; %update to new bins
                    else
                        EPtopos.freqs(rowCounter)=newBin; %update to new bins
                    end
                end
            end
        end
    end
end

ep_showTopos(EPtopos.firstHz,EPtopos.lastHz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeChannel(src,eventdata,theRow)
global EPtopos

if EPtopos.complexData
    rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(theRow/2);
else
    rowCounter=(EPtopos.page-1)*EPtopos.perPage+theRow;
end
theChannel=EPtopos.nonRegChans(get(EPtopos.handles.channel(theRow),'value'));
if EPtopos.synch
    EPtopos.chans(:)=theChannel; %update all rows to new chans
    EPtopos.iChans(:)=theChannel; %update all rows to new chans
else
    if EPtopos.complexData && rem(theRow,2)
        EPtopos.iChans(rowCounter)=theChannel;
    else
        EPtopos.chans(rowCounter)=theChannel;
    end
end

ep_showTopos(EPtopos.firstTime,EPtopos.lastTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D3head(src,eventdata)
global EPtopos EPmain

disp('Using EEGlab function headplot to perform 3D head display.');

[row,col,~] = ind2sub(size(EPtopos.handles.topos.line),find(EPtopos.handles.topos.line==gco));

if isempty(row)
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==gco));
end

if isempty(row)
    tempVar=get(gco);
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==tempVar.Parent));
end

if EPtopos.complexData
    rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(row/2);
else
    rowCounter=(EPtopos.page-1)*EPtopos.perPage+row;
end

%check to see if spline file needs to be generated
if ~isempty(EPtopos.ced) && ~any(strcmp(EPtopos.ced,{'none','internal'}))
    [pathstr, name, ext] = fileparts(EPtopos.ced);
    CEDloc=which(EPtopos.ced);
    if isempty(CEDloc)
        name=EPtopos.dataName;
        CEDloc=[pwd filesep 'temp'];
    end
else
    name=EPtopos.dataName;
    CEDloc=[pwd filesep 'temp'];
end
[eeglabPath,a,b] = fileparts(which('eeglab'));

theMesh=which('mheadnew.mat');
if isempty(theMesh)
    theMesh=[eeglabPath filesep 'functions' filesep 'resources' filesep 'mheadnew.mat'];
    if ~exist(theMesh,'file')
        disp('Sorry, no mesh file available.  You may need to install EEGlab.')
        return;
    end
end

if isempty(which([name '.spl']))
    [pathstr2, name2, ext2] = fileparts(CEDloc);
    dataEloc=EPtopos.eloc(EPtopos.nonRegChans);
    chanList=zeros(length(dataEloc),1);
    for iChan=1:length(dataEloc)
        if ~isempty(dataEloc(iChan).X)
            chanList(iChan)=1;
            dataEloc(iChan).X=dataEloc(iChan).cX;
            dataEloc(iChan).Y=dataEloc(iChan).cY;
            dataEloc(iChan).Z=dataEloc(iChan).cZ;
        end
    end
    chanList=find(chanList);
    dataEloc=dataEloc(chanList);    
    ep_headplot('setup', dataEloc, [pathstr2 filesep name '.spl'],'meshfile',theMesh); %save spline file in same place as the ced file.
end

if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
    theRow=rowCounter;
else
    theRow=1;
end

if (rowCounter==1) && any(strcmp(EPtopos.type,{'time';'freq';'chan'}))
    pointsRange=EPtopos.firstRow.latency1:EPtopos.firstRow.latency2;
    binsRange=EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2;
    if EPtopos.complexData
        iPointsRange=EPtopos.firstRow.latency1:EPtopos.firstRow.latency2;
        iBinsRange=EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2;
    end
else
    pointsRange=EPtopos.points(rowCounter);
    binsRange=EPtopos.freqs(rowCounter);
    if EPtopos.complexData
        iBinsRange=EPtopos.iFreqs(rowCounter);
        iPointsRange=EPtopos.iPoints(rowCounter);
    end
end

if EPmain.view.rel(col) %if relational data
    if EPtopos.complexData && rem(row,2)
        theData=squeeze(mean(mean(EPtopos.totalImagData(EPtopos.iChans(rowCounter),iPointsRange,theRow,iBinsRange,:,col),2),4));
    else
        theData=squeeze(mean(mean(EPtopos.totalData(EPtopos.chans(rowCounter),pointsRange,theRow,binsRange,:,col),2),4));
    end
else
    if EPtopos.complexData && rem(row,2)
        theData=squeeze(mean(mean(EPtopos.totalImagData(:,iPointsRange,theRow,iBinsRange,1,col),2),4));
    else
        theData=squeeze(mean(mean(EPtopos.totalData(:,pointsRange,theRow,binsRange,1,col),2),4));
    end
end

theData(isnan(theData))=0;

figure
[hdaxis cbaraxis] = ep_headplot(theData,which([name '.spl']),'cbar',0,'maplimits',[EPtopos.plotMVmin EPtopos.plotMVmax],'meshfile',theMesh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sources=computeDipoles(topo,eloc,sampleSize)

ep_tictoc('begin');

disp('Using FieldTrip function ft_dipolefitting to perform dipole analysis.');

data=[];
cfg=[];

dataEloc=eloc;
chanList=zeros(length(dataEloc),1);
for iChan=1:length(dataEloc)
    if ~isempty(dataEloc(iChan).X)
        chanList(iChan)=1;
        dataEloc(iChan).X=dataEloc(iChan).cX;
        dataEloc(iChan).Y=dataEloc(iChan).cY;
        dataEloc(iChan).Z=dataEloc(iChan).cZ;
    end
end
chanList=find(chanList);
dataEloc=dataEloc(chanList);

data.elec.pnt=zeros(length(chanList),3);
for iChan=1:length(chanList)
    theChan=chanList(iChan);
    data.elec.label{iChan}=dataEloc(theChan).labels;
    data.elec.pnt(iChan,1)=dataEloc(theChan).X;
    data.elec.pnt(iChan,2)=dataEloc(theChan).Y;
    data.elec.pnt(iChan,3)=dataEloc(theChan).Z;
end
data.elec.pnt=data.elec.pnt(chanList,:);
topo=topo(chanList);
eloc=eloc(chanList);

data.fsample = 1000/(sampleSize);
data.avg  = topo;
%data.trial{1}  = topo;
data.var  = ones(length(topo),1);
data.time(1) = 1;
data.label = { eloc.labels };
data.dimord='rpt_chan_time';
%data.unmixing = diag(ones(length(data.label),1)); % workaround for FieldTrip bug

cfg.numdipoles  = 2;
cfg.symmetry    = 'x';
cfg.channel     = 'all';
cfg.gridsearch  = 'yes';
cfg.nonlinear   = 'yes';
cfg.headmodel     = which('standard_BESA.mat');
cfg.xgrid  = 'auto';
cfg.ygrid  = 'auto';
cfg.zgrid  = 'auto';
cfg.resolution = 10;
cfg.elec=data.elec;
if ~ft_hastoolbox('OPTIM')
    cfg.dipfit.optimfun = 'fminsearch';
end

[grid, cfg] = ft_prepare_sourcemodel(cfg);

[source] = ft_dipolefitting(cfg, data);

sources.momxyz(1,1:3)=source.dip.mom(1:3);
sources.momxyz(2,1:3)=source.dip.mom(4:6);

if (abs(source.dip.pos(1,1)) < .1) && (abs(source.dip.pos(2,1)) < .1)
    disp('Restarting dipole analysis with one dipole - two dipole solution blew up.');
    cfg.numdipoles  = 1;
    cfg.symmetry    = [];
    [source] = ft_dipolefitting(cfg, data);
    sources.momxyz=[];
    sources.momxyz(1,1:3)=source.dip.mom(1:3);
end

if isfield(source.dip,'rv')
    sources.posxyz=source.dip.pos;
    sources.rv=source.dip.rv;
    fprintf('Residual Variance: %04.2f\n',source.dip.rv*100);
    fprintf('Dipole 1 MNI Position: %04.2f %04.2f %04.2f\n',source.dip.pos(1,1),source.dip.pos(1,2),source.dip.pos(1,3));
    normDist=norm(source.dip.mom(1:3));
    fprintf('Dipole 1 Orientation: %04.2f %04.2f %04.2f\n',source.dip.mom(1)/normDist,source.dip.mom(2)/normDist,source.dip.mom(3)/normDist);
    fprintf('Dipole 1 Amplitude: %06.2f\n',normDist);
    if size(source.dip.pos,1)>1
        fprintf('Dipole 2 MNI Position: %04.2f %04.2f %04.2f\n',source.dip.pos(2,1),source.dip.pos(2,2),source.dip.pos(2,3));
        normDist=norm(source.dip.mom(4:6));
        fprintf('Dipole 2 Orientation: %04.2f %04.2f %04.2f\n',source.dip.mom(4)/normDist,source.dip.mom(5)/normDist,source.dip.mom(6)/normDist);
        fprintf('Dipole 2 Amplitude: %06.2f\n',normDist);
    end
    if exist('DrawDipole','file')
        dipLoc=[source.dip.pos(1,1),source.dip.pos(1,2),source.dip.pos(1,3)];
        dipOri=[source.dip.mom(1)/normDist,source.dip.mom(2)/normDist,source.dip.mom(3)/normDist];
        if size(source.dip.pos,1)>1
            dipLoc(2,:)=[source.dip.pos(2,1),source.dip.pos(2,2),source.dip.pos(2,3)];
            dipOri(2,:)=[source.dip.mom(4)/normDist,source.dip.mom(5)/normDist,source.dip.mom(6)/normDist];
        end
        drawDipole
    end
else
    sources=[];
    msg{1}='Error: Source analysis has failed.';
    if ~ft_hastoolbox('OPTIM')
        msg{2}='Probably due to Optimization Toolbox not being installed.';
    end
    [msg]=ep_errorMsg(msg);
    return
end

ep_tictoc('end');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dipoles(src,eventdata)
global EPtopos EPmain

[row,col,lin] = ind2sub(size(EPtopos.handles.topos.line),find(EPtopos.handles.topos.line==gco));

if isempty(row)
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==gco));
end

if isempty(row)
    tempVar=get(gco);
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==tempVar.Parent));
end

if strcmp('FFT',EPmain.view.dataTransform) && ~EPmain.view.rel(col)
    msg{1}='Error: Dipole analysis only works with voltage data.';
    [msg]=ep_errorMsg(msg);
    return
end

if EPtopos.complexData
    rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(row/2);
else
    rowCounter=(EPtopos.page-1)*EPtopos.perPage+row;
end

if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
    theRow=rowCounter;
else
    theRow=1;
end

if (rowCounter==1) && any(strcmp(EPtopos.type,{'time';'freq';'chan'}))
    pointsRange=EPtopos.firstRow.latency1:EPtopos.firstRow.latency2;
    binsRange=EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2;
    if EPtopos.complexData
        iPointsRange=EPtopos.firstRow.latency1:EPtopos.firstRow.latency2;
        iBinsRange=EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2;
    end
else
    pointsRange=EPtopos.points(rowCounter);
    binsRange=EPtopos.freqs(rowCounter);
    if EPtopos.complexData
        iBinsRange=EPtopos.iFreqs(rowCounter);
        iPointsRange=EPtopos.iPoints(rowCounter);
    end
end

if EPmain.view.rel(col) %if relational data
    if EPtopos.complexData && rem(row,2)
        theData=squeeze(mean(mean(EPtopos.totalImagData(EPtopos.iChans(rowCounter),iPointsRange,theRow,iBinsRange,:,col),4),2));
    else
        theData=squeeze(mean(mean(EPtopos.totalData(EPtopos.chans(rowCounter),pointsRange,theRow,binsRange,:,col),4),2));
    end
else
    if EPtopos.complexData && rem(row,2)
        theData=squeeze(mean(mean(EPtopos.totalImagData(:,iPointsRange,theRow,iBinsRange,1,col),4),2));
    else
        theData=squeeze(mean(mean(EPtopos.totalData(:,pointsRange,theRow,binsRange,1,col),4),2));
    end
end

theData=squeeze(theData);
theData(isnan(theData))=0;

eloc=EPtopos.eloc(EPtopos.nonRegChans);

goodChans=find(~isnan(theData));

sources=computeDipoles(theData(goodChans),eloc(goodChans),EPtopos.sampleSize);

if ~isempty(sources)
    disp('Using EEGlab function dipplot to present dipole solution.');
    try
        dipplot(sources,'mri',which('avg152t1.mat'),'color', { 'g' 'b' },'coordformat','MNI','verbose','off');
    catch
        msg{1}='Error: There is something wrong with your EEGlab installation.';
        [msg]=ep_errorMsg(msg);
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function jackknife(src,eventdata)
global EPtopos EPtictoc

[row,col,lin] = ind2sub(size(EPtopos.handles.topos.line),find(EPtopos.handles.topos.line==gco));

if isempty(row)
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==gco));
end

if isempty(row)
    tempVar=get(gco);
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==tempVar.Parent));
end

if EPtopos.complexData
    rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(row/2);
else
    rowCounter=(EPtopos.page-1)*EPtopos.perPage+row;
end

if ~isfield(EPtopos,'jack')
    msg{1}='Error: Jack-knife dipole analysis only works with PCA dataset.';
    [msg]=ep_errorMsg(msg);
    return
end

if ~strcmp(EPtopos.jack(col).PCAmode,'spat')
    msg{1}='Error: Jack-knife dipole analysis only works with spatial or temporo-spatial PCA dataset.';
    [msg]=ep_errorMsg(msg);
    return
end

ep_tictoc('begin');

if (rowCounter==1) && any(strcmp(EPtopos.type,{'time';'freq';'chan'}))
    pointsRange=EPtopos.firstRow.latency1:EPtopos.firstRow.latency2;
    binsRange=EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2;
    if EPtopos.complexData
        iPointsRange=EPtopos.firstRow.latency1:EPtopos.firstRow.latency2;
        iBinsRange=EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2;
    end
else
    pointsRange=EPtopos.points(rowCounter);
    binsRange=EPtopos.freqs(rowCounter);
    if EPtopos.complexData
        iBinsRange=EPtopos.iFreqs(rowCounter);
        iPointsRange=EPtopos.iPoints(rowCounter);
    end
end

eloc=EPtopos.eloc;
if EPtopos.complexData && rem(row,2)
    topo=squeeze(mean(mean(EPtopos.totalImagData(:,iPointsRange,rowCounter,iBinsRange,1,col),4),2));
else
    topo=squeeze(mean(mean(EPtopos.totalData(:,pointsRange,rowCounter,binsRange,1,col),4),2));
end
sources=computeDipoles(topo,eloc,EPtopos.sampleSize);
dipColors{1}='r';
numJK=size(EPtopos.jack(col).FacPat,3);
for iJK=1:numJK
    ep_tictoc;if EPtictoc.stop;return;end
    disp(['Jacknife ' num2str(iJK) ' of ' num2str(numJK)]);
    topo=EPtopos.jack(col).FacPat(:,rowCounter,iJK);
    sources(end+1)=computeDipoles(topo,eloc,EPtopos.sampleSize);
    dipColors{end+1}='b';
end

EPtopos.jack(col).sources=sources;

if ~isempty(sources)
    disp('Using EEGlab function dipplot to present dipole solution.');
    try
        dipplot(sources,'mri',which('avg152t1.mat'),'color', dipColors,'coordformat','MNI','verbose','off','spheres','on','dipolelength',0,'summary','3d');
    catch
        msg{1}='Error: There is something wrong with your EEGlab installation.';
        [msg]=ep_errorMsg(msg);
        return
    end
else
    msg{1}='Error: No sources to plot.';
    [msg]=ep_errorMsg(msg);
    return
end

%jack-knife t-test for hemispheric main effect
%using statistical test presented by: Miller, J., Patterson, T., & Ulrich, R. (1998). Jackknife-based method for measuring LRP onset latency differences. Psychophysiology, 35(1), 99-115.

numJKsources=length(sources);
DipAmp=zeros(numJK,2);

numJKpairs=0;
for JK=1:numJKsources
    momxyz=EPtopos.jack(col).sources(JK).momxyz;
    if size(momxyz,1) == 2
        numJKpairs=numJKpairs+1;
        DipAmp(numJKpairs,1)=sqrt(momxyz(1,1)^2+momxyz(1,2)^2+momxyz(1,3)^2);
        DipAmp(numJKpairs,2)=sqrt(momxyz(2,1)^2+momxyz(2,2)^2+momxyz(2,3)^2);
    end
end

JD=DipAmp(1:numJKpairs,1)-DipAmp(1:numJKpairs,2);
JM=mean(JD);

SDjk=sqrt(((numJKpairs-1)/numJKpairs)*sum((JD-JM).^2));

if numJKpairs > 1
    disp(['The t-value for the amplitude of the two hemispheric dipoles is ' num2str(sum(JD)/SDjk) ' with ' num2str(numJKpairs-1) ' degrees of freedom.']);
    if sum(JD) > 0
        disp('Right hemisphere larger than left.');
    elseif sum(JD) < 0
        disp('Left hemisphere larger than right.');
    else
        disp('No difference between hemispheres.');
    end
    if ft_hastoolbox('STATS', 0, 1)
        disp(['It has a two-sided p-value of ' num2str(2*tcdf(-abs(sum(JD)/SDjk),numJKpairs-1)) '.']);
    end
else
    disp('Too few solutions with hemispheric dipoles to calculate hemispheric t-test.');
end

EPtopos.changed=1;

ep_tictoc('end');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rescaleFigures(src,eventdata)
global EPtopos EPmain

[row,col,lin] = ind2sub(size(EPtopos.handles.topos.line),find(EPtopos.handles.topos.line==gco));

if isempty(row)
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==gco));
end

if isempty(row)
    tempVar=get(gco);
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==tempVar.Parent));
end

if EPtopos.complexData
    rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(row/2);
else
    rowCounter=(EPtopos.page-1)*EPtopos.perPage+row;
end

if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
    theRow=rowCounter;
else
    theRow=1;
end

if (rowCounter==1) && any(strcmp(EPtopos.type,{'time';'freq';'chan'}))
    pointsRange=EPtopos.firstRow.latency1:EPtopos.firstRow.latency2;
    binsRange=EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2;
    if EPtopos.complexData
        iPointsRange=EPtopos.firstRow.latency1:EPtopos.firstRow.latency2;
        iBinsRange=EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2;
    end
else
    pointsRange=EPtopos.points(rowCounter);
    binsRange=EPtopos.freqs(rowCounter);
    if EPtopos.complexData
        iBinsRange=EPtopos.iFreqs(rowCounter);
        iPointsRange=EPtopos.iPoints(rowCounter);
    end
end

if EPmain.view.rel(col) %if relational data
    if EPtopos.complexData && rem(row,2)
        theData=squeeze(mean(mean(EPtopos.totalImagData(EPtopos.iChans(rowCounter),iPointsRange,theRow,iBinsRange,:,col),4),2));
    else
        theData=squeeze(mean(mean(EPtopos.totalData(EPtopos.chans(rowCounter),pointsRange,theRow,binsRange,:,col),4),2));
    end
else
    if EPtopos.complexData && rem(row,2)
        theData=squeeze(mean(mean(EPtopos.totalImagData(:,iPointsRange,theRow,iBinsRange,1,col),4),2));
    else
        theData=squeeze(mean(mean(EPtopos.totalData(:,pointsRange,theRow,binsRange,1,col),4),2));
    end
end

if ~any(theData)
    disp('There is no data to rescale to in this topoplot.')
    return
end

EPtopos.plotMVmin=min(min(min(theData)));
EPtopos.plotMVmax=max(max(max(theData)));

ep_showTopos(EPtopos.firstTime,EPtopos.lastTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%explore effects of reference scheme on the data
function referenceChan(src,eventdata)
global EPtopos EPmain

scrsz = EPmain.scrsz;

[row,col,lin] = ind2sub(size(EPtopos.handles.topos.line),find(EPtopos.handles.topos.line==gco));

if isempty(row)
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==gco));
end

if isempty(row)
    tempVar=get(gco);
    [row,col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==tempVar.Parent));
end

if EPtopos.complexData
    rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(row/2);
else
    rowCounter=(EPtopos.page-1)*EPtopos.perPage+row;
end

if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
    theRow=rowCounter;
else
    theRow=1;
end

if (rowCounter==1) && any(strcmp(EPtopos.type,{'time';'freq';'chan'}))
    pointsRange=EPtopos.firstRow.latency1:EPtopos.firstRow.latency2;
    binsRange=EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2;
    if EPtopos.complexData
        iPointsRange=EPtopos.firstRow.latency1:EPtopos.firstRow.latency2;
        iBinsRange=EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2;
    end
else
    pointsRange=EPtopos.points(rowCounter);
    binsRange=EPtopos.freqs(rowCounter);
    if EPtopos.complexData
        iBinsRange=EPtopos.iFreqs(rowCounter);
        iPointsRange=EPtopos.iPoints(rowCounter);
    end
end

if EPmain.view.rel(col) %if relational data
    if EPtopos.complexData && rem(row,2)
        theData=squeeze(mean(mean(EPtopos.totalImagData(EPtopos.iChans(rowCounter),iPointsRange,theRow,iBinsRange,:,col),4),2));
    else
        theData=squeeze(mean(mean(EPtopos.totalData(EPtopos.chans(rowCounter),pointsRange,theRow,binsRange,:,col),4),2));
    end
else
    if EPtopos.complexData && rem(row,2)
        theData=squeeze(mean(mean(EPtopos.totalImagData(:,iPointsRange,theRow,iBinsRange,1,col),4),2));
    else
        theData=squeeze(mean(mean(EPtopos.totalData(:,pointsRange,theRow,binsRange,1,col),4),2));
    end
end

if ~any(theData)
    disp('There is no data to depict to in this topoplot.')
    return
end

EPtopos.referenceFigure.refchan1=1;
EPtopos.referenceFigure.refchan2=1;
switch EPtopos.reference(col).type
    case 'REG'
        switch length(EPtopos.reference(col).current)
            case 0
                EPtopos.referenceFigure.reference=5;
            case 1
                EPtopos.referenceFigure.reference=1;
                EPtopos.referenceFigure.refchan1=EPtopos.reference(col).current(1);
                EPtopos.referenceFigure.refchan2=1;
            case 2
                EPtopos.referenceFigure.reference=2;
                EPtopos.referenceFigure.refchan1=EPtopos.reference(col).current(1);
                EPtopos.referenceFigure.refchan2=EPtopos.reference(col).current(2);
            otherwise
                EPtopos.referenceFigure.reference=5;
        end
    case 'AVG'
        EPtopos.referenceFigure.reference=3;
    case 'PAR'
        EPtopos.referenceFigure.reference=4;
    otherwise
        EPtopos.referenceFigure.reference=5;
end

EPtopos.handles.referenceFigure.figure = figure('Name', 'Explore effect of rereferencing', 'NumberTitle', 'off', 'MenuBar', 'none', 'Toolbar', 'figure', 'Position',[scrsz(1)+scrsz(3)/2 scrsz(2)+scrsz(4)/2 600 400]);
colormap jet;

% disableDefaultInteractivity(gca)

uicontrol('Style','text',...
    'String','Reference Scheme','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[5 70 125 20]);

EPtopos.handles.referenceFigure.reference= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
    'String',{'One Chan','Two Chans','Average','PARE','unknown'},...
    'CallBack',@toposRereference,...
    'Value',EPtopos.referenceFigure.reference,'Position',[10 50 100 20]);

EPtopos.handles.referenceFigure.refchan1= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
    'String',EPtopos.chanNames(EPtopos.nonRegChans),...
    'CallBack',@toposRereference,...
    'Value',EPtopos.referenceFigure.refchan1,'Position',[10 30 100 20]);

EPtopos.handles.referenceFigure.refchan2= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
    'String',EPtopos.chanNames(EPtopos.nonRegChans),...
    'CallBack',@toposRereference,...
    'Value',EPtopos.referenceFigure.refchan2,'Position',[10 10 100 20]);

if ~any(ismember([1 2],EPtopos.referenceFigure.reference))
    set(EPtopos.handles.referenceFigure.refchan1,'enable','off');
end
if ~any(ismember([2],EPtopos.referenceFigure.reference))
    set(EPtopos.handles.referenceFigure.refchan2,'enable','off');
end

EPtopos.handles.referenceFigure.axes=axes('Units','pixel','position',[200 20 360 360]);

nSides=20;
[sphereCoords, elecInSphere]=ep_sphereHead(nSides, EPtopos.eloc);
sphereValues=ep_interpolateHead(theData, elecInSphere, sphereCoords);
X=reshape(sphereCoords(:,1),nSides+1,nSides+1);
Y=reshape(sphereCoords(:,2),nSides+1,nSides+1);
Z=reshape(sphereCoords(:,3),nSides+1,nSides+1);
sphereValues=reshape(sphereValues,nSides+1,nSides+1);
EPtopos.handles.referenceFigure.sphere=surf(X,Y,Z,sphereValues);
set(EPtopos.handles.referenceFigure.sphere,'CDataSource','EPtopos.referenceFigure.sphereValues');
set(EPtopos.handles.referenceFigure.sphere,'FaceColor','interp');
theMax=max(abs([min(theData) max(theData)]));
set(EPtopos.handles.referenceFigure.axes,'CLim',[-theMax theMax]);
cmap=colormap(EPtopos.handles.referenceFigure.figure);
cmap(floor(size(cmap,1)/2)-1:floor(size(cmap,1)/2)+2,:)=1;
set(EPtopos.handles.referenceFigure.figure,'colormap',cmap);

EPtopos.referenceFigure.elecValues=theData;
EPtopos.referenceFigure.sphereValues=sphereValues;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%change the reference on the Topos rereference panel
function toposRereference(src,eventdata)
global EPtopos

EPtopos.referenceFigure.reference=get(EPtopos.handles.referenceFigure.reference,'Value');
EPtopos.referenceFigure.refchan1=get(EPtopos.handles.referenceFigure.refchan1,'Value');
EPtopos.referenceFigure.refchan2=get(EPtopos.handles.referenceFigure.refchan2,'Value');

if ~any(ismember([1 2],EPtopos.referenceFigure.reference))
    set(EPtopos.handles.referenceFigure.refchan1,'enable','off');
else
    set(EPtopos.handles.referenceFigure.refchan1,'enable','on');
end
if ~any(ismember([2],EPtopos.referenceFigure.reference))
    set(EPtopos.handles.referenceFigure.refchan2,'enable','off');
else
    set(EPtopos.handles.referenceFigure.refchan2,'enable','on');
end

switch EPtopos.referenceFigure.reference
    case 1
        refValue=EPtopos.referenceFigure.elecValues(EPtopos.referenceFigure.refchan1);
    case 2
        refValue=(EPtopos.referenceFigure.elecValues(EPtopos.referenceFigure.refchan1)+EPtopos.referenceFigure.elecValues(EPtopos.referenceFigure.refchan2))/2;
    case 3
        refValue=mean(EPtopos.referenceFigure.elecValues);
    case 4
        %first average reference
        refValue=mean(EPtopos.referenceFigure.elecValues);
        EPtopos.referenceFigure.elecValues=EPtopos.referenceFigure.elecValues-refValue;
        
        %then interpolate surface
        theData=EPtopos.referenceFigure.elecValues;
        nSides=20;
        [sphereCoords, elecInSphere]=ep_sphereHead(nSides, EPtopos.eloc);
        sphereValues=ep_interpolateHead(theData, elecInSphere, sphereCoords);
        X=reshape(sphereCoords(:,1),nSides+1,nSides+1);
        Y=reshape(sphereCoords(:,2),nSides+1,nSides+1);
        Z=reshape(sphereCoords(:,3),nSides+1,nSides+1);
        sphereValues=reshape(sphereValues,nSides+1,nSides+1);
        EPtopos.handles.referenceFigure.sphere=surf(X,Y,Z,sphereValues);
        set(EPtopos.handles.referenceFigure.sphere,'CDataSource','EPtopos.referenceFigure.sphereValues');
        set(EPtopos.handles.referenceFigure.sphere,'FaceColor','interp');
        theMax=max(abs([min(theData) max(theData)]));
        set(EPtopos.handles.referenceFigure.axes,'CLim',[-theMax theMax]);
        cmap=colormap(EPtopos.handles.referenceFigure.figure);
        % cmap(floor(size(cmap,1)/2)-1:floor(size(cmap,1)/2)+2,:)=1;
        % set(EPtopos.handles.referenceFigure.figure,'colormap',cmap);
        EPtopos.referenceFigure.sphereValues=sphereValues;
        
        %then sum up the surface values
        sphereValues=EPtopos.referenceFigure.sphereValues(2:end-1,2:end-1);
        sphereValues=[sphereValues(:); EPtopos.referenceFigure.sphereValues(1,1); EPtopos.referenceFigure.sphereValues(end,end)];
        refValue=mean(sphereValues);
    case 5
        refValue=0;
    otherwise
        refValue=0;
end
EPtopos.referenceFigure.elecValues=EPtopos.referenceFigure.elecValues-refValue;
EPtopos.referenceFigure.sphereValues=EPtopos.referenceFigure.sphereValues-refValue;
refreshdata(EPtopos.handles.referenceFigure.sphere,'caller');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Display two channels and plot sample-by-sample t-test comparison
function twoChan(src,eventdata)
global EPtopos EPmain EPdataset EPtictoc

if isempty(EPtopos.handles.twoChan.figure) || ~isgraphics(EPtopos.handles.twoChan.figure)
    scrsz = EPmain.scrsz;
    
    EPtopos.twoChan.sigTest=0;
    
    [EPtopos.twoChan.row,EPtopos.twoChan.col,lin] = ind2sub(size(EPtopos.handles.topos.line),find(EPtopos.handles.topos.line==gco));
    
    if isempty(EPtopos.twoChan.row)
        [EPtopos.twoChan.row,EPtopos.twoChan.col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==gco));
    end
    
    if isempty(EPtopos.twoChan.row)
        tempVar=get(gco);
        [EPtopos.twoChan.row,EPtopos.twoChan.col] = ind2sub(size(EPtopos.handles.topos.topo),find(EPtopos.handles.topos.topo==tempVar.Parent));
    end
    
    if EPtopos.complexData
        EPtopos.twoChan.rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(EPtopos.twoChan.row/2);
    else
        EPtopos.twoChan.rowCounter=(EPtopos.page-1)*EPtopos.perPage+EPtopos.twoChan.row;
    end
    
    if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
        theRow=EPtopos.twoChan.rowCounter;
    else
        theRow=1;
    end
    
    if EPmain.view.rel(EPtopos.twoChan.col) %if relational data
        if EPtopos.complexData && rem(EPtopos.twoChan.row,2)
            EPtopos.twoChan.theData=EPtopos.totalImagData(EPtopos.iChans(EPtopos.twoChan.rowCounter),EPtopos.iPoints(EPtopos.twoChan.rowCounter),theRow,EPtopos.iFreqs(EPtopos.twoChan.rowCounter),:,EPtopos.twoChan.col);
        else
            EPtopos.twoChan.theData=EPtopos.totalData(EPtopos.chans(EPtopos.twoChan.rowCounter),EPtopos.points(EPtopos.twoChan.rowCounter),theRow,EPtopos.freqs(EPtopos.twoChan.rowCounter),:,EPtopos.twoChan.col);
        end
    else
        if EPtopos.complexData && rem(EPtopos.twoChan.row,2)
            if isempty(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).timeNames) && ~isempty(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).freqNames)
                EPtopos.twoChan.theData=EPtopos.totalImagData(:,EPtopos.iPoints(EPtopos.twoChan.rowCounter),theRow,:,1,EPtopos.twoChan.col);
            else
                EPtopos.twoChan.theData=EPtopos.totalImagData(:,EPtopos.iPoints(EPtopos.twoChan.rowCounter),theRow,EPtopos.iFreqs(EPtopos.twoChan.rowCounter),1,EPtopos.twoChan.col);
            end
        else
            if isempty(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).timeNames) && ~isempty(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).freqNames)
                EPtopos.twoChan.theData=EPtopos.totalData(:,:,theRow,:,1,EPtopos.twoChan.col);
            else
                EPtopos.twoChan.theData=EPtopos.totalData(:,:,theRow,EPtopos.freqs(EPtopos.twoChan.rowCounter),1,EPtopos.twoChan.col);
            end
        end
    end
    
    if ~any(EPtopos.twoChan.theData)
        disp('There is no data to depict to in this topoplot.')
        return
    end
    
    EPtopos.handles.twoChan.figure = figure('Name', 'Compare two channels', 'NumberTitle', 'off', 'MenuBar', 'none', 'Toolbar', 'figure', 'Position',[scrsz(1)+scrsz(3)/2 scrsz(2)+scrsz(4)/2 600 400]);
    
    EPtopos.twoChan.chan1=EPtopos.chans(EPtopos.twoChan.rowCounter);
    eDist=zeros(length(EPtopos.nonRegChans),1);
    for iChan=1:length(EPtopos.nonRegChans)
        theChan=EPtopos.nonRegChans(iChan);
        if theChan==EPtopos.twoChan.chan1
            eDist(iChan)=inf;
        else
            %Y is left-right
            eDist(iChan)=norm([EPtopos.eloc(EPtopos.twoChan.chan1).X-EPtopos.eloc(theChan).X -EPtopos.eloc(EPtopos.twoChan.chan1).Y-EPtopos.eloc(theChan).Y EPtopos.eloc(EPtopos.twoChan.chan1).Z-EPtopos.eloc(theChan).Z]);
        end
        
    end
    [X I]=min(eDist);
    EPtopos.twoChan.chan2=EPtopos.nonRegChans(I);
    
    if isempty(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).timeNames) && ~isempty(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).freqNames)
        EPtopos.twoChan.plotForm='FFT';
    else
        EPtopos.twoChan.plotForm='VLT';
    end
    
    uicontrol('Style','text',...
        'String','Channel 1','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
        'ForegroundColor',EPmain.preferences.view.color(1).RGB,'Position',[15 380 60 20]);
    
    EPtopos.handles.twoChan.chan1=uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',EPtopos.chanNames(EPtopos.nonRegChans),...
        'Value',EPtopos.twoChan.chan1,...
        'Position',[80 380 100 20],...
        'CallBack',@twoChan);

    
    uicontrol('Style','text',...
        'String','Channel 2','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
        'ForegroundColor',EPmain.preferences.view.color(3).RGB,'Position',[185 380 60 20]);
    
    EPtopos.handles.twoChan.chan2=uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
        'String',EPtopos.chanNames(EPtopos.nonRegChans),...
        'Value',EPtopos.twoChan.chan2,...
        'Position',[250 380 100 20],...
        'CallBack',@twoChan);
    
    EPtopos.handles.twoChan.sigTest=uicontrol('Style','checkbox',...
        'String','Sample-by-sample t-test','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
        'Value',EPtopos.twoChan.sigTest,'Position',[355 380 200 20],'CallBack',@twoChan);
    
    if strcmp(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).dataType,'continuous') ||...
        (strcmp(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).dataType,'single_trial') && length(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).trialNames)==1) ||...
        (strcmp(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).dataType,'average') && length(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).subNames)==1)
            set(EPtopos.handles.twoChan.sigTest,'enable','off');
    end
end

EPtopos.twoChan.chan1=EPtopos.nonRegChans(get(EPtopos.handles.twoChan.chan1,'Value'));
EPtopos.twoChan.chan2=EPtopos.nonRegChans(get(EPtopos.handles.twoChan.chan2,'Value'));

EPtopos.twoChan.sigTest=get(EPtopos.handles.twoChan.sigTest,'Value');
if EPtopos.twoChan.sigTest
    
    disp('Will need to load the dataset to perform the sample-by-sample statistical test so there will be a delay.');
    
    cellList=find(strcmp(EPmain.view.theCells{EPmain.view.dataset(EPtopos.twoChan.col)}(EPmain.view.cell(EPtopos.twoChan.col)),EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).cellNames));
    if strcmp(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).dataType,'average')
        subList=find(strcmp(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).subTypes,'AVG'));
    else
        subList=EPmain.view.subject(EPtopos.twoChan.col);
    end
    if strcmp(EPtopos.plotForm,'FFT')
        freqList=[];
    elseif ~isempty(EPdataset.dataset(EPmain.view.dataset(EPtopos.twoChan.col)).freqNames)
        freqList=EPtopos.freqs(EPtopos.twoChan.col);
    else
        freqList=[];
    end
    if strcmp(EPtopos.plotForm,'factor')
        facList=EPtopos.twoChan.rowCounter;
    else
        facList=EPmain.view.factor(EPtopos.twoChan.col);
    end
    
    EPdata=ep_loadEPdataset(EPmain.view.dataset(EPtopos.twoChan.col));
    %convert virtual GAVEs to normal form so subject specs etc are available.
    EPdata=ep_combineData(EPdata,'convert',{[],[],[],[],[],[]},[],[],[]);
    if EPtictoc.stop;EPtictoc.stop=0;end
    if isempty(EPdata)
        return
    end
    EPdataIn=ep_selectData(EPdata,{[EPtopos.twoChan.chan1 EPtopos.twoChan.chan2],[],cellList,[],facList,freqList});
    if isempty(EPdataIn)
        disp('Aborting sample test.')
        return
    end
    [outputData, sampLat1, sampLat2, sampAmp1, sampAmp2]=ep_sampleTest(EPdataIn,'channels',cellList,cellList,subList,subList,'sample',.05,4,[],[],[],'t-test',[],[]);

    outputData=squeeze(outputData);
    outputData=outputData(:);
end

%single-trial, average, FFT, TFT

waveSize=.8;

switch EPtopos.twoChan.plotForm
    case 'VLT'
        numImages=length(find(EPmain.view.allTrials(1:EPmain.numColors) == 2));
        if (length(EPtopos.plotColors)-numImages) > 0
            imageSpace=EPmain.numColors;
        else
            imageSpace=numImages;
        end
        EPtopos.twoChan.handles.waves.hExpandedAxes=[];
        if numImages %if any erpimage
            imageCount=0;
            for iColor=1:EPmain.numColors
                if (EPmain.view.dataset(iColor) <= length(EPdataset.dataset)) && (EPmain.view.allTrials(iColor) == 2)
                    imageCount=imageCount+1;
                    trialList=find(EPtopos.colorIndex==iColor);
                    EPtopos.twoChan.handles.waves.hExpandedAxes(end+1,1) = axes('position',[.05 .10+(waveSize/imageSpace)*(imageSpace-imageCount) waveSize (waveSize/imageSpace)]);
                    EPtopos.twoChan.handles.waves.hExpandedAxes(end+1,1) = imagesc(EPtopos.firstTime:EPtopos.lastTime,1:length(trialList),squeeze(EPtopos.twoChan.theData([EPtopos.twoChan.chan1 EPtopos.twoChan.chan2],:,trialList,:))',[EPtopos.plotMVmin, EPtopos.plotMVmax]);
                    axis([EPtopos.firstTime EPtopos.lastTime 1 length(trialList)]);
                    line([0 0],[1 length(trialList)],'Color','black','LineWidth',EPmain.preferences.view.lineSize) %stimulus onset
                    if ~isempty(EPtopos.marker1)
                        line(repmat(EPtopos.marker1,2),[1 length(trialList)],'Color',[.5 .5 .5],'LineWidth',EPmain.preferences.view.lineSize);
                    end
                    if ~isempty(EPtopos.marker2)
                        line(repmat(EPtopos.marker2,2),[1 length(trialList)],'Color',[.5 .5 .5],'LineWidth',EPmain.preferences.view.lineSize);
                    end
                end
            end
        end
        if (length(EPtopos.plotColors)-numImages) > 0 %if there will be waveforms
            EPtopos.twoChan.handles.waves.hExpandedAxes(end+1,1) = axes('position',[.05 .10 waveSize waveSize*((EPmain.numColors-numImages)/EPmain.numColors)]);
            hold on
            if EPtopos.twoChan.sigTest
                breakList=sort([find(diff([0; (outputData>0); 0])<0)-1; find(diff([0; (outputData>0); 0])>0)]);
                if ~isempty(breakList)
                    theData1=squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan1,:,:,:));
                    theData2=squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan2,:,:,:));
                    theData1=theData1(:);
                    theData2=theData2(:);
                    for iSigArea=1:length(breakList)/2
                        theTimePoints=[breakList((iSigArea-1)*2+1):breakList(iSigArea*2)]';
                        if isscalar(theTimePoints)
                            EPtopos.twoChan.handles.waves.hLines=line(([theTimePoints theTimePoints]*EPtopos.spacing)+(EPtopos.firstTime-EPtopos.spacing),[theData1(theTimePoints) theData2(theTimePoints)],'LineWidth',EPmain.preferences.view.lineSize,'Color',[1 .5 .5]);
                        else
                            EPtopos.twoChan.handles.waves.hLines=patch(([theTimePoints; flip(theTimePoints)]*EPtopos.spacing)+(EPtopos.firstTime-EPtopos.spacing),[theData1(theTimePoints); theData2(flip(theTimePoints))],EPtopos.plotColorIndex(EPtopos.twoChan.col,:),'FaceColor','red','EdgeColor','none','FaceAlpha',.25);
                        end
                    end
                end
            end
            plot([EPtopos.firstTime:EPtopos.spacing:EPtopos.firstTime+(EPtopos.spacing*(EPtopos.numPoints-1))],squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan1,:,1,:)),'LineStyle',EPtopos.plotLineIndex{EPtopos.twoChan.col},'color',EPmain.preferences.view.color(1).RGB);
            plot([EPtopos.firstTime:EPtopos.spacing:EPtopos.firstTime+(EPtopos.spacing*(EPtopos.numPoints-1))],squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan2,:,1,:)),'LineStyle',EPtopos.plotLineIndex{EPtopos.twoChan.col},'color',EPmain.preferences.view.color(3).RGB);
            hold off
            axis([EPtopos.firstTime EPtopos.lastTime EPtopos.plotMVmin EPtopos.plotMVmax]);
            if EPtopos.direction ==2
                set(EPtopos.twoChan.handles.waves.hWave([EPtopos.twoChan.chan1 EPtopos.twoChan.chan2]),'YDir','reverse')
            end
            line([EPtopos.firstTime EPtopos.lastTime-EPtopos.sampleSize],[0 0],'Color','black','LineWidth',EPmain.preferences.view.lineSize) % zero line
            line([0 0],[0 EPtopos.plotMVmax],'Color','black','LineWidth',EPmain.preferences.view.lineSize) %stimulus onset
        end
        
        if ~isempty(EPtopos.marker1)
            try
                eval('line(repmat(EPtopos.marker1,2),[EPtopos.plotMVmin EPtopos.plotMVmax],''Color'',[.5 .5 .5],''LineWidth'',EPmain.preferences.view.lineSize)');
            catch
            end
        end
        if ~isempty(EPtopos.marker2)
            try
                eval('line(repmat(EPtopos.marker2,2),[EPtopos.plotMVmin EPtopos.plotMVmax],''Color'',[.5 .5 .5],''LineWidth'',EPmain.preferences.view.lineSize);');
            catch
            end
        end
    case 'FFT'
        EPtopos.twoChan.handles.waves.hExpandedAxes=[];
        EPtopos.twoChan.handles.waves.hExpandedAxes(end+1,1) = axes('position',[.05 .05 waveSize waveSize]);
        hold on
        if EPtopos.twoChan.sigTest
            breakList=sort([find(diff([0; (outputData>0); 0])<0)-1; find(diff([0; (outputData>0); 0])>0)]);
            if ~isempty(breakList)
                theData1=squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan1,:,:,:));
                theData2=squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan2,:,:,:));
                theData1=theData1(:);
                theData2=theData2(:);
                for iSigArea=1:length(breakList)/2
                    theTimePoints=[breakList((iSigArea-1)*2+1):breakList(iSigArea*2)]';
                    if length(theTimePoints) == 1
                        EPtopos.twoChan.handles.waves.hLines=line(([theTimePoints theTimePoints]*EPtopos.spacing)+(EPtopos.firstHz-EPtopos.spacing),[theData1(theTimePoints) theData2(theTimePoints)],'LineWidth',EPmain.preferences.view.lineSize,'Color',[1 .5 .5]);
                    else
                        EPtopos.twoChan.handles.waves.hLines=patch(([theTimePoints; flip(theTimePoints)]*EPtopos.spacing)+(EPtopos.firstHz-EPtopos.spacing),[theData1(theTimePoints); theData2(flip(theTimePoints))],EPtopos.plotColorIndex(EPtopos.twoChan.col,:),'FaceColor','red','EdgeColor','none','FaceAlpha',.25);
                    end
                end
            end
        end
        plot([EPtopos.firstHz+EPtopos.sampleSize:EPtopos.spacing:EPtopos.lastHz],squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan1,:,1,:)),'LineStyle',EPtopos.plotLineIndex{EPtopos.twoChan.col},'color',EPmain.preferences.view.color(1).RGB);
        plot([EPtopos.firstHz+EPtopos.sampleSize:EPtopos.spacing:EPtopos.lastHz],squeeze(EPtopos.twoChan.theData(EPtopos.twoChan.chan2,:,1,:)),'LineStyle',EPtopos.plotLineIndex{EPtopos.twoChan.col},'color',EPmain.preferences.view.color(3).RGB);
        hold off
        axis([EPtopos.firstHz EPtopos.lastHz EPtopos.plotMVmin EPtopos.plotMVmax]);
        line([EPtopos.firstHz+EPtopos.sampleSize EPtopos.lastHz],[0 0],'Color','black','LineWidth',EPmain.preferences.view.lineSize) % zero line
        if ~isempty(EPtopos.marker1)
            try
                eval('line(repmat(EPtopos.marker1,2),[EPtopos.plotMVmin EPtopos.plotMVmax],''Color'',[.5 .5 .5],''LineWidth'',EPmain.preferences.view.lineSize)');
            catch
            end
        end
        if ~isempty(EPtopos.marker2)
            try
                eval('line(repmat(EPtopos.marker2,2),[EPtopos.plotMVmin EPtopos.plotMVmax],''Color'',[.5 .5 .5],''LineWidth'',EPmain.preferences.view.lineSize);');
            catch
            end
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates figure window for topo grid
function topoMatrix(src, eventdata, topoRow, topoCol)
global EPtopos EPmain

EPtopos.topoGrid.type=1;
EPtopos.topoGrid.col=topoCol;

scrsz = EPmain.scrsz;
windowSize=min(scrsz(3:4));
mainPane=get(EPmain.handles.hMainWindow,'Position');

if EPtopos.complexData
    rowCounter=(EPtopos.page-1)*ceil(EPtopos.perPage/2)+ceil(topoRow/2);
else
    rowCounter=(EPtopos.page-1)*EPtopos.perPage+topoRow;
end

if any(strcmp(EPtopos.type,{'subject','trial','cell','factor'}))
    EPtopos.topoGrid.row=rowCounter;
else
    EPtopos.topoGrid.row=1;
end

if (rowCounter==1) && any(strcmp(EPtopos.type,{'time';'freq';'chan'}))
    pointsRange=EPtopos.firstRow.latency1:EPtopos.firstRow.latency2;
    binsRange=EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2;
    if EPtopos.complexData
        iPointsRange=EPtopos.firstRow.latency1:EPtopos.firstRow.latency2;
        iBinsRange=EPtopos.firstRow.Hz1:EPtopos.firstRow.Hz2;
    end
else
    pointsRange=EPtopos.points(rowCounter);
    binsRange=EPtopos.freqs(rowCounter);
    if EPtopos.complexData
        iBinsRange=EPtopos.iFreqs(rowCounter);
        iPointsRange=EPtopos.iPoints(rowCounter);
    end
end

if EPmain.view.rel(EPtopos.topoGrid.col) %if relational data
    if EPtopos.complexData && rem(EPtopos.topoGrid.row,2)
        theData=squeeze(mean(mean(EPtopos.totalImagData(EPtopos.iChans(rowCounter),iPointsRange,EPtopos.topoGrid.row,iBinsRange,:,EPtopos.topoGrid.col),4),2));
    else
        theData=squeeze(mean(mean(EPtopos.totalData(EPtopos.chans(rowCounter),pointsRange,EPtopos.topoGrid.row,binsRange,:,EPtopos.topoGrid.col),4),2));
    end
else
    if EPtopos.complexData && rem(EPtopos.topoGrid.row,2)
        theData=squeeze(mean(mean(EPtopos.totalImagData(:,iPointsRange,EPtopos.topoGrid.row,iBinsRange,1,EPtopos.topoGrid.col),4),2));
    else
        theData=squeeze(mean(mean(EPtopos.totalData(:,pointsRange,EPtopos.topoGrid.row,binsRange,1,EPtopos.topoGrid.col),4),2));
    end
end

EPtopos.topoGrid.show=[1:rowCounter];
EPtopos.topoGrid.matCols=round(sqrt(rowCounter));
EPtopos.topoGrid.matRows=ceil(rowCounter/EPtopos.topoGrid.matCols);

if ~any(theData)
    disp('There is no data to depict to in this topoplot.')
    return
end

EPtopos.handles.topoGrid.figure = figure('Name', 'Topo Grid', 'NumberTitle', 'off', 'MenuBar', 'none', 'Toolbar', 'figure',...
    'Position',[mainPane(1)+mainPane(3)+1 101 windowSize windowSize]);
colormap jet;

drawTopoGrid
set(EPtopos.handles.topoGrid.figure,'SizeChangedFcn',@resizeTopoGrid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%populates topo grid figure
function drawTopoGrid(src, eventdata)
global EPtopos EPmain EPdataset

topoGridPane=get(EPtopos.handles.topoGrid.figure,'Position');
if min(topoGridPane(3:4))<200
    disp('Figure window too small')
    return
end
windowHeight=topoGridPane(4);

EPtopos.topoGrid.matCols=round(sqrt((topoGridPane(3)/topoGridPane(4))*length(EPtopos.topoGrid.show)));
EPtopos.topoGrid.matRows=ceil(length(EPtopos.topoGrid.show)/EPtopos.topoGrid.matCols);

if EPtopos.topoGrid.type==1
    EPtopos.topoGrid.width=floor(min((topoGridPane(3)/EPtopos.topoGrid.matCols),(topoGridPane(4)/EPtopos.topoGrid.matRows)));
    EPtopos.topoGrid.height=EPtopos.topoGrid.width;
else
    EPtopos.topoGrid.height=floor((min(topoGridPane(3),topoGridPane(4))-20)/max([EPtopos.topoGrid.matCols EPtopos.topoGrid.matRows*1.4]));
    EPtopos.topoGrid.width=EPtopos.topoGrid.height*1.4;    
end

for iTopo=1:length(EPtopos.topoGrid.show)
    theTopo=EPtopos.topoGrid.show(iTopo);
    iRow=floor((iTopo-1)/EPtopos.topoGrid.matCols)+1;
    iCol=mod(iTopo-1,EPtopos.topoGrid.matCols)+1;
    
    EPtopos.handles.topos.topoGrid(iRow,iCol) = axes('Units','pixel','position',[20+EPtopos.topoGrid.width*(iCol-1) windowHeight-(EPtopos.topoGrid.height*iRow) EPtopos.topoGrid.width EPtopos.topoGrid.height]);
    if EPmain.view.rel(EPtopos.topoGrid.col) %if relational data
        if EPtopos.complexData && rem(iRow,2)
            theData=EPtopos.totalImagData(:,EPtopos.iPoints(theTopo),theTopo,EPtopos.iFreqs(theTopo),EPtopos.iChans(theTopo),EPtopos.topoGrid.col);
        else
            theData=EPtopos.totalData(:,EPtopos.points(theTopo),theTopo,EPtopos.freqs(theTopo),EPtopos.chans(theTopo),EPtopos.topoGrid.col);
        end
    else
        if EPtopos.complexData && rem(iRow,2)
            theData=EPtopos.totalImagData(:,EPtopos.iPoints(theTopo),theTopo,EPtopos.iFreqs(theTopo),1,EPtopos.topoGrid.col);
        else
            theData=EPtopos.totalData(:,EPtopos.points(theTopo),theTopo,EPtopos.freqs(theTopo),1,EPtopos.topoGrid.col);
        end
    end
    theData=squeeze(theData);
    theData(isnan(theData))=0;
    el_topoplot(theData, EPtopos.eloc(EPtopos.nonRegChans),'maplimits',[EPtopos.plotMVmin EPtopos.plotMVmax]);
    
    hcmenu = uicontextmenu;
    % Define the context menu items and install their callbacks
    item1 = uimenu(hcmenu, 'Label', 'delete', 'Callback', {@deleteTopo,theTopo});
    % Attach the context
    set(EPtopos.handles.topos.topoGrid(iRow,iCol),'UIContextMenu',hcmenu)
    theChildren=get(EPtopos.handles.topos.topoGrid(iRow,iCol),'Children');
    for iChild=1:length(theChildren)
        set(theChildren(iChild),'UIContextMenu',hcmenu)
    end
    
    switch EPtopos.type
        case 'subject'
            theLabel=EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).subNames{min(theTopo,length(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).subNames))};
        case 'trial'
            theLabel=num2str(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).trialNames(min(theTopo,length(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).trialNames))));
        case 'factor'
            theLabel=EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).facNames{min(theTopo,length(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).facNames))};
        case 'cell'
            theLabel=EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).cellNames{min(theTopo,length(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).cellNames))};
        case 'time'
            theLabel=EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).timeNames(min(theTopo,length(EPdataset.dataset(EPmain.view.dataset(EPtopos.theFirstColor)).timeNames)));
    end
    
    uicontrol('Style','text',...
    'String',theLabel,'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
    'Position',[20+EPtopos.topoGrid.width*(iCol-1)+(EPtopos.topoGrid.width/2)-25 windowHeight-(EPtopos.topoGrid.height*iRow)+EPtopos.topoGrid.height-20 50 20]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%redraws topo grid figure after size change
function resizeTopoGrid(src, eventdata)
global EPtopos

set(EPtopos.handles.topoGrid.figure,'Resize',0)
set(EPtopos.handles.topoGrid.figure,'SizeChangedFcn',[])

clf(EPtopos.handles.topoGrid.figure)

drawTopoGrid

set(EPtopos.handles.topoGrid.figure,'Resize',1)
set(EPtopos.handles.topoGrid.figure,'SizeChangedFcn',@resizeTopoGrid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%deletes topo from Topo Grid
function deleteTopo(src, eventdata, theTopo)
global EPtopos

if isscalar(EPtopos.topoGrid.show)
    return
end

EPtopos.topoGrid.show(find(EPtopos.topoGrid.show==theTopo))=[];

clf(EPtopos.handles.topoGrid.figure)
drawTopoGrid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Confirms closing of the Topos pane.
function closeTopo(src,event)
global EPmain EPtopos

if EPtopos.changed
    selection = questdlg('Close The Topos Window Without Saving Changes?', ...
        'Confirm Closure', ...
        'Yes','No','Yes');
    switch selection
        case 'Yes'
            delete(EPtopos.handles.topos.topoWindow)
        case 'No'
            return
    end
    EPtopos.done=1;
else
    %if not the current waves figure, let it close.
    if ishandle(EPtopos.handles.topos.topoWindow)
        delete(EPtopos.handles.topos.topoWindow)
    else
        delete(gcf)
    end
end

%enable the main pane
ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
ep_tictoc('done');

ep('startView')
