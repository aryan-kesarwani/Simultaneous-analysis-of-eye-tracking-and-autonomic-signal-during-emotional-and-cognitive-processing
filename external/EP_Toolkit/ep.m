function ep(varargin)
% ep - function ep(varargin);
% GUI interface for the EP Toolkit.
%

%History
%  by Joseph Dien (7/30/09)
%  jdien07@mac.com
%
%  bugfix 8/2/09 JD
%  EPdataset.dataset initialized as an empty string rather than as a cell, which was causing a crash when adding the
%  first dataset.
%
%  modified 8/3/09 JD
%  Consolidated EGI EGIS average and session entries in the file format list.
%
%  bugfix 8/11/09 JD
%  For artifact correction, bug resulted in EGI EGIS session format always being used as input format.  Thanks
%  to Grega Repov.
%
%  bugfix 8/27/09 JD
%  Location of windows appearing partly off screen on some systems.
%  ANOVA preferences pane crash fixed and crash when doing ANOVA with only one between group and add option on.  Thanks to Tim Curran.
%  File title missing from first ANOVA in output files.  Thanks to Tim Curran.
%
%  modified 8/30/09 JD
%  Added preferences for importing text files.
%
%  bugfix and modified 9/5/09 JD
%  Crash when generating spatial PCA.  Thanks to Tim Curran.
%  Move factors preference were not used for preprocessing. 
%  Preprocess controls now remember settings.
%  Single cell files option for preprocessing were not functional.
%  Added reference channel field to preprocessing pane to deal with ICA problems when processing mean mastoid data.
%  Blink file was not passed on to blink correction function for fileTemplate option.
%  Crash when saving preferences file on a non-Mac.
%  Dataset information like the cells were not updating when changing datasets in the Window Data Pane.
%
%  bugfix 9/10/09 JD
%  Fixed not saving file when averaging session files and outputting as EGIS format.
%  Fixed crash when preprocessing data and using autoTemplate for blinks.
%
%  bugfix 9/23/09 JD
%  Not finding EPwork directory when in the search path but not in the active directory.
%
%  bugfix 10/12/09 JD
%  Fixed crash when using single file mode with Read function.
%
%  modified & bugfix 10/28/09 JD
%  Added option to disable preprocessing figure for low memory situations.
%  Movement correction preferences not registering.
%
%  modified 10/31/09 JD
%  When reading multiple files, assume all have the same electrode location and montage information rather than
%  repeatedly asking for them.
%  Cached more information in EPdataset to avoid slowdowns in main pane and preprocess pane.
%
%  modified 10/31/09 JD
%  checks to see if EPwork cache is current and regenerates it if not.
%
%  bugfix 11/2/09 JD
%  Fixed crash when using single file mode with Read function.
%
%  bugfix 11/4/09 JD
%  Fixed crash when running PCA scree and no noise PCA was conducted.  Thanks to Jia Wu.
%
%  bugfix 11/12/09 JD
%  Fixed bug where Toolkit was saving all variables into EPwork cache, not just EPdataset, resulting in erratic problems
%  when old values for variables were loaded back into the workspace.
%
%  bugfix 12/3/09 JD
%  Added workaround for sporadic Matlab menu bug that was causing menu to disappear or crash when menu item that is selected is
%  same as before.
%
%  bugfix 12/10/09 JD
%  Fixed crash when viewing Window pane for a dataset with less than four cells under Matlab versions prior to 2008.
%
%  bugfix 1/8/10 JD
%  Fixed crash when entering between group name that is not three letters long.
%  Fixed failure due to Matlab bug to redraw View pane when changing the dataset for one of the colors.
%
%  bugfix & modified 1/17/10 JD
%  Fixed crash when examining a file in the View pane with less than four cells.
%  Modified ANOVA module to accommodate non-ERP data by adding support for "behavioral" keyword.
%
%  modified 1/28/10 JD
%  Can leave out subject and cell fields when reading files in single file mode
%  (will use default values instead and will assume all the subjects/cells are the same).
%  Can accomodate file formats such as Neuroscan .avg where channels can be deleted, not just zeroed out, if bad and
%  channels can appear in an order different than that of the .ced file.  When importing a file with missing channels,
%  will now add them in as zeroed out channels and will mark them as bad.  Will also reorder channels to match that of
%  the ced file.
%  The type field in the .ced files is now used since the EEGlab bug was apparently fixed and it is now functional.  The
%  type field must now be present in the .ced file and assumptions will no longer be made about which ones are REF or FID types.
%  Added ability to plot data from continuous data files by showing only one second at a time.
%  Epoch ms now displayed as from beginning of first sample to end of last sample.
%  Fixed bug where checkboxes and subject and cell fields on Read and Preprocess panes impossible to deselect.
%  Refreshes Edit pane when Done button is pressed so when Data Name is changed, the name in the list of datasets is updated.
%
%  bugfix and modified 2/28/10 JD
%  Average function now correctly checks to see if output file name already exists for EP format files.
% analysis fields no longer optional.
% Now has option to use file's default reference settings in the preprocess pane.
% When running ANOVAs with between subject factors, fixed the addition of grand averages corresponding to the levels of the factors to the dataset.
% Made autoPCA more memory efficient so it wouldn't run out of memory.
% Fixed crash in ANOVA function when level names of between group factors were of different lengths.
% Close files after running batch of ANOVAs to avoid "too many files open" error.
% When running ANOVAs with between subject factors using the output of the autoPCA option, fixed the addition of grand averages corresponding to the levels of the factors to the dataset.
% Corrected the peak channel and peak time point identification of factors made by the autoPCA option of the window function.
% Fixed bug that prevented one from fully deleting the contents of edit fields, such as "baseline" on the preprocessing pane.
% Fixed bug that prevented one from changing the contents of the edit fields on the postprocessing pane.
% Fixed inability to set minimum variance setting for autoPCA preference to anything other than zero.
% Added option to turn off adding montage information to EGIS file format files due to incompatibility issues with some
% versions of NetStation.
%
% bugfix 3/6/10 JD
% When running ANOVAs with between subject factors, fixed having only the first between group waveform being added to
% the datafile.
%
% modified 3/27/10 JD
% Added ability to save data using .set file format.
%
% bugfix 4/26/10 JD
% Fixed wrong peak chans and time points for AutoPCA when not two-step PCA.
% Added topos function to the View Pane.
%
% bugfix 5/5/10 JD
% Fixed crash when in Window Data pane using noTable mode (older versions of Matlab) and there are less than five
% subject specs.
% Fixed not changing table of cells and table of specs when in Window Data pane using noTable mode (older versions of
% Matlab).
% Fixed crash when using Topo button of View EEG pane and not all four colors are being used.
% Added option to use unrotated solution for PCAs.
%
% modified 5/22/09 JD
% Number of factors and title of PCA no longer being changed to blank when changing other PCA settings.
% Added support for EEGlab .study files.
%
% bugfix 6/5/10 JD
% Fixed crash when postprocessing.
% View function no longer crashes when average file erroneously has non-empty trial names field.
% Fixed files sometimes not being recognized as being selected when names are in uppercase.
%
% modified & bugfix 6/17/10 JD
% Marks when a file isn't saved yet.
% Changed the suffix of EP files to ".ept".
% Fixed crash when windowing adds regional channel and there is no electrode coordinate information (eloc).
% Fixed crash when conducting robust ANOVA and there are no spec columns in the ANOVA data file.
%
% bugfix 6/28/10 JD
% Fixed when no within factors for ANOVAs, pane indicated needed one within cell rather than zero wthin cells.
%
% bugfix 7/2/10 JD
% In order fix crash due to Matlab bug on Windows, put in brief pause after spec name requestor diaglog box.
%
% bugfix 7/23/10 JD
% Fixed view pane listing all the trial names of the dataset rather than just those for a specific cell, resulting in
% crashes when they were selected.
% Fixed crash when viewing Waves and data are all zero.
% Fixed grand averages, which are added when computing ANOVAs with between factors, being computed incorrectly (only
% affected waveforms, not the ANOVA results).
%
% bugfix 8/27/10 JD
% Fixed unable to output averages in EGIS format.
% 
%  modified 10/17/10 JD
%  Added support for saccadeTrials and saccadeOnset fields.
% 
%  modified 10/28/10 JD
%  Added support for saccade preprocessing preference settings.  Removed the timepoints option from the preprocessing
%  pane.
%
% bugfix 12/7/10 JD
% Fixed crash in Window function after switching to a dataset with fewer cells.
% Fixed crash in ANOVA function when performing ANOVA with no between group factors under Matlab 2008b.
%
% bugfix 1/5/11 JD
% Fixed crash in ANOVA function when used prior to data being loaded into the work set yet (due to badly initialized variable).
% Now putting preferences file and work directory at default user directory if old one cannot be found.
% 
% modified 1/18/11 JD
% Added support for selecting timepoints in preprocessing (previous implementation was non-functional).
% 
% modified 1/20/11 JD
% Added support for manually specifying EOG channels in the preferences.
%
% bugfix 2/3/11 JD
% Fixed crash when starting program on a computer with spaces in the default user path, as in "Documents and Settings"
%
% bugfix 2/4/11 JD
% Hopefully fixed crash when starting program with certain configurations of EPwork and EPprefs.
%
% bugfix 2/8/11 JD
% Hopefully fixed crashes from buggy reset of preferences.
%
% bugfix 2/10/11 JD
% Fixed crash when performing eyeblink correction due to change in preferences file.
% Fixed not using saccade preference settings.
% Fixed crashing when loading EPdataset files.
% 
% modified 1/20/11 JD
% Added commands to change or create EPwork directory from the EP menu.
% Only stores EPpref files in EPwork directories.
%
% modified 2/26/11 JD
% Added support for EGI Matlab files
%
% bugfix 4/18/11 JD
% Fixed output file from the Window function having the suffix ".txt.txt"
%
% bugfix 4/22/11 JD
% Fixed crash when seeking to form manual blink templates and the only available files have missing channel coordinates.
%
% bugfix 5/5/11 JD
%Fixed not saving reset preference values when preference file found to be out-of-date.
%
% modified 5/5/11 JD
% Added support for FFT analyses.
%
% bugfix 6/2/11 JD
% Fixed crash when userpath is blank and Matlab is version 2007 or presumably earlier.
%
% modified 6/7/11 JD
% Added support for specification of current reference in preprocessing pane.
%
% bugfix 6/20/11 JD
% Fixed not saving preference changes and now looks for the preferences file in the current working directory
%
% bugfix 7/8/11 JD
% Changed separation character for added grand average when computing robust ANOVAs from "|" to "_" as the former was
% causing crash on PCs when saving files to text format.
% 
% modified 1/24/12 JD
% Eliminated REF channel type and added reference field to keep track of original and current reference scheme.
%
% bugfix 1/27/12 JD
% Fixed Window pane not allowing channel groups to be chosen using PCA results when the data file has a regional
% channel.
% 
% modified 2/22/12 JD
% Noise and std fields no longer optional.  Now set to empty if not used.
%
% modified 4/11/12 JD
% Added support for 6th dimension of frequency.  Unlogs frequency data for computations.
%
% modified 4/23/12 JD
% Changed default for channel group windowing to collapsing channels first
% then taking measures.  Added preference option so can use original
% approach too of measuring first.
%
% modified 5/24/12 JD
% Added support for missing data when adding channels, cells, or subjects.
%
% modified 5/24/12 JD
% Eliminated "no ref" option for preprocessing as no longer needed and confusing.
%
% bugfix 6/4/12 JD
% Fixed regional channel waveform added to dataset when windowing not computed correctly 
% when there are more than one channel group and the last one to be edited by the channel group function is not the one being used. 
%
% modified 8/6/12 JD
% Added suport for coherence and phase locking measures.
%
% bugfix 9/3/12 JD
% Fixed crash when trying to add subject ANOVA waveform after computing robust ANOVA and subjects have been trimmed.
% Fixed crash when trying to add subject ANOVA waveform after computing robust ANOVA and data is not from PCA output.
%
% modified 10/16/12 JD
% Added option to set last row to be imported in text files.
%
% modified 1/11/13 JD
% Added option to do internal calculations of frequency data in either amplitude or power form.
%
% bugfix 1/11/13 JD
% Fixed erroneous error message and crash when trying to PCA average file where bad channels were dropped rather than replaced.
%
% modified 1/17/13 JD
% Added support for ERPlab .erp file format.
%
% bugfix 1/18/13 JD
% Fixed erroneous "labels" error message when trying to load .study file.
%
% bugfix 1/23/13 JD
% Fixed problem where information for expanding channels in waveform plot is lost under some circumstances.
%
% modified 1/28/13 JD
% Added manual editing option for artifact correction.
%
% modified 1/29/13 JD
% Added frequency-domain peak measures to windowing function.
%
% modified 1/30/13 JD
% When saving a dataset, name of dataset changes if new name is chosen for saved dataset.
%
% bugfix 2/1/13 JD
% Fixed crash in View pane under when changing a dataset under certain circumstances.
%
% modified and bugfix 2/1/13 JD
% Clearing volt, hz, and sample parameters in View pane no longer crashes.  If value is manually set, will not change
% until a new dataset is chosen or until the value is cleared (in which case it will be replaced with automatic value).
% Added support for reading/writing ERPlab datasets.
%
% bugfix 2/4/13 JD
% Fixed error when setting file format preferences to ERPlab files.
% Fixed failure to save preferences when save preferences button clicked.
% Fixed ms window information on Window pane so that it ranges from onset of sample to offset of sample rather than
% offset to offset (e.g., 4-4 ms rather than 0-4 ms for first sample) after changes to the windowing settings.
%
% modified 2/6/13 JD
% Allow View Scan function to operate on average files.
%
% modified 2/21/13 JD
% Added version number to preferences to allow for automatized updating processes.
%
% modified 4/1/13 JD
% Improved the controls for the Windowing pane.
%
% bugfix 4/2/13 JD
% Markers in waveform plots can be set at zero ms.
%
% modified 4/5/13 JD
% Addressed change to userpath function in Matlab 2013a.
%
% bugfix 4/12/13 JD
% Scaling of topos now obeys the values on the View pane.
%
% bugfix 4/24/13 JD
% Fixed choosing "auto" as Edit Mode setting in Preprocessing pane resulting in error message.
%
% modified 5/6/13 JD
% Implemented adaptive mean in Windows pane by adding width setting for peak measure.
%
% bugfix 5/8/13 JD
% Fixed detrend option in Preprocess Data pane not working.
%
% bugfix 5/9/13 JD
% Fixed Scan allowed to be activated if the scaling is from zero to zero, as with a bad trial, causing a crash.
%
% bugfix 5/18/13 JD
% Fixed Single File Mode in Preprocessing pane crashing when values entered.
%
% bugfix 5/20/13 JD
% Fixed Single File Mode in Preprocessing pane crashing.
%
% modified 5/20/13 JD
% Single-Trial Files from multiple subjects can now be selected using the Read pane's Single File Mode and read in as separate files.
%
% bugfix 6/18/13 JD
% Fixed windowing pane crashing when Herz bins changed.
%
% bugfix 8/20/13 JD
% Fixed Change Work Directory menu function not working.
%
%  modified 9/14/13 JD
%  Initializes mff at the outset so that global variables bug doesn't crash the toolkit down the line.
%
% modified 9/16/13 JD
% Added support for Biosemi bdf files
%
% bugfix 9/17/13 JD
% Fixed not making factor data available for setting channel groups when first in analysis set had different number of channels.
% 
% modified 10/10/13 JD
% Supports presence of non-EEG channels.
%
%  bugfix 11/1/13 JD
%  Fixes font sizes on Windows.
%  Fixed About EP Toolkit spawned new main window.
%
%  modified 11/6/13 JD
%  Scan and Waves functions can now present event markings.
%
%  modified 11/13/13 JD
%  Added Segment Data function.
%
%  modified 11/27/13 JD
%  Added fMRI artifact correction option to Preprocess data function using both EEGlab fMRIb plugin and AMRI EEG fMRI Toolbox.
%
%  modified 12/23/13 JD
%  Added windowing function option to window 'all' the channels.
%
%  bugfix 12/24/13 JD
%  Fixed between group ANOVAs not being calculated correctly when the
%  between group variable is not sorted alphabetically.
%
%  bugfix 1/12/14 JD
%  Workaround for Matlab bug periodically causing screen size to register as having zero size.
%
%  bugfix 1/22/14 JD
%  Added Trim Data option to the Segment Data function.
%
% modified 2/26/14 JD
% Added View function option to plot or erpimage all trials and all subjects.
%
%  bugfix 2/27/14 JD
%  Fixed PCA not correctly recognizing that a file has already undergone two-step PCA.
%
% bugfix 3/5/14 JD
% Fixed when using single file mode to read in single-trial files, all the resulting files are identical to the very first subject.
%
% modified 3/12/14 JD
% Changed MEG channel type to MGA (axial gradiometer) and MGP (planar gradiometer) and MGM (magnetometer)
%
% modified 3/16/14 JD
% Uses internal electrode coordinates provided by MFF and FIFF files.  Added elecPrefs.
%
% modified 3/19/14 JD
% Added support for saving data in Neuromag FIFF file format.
% Changed uses of "temp" as a variable name to "tempVar" due to other Matlab programmers often using it as a function
% name, resulting in collisions.
% Added option to window single-trial data.
% Windowing outputs cells in actual order rather than alphabetical order.
% Eliminated noTable option for old versions of Matlab.
%
% bugfix 3/19/14 JD
% Fixed segmenting table + button mirroring the first line of the table rather than the current settings above the table.
% Fixed segmenting table - button deleting all but second to last line rather than just the last line.
%
% bugfix 3/23/14 JD
% Fixed waveforms added during ANOVAs to correspond to trimmed cell means being all flat.
%
% bugfix 3/26/14 JD
% Fixed crash when running an ANOVA on a PCA dataset results in trimmed cell means being added to it.
%
% modified 4/7/14 JD
% Added "all" and "erpimage" options to the Factors list in View.
%
% bugfix 4/8/14 JD
% Fixed AutoPCA generating nothing but missing data values when maxcentroid and mincentroid measures chosen.
%
% modified 4/14/14 JD
% Added .covNum field.
%
% modified 4/24/14 JD
% Added coherence and phase-locking options and allows choosing either hanning or multi-taper methods for spectral measures.
%
% bugfix 4/24/14 JD
% Fixed 'all' option for Window Data leaving out last channel.
%
% bugfix 4/29/14 JD
% Fixed when spectral range changed in View function and dB or psd options are on, the  values are immediately further transformed.
%
% bugfix 5/13/14 JD
% Fixed crash when running an ANOVA on a windowed text file generated by autoPCA and the adds option is on. 
%
% bugfix 5/20/14 JD
% Fixed minimum and maximum voltages in View pane not reflecting correct values for single trial data.
%
% modified 5/20/14 JD
% Added support for BrainVision EEG files
%
% bugfix 5/29/14 JD
% Allow for manual windowing of PCA files rather than just autoPCA.
%
% bugfix 6/1/14 JD
% Fixed list of trials for single-trial data in View pane not correct.
%
% modified 6/3/14 JD
% Added support for SMI eye tracking files.
%
% modified 6/18/14 JD
% Added starts, ends, contains, and -follows- keywords to the Segment function.
%
% modified 6/19/14 JD
% Added support for sample-by-sample t-tests, including STS chanType.
%
% modified 7/16/14 JD
% Simplified keys code field structure.
%
% modified 8/12/14 JD
% Added support for PARE-corrected average reference.
%
% bugfix 8/27/14 JD
% Fixed crash for files where event values are all numbers rather than
% strings.
%
% modified 8/31/14 JD
% Added field to Transform panel indicating how much of a delay has been added to ERP latencies by one-pass filters.
% Added support for adding additional keys information to events in continuous data and trial specs in single-trial
% data.
%
% modified 9/4/14 JD
% Added delay field to segment function.
%
% modified 9/15/14 JD
% Added contents of type field to time-lock events list in segment function.
%
% modified 10/16/14 JD
% Passes screen size to ep_readData call.
%
% bugfix 10/21/14 JD
% Fixed cannot update the preferences for text event and SMI file suffixes.
%
% bugfix 10/24/14 JD
% Fixed crash when using Transform pane and rereference set to "none."
%
% bugfix 12/8/14 JD
% Fixed crash when selecting mff files on a PC.
%
% bugfix 12/14/14 JD
% Fixed crash in sampTest function when there are bad trials present in single-trial data.
%
% modified 12/23/14 JD
% Added convert option to Save function.
%
% bugfix 1/24/15 JD
% Fixed crash in Save function when name of dataset to be saved is changed.
%
% bugfix 3/25/15 JD
% Fixed crash when reading in multiple files.
%
% bugfix 5/29/15 JD
% Fixed crash in Segment pane when dataset has no events.
%
% modified 5/29/15 JD
% Changed min and max TopoPlot scale to be set by plotted data unless overriden.
%
% bugfix & modified 9/25/15 JD
% Fixed Window button on main pane not becoming active for single-trial
% data being present in the working set.
% Fixed crash when working set newly initialized and running an ANOVA on
% behavioral data with adds option activated.
% Fixed crash when ANOVA conducted on data with characters instead of
% numbers.  Now treated as missing data.
% Added capability to resegment single-trial data to reassign cell
% conditions.
% Added capability to specify OR criteria for a condition by just giving
% them the same name.
% Added trial specs for average files.
% Changed default head rotation for electrode coordinates of mff and fiff files to 90 degrees.
%
% modified 10/13/15 JD
% On the assumption that EGI users have largely migrated from GSN200 to
% Hydrocel nets, the default right mastoid channel has been changed from
% 101 to 100.
%
% bugfix 11/19/15 JD
% Corrected erroneous error message that degrees of freedom of robust ANOVA
% is calculated per multivariate rather than univariate approach.
%
% modified 11/25/14 JD
% Added cross-validation option to PCA.
%
% modified 12/10/15 JD
% Checks if p-value variability exceeds two standard deviations over the threshold.
%
% modified 12/18/15 JD
% On further thought, changed cross-validation option to apply FacPat rather than FacCof.
%
% bugfix 1/15/16 JD
% Fixed when clicking on Views channels, can only get expanded channel window if one clicks along the very top edge for FFT data.
% Fixed upper and lower amplitude/power being changed immediately after new values manually entered into the View panel for FFT data.
% Fixed min and max values displayed in Views pane divided by Hz bin rather than sqrt(Hz bin) when set to amplitudes.
% Now allows power scaled data to be displayed as amplitudes.
% Now handles complex FFT numbers.
%
% modified 1/22/16 JD
% Consolidated spectral unit controls so just four options (cmp, asd, psd, dB).
% If minimum value in View pane is zero, in dB scaling will set to -2 instead of ?inf to maintain useful range.
% Eliminated option to transform data in power rather than amplitude form to avoid potential confusion.
%
% modified 2/19/16 JD
% Fixed jack-knife not calculated for first step of PCAs.
%
% modified 3/8/16 JD
% Eliminated support for .power field.
%
% bugfix 3/21/16 JD
% Fixed PCA Woody option of SampleTest function not saving latencies and amplitudes to the correct trials.
%
% bugfix 5/22/16 JD
% Fixed crash when using average function and the working set is empty.
%
% bugfix 8/23/16 JD
% Fixed crash when changing the Trial Spec Model ofthe Latency-Lock procedure of the Average function.
%
% modified 9/16/16 JD
% Removed saccade and saccademin preferences.
%
% bugfix & modified 11/8/16 JD
% Added support for template Woody and for continuous data to Sample Test function.
% Added support for reading and writing subject spec text files.
% Fixed crash in segment function after changing line after -follows- to 'none'.
%
% modified 1/3/17 JD
% Added eyeTracker option to blink and saccade correction routines.
%
% bugfix 3/2/17 JD
% Fixed crash when switching to jitter-correct in Average Pane and there are no single-trial datasets in the working set.
%
% modified 4/19/17 JD
% When quiting, clear EPwork option leaves both the directory and the preferences file intact so there is no need to redo changes to the preference settings.
%
% bugfix 5/4/17 JD
% Fixed crash when setting the Mark fields in the View function.
%
% bugfix 5/19/17 JD
% Fixed crash when merging multiple subjects with single cell mode in Read function.
% Added option to flip the electrode locations for .mff and .fiff.
%
% bugfix & modified 6/19/17 JD
% Flexible segments implemented in Segmentation function.
% Added option to clear working set to File menu.
% Added support for averaging multi-file subjects.
% Now allows just one Hz bin of TFT data to be chosen for display using the View function.
% PCA now drops non-EEG channels and regional EEG channels.
% Added support for using sampleTest function with TFT data.
% Fixed error when applying sample method to average data in sampleTest function.
% Fixed Window pane crashing when first dataset in working set is not an average file.
% Fixed Window function saves both imaginary and real files for spectral data even when units are not set as complex.
% Scaling input boxes in View pane now allow up to four decimals to better accommodate power units for spectral data.
% Fixed spectral density conversion incorrectly applied to View pane controls for complex and amplitude scaling, resulting in the numbers being changed from what was input.
% Now allows just one sample of TFT data to be chosen for display using the View function.
%
% bugfix & modified 10/4/17 JD
% Added EMG correction option to Preprocessing function.
% Segment function user interface now providing correct default sample range for sampling rates other than 250 Hz.
% Segment function user interface now provides only unique set of values for integer trials specs.
%
% bugfix & modified 11/26/17 JD
% Added -precedes- crit option to the Segment function.
% Fixed crash when changing value of a -precedes- or -follows- crit to an event name that has no keys.
% Segment function user interface only includes events with identical Type and Value once in Event Lock list.
% Segment function user interface does not include TRSP in Event Lock list.
% Segment function user interface always provides "100" as an option in a list of TS-RT values.
% Fixed cross-validation option giving wrong results!  Also changed name to "cross-verification."
%
% modified 2/9/18 JD
% Added support for GFP plots and error bands.
% Added support for stdCM field.
%
% modified 2/23/18 JD
% When computing an add's std, no longer tries to combine std values.  Instead sets to zero.
% Added -all- and -erpimage- options for cells in View Waves function.
%
% bugfix 3/9/18 JD
% Fixed cells popupmenu listing all trials rather than just one of each name.
%
% bugfix 3/27/18 JD
% Fixed crash when performing scree on second step of a two-step PCA.
%
% modified 4/29/18 JD
% Added support for calculation of RT and accuracy summary scores during averaging.
% Added support for outputting behavioral data in the Window function.
% Added support for View Topos listing all trials or cells or subjects.
% Cell box for Window Data function now allows blanks so may edit freely to change order or to drop cells.  'none' no longer recognized as a keyword.
% Cell box for Window Data function no longer sorting cell names alphabetically for average files.
% Window Data pane settings no longer resetting every time one leaves the pane.
%
% bugfix 5/13/18 JD
% Fixed crash sometimes when starting Window Data pane.
%
% bugfix 5/18/18 JD
% Fixed failing to perform PCA scree test when noise field consists of zeros.
%
% bugfix & modified 6/5/18 JD
% Added preference option to turn off BV header encoding.
% Sorts batch files by name prior to running ANOVAs, averages, and preprocessing.
% Fixed not adding new cells as adds (if enabled) if new cell names specified during windowing.
% Fixed crash during ANOVA if a cell name is specified that is not present in the original data if adds are enabled.
% Fixed crash in Save function when clicking on table rather than using Convert mode.
% Don't remove star for unsaved data when save is cancelled for whatever reason.
% Added option to add marks for RT and selected events to View figures.
%
% bugfix & modified 8/6/18 JD
% Modified so that sample test function also works with combined (CMB) cells.
% Fixed crash in View function for frequency-domain data.
% Fixed missing Factor popupmenus for TFT factor data.
%
% modified 8/27/18 JD
% Added RIDE option to the Averaging function.
%
% bugfix 9/26/18 JD
% Fixed crash in ANOVA when there are missing spec values for between group factors.
%
% modified & bugfix 10/22/18 JD
% Added item-averaging option.
% Added support for reading E-Prime text output to add trial spec information during segmentation.
% Fixed not being able to batch run ANOVAs from an autoPCA when the peak channels changed from factor to factor.
% Fixed windowing function crashing with single-trial data.
% Fixed crash in SampleTest function when a dataset in the working set lacks electrode coordinates.
% Fixed crash in ANOVA when there are no between group factors.
% Made EP Toolkit robust to variations in the order of the .ced fields as EEGlab will generate ced files with different field orders 
% under some circumstances and this would previously result in scrambling of the electrode coordinate information.
%
% bugfix 11/9/18 JD
% Fixed file type control in the Average function not working properly.
%
% modified 12/13/18 JD
% Added support for trial names to single file mode.
% Added support for specifying sample rate for text data.
% No longer resets all the preference settings to default when a new preference field is added by an EP version.
%
% bugfix 12/17/18 JD
% Fixed error message in Preprocessing function when clicking on Run and Single File mode note checked.
%
% bugfix & modified 12/25/18 JD
% Fixed crash when starting up the EP Toolkit for the first time on a computer.
% Made fontsize smaller on linux.
%
% bugfix & modified 1/14/19 JD
% Fixed not adding between-subject ANOVA adds.  Also, now only provides waveform for the ANOVA channel.
% Optimized ANOVAadd in batch runs by only loading dataset when it changes.
% Added limited support for writing out data in mff file format.
%
% bugfix 2/12/19 JD
% Fixed bug wherein ANOVA adds waveforms not calculated correctly when grand average data is listed prior to any of the subject average data or combined cells before any of the single cells.
% Fixed crash in Read function when using Single File Mode with averaged non-EP file format data.
%
% bugfix 3/8/19 JD
% Fixed Segment function user interface including a blank when populating lists of potential spec names for datasets with no keys.
%
% modified 4/10/19 JD
% Added support for task level performance measures.
% Now takes takes DisplayScaleFactor setting into account when determining font size.
% Fixed crash in Window function when adding regional Adds to the dataset.
%
% bugfix 4/22/19 JD
% Fixed crash on startup on Matlab earlier than 2018a.
%
% modified 5/01/19 JD
% Added option to fix truncated EP header in continuous BV files.
%
% bugfix 5/9/19 JD
% Fixed crash when loading spec table in Segment function and there are no continuous data in the working set.
% Fixed dropping adds such as the grand average from an existing dataset in working set when adding ANOVAadds.
% Fixed crash when performing three-step PCA.
%
% bugfix 5/10/19 JD
% Fixed crash when the -all- option is used and then the dataset is changed and is then viewed with View Topos.
%
% bugfix & modified 7/7/19 JD
% Added EP reset button to upper right corner.
% Fixed crash when conducting sample test and there are cells listed after sample test cells.
% Improved default setting behavior of the View controls.
%
% bugfix 7/22/19 JD
% Extra carriage returns at end of ANOVA data files now ignored.
% Fixed not reading ANOVA text files properly if using Linux or OS X end-of-line convention (e.g., linefeeds).
% Fixed crash when running PCA.
% Fixed crash when running Window function on a non-PCA dataset after having run it on a PCA dataset without using AutoPCA.
% Fixed crash when using the same Window settings on an ANOVA file with more columns.
% Fixed crash in ANOVA function when Bonferroni field left empty.
% Fixed crash in ANOVA function when applied to non-factor data and adds option is being used.
% Fixed existing adds being deleted by Window adds option when Windowing.
%
% bugfix & modified 8/20/19 JD
% Fixed reset button not working properly for PCA function.
% Added t-map option to sampleTest.
% Added preferences for waveform figure linewidth and channel font size.
%
% modified 9/8/19 JD
% When reading or converting a batch of files, first sorted alphabetically.
%
% bugfix 9/29/19 JD
% Fixed error message during startup that EPmain.preferences.view.lineSize was in error and resetting preferences to default values.
% Fixed Segment function crashing upon entering pane.
% Fixed Template function crashing upon entering pane.
%
% bugfix 10/25/19 JD
% Fixed crash when applying autoPCA in Window function to single-trial PCA data.
%
% bugfix & modified 11/1/19 JD
% Added -rangesamps- functionality.
% eventlock field is now always only based on the .type field, with additional specs needed to utilize the .value field.
%
% modified 11/5/19 JD
% Added sessNums sessNames fields.
% Added option to compare between different sessions.
%
% modified 11/22/19 JD
% Fixed crash in Window function when leaving the pane, changing the cellnames in the chosen dataset, and then returning to the Window pane.
%
% modified 12/24/19 JD
% Upgraded support of std information by adding .covAVE and .GAVsubs fields and eliminating .std and .stdCM fields.
% May now drop SD and/or noise information to reduce memory and time requirements.
%
% bugfix & modified 1/13/20 JD
% Added support for multi-session data.
% Fixed crash when between group level names have more than one character.
% When loading E-Prime spec files, stopped assigning EPoffset to .onsetdelay field.
% Fixed Scan not displaying event line if the initial epoch did not have an event or leaving it unchanged if changing from an epoch with an event to one without.
% Fixed crash when using low-pass filter with sampleTest.
%
% bugfix & modified 3/23/20 JD
% Fixed crash in Read function when using Single File Mode.
% Fixed scree of random data for three-step PCA.
% Fixed crash when deleting dataset in working set that is cached in EPeeg.
% Blink and saccade file templates now default to canonical templates provided in the templates directory.
% Fixed crash in sample test pane when one of the datasets has electrode coordinates with missing values.
% Added support for BOSC analysis option to Transform, View, and Window functions.
% Added support for decomp option to PCA function.
% Correlation matrix type now allowed for ICA.
% Fixed crash in autoPCA with frequency PCA data.
% Fixed crash when switching to EGIS data in sampleTest.
%
% bugfix & modified 4/13/20 JD
% Applies sample tests to all channels, not just EEG channels.
% Added extended Infomax, JADE, SOBI, and fastICA rotations.
% Added preferences setting to set rotation type for blink correction.
% Fixed not saving EPdataset when file is large.
% Automatically regenerates corrupted EPdataset files.
% ep_saveEPdataset now handles replacing existing EPdataset entries.
%
% bugfix & modified 4/24/20 JD
% Fixed BV header fix option no longer appearing.
% Added support for up to eight colors in waveform figures.
% Added preference settings to control the colors of the waveforms.
% Added preference buttons to panes.
%
% bugfix & modified 4/28/20 JD
% Fixed crash when starting this version of EP Toolkit for the first time.
% Added mains noise correction option to Preprocessing function.
%
% bugfix & modified 9/22/20 JD
% Fixed only reading in the first file of a list of mff files.
% Added ICA options for CRD saccade and spike potential artifact corrections.
% Fixed save function adding new file with "-1" added to name rather than offering to replace an existing file.
% When changing name of dataset in working set to match that of saved file, no longer including file suffix.
% Added preference to turn off use of Matlab Distributed Computing Toolbox in EP Toolbox functions.
% Added eye artifact control to specify EMCP, MAAC, or single-step ICA options.
% Fixed EP Toolkit forgetting data in working set had been saved (and therefore no longer displaying asterisk in file lists).
% Added -responseCorr- and -responseErr- options.
% Removed preference option to rotate electrode coordinates when importing mff or eeglab files that have internal coordinates as no longer needed.
% Added control to Segment function to specify the ACC and RT trial specs.
% Added option to salvage NS4 datasets via simple binary files and event text files.
% Fixed EPeeg cache being corrupted after deleting a file from it, resulting in it needing to be reinitialized.
%
% modified 10/19/20 JD
% Added support for preferences to manually specify screen size and RAM information.
%
% bugfix & modified 11/8/20 JD
% Fixed cannot change settings in ANOVA window when the cells of the currently selected dataset are not listed in alphabetical order.
% Moved spike potential file template to separate file.
% Added support for preferences to manually specify screen position.
%
% bugfix & modified 12/25/20 JD
% For Segment function, if there are no trial specs, will instead offer the key names, if any, of the event lock event for the ACC and RT popup menus.
% Added cell to the single-file options in the Average function.
%
% modified 1/28/21 JD
% No longer adding regional channel add in Window function if only one channel was used.
% Fixed crash in Save function when in Convert mode.
%
% bugfix 2/7/21 JD
% Fixed not having separate preference settings for figure markings and the line widths settings were not having any effect.
%
% modified & bugfix 3/5/21 JD
% Moved option to fix mains noise from Preprocess function to Transform function and added preferences setting.
% Eliminated legacy 'factors' file type from Read function.
% PCA data can now be represented without factor vector compression.
% Fixed crashing on start up during certain cases.
% Fixed problems with Topos and Scan buttons after using all eight colors to display factor results.
% Added additional safeguards for initial EPwork setup.
%
% bugfix 7/22/21 JD
% Fixed EOG channels preference control not displaying the channel numbers with correct format and hence not able to read them when Done was hit again on a later occasion.
%
% modified & bugfix 8/11/21 JD
% Fixed error message when reading in a data file where the file type has been designated as grand_average.
% Average function now allows specs for RT and ACC to be directly specified.
% ACC codes can now be set in the Average preferences.
% Added expansion button to the Segment pane to allow for a fuller examination of the segmentation table.
% Added support for canonical 10-05 coordinates in eloc structure.
% Added montage popupmenu to enable specification of montages. 
% 
% modified & bugfix 10/18/21 JD
% Upgraded history field to provide more information on changes.
% Fixed window add for spectral data so that both real and imaginary components generated correctly.
%
% modified & bugfix 11/22/21 JD
% Added meanSME and meanNoise options to the Window function.
% Fixed crash in Segment function pane when there is a single-trial dataset in the working set that has trials with empty events cells.
%
% bugfix & modified 12/23/21 JD
% Fixed crash in ANOVA function when trying to handle no-header file with between-group variables.
% Fixed crash in Segment function when using the specs tab and the specs table is empty.
% Average function now disables the ACC and RT settings if the model has no such trial specs 
% and also the latest dataset is now chosen as the default model rather than the first.
% Added meanERA option to Window function.
%
% bugfix 3/9/22 JD
% Fixed crash with hi-pass filter option in Transform function.
%
% bugfix 5/6/22 JD
% Fixed problems with handling session labels with Read function in Single File Mode.
%
% bugfix 5/11/22 JD
% Fixed crash when running a batch of ANOVA files and the adds option is activated.
%
% modified 7/24/22 JD
% Added support for reading Matlab .mat files.
%
% modified 8/8/22 JD
% Added option to run either SAS or SPSS version of the Promax algorithm.
%
% modified 9/3/22 JD
% Fixed changes to bad trials spec in Preprocessing Preferences not being registered.
%
% bugfix 10/13/22 JD
% Fixed crash when batch reading or converting continuous files.
%
% bugfix 11/24/22 JD
% Fixed not allowing for smoothing to be adjusted for wavelet multiplication and wavelet convolution options for FFT and TFT transforms.
%
% bugfix 1/11/23 JD
% Fixed crash in PCA pane when Deep Learning Toolbox not installed.
%
% bugfix 10/1/23 JD
% Disables Trim button for data other than continuous to avoid crash.
%
% bugfix 11/3/23 JD
% Fixed crash when batching files via Transform function.
%
% bugfix 11/30/23 JD
% Fixed not finding files under some Matlab configurations when batching files via Transform function.
%
% bugfix 5/16/24 JD
% Fixed ignoring EGI montage values when preprocessing ept files with EGI data, resulting in incorrect identification of Cz and EOG channels and thus suboptimal artifact correction.
%
% modified 6/28/24 JD
% Added buttons for loading and saving the Window cells table.
% Added checkboxes to Save function to allow for multiple saves.
%
% bugfix 7/5/24 JD
% Fixed not clearing out old column names in ANOVA interface when loading in a newer dataset with fewer columns.
%
% bugfix 7/5/24 JD
% Fixed Save function checkboxes crashing and also further streamlining functionality.
%
% modified 8/14/24 JD
% Added buttons to ANOVA output html for selecting between cell means tables with SD, Cousineau-Moray SE, robust Cousineau-Moray CI, or plots with CI error bars.
%
% bugfix & modified 11/15/24 JD
% Fixed montage submenu not updating when import file format changed.
% May now specify suffix for output segmented files via Preferences setting and eliminated SMI suffix preference.
%
% bugfix 11/17/24 JD
% Fixed crash in Window pane when the selected dataset has no history field.
%
% bugfix & modified 3/29/25 JD
% Added lineStyle control for the waveform lines.
% Fixed sampTest cell name.
% Added color picker to View color preferences.
% Fixed crash when switching Save function to convert mode.
% ACC and RT default to TS-ACC and TS-RT if present.
% Popup menus of spec values are now alphabetically sorted.
% Fixed sampTest crash when selected dataset has numeric subject specs.
% Fixed sampTest crash with average file when grand averages are not at the end of the list.
% Added support for lastChange field in EPdataset.
% Fixed ANOVAadds not excluding gave subjects if not last and not generating trimmed average waveforms if no between group factors.
% Fixed ANOVAadds not computing which subjects to include correctly.
% Fixed Save not including sample-by-sample waveform results.
% Added support for virtual grand averages.
% Fixed WindowAdds addition of channels to PCA files.
% Fixed spatial PCA channels not being combined correctly for regional channels.
% Added support for Version 3 ANOVA text files, to ensure ANOVAadds are matched to the correct cells.
% Fixed crash when View Waves with TFT data and no events selected.
%
% bugfix & modified 5/9/25 JD
% Fixed crash when abort button hit for certain functions.
% Fixed EP version not updated in preferences, resulting in History record not being accurate.
% Added preference for electrode markings in topo maps.
%
% bugfix & modified 5/9/25 JD
% Fixed crash on startup of EP Toolkit under some conditions.
%
% bugfix 6/14/25 JD
% Fixed sampleTest dataset popupmenu sometimes not appearing if the working set includes datasets with which it cannot be used.
% Fixed crash when running sampleTest on virtual cells.
% Fixed dataset popupmenu in SampleTest function not allowing dataset to be changed.
% Fixed reset button not working on Matlab 2025a.
% Fixed crash at start if EPdataset does not have lastChange field and dataset does not have .time field.
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

global EPdataset EPmain EPoverview EPwaves EPscree EPchanGrp EPtopos EPmanualEdit EPtictoc EP_VER

refList={'ave ref','1-2 refs','CSD','default'};

if isempty(varargin)
    
    if ~exist('ep_initEP','file')
        disp(' ');
        disp('**************************************************************');
        disp('Error: The rest of the EP Toolkit is not on the path.  Make sure to use');
        disp('the ''Add with Subfolders'' button rather than the ''Add Folder'' button.');
        disp('**************************************************************');
        disp(' ');
        warndlg('See command window for warning message.','!! Warning !!')
        return
    end
    
    if ~isempty(EPmain)
        if isfield(EPmain,'mode') %second time through during initialization, after splash page.
            if ~any(strcmp(EPmain.mode,{'agree','disagree'})) ||  ~isfield(EPmain,'handles') || ~isfield(EPmain.handles,'hSPlash') || ~ishandle(EPmain.handles.hSPlash)
                EPmain=[];
            end
        else
            EPmain=[];
        end
    end
    
    if isempty(EPmain)
        if ismac
            EPmain.fontsize=10;
            EPmain.panesize=[200 500];
        elseif ispc
            EPmain.fontsize=8;
            EPmain.panesize=[200 500];
        else
            EPmain.fontsize=8; %linux
            EPmain.panesize=[205 500];
        end
        try
            s = settings;
            EPmain.fontsize=EPmain.fontsize*s.matlab.desktop.DisplayScaleFactor.PersonalValue; %adjust if this parameter has been set.
        catch
        end
        
%         return
    elseif strcmp(EPmain.mode,'disagree')
        close('About EP Toolkit')
        EPmain=[];
        return
    else
        close('About EP Toolkit')
        abortEP=ep_initEP;
        if abortEP
            return
        end
    end
    
    eval(['global EPdataset EPmain EPoverview EPwaves EPscree EPchanGrp EPtopos EPmanualEdit']); %restore the global variable linkages after the mff bug has eliminated them.
    
    EPmain.fileFormatReadList=ep_fileFormats;
    EPmain.fileFormatSaveList=ep_writeData;
    EPmain.montageList=ep_montage;
    EPmain.lineStyleList={'-';'--';':';'-.'};
    
    PathName=userpath;
    if ~isempty(PathName)
        if any(strcmp(PathName(end),{';',':'}))
            PathName=PathName(1:end-1);
        end
        if ismac
            colonPos=findstr(PathName,':');
            if ~isempty(colonPos)
                PathName=PathName(1:colonPos-1);
            end
        end
    end
    %set up work directory
    if isfield(EPmain,'mode') %not first time through
        disp('Loading summary information about datasets in working set (The more datasets in the working set, the longer this will take).');
    end
    if exist([pwd filesep 'EPwork' filesep 'EPdataset.mat'],'file')
         %if there is already a work directory in the active directory, use it.
        if isfield(EPmain,'mode')
            disp(['Loading from EPwork in the active directory: ' pwd])
            try
                eval(['tempVar=load(''EPwork' filesep 'EPdataset.mat'');']);
            catch ME
                if any(strcmp(ME.identifier,{'MATLAB:load:unableToReadMatFile','MATLAB:load:notBinaryFile'}))
                    tempVar=regenerateEPdataset([pwd filesep 'EPwork' filesep 'EPdataset.mat']);
                else
                    msg{1}='Error: Working directory has been corrupted.  Delete the contents of the EPwork directory, except for the EPprefs file, and restart the EP Toolkit.';
                    [msg]=ep_errorMsg(msg);
                    beep
                    return
                end
            end
            if isfield(tempVar,'EPdataset')
                EPdataset=tempVar.EPdataset;
            else
                EPdataset=tempVar;
            end
            clear tempVar;
        end
        EPdataset.EPwork=pwd;
    elseif exist(['~' filesep 'Documents' filesep 'EPwork' filesep 'EPdataset.mat'],'file')
         %on a Mac, work directories kept in this location
        if isfield(EPmain,'mode')
            disp(['Loading from EPwork in the home Documents directory.'])
            try
                eval(['tempVar=load(''~' filesep 'Documents' filesep 'EPwork' filesep 'EPdataset.mat'');']); %if there is already a work directory, use it.
            catch ME
                if any(strcmp(ME.identifier,{'MATLAB:load:unableToReadMatFile','MATLAB:load:notBinaryFile'}))
                    tempVar=regenerateEPdataset(['~' filesep 'Documents' filesep 'EPwork' filesep 'EPdataset.mat']);
                else
                    msg{1}='Error: Working directory has been corrupted.  Delete the contents of the EPwork directory, except for the EPprefs file, and restart the EP Toolkit.';
                    [msg]=ep_errorMsg(msg);
                    beep
                    return
                end
            end
            if isfield(tempVar,'EPdataset')
                EPdataset=tempVar.EPdataset;
            else
                EPdataset=tempVar;
            end
            clear tempVar;
        end
        EPdataset.EPwork=['~' filesep 'Documents'];
    elseif exist([PathName filesep 'EPwork' filesep 'EPdataset.mat'],'file')
         %if there is already a work directory in the home directory, use it.
        if isfield(EPmain,'mode')
            disp(['Loading from EPwork in the home directory: ' PathName])
            try
                eval(['tempVar=load(''' PathName filesep 'EPwork' filesep 'EPdataset.mat'');']);
            catch ME
                if any(strcmp(ME.identifier,{'MATLAB:load:unableToReadMatFile','MATLAB:load:notBinaryFile'}))
                    tempVar=regenerateEPdataset([PathName filesep 'EPwork' filesep 'EPdataset.mat']);
                else
                    msg{1}='Error: Working directory has been corrupted.  Delete the contents of the EPwork directory, except for the EPprefs file, and restart the EP Toolkit.';
                    [msg]=ep_errorMsg(msg);
                    beep
                    return
                end
            end
            if isfield(tempVar,'EPdataset')
                EPdataset=tempVar.EPdataset;
            else
                EPdataset=tempVar;
            end
            clear tempVar;
        end
        EPdataset.EPwork=PathName;
    else
        %otherwise, create a new EPwork at the userpath after reseting in case it was set to root, otherwise try in the active directory.
        if isempty(PathName)
            try
                userpath('reset')
                disp('Resetting userpath as it was set to root');
                PathName=userpath;
                if any(strcmp(PathName(end),{';',':'}))
                    PathName=PathName(1:end-1);
                end
            catch
                %userpath('reset') doesn't work for Matlab 2007 and presumably earlier
                disp('Setting userpath to the active directory.');
                PathName=pwd;
            end
        end
        if isfield(EPmain,'mode')
            EPdataset=[];
            EPdataset.EPwork=PathName;
            EPdataset.dataset=cell(0);
            eval(['save ''' PathName filesep 'EPwork' filesep 'EPdataset'' EPdataset']);
        else
            disp(['Creating EPwork in directory: ' PathName])
            mkdir([PathName filesep 'EPwork']);
        end
        EPdataset.EPwork=PathName;
    end
    if isfield(EPmain,'mode') && ~exist([EPdataset.EPwork filesep 'EPwork'],'dir')
        %try one last time to set up the EPwork file by putting it wherever was specified by the current EPdataset.EPwork
        disp(['The set up of the EPwork folder did not work somehow.  Trying one last time at: ' EPdataset.EPwork]);
        [status,msg,msgID]=mkdir([EPdataset.EPwork filesep 'EPwork']);
        if ~status
            disp(['Error: Failed to create EPwork directory due to ' msgID]);
            disp(msg)
            EPdataset=[];
        else
            EPdataset=[];
            EPdataset.dataset=cell(0);
            eval(['save ''' PathName filesep 'EPwork' filesep 'EPdataset'' EPdataset']);
            if ~exist([PathName filesep 'EPwork'],'dir')
                disp('Error: Failed to create EPdataset file inside EPwork folder.');
                EPdataset=[];
            end
        end
    end
    if isfield(EPmain,'mode') && ~isfield(EPdataset,'dataset')
        disp('Error: Unable to set up EPwork folder.  There is something funky about how your computer is set up.')
        disp('Maybe permissions problem with home directory or disk full?  Cannot run EP Toolkit until this problem is fixed.')
        return
    end
    if isfield(EPdataset,'dataset') && ~isempty(EPdataset.dataset) && (~isfield(EPdataset.dataset,'lastChange'))
        disp('Updating EPdataset records to include lastChange field.')
        for iDataset=1:length(EPdataset.dataset)
            eval(['load(''' PathName filesep 'EPwork' filesep EPdataset.dataset(iDataset).dataName '.mat'')']);
            if ~isempty(EPdata.history) && isfield(EPdata.history,'time')
                EPdataset.dataset(iDataset).lastChange=EPdata.history{end,1}.time;
            else
                EPdataset.dataset(iDataset).lastChange='';
            end
        end
        eval(['save ''' PathName filesep 'EPwork' filesep 'EPdataset'' EPdataset']);
    end
    if isfield(EPmain,'mode')
        disp('Done loading summary information.');
    end
    
    EPmain.maxColors=8;
    EPmain.numColors=4;
    EPmain.defaultColor(1).RGB=[0 0 1];
    EPmain.defaultColor(2).RGB=[1 0 0];
    EPmain.defaultColor(3).RGB=[0 1 0];
    EPmain.defaultColor(4).RGB=[0 0 0];
    EPmain.defaultColor(5).RGB=[0 .8 .9];
    EPmain.defaultColor(6).RGB=[.4 .1 1];
    EPmain.defaultColor(7).RGB=[1 .4 .1];
    EPmain.defaultColor(8).RGB=[1 0 1];

    %set up preferences
    if exist([EPdataset.EPwork filesep 'EPwork' filesep 'EPprefs.mat'],'file')
        checkEPdataset=EPdataset;
        tempVar=load([EPdataset.EPwork filesep 'EPwork' filesep 'EPprefs.mat']);
        if isfield(tempVar,'prefs')
            prefs=tempVar.prefs;
        end
        clear tempVar;
        if exist('EPdataset','var') && isfield(EPdataset,'dataset') && length(EPdataset.dataset) ~= length(checkEPdataset.dataset) %dealing with case where EPdataset was saved as part of prefs file due to bug.
            EPmain.preferences=[];
            resetPrefs;
            prefs=EPmain.preferences;
            eval(['save ''' EPdataset.EPwork filesep 'EPwork' filesep 'EPprefs'' prefs']);
            disp('EPprefs file is corrupted.  Resetting preferences file.');
            EPdataset=checkEPdataset;
        end
        EPmain.preferences=prefs;
    else
        resetPrefs;
        prefs=EPmain.preferences;
        eval(['save ''' EPdataset.EPwork filesep 'EPwork' filesep 'EPprefs'' prefs']);
    end

    [prefErr prefMissing]=checkPrefs;

    if ~isempty(prefErr) || ~isempty(prefMissing)
        if prefErr
            disp('Preferences are corrupt or out-of-date.  Resetting to default values.');
            EPmain.preferences=[];
        elseif prefMissing
            disp('Adding new preference values.');
        end
        resetPrefs;
        prefs=EPmain.preferences;
        eval(['save ''' EPdataset.EPwork filesep 'EPwork' filesep 'EPprefs'' prefs']);
    end

    if any(EPmain.preferences.advanced.monitor)
        EPmain.scrsz=EPmain.preferences.advanced.monitor;
    else
        scrsz = get(0,'ScreenSize');
        if max(scrsz) <2
            msg{1}='Error: Matlab is currently unable to determine the size of the monitor.  Please restart Matlab.';
            [msg]=ep_errorMsg(msg);
            return
        else
            EPmain.scrsz=scrsz;
        end
    end

    if ~isfield(EPmain,'mode')
        splash
        return
    end
    
    EPmain.tempColor=EPmain.preferences.view.color;
    EPmain.prefReturn='preferenceMain'; %pane to return to from preferences panes.
    
    varargin{1}='start';
    EPoverview=[];
    EPmain.handles.playback=[];
    EPmain.handles.hMainWindow =[];
    EPmain.mode='main';
    
    EPmain.segment=[];
    EPmain.segment.importFormat=EPmain.preferences.general.sessionImportFormat;
    EPmain.segment.outputFormat=EPmain.preferences.general.outputFormat;
    EPmain.segment.delay=0;
    EPmain.segment.importMontage=EPmain.preferences.general.defaultMontage;
    
    EPmain.average.importFormat=EPmain.preferences.general.sessionImportFormat;
    EPmain.average.fileType=2;
    EPmain.average.outputFormat=EPmain.preferences.general.outputFormat;
    EPmain.average.method=EPmain.preferences.average.method;
    EPmain.average.averageType=1;
    EPmain.average.freqMethod=1;
    EPmain.average.smoothing=1;
    EPmain.average.minLatency=[];
    EPmain.average.maxLatency=[];
    EPmain.average.RTmethod=1;
    EPmain.average.minRT=100;
    EPmain.average.maxRT=2;
    EPmain.average.dropBad=1;
    EPmain.average.dropError=1;
    EPmain.average.dropTimeout=1;
    EPmain.average.sLatency=0;
    EPmain.average.sMin=0;
    EPmain.average.sMax=500;
    EPmain.average.c1Min=100;
    EPmain.average.c1Max=900;
    EPmain.average.c2Min=100;
    EPmain.average.c2Max=900;
    EPmain.average.rMin=-400;
    EPmain.average.rMax=400;
    EPmain.average.cCutoff=6;
    EPmain.average.cLag=300;
    EPmain.average.peakPolarity=1;
    EPmain.average.dropEvents=1;
    EPmain.average.dropSD=3;
    EPmain.average.dropNoise=0;
    EPmain.average.subject='';
    EPmain.average.session='';
    EPmain.average.cell='';
    EPmain.average.codeCorrect=EPmain.preferences.average.codeCorrect;
    EPmain.average.codeError=EPmain.preferences.average.codeError;
    EPmain.average.codeTimeout=EPmain.preferences.average.codeTimeout;
    EPmain.average.importMontage=EPmain.preferences.general.defaultMontage;
    
    EPmain.transform.importFormat=EPmain.preferences.general.importFormat;
    EPmain.transform.fileType=3;
    EPmain.transform.outputFormat=EPmain.preferences.general.outputFormat;
    EPmain.transform.reference=EPmain.preferences.transform.reference;
    EPmain.transform.refChan1=EPmain.preferences.transform.refChan1;
    EPmain.transform.refChan2=EPmain.preferences.transform.refChan2;
    EPmain.transform.baselineStart=EPmain.preferences.transform.baselineStart;
    EPmain.transform.baselineEnd=EPmain.preferences.transform.baselineEnd;
    EPmain.transform.preStim=-EPmain.preferences.transform.baselineStart;
    EPmain.transform.delay=0;
    EPmain.transform.domain=1;
    EPmain.transform.method=1;
    EPmain.transform.smoothing=1;
    EPmain.transform.waveletwidth=7;
    EPmain.transform.filterPass=1;
    EPmain.transform.filterType=2;
    EPmain.transform.filter1=[];
    EPmain.transform.filter2=[];
    EPmain.transform.filterOrder=6;
    EPmain.transform.detrend=0;
    EPmain.transform.dataMode=4;
    EPmain.transform.whiteNoise=rand(1,250);
    EPmain.transform.BOSC.width=6;
    EPmain.transform.BOSC.threshold=.95;
    EPmain.transform.BOSC.duration=3;
    EPmain.transform.BOSCdataset=1;
    EPmain.transform.BOSC.cell=1;
    EPmain.transform.mainsFix=EPmain.preferences.transform.mainsFix;
    EPmain.transform.importMontage=EPmain.preferences.general.defaultMontage;
    
    EPmain.read.check=0;
%     EPmain.read.ced=0;
    EPmain.read.subject='';
    EPmain.read.cell='';
    EPmain.read.freq='';
    EPmain.read.trial='';
    EPmain.read.session='';
    EPmain.read.format=EPmain.preferences.general.importFormat;
    EPmain.read.fileType=3;   
    EPmain.read.numEEG=EPmain.preferences.general.numEEG;
    EPmain.read.importMontage=EPmain.preferences.general.defaultMontage;
    
    EPmain.preprocess.importFormat=EPmain.preferences.general.sessionImportFormat;
    EPmain.preprocess.outputFormat=EPmain.preferences.general.sessionOutputFormat;
    EPmain.preprocess.importMontage=EPmain.preferences.general.defaultMontage;
    EPmain.preprocess.timepoints='';
    EPmain.preprocess.baseline='';
    EPmain.preprocess.detrend=0;
    EPmain.preprocess.EMG=1;
    EPmain.preprocess.alpha=0;
    EPmain.preprocess.blinkTemplate=3;
    EPmain.preprocess.saccadeTemplate=1;
    EPmain.preprocess.channelMode=1;
    EPmain.preprocess.trialMode=1;
    EPmain.preprocess.check=0;
    EPmain.preprocess.subject='';
    EPmain.preprocess.cell='';
    EPmain.preprocess.trial='';
    EPmain.preprocess.session='';
    EPmain.preprocess.fileType=2;
    EPmain.preprocess.origReference='';
    EPmain.preprocess.origRefType=4;
    EPmain.preprocess.currReference='';
    EPmain.preprocess.currRefType=4;
    EPmain.preprocess.editMode=3;
    EPmain.preprocess.fMRI=0;
    EPmain.preprocess.SP=1;
    EPmain.preprocess.eogMode='MAAC';
    
    EPmain.view=[];
    
    EPmain.pca.facNum=0;
    EPmain.pca.name='';
    EPmain.pca.mode=EPmain.preferences.pca.mode;
    EPmain.pca.rotation=EPmain.preferences.pca.rotation;
    EPmain.pca.rotopt=EPmain.preferences.pca.rotopt;
    EPmain.pca.decomp=1;
    EPmain.pca.rel=EPmain.preferences.pca.rel;
    EPmain.pca.loadings=EPmain.preferences.pca.loadings;
    EPmain.pca.parametric=false;
    EPmain.pca.crossVerifyPCA=[];
    EPmain.pca.rotFlag=0;
    EPmain.pca.theAlgorithm='SAS';
    EPscree=[];
    
    EPmain.sampleTest=[];
    
    EPmain.window.minFacVar=[];
    EPmain.window.FFTunits=4;
    EPmain.window.sampAdapt=0;
    EPmain.window.datasetName='';
    EPmain.window.undo=[];
    EPmain.window.lastChange=cell(0);
    
    EPmain.anova.numComps=0;
    EPmain.anova.data=[];
    EPmain.anova.method=1;
    EPmain.anova.undo=[];
%     EPmain.anova.epsilon=EPmain.preferences.anova.epsilon;
%     EPmain.anova.posthoc=EPmain.preferences.anova.posthoc;
    
    EPmain.save.format=EPmain.preferences.general.outputFormat;
    EPmain.save.SGLchan=1;
    EPmain.save.REGchan=1;
    EPmain.save.SGLcell=1;
    EPmain.save.CMBcell=1;
    EPmain.save.RAW=1;
    EPmain.save.AVG=1;
    EPmain.save.GAV=1;
    EPmain.save.SGLfac=1;
    EPmain.save.CMBfac=1;
    EPmain.save.batch=1;
    EPmain.save.EPheaderFix=0;
    EPmain.save.NS4fix=0;
    
    EPmain.save.check=0;
    EPmain.save.subject='';
    EPmain.save.cell='';
    EPmain.save.freq='';
    EPmain.save.trial='';
    EPmain.save.session='';
    EPmain.save.readFormat=EPmain.preferences.general.importFormat;
    EPmain.save.fileType=3;   
    
    EPchanGrp=[];
    
    EPtictoc.start=[];
    EPtictoc.step=1;
    EPtictoc.stop=0;
    
    ep_checkEPworkCache;
    
    if prefErr
        ep('savePrefs'); %save reset preference settings
    end
end

[EPpath, name, ext]=fileparts(which('ep.m'));
[a, b]=strtok(EPpath,filesep);
b=b(2:end);
EPver=ver(b);
if ~strcmp(EPver.Version,EP_VER)
    warndlg('It looks like you have updated your copy of the EP Toolkit without restarting it.','!! ERROR !!')
    return
end
EPmain.preferences.EPver=EPver.Version;

scrsz = EPmain.scrsz;

switch varargin{1}
    case 'start'
        if ~isempty(EPmain.handles.hMainWindow)
            if ishandle(EPmain.handles.hMainWindow)
                clf(EPmain.handles.hMainWindow)
                figure(EPmain.handles.hMainWindow)
            else
                disp('Sorry, the main window is missing.  You will need to restart the program.')
                return
            end
        else
            EPmain.handles.hMainWindow = figure('Name', 'EP Toolkit', 'NumberTitle', 'off', 'Position',[scrsz(1) scrsz(4)-550 EPmain.panesize(1) EPmain.panesize(2)], 'MenuBar', 'none');
            colormap jet;
            drawnow
        end
        
        figure(EPmain.handles.hMainWindow)
        EPtictoc.handles.T(1) = uicontrol('Style', 'pushbutton', 'FontSize',EPmain.fontsize,...
            'Tag','T1','Position', [175 495 20 5]);
        EPtictoc.handles.T(2) = uicontrol('Style', 'pushbutton', 'FontSize',EPmain.fontsize,...
            'Tag','T2','Position', [195 475 5 20]);
        EPtictoc.handles.T(3) = uicontrol('Style', 'pushbutton', 'FontSize',EPmain.fontsize,...
            'Tag','T3','Position', [175 470 20 5]);
        EPtictoc.handles.T(4) = uicontrol('Style', 'pushbutton', 'FontSize',EPmain.fontsize,...
            'Tag','T4','Position', [170 475 5 20]);
        % EPtictoc.handles.reset = uicontrol('Style', 'pushbutton', 'String', 'EP','FontSize',EPmain.fontsize,...
        %     'Tag','Reset','Position', [175 475 20 20], 'Callback', ['global EPtictoc;','if isempty(dbstack(4)),','EPtictoc.stop=0;','ep(''start'');','else EPtictoc.stop=1;end;',]);
         EPtictoc.handles.reset = uicontrol('Style', 'pushbutton', 'String', 'EP','FontSize',EPmain.fontsize,...
            'Tag','Reset','Position', [175 475 20 20], 'Callback', @resetEP);
       
        drawnow
        
        switch EPmain.mode
            case 'main'
                ep('startMain');
                
            case 'segment'
                ep('startSegment');
                
            case 'preprocess'
                ep('startPreprocess');
                
            case 'average'
                ep('startAverage');
                
            case 'transform'
                ep('startTransform');
                
            case 'read'
                ep('startRead');
                
            case 'edit'
                ep('startEdit');
                
            case 'view'
                ep('startView');
                
            case 'sampleTest'
                ep('startSampleTest');
                
            case 'PCA'
                ep('startPCA');
                
            case 'window'
                ep('startWindow');
                
            case 'ANOVA'
                ep('startANOVA');
                
            case 'save'
                ep('startSave');
                
            case 'preferenceMain'
                ep('startPreferenceMain');
                
            case 'preferenceGeneral'
                ep('startPreferenceGeneral');
                
            case 'preferencePreprocess'
                ep('startPreferencePreprocess');
                
            case 'preferenceAverage'
                ep('startPreferenceAverage');
                
            case 'preferenceTransform'
                ep('startPreferenceTransform');
                
            case 'preferenceView'
                ep('startPreferenceView');
                
            case 'preferencePCA'
                ep('startPreferencePCA');
                
            case 'preferenceWindow'
                ep('startPreferenceWindow');
                
            case 'preferenceANOVA'
                ep('startPreferenceANOVA');
                
            case 'preferenceAdvanced'
                ep('startPreferenceAdvanced');
                
            case 'preferenceRecords'
                ep('startPreferenceRecords');
            otherwise
                error([EPmain.mode ' is not a valid mode.']);
        end
        
    case 'startMain'
        
        uicontrol('Style','text',...
            'String',[EPmain.preferences.records.user '-' EPmain.preferences.records.lab '-' EPmain.preferences.records.institution '-' EPmain.preferences.records.project '-' EPmain.preferences.records.experiment],'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Tag','Records','Position',[20 475 150 20]);
        
        screeFigure=findobj('Name', 'ScreeWindow');
        if ~isempty(screeFigure)
            close(screeFigure)
        end
        
        chanGrpFigure=findobj('Name', 'Channel Group Window');
        if ~isempty(chanGrpFigure)
            close(chanGrpFigure)
        end
        
        editDataFigure=findobj('Name', 'EditData');
        if ~isempty(editDataFigure)
            close(editDataFigure)
        end
        
        set(EPmain.handles.hMainWindow,'Name', 'EP Toolkit');
        
        figure(EPmain.handles.hMainWindow)
        
        EPmain.handles.hMenu = uimenu('Label','File');
        uimenu(EPmain.handles.hMenu,'Label','About EP Toolkit','Callback',@splash);
        uimenu(EPmain.handles.hMenu,'Label','Preferences','Callback',['global EPmain;','EPmain.mode=''preferenceMain'';','ep(''start'');']);
        uimenu(EPmain.handles.hMenu,'Label','Change Work Directory','Callback',@changeWork);
        uimenu(EPmain.handles.hMenu,'Label','Create Work Directory','Callback',@createWork);
        uimenu(EPmain.handles.hMenu,'Label','Clear Work Directory','Callback',@clearWorkingSet);
        uimenu(EPmain.handles.hMenu,'Label','Quit','Callback',@quit);
        
        EPmain.handles.main.segment = uicontrol('Style', 'pushbutton', 'String', 'Segment','FontSize',EPmain.fontsize,...
            'Tag','Segment','Position', [20 450 100 30], 'Callback', ['global EPmain;','EPmain.mode=''segment'';','ep(''start'');']);
        
        EPmain.handles.main.preprocess = uicontrol('Style', 'pushbutton', 'String', 'Preprocess','FontSize',EPmain.fontsize,...
            'Tag','Preprocess','Position', [20 410 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preprocess'';','ep(''start'');']);
        
        EPmain.handles.main.transform = uicontrol('Style', 'pushbutton', 'String', 'Transform','FontSize',EPmain.fontsize,...
            'Tag','Transform','Position', [20 370 100 30], 'Callback', ['global EPmain;','EPmain.mode=''transform'';','ep(''start'');']);
        
        EPmain.handles.main.average = uicontrol('Style', 'pushbutton', 'String', 'Average','FontSize',EPmain.fontsize,...
            'Tag','Average','Position', [20 330 100 30], 'Callback', ['global EPmain;','EPmain.mode=''average'';','ep(''start'');']);
        
        EPmain.handles.main.read = uicontrol('Style', 'pushbutton', 'String', 'Read','FontSize',EPmain.fontsize,...
            'Tag','Read','Position', [20 290 100 30], 'Callback', ['global EPmain;','EPmain.mode=''read'';','ep(''start'');']);
        
        EPmain.handles.main.edit = uicontrol('Style', 'pushbutton', 'String', 'Edit','FontSize',EPmain.fontsize,...
            'Tag','Edit','Position', [20 250 100 30], 'Callback', ['global EPmain;','EPmain.mode=''edit'';','ep(''start'');']);
        
        if isempty(EPdataset.dataset)
            set(EPmain.handles.main.edit,'enable','off');
        end
        
        EPmain.handles.main.view = uicontrol('Style', 'pushbutton', 'String', 'View','FontSize',EPmain.fontsize,...
            'Tag','View','Position', [20 210 100 30], 'Callback', ['global EPmain;','EPmain.mode=''view'';','ep(''start'');']);
        
        if isempty(EPdataset.dataset)
            set(EPmain.handles.main.view,'enable','off');
        end
        
        EPmain.handles.main.PCA = uicontrol('Style', 'pushbutton', 'String', 'PCA','FontSize',EPmain.fontsize,...
            'Tag','PCA','Position', [20 170 100 30], 'Callback', ['global EPmain;','EPmain.mode=''PCA'';','ep(''start'');']);
        
        if isempty(EPdataset.dataset)
            set(EPmain.handles.main.PCA,'enable','off');
        end
        
        EPmain.handles.main.sampleTest = uicontrol('Style', 'pushbutton', 'String', 'Sample Test','FontSize',EPmain.fontsize,...
            'Tag','Sample Test','Position', [20 130 100 30], 'Callback', ['global EPmain;','EPmain.mode=''sampleTest'';','ep(''start'');']);
        
        if ~ft_hastoolbox('STATS', 0, 1)
            set(EPmain.handles.main.sampleTest,'enable','off');
        end
        
        EPmain.handles.main.window = uicontrol('Style', 'pushbutton', 'String', 'Window','FontSize',EPmain.fontsize,...
            'Tag','Window','Position', [20 90 100 30], 'Callback', ['global EPmain;','EPmain.mode=''window'';','ep(''start'');']);
        
        if ~isempty(EPdataset.dataset)
            isAve=false;
            for i=1:length(EPdataset.dataset)
                if any(strcmp(EPdataset.dataset(i).dataType,{'average','single_trial'})) && any(ismember(EPdataset.dataset(i).chanTypes,{'EEG','MGM','MGA','MGP'})) && any(strcmp(EPdataset.dataset(i).cellTypes,'SGL')) && (any(strcmp(EPdataset.dataset(i).facTypes,'SGL')) || isempty(EPdataset.dataset(i).facNames))
                    isAve=true;
                end
            end
            if ~isAve
                set(EPmain.handles.main.window,'enable','off');
            end
        else
            set(EPmain.handles.main.window,'enable','off');
        end
        
        EPmain.handles.main.ANOVA = uicontrol('Style', 'pushbutton', 'String', 'ANOVA','FontSize',EPmain.fontsize,...
            'Tag','ANOVA','Position', [20 50 100 30], 'Callback', ['global EPmain;','EPmain.mode=''ANOVA'';','ep(''start'');']);
        
        EPmain.handles.main.save = uicontrol('Style', 'pushbutton', 'String', 'Save','FontSize',EPmain.fontsize,...
            'Tag','Save','Position', [20 10 100 30], 'Callback', ['global EPmain;','EPmain.mode=''save'';','ep(''start'');']);
        
%         if isempty(EPdataset.dataset)
%             set(EPmain.handles.main.save,'enable','off');
%         end
        
        figure(EPmain.handles.hMainWindow)
        drawnow
        
    case 'startSegment'
        
        set(EPmain.handles.hMainWindow,'Name', 'Segment Data');
        figure(EPmain.handles.hMainWindow)
        
        contData=[];
        if isfield(EPmain.segment,'cellTable')
            contData=[];
            for iFile=1:length(EPdataset.dataset)
                if any(strcmp(EPdataset.dataset(iFile).dataType,{'continuous','single_trial'})) && ~isempty(EPdataset.dataset(iFile).events) && ~isempty(EPdataset.dataset(iFile).events{1})
                    contData(end+1)=iFile;
                end
            end
        end
        
        if isfield(EPmain.segment,'contData') && ~isequal(contData,EPmain.segment.contData)
            EPmain.segment=[];
        end
        
        if ~isfield(EPmain.segment,'cellTable') %just entering segment function
            EPmain.segment.contData=[];
            for iFile=1:length(EPdataset.dataset)
                if any(strcmp(EPdataset.dataset(iFile).dataType,{'continuous','single_trial'})) && ~isempty(EPdataset.dataset(iFile).events) && ~isempty(EPdataset.dataset(iFile).events{1})
                    EPmain.segment.contData(end+1)=iFile;
                end
            end
            
            EPmain.segment.subPane='cells';
            
            clearSpec(0) %initialize specs subpane
            
            EPmain.segment.numSpecs=6; %number of specs in the specs table.
            EPmain.segment.numFixed=9; %number of fixed columns at the start of the specs table, including the # column that is dropped during the segment function.
            EPmain.segment.relList={'=','~=','<','>','<=','>=','starts','ends','contains'};
            EPmain.segment.delay=0;
            EPmain.segment.cellTable=cell(1,EPmain.segment.numFixed+EPmain.segment.numSpecs*3);
            EPmain.segment.importFormat=EPmain.preferences.general.sessionImportFormat;
            EPmain.segment.outputFormat=EPmain.preferences.general.outputFormat;
            EPmain.segment.importMontage=EPmain.preferences.general.defaultMontage;
            EPmain.segment.flexible=0;
            EPmain.segment.flexEnd=1;
            EPmain.segment.flexLength=20;
            EPmain.segment.task='';
            EPmain.segment.ACC=1;
            EPmain.segment.RT=1;
            EPmain.segment.expand=0;
            
            if ~isempty(EPmain.segment.contData)
                changeSegmentDataset;
            end
        end

        
        EPmain.handles.segment.prefs = uicontrol('Style', 'pushbutton', 'String', '•','FontSize',EPmain.fontsize,...
            'Tag','prefs',...
            'Position', [5 485 10 10], 'Callback', ['global EPmain;','EPmain.prefReturn=''segment'';','EPmain.mode=''preferenceGeneral'';','ep(''start'');']);

        EPmain.handles.segment.cells = uicontrol('Style', 'pushbutton', 'String', 'cells','FontSize',EPmain.fontsize,...
            'Tag','cells',...
            'Position', [20 480 40 25], 'Callback', ['global EPmain;','EPmain.segment.subPane=''cells'';','ep(''start'');']);

        EPmain.handles.segment.specs = uicontrol('Style', 'pushbutton', 'String', 'specs','FontSize',EPmain.fontsize,...
            'Tag','specs',...
            'Position', [60 480 40 25], 'Callback', ['global EPmain;','EPmain.segment.subPane=''specs'';','ep(''start'');']);

        if strcmp(EPmain.segment.subPane,'cells')
            
            set(EPmain.handles.segment.cells,'enable','off');

            if isempty(EPmain.segment.contData)
                uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','No continuous or single-trial dataset with events in working set to serve as template.','Position',[5 400 150 60]);
            else
                EPmain.handles.segment.dataset = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{EPdataset.dataset(EPmain.segment.contData).dataName},...
                    'Value',find(EPmain.segment.dataset==EPmain.segment.contData),'Position',[5 440 150 40],...
                    'Tag','dataset',...
                    'Callback', @changeSegmentDataset);
                
                EPmain.handles.segment.flexible = uicontrol('Style','popupmenu',...
                    'String',{'Time Event','Flex Event'},'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Value',EPmain.segment.flexible+1,'Position',[5 440 90 20],...
                    'Tag','flexible',...
                    'Callback',['global EPmain;','EPmain.segment.flexible=get(EPmain.handles.segment.flexible,''Value'')-1;','ep(''start'')']);
                
                if ~isempty(EPmain.segment.eventLocks)
                    EPmain.handles.segment.eventValues = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.eventLocks,...
                        'Value',EPmain.segment.event,'Position',[105 440 100 20],...
                        'Tag','eventValues',...
                        'Callback', @changeSegmentEvent);
                else
                    EPmain.handles.segment.eventValues = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Tag','eventValues',...
                        'Position', [105 440 80 20],'enable','off');
                end
                
                if ~EPmain.segment.flexible
                    
                    uicontrol('Style','text',...
                        'String','samples','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                        'Position',[5 420 60 20]);
                    
                    uicontrol('Style','text',...
                        'String','ms','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                        'Position',[75 420 50 20]);
                    
                    uicontrol('Style','text',...
                        'String','delay ms','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                        'Position',[145 420 60 20]);
                    
                    EPmain.handles.segment.sampStart= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.sampStart,...
                        'TooltipString','First sample of epoch in relation to the event, where negative is before it.',...
                        'Tag','sampStart',...
                        'Position',[5 400 35 20],'Callback',@segmentSampStart);
                    
                    EPmain.handles.segment.sampEnd= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.sampEnd,...
                        'TooltipString','Last sample of epoch in relation to the event, where negative is before it.',...
                        'Tag','sampEnd',...
                        'Position',[40 400 35 20],'Callback',@segmentSampEnd);
                    
                    EPmain.handles.segment.msStart= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                        'String',round((EPmain.segment.sampStart)*(1000/EPdataset.dataset(EPmain.segment.dataset).Fs)),...
                        'TooltipString','Last ms of epoch in relation to the event, where negative is before it.',...
                        'Tag','msStart',...
                        'Position',[75 400 35 20],'Callback',@segmentSampStart);
                    
                    EPmain.handles.segment.msEnd= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                        'String',round((EPmain.segment.sampEnd)*(1000/EPdataset.dataset(EPmain.segment.dataset).Fs)),...
                        'TooltipString','First ms of epoch in relation to the event, where negative is before it.',...
                        'Tag','msEnd',...
                        'Position',[110 400 35 20],'Callback',@segmentSampEnd);
                else
                    
                    uicontrol('Style','text',...
                        'String','Flex end event','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                        'Position',[5 420 100 20]);
                    
                    if ~isempty(EPmain.segment.eventLocks)
                        EPmain.handles.segment.flexEnd = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                            'String',EPmain.segment.eventLocks,...
                            'Value',EPmain.segment.flexEnd,'Position',[105 420 100 20],...
                            'Tag','flexEnd',...
                            'Callback',['global EPmain;','EPmain.segment.flexEnd=get(EPmain.handles.segment.flexEnd,''Value'');','refresh']);
                    else
                        EPmain.handles.segment.flexEnd = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                            'Tag','flexEnd',...
                            'Position', [5 420 80 20],'enable','off');
                    end
                    
                    uicontrol('Style','text',...
                        'String','Samples','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                        'Position',[5 400 55 20]);
                    
                    EPmain.handles.segment.flexLength= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                        'String',sprintf('%d',EPmain.segment.flexLength),...
                        'TooltipString','Number of samples in the flexible length segment.',...
                        'Position',[60 400 35 20],...
                        'Tag','flexLength',...
                        'Callback',['global EPmain;','EPmain.segment.flexLength=str2num(get(EPmain.handles.segment.flexLength,''String''));','refresh']);
                    
                    uicontrol('Style','text',...
                        'String','delay','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                        'Position',[95 400 50 20]);
                    
                end
                
                EPmain.handles.segment.delay= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                    'String',sprintf('%d',EPmain.segment.delay),...
                    'TooltipString','Ms correction for nominal event timing (positive means delay).',...
                    'Position',[145 400 35 20],...
                    'Tag','delay',...
                    'Callback',['global EPmain;','EPmain.segment.delay=str2num(get(EPmain.handles.segment.delay,''String''));','refresh']);
                
                if ~isempty(EPmain.segment.critSpecNames)
                    EPmain.handles.segment.trialSpecNames(1) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.critSpecItemNames{1},...
                        'Value',EPmain.segment.trialSpec(1),'Position',[5 380 80 20],...
                        'Tag','trialSpecNames1',...
                        'Callback',@changeSegmentCritName);
                    EPmain.handles.segment.trialSpecRel(1) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.relList,...
                        'Value',EPmain.segment.trialSpecRel(1),'Position',[75 380 70 20],...
                        'Tag','trialSpecRel1',...
                        'Callback',['global EPmain;','EPmain.segment.trialSpecRel(1)=get(EPmain.handles.segment.trialSpecRel(1),''Value'');','refresh']);
                    EPmain.handles.segment.trialSpecVal(1) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.critSpecItemValues{1},...
                        'Value',EPmain.segment.trialSpecVal(1),'Position',[135 380 70 20],...
                        'Tag','trialSpecVal1',...
                        'Callback',@changeSegmentCritValue);
                else
                    EPmain.handles.segment.trialSpecNames(1) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [5 380 80 20],'enable','off');
                    EPmain.handles.segment.trialSpecRel(1) = uicontrol('Style', 'popupmenu', 'String', EPmain.segment.relList{EPmain.segment.trialSpecRel(1)},'FontSize',EPmain.fontsize,...
                        'Position', [75 380 70 20],'enable','off');
                    EPmain.handles.segment.trialSpecVal(1) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [135 380 70 20],'enable','off');
                end
                
                if ~isempty(EPmain.segment.critSpecNames)
                    EPmain.handles.segment.trialSpecNames(2) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.critSpecItemNames{2},...
                        'Value',EPmain.segment.trialSpec(2),'Position',[5 360 80 20],...
                        'Tag','trialSpec2',...
                        'Callback',@changeSegmentCritName);
                    EPmain.handles.segment.trialSpecRel(2) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.relList,...
                        'Value',EPmain.segment.trialSpecRel(2),'Position',[75 360 70 20],...
                        'Tag','trialSpecRel2',...
                        'Callback',['global EPmain;','EPmain.segment.trialSpecRel(2)=get(EPmain.handles.segment.trialSpecRel(2),''Value'');','refresh']);
                    EPmain.handles.segment.trialSpecVal(2) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.critSpecItemValues{2},...
                        'Value',EPmain.segment.trialSpecVal(2),'Position',[135 360 70 20],...
                        'Tag','trialSpecVal2',...
                        'Callback',@changeSegmentCritValue);
                else
                    EPmain.handles.segment.trialSpecNames(2) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [5 360 80 20],'enable','off');
                    EPmain.handles.segment.trialSpecRel(2) = uicontrol('Style', 'popupmenu', 'String', EPmain.segment.relList{EPmain.segment.trialSpecRel(2)},'FontSize',EPmain.fontsize,...
                        'Position', [75 360 70 20],'enable','off');
                    EPmain.handles.segment.trialSpecVal(2) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [135 360 70 20],'enable','off');
                end
                
                if ~isempty(EPmain.segment.critSpecNames)
                    EPmain.handles.segment.trialSpecNames(3) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.critSpecItemNames{3},...
                        'Value',EPmain.segment.trialSpec(3),'Position',[5 340 80 20],...
                        'Tag','trialSpecNames3',...
                        'Callback',@changeSegmentCritName);
                    EPmain.handles.segment.trialSpecRel(3) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.relList,...
                        'Value',EPmain.segment.trialSpecRel(3),'Position',[75 340 70 20],...
                        'Tag','trialSpecRel3',...
                        'Callback',['global EPmain;','EPmain.segment.trialSpecRel(3)=get(EPmain.handles.segment.trialSpecRel(3),''Value'');','refresh']);
                    EPmain.handles.segment.trialSpecVal(3) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.critSpecItemValues{3},...
                        'Value',EPmain.segment.trialSpecVal(3),'Position',[135 340 70 20],...
                        'Tag','trialSpecVal3',...
                        'Callback',@changeSegmentCritValue);
                    
                else
                    EPmain.handles.segment.trialSpecNames(3) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [5 340 80 20],'enable','off');
                    EPmain.handles.segment.trialSpecRel(3) = uicontrol('Style', 'popupmenu', 'String', EPmain.segment.relList{EPmain.segment.trialSpecRel(3)},'FontSize',EPmain.fontsize,...
                        'Position', [75 340 70 20],'enable','off');
                    EPmain.handles.segment.trialSpecVal(3) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [135 340 70 20],'enable','off');
                end
                
                if ~isempty(EPmain.segment.critSpecNames)
                    EPmain.handles.segment.trialSpecNames(4) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.critSpecItemNames{4},...
                        'Value',EPmain.segment.trialSpec(4),'Position',[5 320 80 20],...
                        'Tag','trialSpecNames4',...
                        'Callback',@changeSegmentCritName);
                    EPmain.handles.segment.trialSpecRel(4) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.relList,...
                        'Value',EPmain.segment.trialSpecRel(4),'Position',[75 320 70 20],...
                        'Tag','trialSpecRel4',...
                        'Callback',['global EPmain;','EPmain.segment.trialSpecRel(4)=get(EPmain.handles.segment.trialSpecRel(4),''Value'');','refresh']);
                    EPmain.handles.segment.trialSpecVal(4) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.critSpecItemValues{4},...
                        'Value',EPmain.segment.trialSpecVal(4),'Position',[135 320 70 20],...
                        'Tag','trialSpecVal4',...
                        'Callback',@changeSegmentCritValue);
                    
                else
                    EPmain.handles.segment.trialSpecNames(4) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [5 320 80 20],'enable','off');
                    EPmain.handles.segment.trialSpecRel(4) = uicontrol('Style', 'popupmenu', 'String', EPmain.segment.relList{EPmain.segment.trialSpecRel(4)},'FontSize',EPmain.fontsize,...
                        'Position', [75 320 70 20],'enable','off');
                    EPmain.handles.segment.trialSpecVal(4) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [135 320 70 20],'enable','off');
                end
                
                if ~isempty(EPmain.segment.critSpecNames)
                    EPmain.handles.segment.trialSpecNames(5) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.critSpecItemNames{5},...
                        'Value',EPmain.segment.trialSpec(5),'Position',[5 300 80 20],...
                        'Tag','trialSpecNames5',...
                        'Callback',@changeSegmentCritName);
                    EPmain.handles.segment.trialSpecRel(5) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.relList,...
                        'Value',EPmain.segment.trialSpecRel(5),'Position',[75 300 70 20],...
                        'Tag','trialSpecRel5',...
                        'Callback',['global EPmain;','EPmain.segment.trialSpecRel(5)=get(EPmain.handles.segment.trialSpecRel(5),''Value'');','refresh']);
                    EPmain.handles.segment.trialSpecVal(5) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.critSpecItemValues{5},...
                        'Value',EPmain.segment.trialSpecVal(5),'Position',[135 300 70 20],...
                        'Tag','trialSpecVal5',...
                        'Callback',@changeSegmentCritValue);
                    
                else
                    EPmain.handles.segment.trialSpecNames(5) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [5 300 80 20],'enable','off');
                    EPmain.handles.segment.trialSpecRel(5) = uicontrol('Style', 'popupmenu', 'String', EPmain.segment.relList{EPmain.segment.trialSpecRel(5)},'FontSize',EPmain.fontsize,...
                        'Position', [75 300 70 20],'enable','off');
                    EPmain.handles.segment.trialSpecVal(5) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [135 300 70 20],'enable','off');
                end
                
                if ~isempty(EPmain.segment.critSpecNames)
                    EPmain.handles.segment.trialSpecNames(6) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.critSpecItemNames{6},...
                        'Value',EPmain.segment.trialSpec(6),'Position',[5 280 80 20],...
                        'Tag','trialSpecNames6',...
                        'Callback',@changeSegmentCritName);
                    EPmain.handles.segment.trialSpecRel(6) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.relList,...
                        'Value',EPmain.segment.trialSpecRel(6),'Position',[75 280 70 20],...
                        'Tag','trialSpecRel6',...
                        'Callback',['global EPmain;','EPmain.segment.trialSpecRel(6)=get(EPmain.handles.segment.trialSpecRel(6),''Value'');','refresh']);
                    EPmain.handles.segment.trialSpecVal(6) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.critSpecItemValues{6},...
                        'Value',EPmain.segment.trialSpecVal(6),'Position',[135 280 70 20],...
                        'Tag','trialSpecVal6',...
                        'Callback',@changeSegmentCritValue);
                    
                else
                    EPmain.handles.segment.trialSpecNames(6) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [5 280 80 20],'enable','off');
                    EPmain.handles.segment.trialSpecRel(6) = uicontrol('Style', 'popupmenu', 'String', EPmain.segment.relList{EPmain.segment.trialSpecRel(6)},'FontSize',EPmain.fontsize,...
                        'Position', [75 280 70 20],'enable','off');
                    EPmain.handles.segment.trialSpecVal(6) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [135 280 70 20],'enable','off');
                end
                
                uicontrol('Style','text',...
                    'String','Task','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[5 260 60 20]);
                
                EPmain.handles.segment.task= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.task,...
                    'TooltipString','Optional task label.  Must end in + or -.',...
                    'Position',[5 240 60 20],...
                    'Tag','task',...
                    'Callback',['global EPmain;','EPmain.segment.task=get(EPmain.handles.segment.task,''String'');','refresh']);
                
                uicontrol('Style','text',...
                    'String','ACC','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[65 260 40 20]);
                
                if ~isempty(EPmain.segment.trialSpecNames) || ~isempty(EPmain.segment.stimSpecNames{EPmain.segment.event})
                    if isempty(EPmain.segment.trialSpecNames)
                        theList=EPmain.segment.stimSpecNames{EPmain.segment.event};
                    else
                        theList=EPmain.segment.trialSpecNames;
                    end
                    EPmain.handles.segment.ACC = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',theList,...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.segment.ACC,''Value'');','if tempVar ~=0,EPmain.segment.ACC=tempVar;end;','if isempty(tempVar),EPmain.segment.ACC=tempVar;end'],...
                        'Value',EPmain.segment.ACC,'Position',[65 240 65 20]);
                else
                    EPmain.handles.segment.ACC = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String','none',...
                        'Value',1,'Position',[65 240 65 20],'enable','off');
                end
                
                uicontrol('Style','text',...
                    'String','RT','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[125 260 40 20]);
                
                if ~isempty(EPmain.segment.trialSpecNames) || ~isempty(EPmain.segment.stimSpecNames{EPmain.segment.event})
                    if isempty(EPmain.segment.trialSpecNames)
                        theList=EPmain.segment.stimSpecNames{EPmain.segment.event};
                    else
                        theList=EPmain.segment.trialSpecNames;
                    end
                    EPmain.handles.segment.RT = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',theList,...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.segment.RT,''Value'');','if tempVar ~=0,EPmain.segment.RT=tempVar;end;','if isempty(tempVar),EPmain.segment.RT=tempVar;end'],...
                        'Tag','RT',...
                        'Value',EPmain.segment.RT,'Position',[130 240 65 20]);
                else
                    EPmain.handles.segment.RT = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String','none',...
                        'Value',1,'Position',[130 240 65 20],'enable','off');
                end
                
                EPmain.handles.segment.addCell = uicontrol('Style', 'pushbutton', 'String', '+','FontSize',EPmain.fontsize*2,...
                    'Tag','addCell',...
                    'Position', [10 220 40 20], 'Callback', @segmentAddCell);
                
                EPmain.handles.segment.delCell = uicontrol('Style', 'pushbutton', 'String', '-','FontSize',EPmain.fontsize*2,...
                    'Tag','delCell',...
                    'Position', [55 220 40 20], 'Callback', @segmentDelCell);
                
            end
            
            EPmain.handles.segment.load = uicontrol('Style', 'pushbutton', 'String', 'Load','FontSize',EPmain.fontsize,...
                'TooltipString','Load segmentation table.',...
                'Tag','load',...
                'Position', [100 220 40 20], 'Callback', @segmentLoad);
            
            EPmain.handles.segment.save = uicontrol('Style', 'pushbutton', 'String', 'Save','FontSize',EPmain.fontsize,...
                'TooltipString','Save segmentation table.',...
                'Tag','save',...
                'Position', [145 220 40 20], 'Callback', @segmentSave);
            
            tableData=EPmain.segment.cellTable;
            
            tableNames{1}='#';
            tableNames{2}='name';
            tableNames{3}='stim';
            if ~EPmain.segment.flexible
                tableNames{4}='prestim';
                tableNames{5}='poststim';
            else
                tableNames{4}='flexEnd';
                tableNames{5}='length';
            end
            tableNames{6}='delay';
            tableNames{7}='task';
            tableNames{8}='ACC';
            tableNames{9}='RT';
            
            for iSpec=1:EPmain.segment.numSpecs
                tableNames{EPmain.segment.numFixed+1+(iSpec-1)*3}=['spec' num2str(iSpec)];
                tableNames{EPmain.segment.numFixed+2+(iSpec-1)*3}=['rel' num2str(iSpec)];
                tableNames{EPmain.segment.numFixed+3+(iSpec-1)*3}=['value' num2str(iSpec)];
            end
            
            columnEditable =  [false true(1,(EPmain.segment.numFixed-1)+EPmain.segment.numSpecs*3)];
            ColumnFormat{1}='numeric';
            ColumnFormat{2}='char';
            ColumnFormat{3}='char';
            ColumnFormat{4}='char';
            ColumnFormat{5}='char';
            ColumnFormat{6}='numeric';
            ColumnFormat{7}='char';
            ColumnFormat{8}='char';
            ColumnFormat{9}='char';
            for iSpec=1:EPmain.segment.numSpecs
                ColumnFormat{EPmain.segment.numFixed+1+(iSpec-1)*3}='char';
                ColumnFormat{EPmain.segment.numFixed+2+(iSpec-1)*3}='char';
                ColumnFormat{EPmain.segment.numFixed+3+(iSpec-1)*3}='char';
            end
            ColumnWidth=num2cell(repmat(50,1,EPmain.segment.numFixed+EPmain.segment.numSpecs*3));
            
            EPmain.handles.segment.cellTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
                'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,'ColumnWidth',ColumnWidth,...
                'TooltipString','The segmentation table.',...
                'RearrangeableColumns','on',...
                'Tag','cellTable',...
                'CellEditCallback',['global EPmain;','EPmain.segment.cellTable=get(EPmain.handles.segment.cellTable,''Data'');','ep(''start'');'],'Position',[5 100 180+EPmain.segment.expand*180*4 120]);
            
            if ~EPmain.segment.expand
                EPmain.handles.segment.expand = uicontrol('Style', 'pushbutton', 'String', '>','FontSize',EPmain.fontsize,...
                    'Tag','expand',...
                    'Position', [185 220 20 20], 'Callback', ['global EPmain;','EPmain.segment.expand=1;','temp=get(EPmain.handles.hMainWindow,''Position'');','temp(3)=temp(3)*5;','set(EPmain.handles.hMainWindow,''Position'',temp);','ep(''start'');']);
            else
                EPmain.handles.segment.expand = uicontrol('Style', 'pushbutton', 'String', '<','FontSize',EPmain.fontsize,...
                    'Tag','expand',...
                    'Position', [185 220 20 20], 'Callback', ['global EPmain;','EPmain.segment.expand=0;','temp=get(EPmain.handles.hMainWindow,''Position'');','temp(3)=temp(3)/5;','set(EPmain.handles.hMainWindow,''Position'',temp);','ep(''start'');']);
            end
            
            uicontrol('Style','text',...
                'String','In','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[5 80 50 20]);
            
            EPmain.handles.segment.importFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',EPmain.fileFormatReadList,...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.segment.importFormat,''Value'');','if tempVar ~=0,EPmain.segment.importFormat=tempVar;end;','if isempty(tempVar),EPmain.segment.importFormat=tempVar;end;','ep(''start'');'],...
                'TooltipString','The file format of the input data to be segmented.',...
                'Tag','importFormat',...
                'Value',EPmain.segment.importFormat,'Position',[50 80 150 20]);
            
            uicontrol('Style','text',...
                'String','Mont','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[5 60 50 20]);
            
            EPmain.handles.segment.importMontage = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',EPmain.montageList,...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.segment.importMontage,''Value'');','if tempVar ~=0,EPmain.segment.importMontage=tempVar;end;','if isempty(tempVar),EPmain.segment.importMontage=tempVar;end;'],...
                'TooltipString','The montage of the input data to be segmented.',...
                'Tag','importMontage',...
                'Value',find(strcmp(EPmain.segment.importMontage,EPmain.montageList)),'Position',[50 60 150 20]);
            
            [importSuffix,importFormatName,importFormat]=ep_fileFormats('average',EPmain.fileFormatReadList{EPmain.segment.importFormat});
            if strcmp(importFormat,'ep_mat')
                set(EPmain.handles.segment.importMontage,'enable','off');
            end
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Out','HorizontalAlignment','left',...
                'Position',[5 40 50 20]);
            
            EPmain.handles.segment.outputFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',EPmain.fileFormatSaveList,...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.segment.outputFormat,''Value'');','if tempVar ~=0,EPmain.segment.outputFormat=tempVar;end;','if isempty(tempVar),EPmain.segment.outputFormat=tempVar;end'],...
                'TooltipString','The file format of the segmented output data.',...
                'Tag','outputFormat',...
                'Value',EPmain.segment.outputFormat,'Position',[50 40 150 20]);
            
            EPmain.handles.segment.trim = uicontrol('Style', 'pushbutton', 'String', 'Trim','FontSize',EPmain.fontsize,...
                'TooltipString','Function for manually trimming unwanted EEG data.',...
                'Tag','trim',...
                'Position', [2 0 50 35], 'Callback', ['global EPtrimData;','EPtrimData=[];','ep_trimData']);
            
            if isempty(EPmain.segment.contData) || ~isempty(EPdataset.dataset(EPmain.segment.dataset).freqNames) || ~isempty(EPdataset.dataset(EPmain.segment.dataset).facNames) || ~strcmp(EPdataset.dataset(EPmain.segment.dataset).dataType,'continuous')
                set(EPmain.handles.segment.trim,'enable','off');
            end
            
            EPmain.handles.segment.preview = uicontrol('Style', 'pushbutton', 'String', 'Preview','FontSize',EPmain.fontsize,...
                'TooltipString','Previews effect of segmentation table on data without yet performing segmentation.',...
                'Tag','preview',...
                'Position', [52 0 50 35], 'Callback', ['EPmain.segment.preview=1;','ep(''segmentData'')']);
            
            if isempty(EPmain.segment.contData)
                set(EPmain.handles.segment.preview,'enable','off');
            end
            
            EPmain.handles.segment.segment = uicontrol('Style', 'pushbutton', 'String', 'Run','FontSize',EPmain.fontsize,...
                'TooltipString','Initiate segmentation function.',...
                'Tag','segment',...
                'Position', [102 0 50 35], 'Callback', ['EPmain.segment.preview=0;','ep(''segmentData'')']);
            
            if any(any(cellfun(@isempty,EPmain.segment.cellTable(:,2:4))))
                set(EPmain.handles.segment.preview,'enable','off');
                set(EPmain.handles.segment.segment,'enable','off');
            end

        else %specs
            set(EPmain.handles.segment.specs,'enable','off');

            if (length(EPmain.segment.spec.specNames)>1) || ~strcmp('none',EPmain.segment.spec.specNames)

                uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Exclude','HorizontalAlignment','left',...
                    'Position',[5 460 50 20]);

                EPmain.handles.segment.excludeSpec = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.spec.specNames,...
                    'CallBack',@changeSpecName,...
                    'Tag','excludeSpec',...
                    'Value',EPmain.segment.spec.excludeSpec,'Position',[1 440 100 20]);
                
                if ~strcmp('none',EPmain.segment.spec.specNames{EPmain.segment.spec.excludeSpec})
                    EPmain.handles.segment.excludeValue = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.spec.specExcludeValues,...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.segment.excludeValue,''Value'');','if tempVar ~=0,EPmain.segment.spec.excludeValue=tempVar;end;','if isempty(tempVar),EPmain.segment.spec.excludeValue=tempVar;end'],...
                        'Tag','excludeValue',...
                        'Value',EPmain.segment.spec.excludeValue,'Position',[100 440 100 20]);
                else
                    EPmain.handles.segment.excludeValue = uicontrol('Style', 'edit', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [100 440 100 20],'enable','off');
                end
                
                EPmain.handles.segment.excludeSpec2 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.spec.specNames,...
                    'CallBack',@changeSpecName,...
                    'Tag','excludeSpec2',...
                    'Value',EPmain.segment.spec.excludeSpec2,'Position',[1 420 100 20]);
                
                if ~strcmp('none',EPmain.segment.spec.specNames{EPmain.segment.spec.excludeSpec2})
                    EPmain.handles.segment.excludeValue2 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.spec.specExcludeValues2,...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.segment.excludeValue2,''Value'');','if tempVar ~=0,EPmain.segment.spec.excludeValue2=tempVar;end;','if isempty(tempVar),EPmain.segment.spec.excludeValue2=tempVar;end'],...
                        'Tag','excludeValue2',...
                        'Value',EPmain.segment.spec.excludeValue2,'Position',[100 420 100 20]);
                else
                    EPmain.handles.segment.excludeValue2 = uicontrol('Style', 'edit', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [100 420 100 20],'enable','off');
                end
                
                EPmain.handles.segment.excludeSpec3 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPmain.segment.spec.specNames,...
                    'CallBack',@changeSpecName,...
                    'Tag','excludeSpec3',...
                    'Value',EPmain.segment.spec.excludeSpec3,'Position',[1 400 100 20]);
                
                if ~strcmp('none',EPmain.segment.spec.specNames{EPmain.segment.spec.excludeSpec3})
                    EPmain.handles.segment.excludeValue3 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.spec.specExcludeValues3,...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.segment.excludeValue3,''Value'');','if tempVar ~=0,EPmain.segment.spec.excludeValue3=tempVar;end;','if isempty(tempVar),EPmain.segment.spec.excludeValue3=tempVar;end'],...
                        'Tag','excludeValue3',...
                        'Value',EPmain.segment.spec.excludeValue3,'Position',[100 400 100 20]);
                else
                    EPmain.handles.segment.excludeValue3 = uicontrol('Style', 'edit', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [100 400 100 20],'enable','off');
                end

                for iSpec=1:8
                    EPmain.handles.segment.specNames(iSpec) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.segment.spec.specNames,...
                        'CallBack',{@updateSpec,iSpec},...
                        'Tag',['specNames(' num2str(iSpec) ' )'],...
                        'Value',EPmain.segment.spec.specField(iSpec),'Position',[1 400-(iSpec*20) 150 20]);
                    EPmain.handles.segment.specLabels(iSpec) = uicontrol('Style', 'edit', 'String', EPmain.segment.spec.specLabels{iSpec},'FontSize',EPmain.fontsize,...
                        'CallBack',{@updateSpec,iSpec},...
                        'Tag',['specLabels(' num2str(iSpec) ' )'],...
                        'Position', [150 400-(iSpec*20) 50 20]);
                end
            else
                uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Exclude','HorizontalAlignment','left',...
                    'Position',[5 460 50 20]);

                EPmain.handles.segment.excludeSpec = uicontrol('Style', 'edit', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [1 440 100 20],'enable','off');

                    EPmain.handles.segment.excludeValue = uicontrol('Style', 'edit', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [100 440 100 20],'enable','off');
                    
                EPmain.handles.segment.excludeSpec2 = uicontrol('Style', 'edit', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [1 420 100 20],'enable','off');

                    EPmain.handles.segment.excludeValue2 = uicontrol('Style', 'edit', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [100 420 100 20],'enable','off');

                EPmain.handles.segment.excludeSpec3 = uicontrol('Style', 'edit', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [1 400 100 20],'enable','off');

                    EPmain.handles.segment.excludeValue3 = uicontrol('Style', 'edit', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [100 400 100 20],'enable','off');

               for iSpec=1:8
                    EPmain.handles.segment.specNames(iSpec) = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [1 400-(iSpec*20) 150 20],'enable','off');
                    EPmain.handles.segment.specLabels(iSpec) = uicontrol('Style', 'edit', 'String', 'none','FontSize',EPmain.fontsize,...
                        'Position', [150 400-(iSpec*20) 50 20],'enable','off');
                end
                
            end
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Format','HorizontalAlignment','left',...
                'Position',[5 220 50 20]);
            
            EPmain.handles.segment.specFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',{'E-Prime Txt','Tab-Text'},...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.segment.specFormat,''Value'');','if tempVar ~=0,EPmain.segment.spec.format=tempVar;end;','if isempty(tempVar),EPmain.segment.spec.format=tempVar;end'],...
                'Tag','specFormat',...
                'Value',EPmain.segment.spec.format,'Position',[40 220 160 20]);
            
            EPmain.handles.segment.clear = uicontrol('Style', 'pushbutton', 'String', 'Clear','FontSize',EPmain.fontsize,...
                'Tag','clear',...
                'Position', [10 190 40 25], 'Callback', {@clearSpec,0});
            
            EPmain.handles.segment.template = uicontrol('Style', 'pushbutton', 'String', 'File','FontSize',EPmain.fontsize,...
                'Tag','template',...
                'Position', [55 190 40 25], 'Callback', @specTemplate);
            
            EPmain.handles.segment.loadSpecs = uicontrol('Style', 'pushbutton', 'String', 'Load','FontSize',EPmain.fontsize,...
                'Tag','loadSpecs',...
                'Position', [100 190 40 25], 'Callback', @specLoad);
            
            EPmain.handles.segment.saveSpecs = uicontrol('Style', 'pushbutton', 'String', 'Save','FontSize',EPmain.fontsize,...
                'Tag','saveSpecs',...
                'Position', [145 190 40 25], 'Callback', @specSave);
            
            EPmain.handles.segment.specTable = uitable('Data',EPmain.segment.spec.specTable,'ColumnName',EPmain.segment.spec.specNames(1:end-1),'FontSize',EPmain.fontsize,...
                'ColumnEditable', EPmain.segment.spec.columnEditable, 'ColumnFormat', EPmain.segment.spec.ColumnFormat,'ColumnWidth',EPmain.segment.spec.ColumnWidth,...
                'RearrangeableColumns','on',...
                'Tag','specTable',...
                'CellEditCallback',['global EPmain;','EPmain.segment.spec.specTable=get(EPmain.handles.segment.specTable,''Data'');','ep(''start'');'],'Position',[5 40 180 150]);
            
        end
        
        EPmain.handles.segment.main = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Tag','main',...
            'Position', [152 0 50 35], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        
        refresh;
        
    case 'segmentData'
        
        ep_tictoc('begin');
        %check cell table
        EPmain.segment.cellTable(:,2)=cellfun(@deblank,EPmain.segment.cellTable(:,2),'UniformOutput',false);
        
        for iCell=1:size(EPmain.segment.cellTable,1)
            if iCell==1
                flexMode=strcmpi(EPmain.segment.cellTable{iCell,5}(1),'F');
            elseif flexMode && ~strcmpi(EPmain.segment.cellTable{iCell,5}(1),'F')
                msg{1}='Error: all cells need to be uniformly either flexmode or not.';
                [msg]=ep_errorMsg(msg);
                return
            end
            if flexMode
                flexLength=str2num(EPmain.segment.cellTable{iCell,5}(2:end));
                if isempty(flexLength) || isnan(flexLength)
                    msg{1}='Error: Flex mode requires a number in addition to the F prefix, as in F20.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if flexLength < 1
                    msg{1}='Error: Flex length needs to be at least one.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if floor(flexLength) ~= flexLength
                    msg{1}='Error: Flex length needs to be an integer.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
            else
                startTime=str2num(EPmain.segment.cellTable{iCell,4});
                endTime=str2num(EPmain.segment.cellTable{iCell,5});
                if isempty(startTime) || isnan(startTime)
                    msg{1}='Error: prestim value is not a number.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if isempty(endTime) || isnan(endTime)
                    msg{1}='Error: poststim value is not a number.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if startTime > endTime
                    msg{1}='Error: cell has starting sample after ending sample.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if (startTime == endTime) && ~EPmain.segment.flexible
                    msg{1}='Error: cell has a zero segment length.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if iCell==1
                    epochLength=endTime-startTime;
                elseif epochLength~=(endTime-startTime)
                    msg{1}='Error: cells do not have same epoch length.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
            end
            
            delayTime=str2num(EPmain.segment.cellTable{iCell,6});
            if isempty(delayTime) || isnan(delayTime)
                msg{1}='Error: delay value is not a number.';
                [msg]=ep_errorMsg(msg);
                return
            end
            
            taskLabel=EPmain.segment.cellTable{iCell,7};
            if ~isempty(taskLabel) && ~any(strcmp(taskLabel(end),{'+','-'}))
                msg{1}='Error: task labels must end in either + or -.';
                [msg]=ep_errorMsg(msg);
                return
            end
            
            for iSpec=1:EPmain.segment.numSpecs
                if strcmp('none',EPmain.segment.cellTable{iCell,EPmain.segment.numFixed+1+(iSpec-1)*3})
                    EPmain.segment.cellTable{iCell,EPmain.segment.numFixed+1+(iSpec-1)*3}='';
                end
                if isempty(strcmp(EPmain.segment.cellTable{iCell,EPmain.segment.numFixed+2+(iSpec-1)*3},EPmain.segment.relList)) && ~isempty(EPmain.segment.cellTable{iCell,EPmain.segment.numFixed+2+(iSpec-1)*3})
                    msg{1}=['Error: spec relationship (' EPmain.segment.cellTable{iCell,EPmain.segment.numFixed+2+(iSpec-1)*3} ') does not match list of options.'];
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if strcmp('none',EPmain.segment.cellTable{iCell,EPmain.segment.numFixed+3+(iSpec-1)*3})
                    EPmain.segment.cellTable{iCell,EPmain.segment.numFixed+3+(iSpec-1)*3}='';
                end
                if any(strcmp(EPmain.segment.cellTable{iCell,EPmain.segment.numFixed+1+(iSpec-1)*3},{'-precedes-','-follows-','-responseCorr-','-responseErr-'}))
                    if  iSpec ~= EPmain.segment.numSpecs
                        if any(strcmp(EPmain.segment.cellTable{iCell,EPmain.segment.numFixed+1+(iSpec)*3},{'-precedes-','-follows-','-responseCorr-','-responseErr-'}))
                            msg{1}='Error: Consecutive ''-precedes-'' and ''-follows-'' and ''-responseCorr-'' and ''-responseErr-'' keywords are not allowed.  The criterion following them must be used to finish their specifications.';
                            [msg]=ep_errorMsg(msg);
                            return
                        end
                    else
                        msg{1}='Error: The final criterion cannot use the ''-precedes-'' or ''-follows-'' or ''-responseCorr-'' or ''-responseErr-'' keywords as a subsequent criterion is needed to finish their specifications.';
                        [msg]=ep_errorMsg(msg);
                        return
                    end
                end
            end
        end
        
        EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);
        drawnow
        
        if EPmain.segment.preview
            sessionFiles=EPdataset.dataset(EPmain.segment.dataset);
            importFormat='ep_mat';
            outputFormat=[];
        else
            [importSuffix,importFormatName,importFormat]=ep_fileFormats('continuous',EPmain.fileFormatReadList{EPmain.segment.importFormat});
            [outputSuffix,outputFormatName,outputFormat]=ep_fileFormats('single_trial',EPmain.fileFormatSaveList{EPmain.segment.outputFormat});
            [sessionFiles, activeDirectory]=ep_getFilesUI(importFormat);
            if isempty(sessionFiles)
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            if (sessionFiles{1}==0) & ~ischar(sessionFiles{1})
                msg{1}='No filenames selected. You have to click on a name.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            for iFile=1:size(sessionFiles,2)
                sessionFiles{iFile}=[activeDirectory sessionFiles{iFile}];
            end
        end
        
        specData=[];
        if ~all(cellfun(@isempty,EPmain.segment.spec.specLabels))
            switch EPmain.segment.spec.format
                case 1
                    specData.specFormat='EPM';
                case 2
                    specData.specFormat='TSP';
                otherwise
                    disp('oops')
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
            end
            specData.excludeSpec=EPmain.segment.spec.specNames{EPmain.segment.spec.excludeSpec};
            if EPmain.segment.spec.excludeSpec==length(EPmain.segment.spec.specNames)
                specData.excludeSpec='';
                specData.excludeValue='';
            else
                specData.excludeValue=EPmain.segment.spec.specExcludeValues{EPmain.segment.spec.excludeValue};
            end
            specData.excludeSpec2=EPmain.segment.spec.specNames{EPmain.segment.spec.excludeSpec2};
            if EPmain.segment.spec.excludeSpec2==length(EPmain.segment.spec.specNames)
                specData.excludeSpec2='';
                specData.excludeValue2='';
            else
                specData.excludeValue2=EPmain.segment.spec.specExcludeValues2{EPmain.segment.spec.excludeValue2};
            end
            specData.excludeSpec3=EPmain.segment.spec.specNames{EPmain.segment.spec.excludeSpec3};
            if EPmain.segment.spec.excludeSpec3==length(EPmain.segment.spec.specNames)
                specData.excludeSpec3='';
                specData.excludeValue3='';
            else
                specData.excludeValue3=EPmain.segment.spec.specExcludeValues3{EPmain.segment.spec.excludeValue3};
            end
            for iSpec=1:length(EPmain.segment.spec.specLabels)
                specData.specField{iSpec}=EPmain.segment.spec.specNames{EPmain.segment.spec.specField(iSpec)};
                if EPmain.segment.spec.specField(iSpec)==length(EPmain.segment.spec.specNames)
                    specData.specField{iSpec}='';
                end
                specData.specLabels{iSpec}=EPmain.segment.spec.specLabels{iSpec};
            end
            specData.specTable=EPmain.segment.spec.specTable;
            specData.specNames=EPmain.segment.spec.specNames;
        end
        
        cellNums=ep_segmentData(EPmain.segment.cellTable(:,2:end),sessionFiles,importFormat,outputFormat,EPmain.segment.preview,specData,EPmain.preferences.general.segSuffix);
        if EPtictoc.stop
            EPtictoc.stop=0;
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            ep('start');
        end
        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        ep_tictoc('end');
        
        if ~isempty(cellNums)
            EPmain.segment.cellTable=get(EPmain.handles.segment.cellTable,'Data');
            if size(cellNums,1) ==1
                theSums=cellNums;
            else
                theSums=sum(cellNums);
            end
            EPmain.segment.cellTable(:,1)=num2cell(theSums);
            set(EPmain.handles.segment.cellTable,'Data',EPmain.segment.cellTable);
            ep('start');
        end
        
    case 'startPreprocess'
        
        set(EPmain.handles.hMainWindow,'Name', 'Preprocess Data');
        
        EPmain.handles.preprocess.prefs = uicontrol('Style', 'pushbutton', 'String', '•','FontSize',EPmain.fontsize,...
            'Tag','prefs',...
            'Position', [5 485 10 10], 'Callback', ['global EPmain;','EPmain.prefReturn=''preprocess'';','EPmain.mode=''preferencePreprocess'';','ep(''start'');']);
        
        uicontrol('Style','text',...
            'String','In','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 430 50 20]);
        
        EPmain.handles.preprocess.importFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatReadList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.importFormat,''Value'');','if tempVar ~=0,EPmain.preprocess.importFormat=tempVar;end;','if isempty(tempVar),EPmain.preprocess.importFormat=tempVar;end;','ep(''start'');'],...
            'Tag','importFormat',...
            'Value',EPmain.preprocess.importFormat,'Position',[50 430 150 20]);
        
        uicontrol('Style','text',...
            'String','Type','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 410 50 20]);
        
        EPmain.handles.preprocess.fileType= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'continuous','single_trial','average','grand_average','factors'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.fileType,''Value'');','if tempVar ~=0,EPmain.preprocess.fileType=tempVar;end;','if isempty(tempVar),EPmain.preprocess.fileType=tempVar;end;','ep(''start'');'],...
            'Tag','fileType',...
            'Value',EPmain.preprocess.fileType,'Position',[50 410 150 20]);
        
        uicontrol('Style','text',...
            'String','Mont','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 390 50 20]);
        
        EPmain.handles.preprocess.importMontage = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.montageList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.importMontage,''Value'');','if tempVar ~=0,EPmain.preprocess.importMontage=EPmain.montageList{tempVar};end;','if isempty(tempVar),EPmain.preprocess.importMontage=EPmain.montageList{tempVar};end'],...
            'Tag','importMontage',...
            'Value',find(strcmp(EPmain.preprocess.importMontage,EPmain.montageList)),'Position',[50 390 150 20]);
        
        [importSuffix,importFormatName,importFormat]=ep_fileFormats('average',EPmain.fileFormatReadList{EPmain.preprocess.importFormat});
        if strcmp(importFormat,'ep_mat')
            set(EPmain.handles.preprocess.fileType,'enable','off');
            set(EPmain.handles.preprocess.importMontage,'enable','off');
        end

        uicontrol('Style','text',...
            'String','Out','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 370 50 20]);
        
        EPmain.handles.preprocess.outputFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatSaveList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.outputFormat,''Value'');','if tempVar ~=0,EPmain.preprocess.outputFormat=tempVar;end;','if isempty(tempVar),EPmain.preprocess.outputFormat=tempVar;end;','ep(''start'');'],...
            'Tag','outputFormat',...
            'Value',EPmain.preprocess.outputFormat,'Position',[50 370 150 20]);        
        
        EPmain.handles.preprocess.timepointsLabel = uicontrol('Style','text','HorizontalAlignment','left','String', 'Points','FontSize',EPmain.fontsize,...
            'Position',[25 350 50 20]);
        
        if isempty(EPmain.preprocess.timepoints)
            set(EPmain.handles.preprocess.timepointsLabel,'enable','off');
        elseif isempty(str2num(EPmain.preprocess.timepoints))
            set(EPmain.handles.preprocess.timepointsLabel,'enable','off');
        end
        
        EPmain.handles.preprocess.timepoints = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.preprocess.timepoints,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.timepoints,''String'');','if tempVar ~=0,EPmain.preprocess.timepoints=tempVar;end;','if isempty(tempVar),EPmain.preprocess.timepoints=tempVar;end;','ep(''start'');'],...
            'Tag','timepoints',...
            'Position',[25 330 50 20],'TooltipString','For retaining a subset of timepoints - example 1:250');
        
        EPmain.handles.preprocess.fMRI= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','fMRI',...
            'CallBack',['global EPmain;','EPmain.preprocess.fMRI=get(EPmain.handles.preprocess.fMRI,''Value'');','ep(''start'');'],...
            'Tag','fMRI',...
            'Value',EPmain.preprocess.fMRI,'Position',[130 350 65 20],'TooltipString','Gradient and BCG artifact correction.  Requires signal processing toolbox.');
        
        EPmain.handles.preprocess.detrend= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Detrend',...
            'CallBack',['global EPmain;','EPmain.preprocess.detrend=get(EPmain.handles.preprocess.detrend,''Value'');','ep(''start'');'],...
            'Tag','detrend',...
            'Value',EPmain.preprocess.detrend,'Position',[130 330 65 20],'TooltipString','Not recommended for segmented ERP trials as the ERP components would be attenuated.');
        
        EPmain.handles.preprocess.EMG= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','EMG',...
            'CallBack',['global EPmain;','EPmain.preprocess.EMG=get(EPmain.handles.preprocess.EMG,''Value'');','ep(''start'');'],...
            'Tag','EMG',...
            'Value',EPmain.preprocess.EMG,'Position',[130 310 65 20],'TooltipString','Removes EMG using BSS-CCA.');
        
        EPmain.handles.preprocess.alpha= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','alpha',...
            'CallBack',['global EPmain;','EPmain.preprocess.alpha=get(EPmain.handles.preprocess.alpha,''Value'');','ep(''start'');'],...
            'Tag','alpha',...
            'Value',EPmain.preprocess.alpha,'Position',[130 290 65 20],'TooltipString','Removes alpha using CWT-PCA.');
                
        EPmain.handles.preprocess.SP= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','SP',...
            'CallBack',['global EPmain;','EPmain.preprocess.SP=get(EPmain.handles.preprocess.SP,''Value'');','ep(''start'');'],...
            'Tag','SP',...
            'Value',EPmain.preprocess.SP,'Position',[130 270 65 20],'TooltipString','Removes saccadic spike potential.');
        
        EPmain.handles.preprocess.baselineLabel=uicontrol('Style','text','HorizontalAlignment','left','String', 'Baseline','FontSize',EPmain.fontsize,...
            'Position',[75 350 60 20]);
        
        if isempty(EPmain.preprocess.baseline)
            set(EPmain.handles.preprocess.baselineLabel,'enable','off');
        elseif isempty(str2num(EPmain.preprocess.baseline))
            set(EPmain.handles.preprocess.baselineLabel,'enable','off');
        end
        
        EPmain.handles.preprocess.baseline = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.preprocess.baseline,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.baseline,''String'');','if tempVar ~=0,EPmain.preprocess.baseline=tempVar;end;','if isempty(tempVar),EPmain.preprocess.baseline=tempVar;end;','ep(''start'');'],...
            'Tag','baseline',...
            'Position',[75 330 50 20],'TooltipString','Recommended - example 1:50 to indicate first fifty timepoints, which would be 200 ms at 250Hz (if timepoints were dropped, 1st retained sample is sample #1)');
        
        uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Edit Mode','HorizontalAlignment','left',...
            'Position',[25 250 100 20]);
        
        EPmain.handles.preprocess.editMode = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'auto','manual','both'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.editMode,''Value'');','if tempVar ~=0,EPmain.preprocess.editMode=tempVar;end;','if isempty(tempVar),EPmain.preprocess.editMode=tempVar;end;','ep(''start'');'],...
            'Tag','editMode',...
            'Value',EPmain.preprocess.editMode,'Position',[20 230 70 20],'TooltipString','Whether bad channel and trial correction is based on manual, automatic, or both types of criteria.');
        
        uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Eye Mode','HorizontalAlignment','left',...
            'Position',[105 250 100 20]);
        
        theList={'MAAC','ICA','EMCP'};
        EPmain.handles.preprocess.eogMode = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',theList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.eogMode,''String'');','EPmain.preprocess.eogMode=tempVar{get(EPmain.handles.preprocess.eogMode,''Value'')};','ep(''start'');'],...
            'Tag','eogMode',...
            'Value',find(strcmp(EPmain.preprocess.eogMode,theList)),'Position',[105 230 80 20],'TooltipString','Procedure for eye movement correction.');
        
        uicontrol('Style','text',...
            'String','Blinks','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 210 70 20]);
        
        EPmain.handles.preprocess.blinkTemplate = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'file','auto','both','eye-track','none'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.blinkTemplate,''Value'');','if tempVar ~=0,EPmain.preprocess.blinkTemplate=tempVar;end;','if isempty(tempVar),EPmain.preprocess.blinkTemplate=tempVar;end;','ep(''start'');'],...
            'Tag','blinkTemplate',...
            'Value',EPmain.preprocess.blinkTemplate,'Position',[20 190 80 20]);
        
        uicontrol('Style','text',...
            'String','Saccades','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[105 210 70 20]);
        
        EPmain.handles.preprocess.saccadeTemplate = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'file','auto','both','eye-track','none'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.saccadeTemplate,''Value'');','if tempVar ~=0,EPmain.preprocess.saccadeTemplate=tempVar;end;','if isempty(tempVar),EPmain.preprocess.saccadeTemplate=tempVar;end;','ep(''start'');'],...
            'Tag','saccadeTemplate',...
            'Value',EPmain.preprocess.saccadeTemplate,'Position',[100 190 80 20]);
        
%         if strcmp(EPmain.preprocess.eogMode,'EMCP')
% %             set(EPmain.handles.preprocess.blinkTemplate,'enable','off');
%             set(EPmain.handles.preprocess.saccadeTemplate,'enable','off');
%         end
        
        uicontrol('Style','text',...
            'String','Bad Chan','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 170 70 20]);
        
        EPmain.handles.preprocess.channelMode = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'replace','mark','none'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.channelMode,''Value'');','if tempVar ~=0,EPmain.preprocess.channelMode=tempVar;end;','if isempty(tempVar),EPmain.preprocess.channelMode=tempVar;end;','ep(''start'');'],...
            'Tag','channelMode',...
            'Value',EPmain.preprocess.channelMode,'Position',[20 150 80 20]);
        
        uicontrol('Style','text',...
            'String','Movement','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[105 170 70 20]);
        
        EPmain.handles.preprocess.trialMode = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'fix','none'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.trialMode,''Value'');','if tempVar ~=0,EPmain.preprocess.trialMode=tempVar;end;','if isempty(tempVar),EPmain.preprocess.trialMode=tempVar;end;','ep(''start'');'],...
            'Tag','trialMode',...
            'Value',EPmain.preprocess.trialMode,'Position',[100 150 80 20]);
        
        uicontrol('Style','text',...
            'String','O Ref','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 130 40 20]);
        
        EPmain.handles.preprocess.origRefType = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',refList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.origRefType,''Value'');','if tempVar ~=0,EPmain.preprocess.origRefType=tempVar;end;','if isempty(tempVar),EPmain.preprocess.origRefType=tempVar;end;','ep(''start'');'],...
            'Tag','origRefType',...
            'Value',EPmain.preprocess.origRefType,'Position',[50 130 100 20],'TooltipString','Original reference channels explicitly represented as a waveform in the data file.');
        
        EPmain.handles.preprocess.origReference = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.preprocess.origReference,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.origReference,''String'');','if tempVar ~=0,EPmain.preprocess.origReference=tempVar;end;','if isempty(tempVar),EPmain.preprocess.origReference=tempVar;end;','ep(''start'');'],...
            'Tag','origReference',...
            'Position',[145 130 45 20],'TooltipString','example: 20 21');
        
        if EPmain.preprocess.origRefType ~= 2
            set(EPmain.handles.preprocess.origReference,'enable','off');
        end
        
        uicontrol('Style','text',...
            'String','C Ref','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 110 40 20]);
        
        EPmain.handles.preprocess.currRefType = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',refList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.currRefType,''Value'');','if tempVar ~=0,EPmain.preprocess.currRefType=tempVar;end;','if isempty(tempVar),EPmain.preprocess.currRefType=tempVar;end;','ep(''start'');'],...
            'Tag','currRefType',...
            'Value',EPmain.preprocess.currRefType,'Position',[50 110 100 20],'TooltipString','Current reference channels explicitly represented as a waveform in the data file.');
        
        EPmain.handles.preprocess.currReference = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.preprocess.currReference,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.preprocess.currReference,''String'');','if tempVar ~=0,EPmain.preprocess.currReference=tempVar;end;','if isempty(tempVar),EPmain.preprocess.currReference=tempVar;end;','ep(''start'');'],...
            'Tag','currReference',...
            'Position',[145 110 45 20],'TooltipString','example: 20 21');
        
        if EPmain.preprocess.currRefType ~= 2
            set(EPmain.handles.preprocess.currReference,'enable','off');
        end
        
        uicontrol('Style','frame',...
            'Position',[5 35 195 70]);
        
        EPmain.handles.preprocess.check= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Single File Mode',...
            'CallBack',['global EPmain;','EPmain.preprocess.check=get(EPmain.handles.preprocess.check,''Value'');','ep(''start'');'],...
            'Tag','check',...
            'Value',EPmain.preprocess.check,'Position',[10 80 120 20]);
        
        if EPmain.preprocess.check
            
            EPmain.handles.preprocess.subject= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.preprocess.subject,...
                'CallBack',['global EPmain;','EPmain.preprocess.subject=get(EPmain.handles.preprocess.subject,''String'');','ep(''start'');'],...
                'Tag','subject',...
                'Position',[10 60 50 20],'TooltipString','example 4:6');
            
            EPmain.handles.preprocess.subjectLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Subject','HorizontalAlignment','left',...
                'Position',[70 60 50 20]);
            
            if isempty(EPmain.preprocess.subject)
                set(EPmain.handles.preprocess.subjectLabel,'enable','off');
            elseif isempty(str2double(EPmain.preprocess.subject))
                set(EPmain.handles.preprocess.subjectLabel,'enable','off');
            end
            
            EPmain.handles.preprocess.trial= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.preprocess.trial,...
                'CallBack',['global EPmain;','EPmain.preprocess.trial=get(EPmain.handles.preprocess.trial,''String'');','ep(''start'');'],...
                'Tag','trial',...
                'Position',[110 60 50 20],'TooltipString','example 10:12');
            
            EPmain.handles.preprocess.trialLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Trial','HorizontalAlignment','left',...
                'Position',[160 60 30 20]);
            
            if isempty(EPmain.preprocess.trial)
                set(EPmain.handles.preprocess.trialLabel,'enable','off');
            elseif isempty(str2double(EPmain.preprocess.trial))
                set(EPmain.handles.preprocess.trialLabel,'enable','off');
            end
            
            EPmain.handles.preprocess.cell= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.preprocess.cell,...
                'CallBack',['global EPmain;','EPmain.preprocess.cell=get(EPmain.handles.preprocess.cell,''String'');','ep(''start'');'],...
                'Tag','cell',...
                'Position',[10 40 50 20],'TooltipString','example 7:9');
            
            EPmain.handles.preprocess.cellLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Cell','HorizontalAlignment','left',...
                'Position',[70 40 50 20]);
            
            if isempty(EPmain.preprocess.cell)
                set(EPmain.handles.preprocess.cellLabel,'enable','off');
            elseif isempty(str2double(EPmain.preprocess.cell))
                set(EPmain.handles.preprocess.cellLabel,'enable','off');
            end
            
            EPmain.handles.preprocess.session= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.preprocess.session,...
                'CallBack',['global EPmain;','EPmain.preprocess.session=get(EPmain.handles.preprocess.session,''String'');','ep(''start'');'],...
                'Tag','session',...
                'Position',[110 40 50 20],'TooltipString','example 7:9');
            
            EPmain.handles.preprocess.sessionLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Session','HorizontalAlignment','left',...
                'Position',[160 40 50 20]);
            
            if isempty(EPmain.preprocess.session)
                set(EPmain.handles.preprocess.sessionLabel,'enable','off');
            elseif isempty(str2double(EPmain.preprocess.session))
                set(EPmain.handles.preprocess.sessionLabel,'enable','off');
            end
        end
        
        EPmain.handles.preprocess.template = uicontrol('Style', 'pushbutton', 'String', 'Template','FontSize',EPmain.fontsize,...
            'Tag','template',...
            'Position', [10 0 60 35], 'Callback', ['global EPtemplate;','EPtemplate=[];','ep_template']);
        
        noData=true;
        for theData=1:length(EPdataset.dataset)
            EEGchans=find(strcmp('EEG',EPdataset.dataset(theData).chanTypes));
            if any(strcmp(EPdataset.dataset(theData).dataType,{'single_trial','continuous'}))
                if ~isempty(EPdataset.dataset(theData).eloc) && (length([EPdataset.dataset(theData).eloc(EEGchans).theta]) == length(EEGchans)) && (length([EPdataset.dataset(theData).eloc(EEGchans).radius]) == length(EEGchans))
                    noData=false;
                end
            end
        end
        if noData
            set(EPmain.handles.preprocess.template,'enable','off');
            if length(EPdataset.dataset) > 1
                set(EPmain.handles.preprocess.template,'string','No Files');
            end
        end
        
        EPmain.handles.preprocess.preprocess = uicontrol('Style', 'pushbutton', 'String', 'Run','FontSize',EPmain.fontsize,...
            'Tag','preprocess',...
            'Position', [70 0 60 35], 'Callback', 'ep(''preprocessData'')');
        
        EPmain.handles.preprocess.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Tag','done',...
            'Position', [130 0 60 35], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        
    case 'preprocessData' %start preprocessing the data
        
        EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);
        
        textPrefs.firstRow=EPmain.preferences.general.firstRow;
        textPrefs.lastRow=EPmain.preferences.general.lastRow;
        textPrefs.firstCol=EPmain.preferences.general.firstCol;
        textPrefs.lastCol=EPmain.preferences.general.lastCol;
        textPrefs.orientation=EPmain.preferences.general.orientation;
        textPrefs.sampleRate=EPmain.preferences.general.sampleRate;
        
        typeNum = get(EPmain.handles.preprocess.fileType,'value');
        switch typeNum
            case 1
                dataType='continuous';
            case 2
                dataType='single_trial';
            case 3
                dataType='average';
            case 4
                dataType='grand_average';
            case 5
                dataType='factors';
        end
        
        importFormatNum = get(EPmain.handles.preprocess.importFormat,'value');
        [importSuffix,importFormatName,importFormat]=ep_fileFormats(dataType,EPmain.fileFormatReadList{importFormatNum});
        
        if strcmp(importFormat,'ep_mat')
            dataType='';
        end
        
        outputFormatNum = get(EPmain.handles.preprocess.outputFormat,'value');
        [outputSuffix,outputFormatName,outputFormat]=ep_fileFormats(dataType,EPmain.fileFormatSaveList{outputFormatNum});
        
        EPmain.preprocess.format=importFormatNum;
        EPmain.preprocess.fileType=typeNum;
        
        timePoints=str2num(get(EPmain.handles.preprocess.timepoints,'string'));
        baseline=str2num(get(EPmain.handles.preprocess.baseline,'string'));
        detrend=get(EPmain.handles.preprocess.detrend,'value');
        fMRI=get(EPmain.handles.preprocess.fMRI,'value');
        if fMRI
            fMRI=EPmain.preferences.preprocess.fMRI;
        end
        
        switch get(EPmain.handles.preprocess.editMode,'value')
            case 1
                editMode='automatic';
            case 2
                editMode='manual';
            case 3
                editMode='both';
        end
        
        switch get(EPmain.handles.preprocess.blinkTemplate,'value')
            case 1
                blinkTemplate='fileTemplate';
            case 2
                blinkTemplate='autoTemplate';
            case 3
                blinkTemplate='bothTemplate';
            case 4
                blinkTemplate='eyeTracker';
            case 5
                blinkTemplate='none';
        end
        
        switch get(EPmain.handles.preprocess.saccadeTemplate,'value')
            case 1
                saccadeTemplate='fileTemplate';
            case 2
                saccadeTemplate='autoTemplate';
            case 3
                saccadeTemplate='bothTemplate';
            case 4
                saccadeTemplate='eyeTracker';
            case 5
                saccadeTemplate='none';
        end
        
        if strcmp(EPmain.preprocess.eogMode,'EMCP')
            if strcmp(saccadeTemplate,'autoTemplate')
                disp('Saccade autotemplate not an option for EMCP')
                saccadeTemplate='none';
            end
            if strcmp(saccadeTemplate,'bothTemplate')
                disp('Saccade autotemplate not an option for EMCP')
                saccadeTemplate='fileTemplate';
            end
            if strcmp(saccadeTemplate,'eyeTracker')
                disp('Saccade eyeTracker not an option for EMCP')
                saccadeTemplate='none';
            end
        end

        if strcmp(EPmain.preprocess.eogMode,'ICA')
            saccMethod=EPmain.preferences.preprocess.blinkRotation;
            SPmethod=EPmain.preferences.preprocess.blinkRotation;
        else
            saccMethod=EPmain.preferences.preprocess.saccRotation;
            SPmethod=EPmain.preferences.preprocess.SProtation;
        end
        switch get(EPmain.handles.preprocess.channelMode,'value')
            case 1
                channelMode='replace';
            case 2
                channelMode='mark';
            case 3
                channelMode='none';
        end
        switch get(EPmain.handles.preprocess.trialMode,'value')
            case 1
                trialMode='fix';
            case 2
                trialMode='none';
        end
        
        [sessionFiles, activeDirectory]=ep_getFilesUI(importFormat);
        if isempty(sessionFiles)
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
        if (sessionFiles{1}==0) & ~ischar(sessionFiles{1})
            msg{1}='No filenames selected. You have to click on a name.';
            [msg]=ep_errorMsg(msg);
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
        
        for iFile=1:size(sessionFiles,2)
            sessionFiles{iFile}=[activeDirectory sessionFiles{iFile}];
        end
        
        sessionFiles=sort(sessionFiles);
        [epDir, ~, ~] = fileparts(which('ep'));
        
        if any(strcmp(blinkTemplate, {'fileTemplate','bothTemplate'}))
            if ismac
                h=msgbox('Please select blinks template file.'); %workaround for Matlab bug
            end
            theDir=[epDir filesep 'templates'];
            if exist([pwd filesep 'blinks.mat'],'file')
                theDir=pwd;
            end
            [FileName,PathName,~] = uigetfile('*.mat','Blink Template',[theDir filesep 'blinks.mat']);
            if ismac
                close(h);
            end
            if FileName == 0
                warndlg('No blink template file selected.');
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return;
            end
            blinkFile=[PathName FileName];
        else
            blinkFile=[];
        end
        
        if any(strcmp(saccadeTemplate, {'fileTemplate','bothTemplate'}))
            if ismac
                h=msgbox('Please select saccade template file.'); %workaround for Matlab bug
            end
            theDir=[epDir filesep 'templates'];
            if exist([pwd filesep 'saccades.mat'],'file')
                theDir=pwd;
            end
            [FileName,PathName,~] = uigetfile('*.mat','Saccade Template',[theDir filesep 'saccades.mat']);
            if ismac
                close(h);
            end
            if FileName == 0
                warndlg('No saccade template file selected.');
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return;
            end
            saccadeFile=[PathName FileName];
        else
            saccadeFile=[];
        end
        
        if EPmain.preprocess.SP
            if ismac
                h=msgbox('Please select spike potential template file.'); %workaround for Matlab bug
            end
            theDir=[epDir filesep 'templates'];
            if exist([pwd filesep 'spikePot.mat'],'file')
                theDir=pwd;
            end
            [FileName,PathName,~] = uigetfile('*.mat','Spike Potential Template',[theDir filesep 'spikePot.mat']);
            if ismac
                close(h);
            end
            if FileName == 0
                warndlg('No spike potential template file selected.');
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return;
            end
            SPfile=[PathName FileName];
            SPtemplate='fileTemplate';
        else
            SPfile=[];
            SPtemplate='none';
        end
        
        inArg=[];
        inArg{1}='files';
        inArg{2}=sessionFiles;
        inArg{3}='format';
        inArg{4}=importFormat;
        inArg{5}='outputFormat';
        inArg{6}=outputFormat;
        inArg{7}='template';
        inArg{8}=blinkTemplate;
        inArg{9}='channelMode';
        inArg{10}=channelMode;
        inArg{11}='saturation';
        inArg{12}=[-EPmain.preferences.preprocess.saturation EPmain.preferences.preprocess.saturation];
        inArg{13}='window';
        inArg{14}=EPmain.preferences.preprocess.window;
        inArg{15}='minmax';
        inArg{16}=EPmain.preferences.preprocess.minmax;
        inArg{17}='badnum';
        inArg{18}=EPmain.preferences.preprocess.badnum;
        inArg{19}='neighbors';
        inArg{20}=EPmain.preferences.preprocess.neighbors;
        inArg{21}='maxneighbor';
        inArg{22}=EPmain.preferences.preprocess.maxneighbor;
        inArg{23}='badchan';
        inArg{24}=EPmain.preferences.preprocess.badchan;
        inArg{25}='blink';
        inArg{26}=EPmain.preferences.preprocess.blink;
        inArg{27}='badtrials';
        inArg{28}=EPmain.preferences.preprocess.badtrials;
        inArg{29}='chunkSize';
        inArg{30}=EPmain.preferences.preprocess.chunkSize;
        inArg{31}='minTrialsPerCell';
        inArg{32}=EPmain.preferences.preprocess.minTrialsPerCell;
        inArg{33}='noadjacent';
        inArg{34}=EPmain.preferences.preprocess.noadjacent;
        inArg{35}='trialMode';
        inArg{36}=trialMode;
        inArg{37}='trialminmax';
        inArg{38}=EPmain.preferences.preprocess.trialminmax;
        inArg{39}='movefacs';
        inArg{40}=EPmain.preferences.preprocess.movefacs;
        inArg{41}='textPrefs';
        inArg{42}=textPrefs;
        inArg{43}='noFigure';
        inArg{44}=EPmain.preferences.preprocess.noFigure;
        inArg{45}='saccTemplate';
        inArg{46}=saccadeTemplate;
        inArg{47}='sacpot';
        inArg{48}=EPmain.preferences.preprocess.sacPot;
        inArg{49}='editMode';
        inArg{50}=editMode;
        inArg{51}='detrend';
        inArg{52}=detrend;
        inArg{53}='fMRI';
        inArg{54}=fMRI;
        inArg{57}='SMIsuffix';
        inArg{58}=EPmain.preferences.general.SMIsuffix;
        inArg{59}='screenSize';
        inArg{60}=EPmain.scrsz;
        inArg{61}='FontSize';
        inArg{62}=EPmain.fontsize;
        inArg{63}='EMG';
        inArg{64}=EPmain.preprocess.EMG;
        inArg{65}=EPmain.preferences.preprocess.EMGratio;
        inArg{66}=EPmain.preferences.preprocess.EMGthresh;
        inArg{67}='alpha';
        inArg{68}=EPmain.preprocess.alpha;
        inArg{71}='blinkMethod';
        inArg{72}=EPmain.preferences.preprocess.blinkRotation;
        inArg{73}='saccMethod';
        inArg{74}=saccMethod;
        inArg{75}='SPmethod';
        inArg{76}=SPmethod;
        inArg{77}='eogMethod';
        inArg{78}=EPmain.preprocess.eogMode;
        inArg{79}='noInternal';
        inArg{80}=EPmain.preferences.general.noInternal;
        inArg{81}='SPtemplate';
        inArg{82}=SPtemplate;

        if ~strcmp(importFormat,'ep_mat')
            inArg{end+1}='montage';
            inArg{end+1}=EPmain.preprocess.importMontage;
        end

        if any(strcmp(blinkTemplate,{'bothTemplate','fileTemplate'}))
            inArg{end+1}='blinkFile';
            inArg{end+1}=blinkFile;
        end
        if any(strcmp(saccadeTemplate,{'bothTemplate','fileTemplate'}))
            inArg{end+1}='saccadeFile';
            inArg{end+1}=saccadeFile;
        end
        if any(strcmp(SPtemplate,{'bothTemplate','fileTemplate'}))
            inArg{end+1}='SPfile';
            inArg{end+1}=SPfile;
        end
        if ~isempty(EPmain.preferences.preprocess.EOGchans)
            inArg{end+1}='eog';
            inArg{end+1}=EPmain.preferences.preprocess.EOGchans;
        end
        SMIsuffix=EPmain.preferences.general.SMIsuffix;
        if ~isempty(SMIsuffix)
            inArg{end+1}='SMIsuffix';
            inArg{end+1}=SMIsuffix;
        end
        specSuffix=EPmain.preferences.general.specSuffix;
        if ~isempty(specSuffix)
            inArg{end+1}='specSuffix';
            inArg{end+1}=specSuffix;
        end
        subjectSpecSuffix=EPmain.preferences.general.subjectSpecSuffix;
        if ~isempty(subjectSpecSuffix)
            inArg{end+1}='subjectSpecSuffix';
            inArg{end+1}=subjectSpecSuffix;
        end
        if ~isempty(dataType)
            inArg{end+1}='type';
            inArg{end+1}=dataType;
        end
        
        errorFlag=0;
        msg=cell(0);
        switch EPmain.preprocess.origRefType
            case 1
                msg{end+1}='Original recording reference could not have been average reference.';
                errorFlag=1;
            case 2
                refChan=str2num(EPmain.preprocess.origReference);
                if isempty(EPmain.preprocess.origReference)
                    msg{end+1}='Please specify explicit recording reference channel(s).';
                    errorFlag=1;
                elseif isempty(refChan)
                    msg{end+1}='Recording reference channel(s) need to be numbers.';
                    errorFlag=1;
                elseif length(refChan) > 2
                    msg{end+1}='No more than two recording reference channels can be indicated.';
                    errorFlag=1;
                end
                inArg{end+1}='origReference';
                inArg{end+1}=refChan;
            case 3
                msg{end+1}='Cannot preprocess CSD data.';
                [msg]=ep_errorMsg(msg);
                return
            case 4
                if ~isempty(EPmain.preprocess.origReference)
                    msg{end+1}='Recording reference channel field should be empty.';
                    errorFlag=1;
                end
        end
        
        switch EPmain.preprocess.currRefType
            case 1
                if ~isempty(EPmain.preprocess.currReference)
                    msg{end+1}='Current reference channel field should be empty.';
                    errorFlag=1;
                else
                    inArg{end+1}='currReference';
                    inArg{end+1}='AVG';
                end
                %             case 2
                %                 if ~isempty(EPmain.preprocess.currReference)
                %                     msg{1}='Current reference channel field should be empty.';
                %                     [msg]=ep_errorMsg(msg);
                %                     return
                %                 else
                %                     inArg{end+1}='currReference';
                %                     inArg{end+1}='none';
                %                 end
            case 2
                refChan=str2num(EPmain.preprocess.currReference);
                if isempty(EPmain.preprocess.currReference)
                    msg{end+1}='Please specify explicit current reference channel(s).';
                    errorFlag=1;
                elseif isempty(refChan)
                    msg{end+1}='Current reference channel(s) need to be numbers.';
                    errorFlag=1;
                elseif length(refChan) > 2
                    msg{end+1}='No more than two current reference channels can be indicated.';
                    errorFlag=1;
                end
                inArg{end+1}='currReference';
                inArg{end+1}=refChan;
            case 3
                msg{end+1}='Cannot preprocess CSD data.';
                [msg]=ep_errorMsg(msg);
                return
            case 4
                if ~isempty(EPmain.preprocess.currReference)
                    msg{1}='Current reference channel field should be empty.';
                    errorFlag=1;
                end
        end
        
        if ~isempty(timePoints)
            if timePoints
                inArg{end+1}='timePoints';
                inArg{end+1}=timePoints;
            end
        end
        
        if ~isempty(baseline)
            if baseline
                inArg{end+1}='baseline';
                inArg{end+1}=baseline;
            end
        end
        
        if get(EPmain.handles.preprocess.check,'value') %single file mode
            subPos=str2num(get(EPmain.handles.preprocess.subject,'string'));
            cellPos=str2num(get(EPmain.handles.preprocess.cell,'string'));
            trialPos=str2num(get(EPmain.handles.preprocess.trial,'string'));
            sessPos=str2num(get(EPmain.handles.preprocess.session,'string'));
            
            if min(subPos) < 1
                beep
                set(EPmain.handles.preprocess.subject,'ForegroundColor','red');
                drawnow
                return;
            end
            if min(cellPos) < 1
                beep
                set(EPmain.handles.preprocess.cell,'ForegroundColor','red');
                drawnow
                return;
            end
            if min(trialPos) < 1
                beep
                set(EPmain.handles.preprocess.trial,'ForegroundColor','red');
                drawnow
                return;
            end
            if min(sessPos) < 1
                beep
                set(EPmain.handles.preprocess.session,'ForegroundColor','red');
                drawnow
                return;
            end
            if isempty(subPos) && isempty(cellPos) && isempty(trialPos) && isempty(sessPos)
                beep
                disp('None of the single file mode fields were set.')
                drawnow
                return;
            end
            
            for iFile=1:size(sessionFiles,2)
                [pathstr, fileName, fileSuffix] = fileparts(sessionFiles{iFile});
                if ~isempty(subPos)
                    if max(subPos) > length(fileName)
                        msg{end+1}=['The file name ' [fileName fileSuffix] ' does not have ' num2str(max(subPos)) ' characters for the subject label.'];
                        errorFlag=1;
                    else
                        theSubs{iFile}=fileName(subPos);
                    end
                else
                    theSubs{iFile}='S01';
                end
                if ~isempty(cellPos)
                    if max(cellPos) > length(fileName)
                        msg{end+1}=['The file name ' [fileName fileSuffix] ' does not have ' num2str(max(cellPos)) ' characters for the cell label.'];
                        errorFlag=1;
                    else
                        theCells{iFile}=fileName(cellPos);
                    end
                else
                    theCells{iFile}='';
                end
                if ~isempty(trialPos)
                    if max(trialPos) > length(fileName)
                        msg{end+1}=['The file name ' [fileName fileSuffix] ' does not have ' num2str(max(trialPos)) ' characters for the trial label.'];
                        errorFlag=1;
                    else
                        theTrials{iFile}=fileName(trialPos);
                    end
                else
                    theTrials{iFile}=[];
                end
                if ~isempty(sessPos)
                    if max(sessPos) > length(fileName)
                        msg{end+1}=['The file name ' [fileName fileSuffix] ' does not have ' num2str(max(sessPos)) ' characters for the session label.'];
                        errorFlag=1;
                    else
                        theSess{iFile}=fileName(sessPos);
                    end
                else
                    theSess{iFile}='';
                end
                if (iFile==1) && ~errorFlag
                    disp(['According to your settings, the first file(' fileName ') is for:']);
                    if ~isempty(subPos)
                        disp(['Subject: ' fileName(subPos)]);
                    end
                    if ~isempty(cellPos)
                        disp(['Cell: ' fileName(cellPos)]);
                    end
                    if ~isempty(trialPos)
                        disp(['Trial: ' fileName(trialPos)]);
                    end
                    if ~isempty(sessPos)
                        disp(['Session: ' fileName(sessPos)]);
                    end
                end
            end
            
            if errorFlag
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            
            uniqueSubs=unique(theSubs);
            
            mergeName = char(inputdlg('Name of new merged dataset?','Dataset name'));
            pause(1); %pause to avoid crash on Windows due to bug in Matlab 2009
            
            inArg{2}=[];
            inArg{4}='ep_mat';
            eloc=[];
            ced=[];
            
            for iSub=1:length(uniqueSubs)
                mergeArg=[];
                theFiles=find(strcmp(uniqueSubs(iSub),theSubs));
                for file=1:length(theFiles)
                    mergeArg{file,1}=sessionFiles{theFiles(file)};
                    mergeArg{file,2}='format';
                    mergeArg{file,3}=importFormat;
                    mergeArg{file,4}='labels';
                    mergeArg{file,5}={theCells(theFiles(file)) theSubs(theFiles(file)) theSess(theFiles(file)) cell(0) theTrials(theFiles(file))};
                    mergeArg{file,6}='textPrefs';
                    mergeArg{file,7}=textPrefs;
                    mergeArg{file,10}='screenSize';
                    mergeArg{file,11}=EPmain.scrsz;
                    mergeArg{file,12}='FontSize';
                    mergeArg{file,13}=EPmain.fontsize;
                    if ~isempty(eloc)
                        mergeArg{file,16}='eloc';
                        mergeArg{file,17}=eloc;
                    end
                    if ~isempty(ced)
                        mergeArg{file,18}='ced';
                        mergeArg{file,19}=ced;
                    end
                    if ~isempty(EPmain.preferences.general.SMIsuffix)
                        mergeArg{file,20}='SMIsuffix';
                        mergeArg{file,21}=EPmain.preferences.general.SMIsuffix;
                    end
                    if ~isempty(EPmain.preferences.general.specSuffix)
                        mergeArg{file,22}='specSuffix';
                        mergeArg{file,23}=EPmain.preferences.general.specSuffix;
                    end
                    if ~isempty(EPmain.preferences.general.subjectSpecSuffix)
                        mergeArg{file,22}='subjectSpecSuffix';
                        mergeArg{file,23}=EPmain.preferences.general.subjectSpecSuffix;
                    end
                    if ~isempty(dataType)
                        mergeArg{file,4}='type';
                        mergeArg{file,5}=dataType;
                    end
                    
                end
                [EPdata eloc]=ep_mergeEPfiles(mergeArg,mergeName);
                if isempty(EPdata)
                    beep();
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    if EPtictoc.stop
                        EPtictoc.stop=0;
                    end
                    return
                end
                ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);ep('start');return;end
                ced=EPdata.ced;
                if ~strcmp(importFormat,'ep_mat')
                    theDescription=['Merged imported single mode files.'];
                    EPdata.history=ep_addHistory(EPdata.history, EPmain.preferences.records, theDescription, EPmain.handles.hMainWindow, EPmain.preferences, EPmain.preferences.EPver,sessionFiles);
                end
                if length(uniqueSubs)==1
                    mergeSubName=mergeName;
                else
                    mergeSubName=[mergeName '-sub' num2str(iSub)];
                end
                [err]=ep_writeData(EPdata,mergeSubName,EPmain.preferences.general.specSuffix,EPmain.preferences.general.subjectSpecSuffix,'ep_mat');
                ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);ep('start');return;end
                if ~isempty(err)
                    inArg{2}{end+1}=[mergeSubName '.ept'];
                    disp(['Generating temporary work file: ' inArg{2}{end} '.']);
                else
                    msg{1}='Writing temporary work file during file merging failed.  Aborting attempt.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
            end
        end
        
        ep_artifactCorrection(inArg);
        if EPtictoc.stop
            EPtictoc.stop=0;
            if isfield(EPmain.preprocess,'handles')
                for iChunk=1:length(EPmain.preprocess.handles.butterflyFig)
                    if ishandle(EPmain.preprocess.handles.butterflyFig{iChunk})
                        EPmain.preprocess.autoClose=1;
                        close(EPmain.preprocess.handles.butterflyFig{iChunk});
                    end
                end
            end
            ep('start')
        end
        
        if get(EPmain.handles.preprocess.check,'value') %single cell mode
            disp('Deleting temporary work files.');
            for i=1:length(inArg{2})
                delete(inArg{2}{i}); %delete temporary merged files
            end
        end
        
        try
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        catch
            ep('start')
        end
        drawnow
        
        ep('start');
        
    case 'startTransform'
        
        set(EPmain.handles.hMainWindow,'Name', 'Transform Data');
        
        EPmain.handles.transform.prefs = uicontrol('Style', 'pushbutton', 'String', '•','FontSize',EPmain.fontsize,...
            'Tag','prefs',...
            'Position', [5 485 10 10], 'Callback', ['global EPmain;','EPmain.prefReturn=''transform'';','EPmain.mode=''preferenceTransform'';','ep(''start'');']);

        uicontrol('Style','text',...
            'String','In','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 430 50 20]);
        
        EPmain.handles.transform.importFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatReadList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.importFormat,''Value'');','if tempVar ~=0,EPmain.transform.importFormat=tempVar;end;','if isempty(tempVar),EPmain.transform.importFormat=tempVar;end;','ep(''start'');'],...
            'Tag','importFormat',...
            'Value',EPmain.transform.importFormat,'Position',[50 430 150 20]);
        
        uicontrol('Style','text',...
            'String','Type','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 410 50 20]);
        
        EPmain.handles.transform.fileType= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'continuous','single_trial','average','grand_average','factors'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.fileType,''Value'');','if tempVar ~=0,EPmain.transform.fileType=tempVar;end;','if isempty(tempVar),EPmain.transform.fileType=tempVar;end;','ep(''start'');'],...
            'Tag','fileType',...
            'Value',EPmain.transform.fileType,'Position',[50 410 150 20]);
        
        uicontrol('Style','text',...
            'String','Mont','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 390 50 20]);
        
        EPmain.handles.transform.importMontage = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.montageList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.importMontage,''Value'');','if tempVar ~=0,EPmain.transform.importMontage=EPmain.montageList{tempVar};end;','if isempty(tempVar),EPmain.transform.importMontage=EPmain.montageList{tempVar};end'],...
            'Tag','importMontage',...
            'Value',find(strcmp(EPmain.transform.importMontage,EPmain.montageList)),'Position',[50 390 150 20]);
        
        [importSuffix,importFormatName,importFormat]=ep_fileFormats('average',EPmain.fileFormatReadList{EPmain.transform.importFormat});
        if strcmp(importFormat,'ep_mat')
            set(EPmain.handles.transform.fileType,'enable','off');
            set(EPmain.handles.transform.importMontage,'enable','off');
        end

        uicontrol('Style','text',...
            'String','Out','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 370 50 20]);
        
        EPmain.handles.transform.outputFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatSaveList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.outputFormat,''Value'');','if tempVar ~=0,EPmain.transform.outputFormat=tempVar;end;','if isempty(tempVar),EPmain.transform.outputFormat=tempVar;end;','ep(''start'');'],...
            'Tag','outputFormat',...
            'Value',EPmain.transform.outputFormat,'Position',[50 370 150 20]);        
        
        uicontrol('Style','text',...
            'String','Data mode:','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 350 105 20]);
        
        [~, chanModes, ~]=ep_chanTypes;
        chanModes=unique(chanModes);
        EPmain.handles.transform.dataMode= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',chanModes,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.dataMode,''Value'');','if tempVar ~=0,EPmain.transform.dataMode=tempVar;end;','if isempty(tempVar),EPmain.transform.dataMode=tempVar;end;','ep(''start'');'],...
            'Tag','dataMode',...
            'Value',EPmain.transform.dataMode,'Position',[105 350 100 20]);
        
        EPmain.handles.transform.referenceLabel= uicontrol('Style','text',...
            'String','Rereference EEG','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 330 105 20]);
        
        EPmain.handles.transform.reference= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Average','Traditional','CSD','PARE','none'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.reference,''Value'');','if tempVar ~=0,EPmain.transform.reference=tempVar;end;','if isempty(tempVar),EPmain.transform.reference=tempVar;end;','ep(''start'');'],...
            'Tag','reference',...
            'Value',EPmain.transform.reference,'Position',[105 330 100 20]);
        
        if ~strcmp(chanModes{EPmain.transform.dataMode},'EEG')
            set(EPmain.handles.transform.referenceLabel,'enable','off');
            set(EPmain.handles.transform.reference,'enable','off');
        end
        
        EPmain.handles.transform.referenceLabel=uicontrol('Style','text','HorizontalAlignment','left','String', 'Reference Channel(s)','FontSize',EPmain.fontsize,...
            'Position',[25 310 150 20]);
        
        EPmain.handles.transform.refChan1 = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.transform.refChan1,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.refChan1,''String''),','if tempVar ~=0,EPmain.transform.refChan1=tempVar;end;','if isempty(tempVar),EPmain.transform.refChan1=tempVar;end;','ep(''start'');'],...
            'Tag','refChan1',...
            'Position',[25 290 50 20]);
        
        EPmain.handles.transform.refChan2 = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.transform.refChan2,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.refChan2,''String'');','if tempVar ~=0,EPmain.transform.refChan2=tempVar;end;','if isempty(tempVar),EPmain.transform.refChan2=tempVar;end;','ep(''start'');'],...
            'Tag','refChan2',...
            'Position',[85 290 50 20]);
        
        if (EPmain.transform.reference ~= 2) || ~strcmp(chanModes{EPmain.transform.dataMode},'EEG')
            set(EPmain.handles.transform.referenceLabel,'enable','off');
            set(EPmain.handles.transform.refChan1,'enable','off');
            set(EPmain.handles.transform.refChan2,'enable','off');
        end
        
        EPmain.handles.transform.detrend= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Detrend',...
            'CallBack',['global EPmain;','EPmain.transform.detrend=get(EPmain.handles.transform.detrend,''Value'');','ep(''start'');'],...
            'Tag','detrend',...
            'Value',EPmain.transform.detrend,'Position',[20 270 160 20],'TooltipString','Detrend data (recommended only for continuous data as it will also detrend ERP components).');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Prestim','FontSize',EPmain.fontsize,...
            'Position',[5 250 50 20]);
        
        EPmain.handles.transform.baselineLabel=uicontrol('Style','text','HorizontalAlignment','left','String', 'Baseline','FontSize',EPmain.fontsize,...
            'Position',[55 250 55 20]);
        
        EPmain.handles.transform.preStim = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.transform.preStim,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.preStim,''String'');','if tempVar ~=0,EPmain.transform.preStim=tempVar;end;','if isempty(tempVar),EPmain.transform.preStim=tempVar;end;','ep(''start'');'],...
            'Tag','preStim',...
            'Position',[5 230 40 20],'TooltipString','Msec start of epoch with respect to stimulus with positive being prior to stimulus.');
        
        EPmain.handles.transform.baselineStart = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.transform.baselineStart,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.baselineStart,''String'');','if tempVar ~=0,EPmain.transform.baselineStart=tempVar;end;','if isempty(tempVar),EPmain.transform.baselineStart=tempVar;end;','ep(''start'');'],...
            'Tag','baselineStart',...
            'Position',[55 230 40 20],'TooltipString','Msec start (left side of sample) of baseline period.');
        
        EPmain.handles.transform.baselineEnd = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.transform.baselineEnd,'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.baselineEnd,''String'');','if tempVar ~=0,EPmain.transform.baselineEnd=tempVar;end;','if isempty(tempVar),EPmain.transform.baselineEnd=tempVar;end;','ep(''start'');'],...
            'Tag','baselineEnd',...
            'Position',[95 230 40 20],'TooltipString','Msec end (right side of sample) of baseline period.');
        
        EPmain.handles.transform.mainsLabel=uicontrol('Style','text','HorizontalAlignment','left','String', 'Mains','FontSize',EPmain.fontsize,...
            'Position',[135 250 35 20]);
        
        theString={'None','50','60'};
        EPmain.handles.transform.mainsFix = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',theString,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.mainsFix,''Value'');','tempVar2=get(EPmain.handles.transform.mainsFix,''String'');','if tempVar ~=0,EPmain.transform.mainsFix=tempVar2{tempVar};end;','if isempty(tempVar),EPmain.transform.mainsFix=tempVar2{tempVar};end;','ep(''start'');'],...
            'Tag','mainsFix',...
            'Value',find(strcmp(EPmain.transform.mainsFix,theString)),'Position',[135 230 75 20],'TooltipString','Fix EMF mains noise.');
                
        if ~ft_hastoolbox('SIGNAL')
            set(EPmain.handles.transform.mainsFix,'enable','off');
            EPmain.transform.mainsFix='None';
        end
        
        EPmain.handles.transform.filterPass= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'none','Low Pass','High Pass','Band Pass','Band Stop','Notch'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.filterPass,''Value'');','if tempVar ~=0,EPmain.transform.filterPass=tempVar;end;','if isempty(tempVar),EPmain.transform.filterPass=tempVar;end;','ep(''start'');'],...
            'Tag','filterPass',...
            'Value',EPmain.transform.filterPass,'Position',[5 210 80 20]);
        
        EPmain.handles.transform.filter1 = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.transform.filter1),'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.transform.filter1,''String''));','if tempVar ~=0,EPmain.transform.filter1=tempVar;end;','if isempty(tempVar),EPmain.transform.filter1=tempVar;end;','ep(''start'');'],...
            'Tag','filter1',...
            'Position',[85 210 30 20],'TooltipString','Lower frequency limit.');
        
        EPmain.handles.transform.filter2 = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.transform.filter2),'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.transform.filter2,''String''));','if tempVar ~=0,EPmain.transform.filter2=tempVar;end;','if isempty(tempVar),EPmain.transform.filter2=tempVar;end;','ep(''start'');'],...
            'Tag','filter2',...
            'Position',[115 210 30 20],'TooltipString','Upper frequency limit.');
        
        EPmain.handles.transform.delay = uicontrol('Style','text','HorizontalAlignment','left','String', [num2str(EPmain.transform.delay) ' delay'],'FontSize',EPmain.fontsize,...
            'Position',[145 210 200 20],'TooltipString','Msec latency delay from use of filter.  Correct during segmentation.');
        
        EPmain.handles.transform.filterType= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'One-Pass Butterworth','Two-Pass Butterworth','One-Pass FIR','Two-Pass FIR','One-Pass FIRLS','Two-Pass FIRLS'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.filterType,''Value'');','if tempVar ~=0,EPmain.transform.filterType=tempVar;end;','if isempty(tempVar),EPmain.transform.filterType=tempVar;end;','ep(''start'');'],...
            'Tag','filterType',...
            'Value',EPmain.transform.filterType,'Position',[5 190 160 20]);
        
        EPmain.handles.transform.filterOrder = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.transform.filterOrder),'FontSize',EPmain.fontsize,...
            'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.transform.filterOrder,''String''));','if tempVar ~=0,EPmain.transform.filterOrder=tempVar;end;','if isempty(tempVar),EPmain.transform.filterOrder=tempVar;end;','ep(''start'');'],...
            'Tag','filterOrder',...
            'Position',[160 190 30 20],'TooltipString','Order of the filter.');
        
        if EPmain.transform.filterPass==1
            set(EPmain.handles.transform.filter1,'enable','off');
            set(EPmain.handles.transform.filter2,'enable','off');
            set(EPmain.handles.transform.filterType,'enable','off');
            set(EPmain.handles.transform.filterOrder,'enable','off');
        end
        
        freqFilter=[];
        freqFilter.trial{1}=EPmain.transform.whiteNoise;
        freqFilter.label{1}='e1';
        freqFilter.time{1}=[-50:199]*.004;
        freqFilter.sampleinfo=[1 250];
        freqFilter.fsample=250;
        
        timeFilter=[];
        timeFilter.trial{1}=[repmat(0,1,100) 1 repmat(0,1,149)];
        timeFilter.label{1}='e1';
        timeFilter.time{1}=[-50:199]*.004;
        timeFilter.sampleinfo=[1 250];
        timeFilter.fsample=250;
        
        cfg=[];
        cfg.pad = 'nextpow2';
        switch EPmain.transform.filterType
            case 1
                filterDirection='onepass';
                filterType='but';
            case 2
                filterDirection='twopass';
                filterType='but';
            case 3
                filterDirection='onepass';
                filterType='fir';
            case 4
                filterDirection='twopass';
                filterType='fir';
            case 5
                filterDirection='onepass';
                filterType='firls';
            case 6
                filterDirection='twopass';
                filterType='firls';
        end
        
        switch EPmain.transform.filterPass
            case 1 %none
                cfg=[];
            case 2 %low pass
                cfg.lpfilter='yes';
                cfg.lpfreq=EPmain.transform.filter1;
                cfg.lpfiltdir=filterDirection;
                cfg.lpfilttype=filterType;
                cfg.lpfiltord=EPmain.transform.filterOrder;
                set(EPmain.handles.transform.filter2,'enable','off');
            case 3 %high pass
                cfg.hpfilter='yes';
                cfg.hpfreq=EPmain.transform.filter1;
                cfg.hpfiltdir=filterDirection;
                cfg.hpfilttype=filterType;
                cfg.hpfiltord=EPmain.transform.filterOrder;
                set(EPmain.handles.transform.filter2,'enable','off');
            case 4 %band pass
                cfg.bpfilter='yes';
                cfg.bpfreq=[EPmain.transform.filter1 EPmain.transform.filter2];
                cfg.bpfiltdir=filterDirection;
                cfg.bpfilttype=filterType;
                cfg.bpfiltord=EPmain.transform.filterOrder;
            case 5 %band stop
                cfg.bsfilter='yes';
                cfg.bsfreq=[EPmain.transform.filter1 EPmain.transform.filter2];
                cfg.bsfiltdir=filterDirection;
                cfg.bsfilttype=filterType;
                cfg.bsfiltord=EPmain.transform.filterOrder;
            case 6 %notch
                cfg.dftfilter='yes';
                cfg.dftfreq=EPmain.transform.filter1;
                cfg.dftfiltdir=filterDirection;
                cfg.dftfilttype=filterType;
                cfg.dftfiltord=EPmain.transform.filterOrder;
                set(EPmain.handles.transform.filter2,'enable','off');
        end
        
        EPmain.transform.cfgFilter=cfg;
        cfg2=[];
        cfg2.taper='hanning';
        cfg2.method='mtmfft';
        cfg2.output='fourier';
        cfg2.pad = 'nextpow2';
        filterProb=0;
        
        if (EPmain.transform.filterPass==1) || ~isempty(EPmain.transform.filter1) && ~(isempty(EPmain.transform.filter2) && ismember(EPmain.transform.filterPass,[4 5]))
            try
                evalc('[timeFiltered] = ft_preprocessing(cfg, timeFilter);');
                evalc('[freqFiltered] = ft_preprocessing(cfg, freqFilter);');
                evalc('[freqFilterFFT] = ft_freqanalysis(cfg2, freqFilter);');
                evalc('[freqFilteredFFT] = ft_freqanalysis(cfg2, freqFiltered);');
            catch
                disp('Filter settings did not work.  Try different settings.  On a Mac, it may be that the security settings are preventing the use of FieldTrip mex files.  See tutorial.')
                lasterr
                filterProb=1;
                freqFilterFFT.freq=[0:125];
                freqFilterFFT.fourierspctrm=zeros(126,1);
                freqFilteredFFT.fourierspctrm=zeros(126,1);
                timeFiltered.trial{1}=timeFilter.trial{1};
            end
            [A originalLatency]=max(timeFilter.trial{1});
            [A newLatency]=max(timeFiltered.trial{1});
            EPmain.transform.delay=(newLatency-originalLatency)*(1000/timeFilter.fsample);
            set(EPmain.handles.transform.delay,'string',[num2str(EPmain.transform.delay) ' delay']);
        else
            freqFilterFFT.freq=[0:125];
            freqFilterFFT.fourierspctrm=zeros(126,1);
            freqFilteredFFT.fourierspctrm=zeros(126,1);
            timeFiltered.trial{1}=timeFilter.trial{1};
        end
        
        %Frequency domain plot of filter effects
        EPmain.transform.freqData=[abs(squeeze(freqFilteredFFT.fourierspctrm))';abs(squeeze(freqFilterFFT.fourierspctrm))'];
        EPmain.transform.freqScale=freqFilterFFT.freq;
        EPmain.transform.freqAxis=[freqFilterFFT.freq(1) freqFilterFFT.freq(end) -.1 max([.1; ([abs(squeeze(freqFilterFFT.fourierspctrm))])])];
        EPmain.handles.transform.frequencyFilter = axes('units','pixels','position',[15 130 75 50]);
        EPmain.handles.transform.freqWaves = plot(EPmain.transform.freqScale,EPmain.transform.freqData);
        axis(EPmain.transform.freqAxis);
        set(EPmain.handles.transform.frequencyFilter,'ButtonDownFcn',@expandTransformChan);
        set(EPmain.handles.transform.freqWaves(1),'ButtonDownFcn',@expandTransformChan);
        set(EPmain.handles.transform.freqWaves(2),'ButtonDownFcn',@expandTransformChan);
        
        %Time domain plot of filter effects
        EPmain.transform.timeData=[timeFiltered.trial{1};timeFilter.trial{1}];
        EPmain.transform.timeScale=timeFilter.time{1}*1000;
        EPmain.transform.timeAxis=[EPmain.transform.timeScale(1) EPmain.transform.timeScale(end) min([-1.2, timeFilter.time{1}]) max([1.2, timeFilter.time{1}])];
        EPmain.handles.transform.timeFilter = axes('units','pixels','position',[115 130 75 50]);
        EPmain.handles.transform.timeWaves = plot(EPmain.transform.timeScale,EPmain.transform.timeData);
        axis(EPmain.transform.timeAxis);
        set(EPmain.handles.transform.timeFilter,'ButtonDownFcn',@expandTransformChan);
        set(EPmain.handles.transform.timeWaves(1),'ButtonDownFcn',@expandTransformChan);
        set(EPmain.handles.transform.timeWaves(2),'ButtonDownFcn',@expandTransformChan);
        set(EPmain.handles.hMainWindow,'DefaultAxesColorOrder',[[1 0 0;0 0 1]]);
        
        EPmain.handles.transform.domainLabel=uicontrol('Style','text',...
            'String','FFT','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 90 40 20]);
        
        EPmain.handles.transform.domain= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'no change','Frequency','Time-Frequency','BOSC','eBOSC'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.domain,''Value'');','if tempVar ~=0,EPmain.transform.domain=tempVar;EPmain.transform.method=1;end;','if isempty(tempVar),EPmain.transform.domain=tempVar;end;','ep(''start'');'],...
            'Tag','domain',...
            'Value',EPmain.transform.domain,'Position',[60 90 140 20]);
        
        if any(ismember(EPmain.transform.domain,[4 5]))
            if isempty(EPdataset.dataset)
                uicontrol('Style','text','String','No files in working set','HorizontalAlignment','left','FontSize',EPmain.fontsize,'Position',[5 70 200 20]);
            else
                EPmain.sampleTest.datasetList=find(~strcmp('continuous',{EPdataset.dataset.dataType}));
                if isempty(EPmain.sampleTest.datasetList)
                    uicontrol('Style','text','HorizontalAlignment','left','String', 'No segmented files in working set','FontSize',EPmain.fontsize,...
                        'Position',[5 70 200 20]);
                else
                    EPmain.handles.transform.BOSCdataset= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',{EPdataset.dataset(EPmain.sampleTest.datasetList).dataName},...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.BOSCdataset,''Value'');','if tempVar ~=0,EPmain.transform.BOSCdataset=tempVar;end;','if isempty(tempVar),EPmain.transform.BOSCdataset=tempVar;end;','ep(''start'');'],...
                        'Tag','BOSCdataset',...
                        'Value',EPmain.transform.BOSCdataset,'Position',[5 70 130 20],'TooltipString','Template file to specify cell for estimate of background activity.');

                    EPmain.transform.BOSC.cellList=unique(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.transform.BOSCdataset)).cellNames);
                    EPmain.handles.transform.BOSCcell= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.transform.BOSC.cellList,...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.BOSC.cell,''Value'');','if tempVar ~=0,EPmain.transform.BOSC.cell=tempVar;;end;','if isempty(tempVar),EPmain.transform.BOSC.cell=tempVar;end;','ep(''start'');'],...
                        'Tag','BOSCcell',...
                        'Value',EPmain.transform.BOSC.cell,'Position',[125 70 80 20],'TooltipString','Cell to use for estimate of background activity.');
                    
                    EPmain.handles.transform.BOSCwidthLabel = uicontrol('Style','text','HorizontalAlignment','left','String', 'Width','FontSize',EPmain.fontsize,...
                        'Position',[5 50 50 20]);
                    
                    EPmain.handles.transform.BOSCwidth = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.transform.BOSC.width),'FontSize',EPmain.fontsize,...
                        'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.transform.BOSC.width,''String''));','if tempVar ~=0,EPmain.transform.BOSC.width=tempVar;end;','if isempty(tempVar),EPmain.transform.BOSC.width=tempVar;end;','ep(''start'');'],...
                        'Tag','BOSCwidth',...
                        'Position',[55 50 50 20],'TooltipString','Width of Morlet wavelet.');
                    
                    EPmain.handles.transform.BOSCthresholdLevel = uicontrol('Style','text','HorizontalAlignment','left','String', 'Thresh','FontSize',EPmain.fontsize,...
                        'Position',[105 50 50 20]);
                    
                    EPmain.handles.transform.BOSCthreshold = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.transform.BOSC.threshold),'FontSize',EPmain.fontsize,...
                        'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.transform.BOSC.threshold,''String''));','if tempVar ~=0,EPmain.transform.BOSC.threshold=tempVar;end;','if isempty(tempVar),EPmain.transform.BOSC.threshold=tempVar;end;','ep(''start'');'],...
                        'Tag','BOSCthreshold',...
                        'Position',[155 50 40 20],'TooltipString','Confidence threshold.');
                    
                    EPmain.handles.transform.BOSCdurationLabel = uicontrol('Style','text','HorizontalAlignment','left','String', 'Duration','FontSize',EPmain.fontsize,...
                        'Position',[5 30 50 20]);
                    
                    EPmain.handles.transform.BOSCduration = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.transform.BOSC.duration),'FontSize',EPmain.fontsize,...
                        'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.transform.BOSC.duration,''String''));','if tempVar ~=0,EPmain.transform.BOSC.duration=tempVar;end;','if isempty(tempVar),EPmain.transform.BOSC.duration=tempVar;end;','ep(''start'');'],...
                        'Tag','BOSCduration',...
                        'Position',[55 30 50 20],'TooltipString','Duration threshold.');
                end
            end
        else
            EPmain.handles.transform.methodLabel=uicontrol('Style','text',...
                'String','Method','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[5 70 55 20]);
            
            switch EPmain.transform.domain
                case 1 %Time
                    menuItems={'none'};
                case 2 %Frequency
                    menuItems={'multi-taper','Hanning'};
                case 3 %Time-Frequency
                    menuItems={'multi-taper','Hanning','wavelet multiplication','wavelet convolution'};
            end
            
            EPmain.handles.transform.method= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',menuItems,...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.method,''Value'');','if tempVar ~=0,EPmain.transform.method=tempVar;end;','if isempty(tempVar),EPmain.transform.method=tempVar;end;','ep(''start'');'],...
                'Tag','method',...
                'Value',EPmain.transform.method,'Position',[60 70 150 20]);
            
            if ~any(strcmp(chanModes{EPmain.transform.dataMode},{'EEG','MEG'}))
                set(EPmain.handles.transform.domainLabel,'enable','off');
                set(EPmain.handles.transform.methodLabel,'enable','off');
                set(EPmain.handles.transform.method,'enable','off');
                set(EPmain.handles.transform.domain,'enable','off');
            end
            
            EPmain.handles.transform.smoothingLabel = uicontrol('Style','text','HorizontalAlignment','left','String', 'Smooth','FontSize',EPmain.fontsize,...
                'Position',[5 50 55 20]);
            EPmain.handles.transform.smoothing = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.transform.smoothing,...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.transform.smoothing,''String'');','if tempVar ~=0,EPmain.transform.smoothing=tempVar;end;','if isempty(tempVar),EPmain.transform.smoothing=tempVar;end;','ep(''start'');'],...
                'Tag','smoothing',...
                'Position',[60 50 50 20]);
            
            if ~ismember(EPmain.transform.domain,[2 3]) || ~ismember(EPmain.transform.method,[1 3 4]) || ~any(strcmp(chanModes{EPmain.transform.dataMode},{'EEG','MEG'}))
                set(EPmain.handles.transform.smoothingLabel,'enable','off');
                set(EPmain.handles.transform.smoothing,'enable','off'); %smoothing not for Hanning
            end
            
            if (EPmain.transform.domain == 2) %Frequency domain analysis
                set(EPmain.handles.transform.baselineLabel,'enable','off');
                set(EPmain.handles.transform.baselineStart,'enable','off');
                set(EPmain.handles.transform.baselineEnd,'enable','off');
            end
            
        end
        
        EPmain.handles.transform.transform = uicontrol('Style', 'pushbutton', 'String', 'Transform','FontSize',EPmain.fontsize,...
            'Tag','transform',...
            'Position', [20 0 80 30], 'Callback', 'ep(''transformData'')');
        
        if filterProb
            set(EPmain.handles.transform.transform,'enable','off');
        end
        
        EPmain.handles.transform.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Tag','done',...
            'Position', [100 0 80 30], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        
    case 'transformData'
        
        switch EPmain.transform.fileType
            case 1
                dataType='continuous';
            case 2
                dataType='single_trial';
            case 3
                dataType='average';
            case 4
                dataType='grand_average';
            case 5
                dataType='factors';
        end
        
        importFormatNum = get(EPmain.handles.transform.importFormat,'value');
        [importSuffix,importFormatName,importFormat]=ep_fileFormats(dataType,EPmain.fileFormatReadList{importFormatNum});
        
        
        outputFormatNum = get(EPmain.handles.transform.outputFormat,'value');
        [outputSuffix,outputFormatName,outputFormat]=ep_fileFormats(dataType,EPmain.fileFormatSaveList{outputFormatNum});
        
        [~, chanModes, ~]=ep_chanTypes;
        chanModes=unique(chanModes);
        dataMode=chanModes{EPmain.transform.dataMode};
        if strcmp(dataMode,'EEG')
            switch EPmain.transform.reference
                case 1
                    referenceMethod='Average';
                case 2
                    referenceMethod='Traditional';
                case 3
                    referenceMethod='CSD';
                case 4
                    referenceMethod='PARE';
                case 5
                    referenceMethod='none';
            end
        else
            referenceMethod='none';
        end
        
        EPmain.transform.detrend=get(EPmain.handles.transform.detrend,'value');
        
%         mainsList=get(EPmain.handles.transform.mainsFix,'string');
%         mainsFix=mainsList{get(EPmain.handles.transform.mainsFix,'value')};
        
        if ~isempty(EPmain.transform.filter1) && ~(isempty(EPmain.transform.filter2) && ismember(EPmain.transform.filterPass,[4 5]))
            if (EPmain.transform.filterPass==4) && (EPmain.transform.filter1 >= EPmain.transform.filter2)
                msg{1}='The start of the bandpass window must come before the end of the bandpass window.';
                [msg]=ep_errorMsg(msg);
                return
            end
            
            if (EPmain.transform.filterPass==5) && (EPmain.transform.filter1 >= EPmain.transform.filter2)
                msg{1}='The start of the bandstop window must come before the end of the bandstop window.';
                [msg]=ep_errorMsg(msg);
                return
            end
            
            if (length(EPmain.transform.filter1) > 1) || (length(EPmain.transform.filter2) > 1) || (length(EPmain.transform.filterOrder) > 1)
                msg{1}='The filter setting in each field needs to be a single number.';
                [msg]=ep_errorMsg(msg);
                return
            end
            
            if (EPmain.transform.filter1 < 0) || (EPmain.transform.filterOrder < 0)
                msg{1}='Filter settings cannot be negative numbers.';
                [msg]=ep_errorMsg(msg);
                return
            end
            
            if ~isempty(EPmain.transform.filter2) && (EPmain.transform.filter2 < 0)
                msg{1}='Filter settings cannot be negative numbers.';
                [msg]=ep_errorMsg(msg);
                return
            end
            
        else
            EPmain.transform.cfgFilter=[];
        end
        
        if any(ismember(EPmain.transform.dataMode,[4 5])) %EEG or MEG
            switch EPmain.transform.domain
                case 1
                    domainName='Time';
                    methodName='None';
                case 2
                    domainName='Frequency';
                    switch EPmain.transform.method
                        case 1
                            methodName='multi-taper';
                        case 2
                            methodName='Hanning';
                    end
                case 3
                    domainName='Time-Frequency';
                    switch EPmain.transform.method
                        case 1
                            methodName='multi-taper';
                        case 2
                            methodName='Hanning';
                        case 3
                            methodName='wavelet multiplication';
                        case 4
                            methodName='wavelet convolution';
                    end
                case 4
                    domainName='BOSC';
                    methodName='None';
                case 5
                    domainName='eBOSC';
                    methodName='None';
            end
        else
            domainName='Time';
            methodName='None';
        end
        
        EPmain.transform.importFormat=importFormatNum;
        EPmain.transform.outputFormat=outputFormatNum;
        EPmain.transform.refChan1=str2num(get(EPmain.handles.transform.refChan1,'String'));
        EPmain.transform.refChan2=str2num(get(EPmain.handles.transform.refChan2,'String'));
        EPmain.transform.baselineStart=str2num(get(EPmain.handles.transform.baselineStart,'String'));
        EPmain.transform.baselineEnd=str2num(get(EPmain.handles.transform.baselineEnd,'String'));
        EPmain.transform.preStim=str2num(get(EPmain.handles.transform.preStim,'String'));
        
        if ~isempty(EPmain.transform.baselineStart) && ~isempty(EPmain.transform.baselineEnd) && (EPmain.transform.baselineStart >= EPmain.transform.baselineEnd) && ~strcmp(domainName,'Frequency')
            msg{1}='The start of the baseline period must come before the end of the baseline period.';
            if (EPmain.transform.baselineStart == 0) && (EPmain.transform.baselineEnd ==0)
                msg{2}='To disable baseline correction, set both fields to blank.';
            end
            [msg]=ep_errorMsg(msg);
            return
        end
        
        if any(ismember(EPmain.transform.domain,[4,5]))
            if (length(EPmain.transform.BOSC.width) > 1) || (length(EPmain.transform.BOSC.duration) > 1) || (length(EPmain.transform.BOSC.threshold) > 1)
                msg{1}='The BOSC settings need to be a single number.';
                [msg]=ep_errorMsg(msg);
                return
            end
            
            if (length(EPmain.transform.BOSC.width) < 0) || (length(EPmain.transform.BOSC.duration) < 0) || (length(EPmain.transform.BOSC.threshold) < 0)
                msg{1}='BOSC settings cannot be negative numbers.';
                [msg]=ep_errorMsg(msg);
                return
            end
            
            if (EPmain.transform.BOSC.threshold < 0) || (EPmain.transform.BOSC.threshold > 1)
                msg{1}='The BOSC threshold must be between 0 and 1, with higher being more conservative.';
                [msg]=ep_errorMsg(msg);
                return
            end
            
            FS=EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.transform.BOSCdataset)).Fs;
            numPoints=length(EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.transform.BOSCdataset)).timeNames);
            disp(['Each shoulder for the background estimate is ' num2str(EPmain.transform.BOSC.width) ' seconds long for 1 Hz.'])
            if numPoints < (2*EPmain.transform.BOSC.width*FS)
                msg{1}=['The data chosen for the background estimation is only ' num2str(numPoints/FS) ' seconds long, so too short.'];
                [msg]=ep_errorMsg(msg);
                return
            end
            disp(['Each shoulder for the data estimate is ' num2str(EPmain.transform.BOSC.width+EPmain.transform.BOSC.duration) ' seconds long for 1 Hz.'])
            if numPoints < (2*EPmain.transform.BOSC.width*EPmain.transform.BOSC.duration*FS)
                msg{1}=['The data are only ' num2str(numPoints/FS) ' seconds long, so too short.'];
                [msg]=ep_errorMsg(msg);
                return
            end
        end

        [inputFiles, activeDirectory]=ep_getFilesUI(importFormat);
        if isempty(inputFiles)
            return
        end
        if inputFiles{1}==0
            return %user hit cancel on file requestor
        end
        for iFile=1:length(inputFiles)
            inputFiles{iFile}=[activeDirectory inputFiles{iFile}];
        end
        
        EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);
        ep_transformData(inputFiles,importFormat,dataType,outputFormat,referenceMethod,EPmain.transform,domainName,methodName,dataMode,EPmain.transform.BOSC);
        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        if EPtictoc.stop
            EPtictoc.stop=0;
        end
        
        ep('start');
        
    case 'startAverage'
        
        set(EPmain.handles.hMainWindow,'Name', 'Average Data');
        
        if ~isfield(EPmain.average,'trialSpecDataset') %just entering average function
            EPmain.average.datasetList=[];
            EPmain.average.channelDataset=1;
            EPmain.average.channel=1;
            EPmain.average.rsChannelDataset=1;
            EPmain.average.rsChannel=cell(5,1);
            for iWindow=1:length(EPmain.average.rsChannel)
                EPmain.average.rsChannel{iWindow}=1;
            end
            EPmain.average.rsMin=cell(5,1);
            EPmain.average.rsMax=cell(5,1);
            EPmain.average.rsMin{1}=300;
            EPmain.average.rsMax{1}=500;
            for iFile=1:length(EPdataset.dataset)
                if ~isempty(EPdataset.dataset(iFile).trialSpecs) && strcmp(EPdataset.dataset(iFile).dataType,'single_trial')
                    EPmain.average.datasetList(end+1)=iFile;
                end
            end
            EPmain.average.trialSpecDataset=0;
            changeAverageDataset;
        end
        
        EPmain.handles.average.prefs = uicontrol('Style', 'pushbutton', 'String', '•','FontSize',EPmain.fontsize,...
            'Position', [5 485 10 10],...
            'Tag','prefs',...
            'Callback', ['global EPmain;','EPmain.prefReturn=''average'';','EPmain.mode=''preferenceAverage'';','ep(''start'');']);
        
        uicontrol('Style','text',...
            'String','In','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 460 50 20]);
        
        EPmain.handles.average.importFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatReadList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.importFormat,''Value'');','if tempVar ~=0,EPmain.average.importFormat=tempVar;end;','if isempty(tempVar),EPmain.average.importFormat=tempVar;end;','ep(''start'');'],...
            'Tag','importFormat',...
            'Value',EPmain.average.importFormat,'Position',[50 460 150 20]);
        
        uicontrol('Style','text',...
            'String','Type','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 440 50 20]);
        
        EPmain.handles.average.fileType= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'continuous','single_trial','average','grand_average','factors'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.fileType,''Value'');','if tempVar ~=0,EPmain.average.fileType=tempVar;end;','if isempty(tempVar),EPmain.average.fileType=tempVar;end;','ep(''start'');'],...
            'Tag','fileType',...
            'Value',EPmain.average.fileType,'Position',[50 440 150 20]);
        
        uicontrol('Style','text',...
            'String','Mont','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 420 50 20]);
        
        EPmain.handles.average.importMontage = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.montageList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.importMontage,''Value'');','if tempVar ~=0,EPmain.average.importMontage=EPmain.montageList{tempVar};end;','if isempty(tempVar),EPmain.average.importMontage=EPmain.montageList{tempVar};end'],...
            'Tag','importMontage',...
            'Value',find(strcmp(EPmain.average.importMontage,EPmain.montageList)),'Position',[50 420 150 20]);
        
        [importSuffix,importFormatName,importFormat]=ep_fileFormats('average',EPmain.fileFormatReadList{EPmain.average.importFormat});
        if strcmp(importFormat,'ep_mat')
            set(EPmain.handles.average.fileType,'enable','off');
            set(EPmain.handles.average.importMontage,'enable','off');
        end

        uicontrol('Style','text',...
            'String','Out','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 400 50 20]);
        
        EPmain.handles.average.outputFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatSaveList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.outputFormat,''Value'');','if tempVar ~=0,EPmain.average.outputFormat=tempVar;end;','if isempty(tempVar),EPmain.average.outputFormat=tempVar;end;','ep(''start'');'],...
            'Tag','outputFormat',...
            'Value',EPmain.average.outputFormat,'Position',[50 400 150 20]);

        uicontrol('Style','text','HorizontalAlignment','left','String', 'Model','FontSize',EPmain.fontsize,...
            'Position',[5 380 40 20]);
        
        if ~isempty(EPmain.average.datasetList)
            EPmain.handles.average.trialSpecDataset = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',{EPdataset.dataset(EPmain.average.datasetList).dataName},...
                'Value',EPmain.average.trialSpecDataset,'Position',[50 380 150 20],'TooltipString','Model single-trial dataset for providing trial spec settings.',...
                'Tag','trialSpecDataset',...
                'CallBack',@changeAverageDataset);
        else
            EPmain.handles.average.trialSpecDataset = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String','none',...
                'Value',1,'Position',[50 380 150 20],'TooltipString','Model single-trial dataset for providing trial specs for jitter correction.');
        end
        
        EPmain.handles.average.dropEvents= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Drop events',...
            'CallBack',['global EPmain;','EPmain.average.dropEvents=get(EPmain.handles.average.dropEvents,''Value'');','ep(''start'');'],'FontSize',EPmain.fontsize,...
            'Tag','dropEvents',...
            'Value',EPmain.average.dropEvents,'Position',[5 360 100 20]);
        
        EPmain.handles.average.dropNoise= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Drop noise',...
            'CallBack',['global EPmain;','EPmain.average.dropNoise=get(EPmain.handles.average.dropNoise,''Value'');','ep(''start'');'],'FontSize',EPmain.fontsize,...
            'Tag','dropNoise',...
            'Value',EPmain.average.dropNoise,'Position',[100 360 100 20],'TooltipString','Whether to keep the noise information in order to reduce file size.');
        
        uicontrol('Style','text',...
            'String','Trial SD','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 340 50 20]);
        
        EPmain.handles.average.dropSD= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'SD only','SD and covariances','none'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.dropSD,''Value'');','if tempVar ~=0,EPmain.average.dropSD=tempVar;end;','if isempty(tempVar),EPmain.average.dropSD=tempVar;end;','ep(''start'');'],...
            'Tag','dropSD',...
            'Value',EPmain.average.dropSD,'Position',[50 340 150 20],'TooltipString','Whether to keep the SD information in order to reduce file size.  Dropping covariance information greatly reduces size and time but loses ability to rereference etc. SD');
        
        uicontrol('Style','text',...
            'String','Type','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 320 50 20]);
        
        EPmain.handles.average.averageType = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'subject','item'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.averageType,''Value'');','if tempVar ~=0,EPmain.average.averageType=tempVar;end;','if isempty(tempVar),EPmain.average.averageType=tempVar;end;','ep(''start'');'],...
            'Tag','averageType',...
            'Value',EPmain.average.averageType,'Position',[50 320 150 20]);
        
        if EPmain.average.averageType==1 %subject average
            
            EPmain.handles.average.subjectLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,'ForegroundColor','red',...
                'String','Subject','HorizontalAlignment','left',...
                'Tag','subjectLabel',...
                'Position',[5 300 50 20]);
            
            EPmain.handles.average.subject= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.average.subject,...
                'CallBack',['global EPmain;','EPmain.average.subject=get(EPmain.handles.average.subject,''String'');','ep(''start'');'],...
                'Tag','subject',...
                'Position',[50 300 50 20],'TooltipString','For when a subject is spread across multiple files and its ID is encoded in the file name, example 4:6 for 4th-6th characters of ''sub001.ept''.');
            
            EPmain.handles.average.sessionLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,'ForegroundColor','black',...
                'String','Session','HorizontalAlignment','left',...
                'Position',[105 300 50 20]);
            
            EPmain.handles.average.session= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.average.session,...
                'CallBack',['global EPmain;','EPmain.average.session=get(EPmain.handles.average.session,''String'');','ep(''start'');'],...
                'Tag','session',...
                'Position',[155 300 45 20],'TooltipString','For when a session is spread across multiple files and its ID is encoded in the file name, example 7:9 for 7th-9th characters of ''sub001pre.ept''.');
            
            EPmain.handles.average.cellLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,'ForegroundColor','black',...
                'String','Cell','HorizontalAlignment','left',...
                'Position',[5 280 50 20]);
            
            EPmain.handles.average.cell= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.average.cell,...
                'CallBack',['global EPmain;','EPmain.average.cell=get(EPmain.handles.average.cell,''String'');','ep(''start'');'],...
                'Tag','cell',...
                'Position',[50 280 50 20],'TooltipString','For when each cell is in a separate file and its ID is encoded in the file name, example 8:10 for 8th-10th characters of ''sub001TAR.ept''.');
        else
            
            uicontrol('Style','text','HorizontalAlignment','left','String', 'Stim','FontSize',EPmain.fontsize,...
                'Position',[5 280 50 20]);
            
            if ~isempty(EPmain.average.datasetList)
                EPmain.handles.average.itemSpec = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPdataset.dataset(EPmain.average.datasetList(EPmain.average.trialSpecDataset)).trialSpecNames,...
                    'Value',EPmain.average.itemSpec,'Position',[50 280 150 20],'TooltipString','Trial spec to be used for item averaging.',...
                    'Tag','itemSpec',...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.itemSpec,''Value'');','if tempVar ~=0,EPmain.average.itemSpec=tempVar;end;','if isempty(tempVar),EPmain.average.itemSpec=tempVar;end;','ep(''start'');']);
            else
                EPmain.handles.average.itemSpec = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String','none',...
                    'Tag','itemSpec',...
                    'Value',1,'Position',[50 280 150 20],'TooltipString','Trial spec to be used for item averaging.');
            end
            
        end
        
        uicontrol('Style','text',...
            'String','Proc','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 260 50 20]);
        
        EPmain.handles.average.method= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Average','Latency-Lock','Jitter-Correct','Frequency-Coherence','Time-Frequency-Coherence','RIDE','ReSync'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.method,''Value'');','if tempVar ~=0,EPmain.average.method=tempVar;end;','if isempty(tempVar),EPmain.average.method=tempVar;end;','EPmain.average.freqMethod=1;','ep(''start'');'],...
            'Tag','method',...
            'Value',EPmain.average.method,'Position',[50 260 150 20]);
        
        uicontrol('Style','text',...
            'String','Method','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 240 50 20]);
        
        switch EPmain.average.method
            case 1 %Average
                menuItems={'Mean','Median','Trimmed Mean'};
            case 2 %Latency Lock
                menuItems={'Mean','Median','Trimmed Mean'};
            case 3 %Jitter-Correct
                menuItems={'Mean','Median','Trimmed Mean'};
            case 4 %Frequency-Coherence
                menuItems={'multi-taper','Hanning'};
            case 5 %Time-Frequency-Coherence
                menuItems={'phase-lock wavelet'};
            case 6 %RIDE
                menuItems={'s-c1-r','s-r','s-c1','s-c1-c2-r','s-c1-c2'};
            case 7 %ReSync
                menuItems={'Mean','Median','Trimmed Mean'};
        end
        
        EPmain.handles.average.freqMethod= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',menuItems,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.freqMethod,''Value'');','if tempVar ~=0,EPmain.average.freqMethod=tempVar;end;','if isempty(tempVar),EPmain.average.freqMethod=tempVar;end;','ep(''start'');'],...
            'Tag','freqMethod',...
            'Value',EPmain.average.freqMethod,'Position',[50 240 150 20]);
        
        switch EPmain.average.method
            case 1
                
            case 2
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Latency','FontSize',EPmain.fontsize,...
                    'Position',[25 160 88 20]);
                
                if ~isempty(EPmain.average.datasetList)
                    EPmain.handles.average.trialSpec = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPdataset.dataset(EPmain.average.datasetList(EPmain.average.trialSpecDataset)).trialSpecNames,...
                        'Value',EPmain.average.latencySpec,'Position',[20 140 160 20],'TooltipString','Trial spec to be used for jitter correction.',...
                        'Tag','trialSpec',...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.trialSpec,''Value'');','if tempVar ~=0,EPmain.average.latencySpec=tempVar;end;','if isempty(tempVar),EPmain.average.latencySpec=tempVar;end;','ep(''start'');']);
                else
                    EPmain.handles.average.trialSpec = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String','none',...
                        'Tag','trialSpec',...
                        'Value',1,'Position',[20 140 160 20],'TooltipString','Trial spec to be used for jitter correction.');
                end
                
                uicontrol('Style','text',...
                    'String','Latency Range','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[25 120 88 20]);
                
                uicontrol('Style','text',...
                    'String','Prestim','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[125 120 88 20]);
                
                EPmain.handles.average.minLatency = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.minLatency,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.minLatency=str2num(get(EPmain.handles.average.minLatency,''String''));','ep(''start'');'],...
                    'Tag','minLatency',...
                    'Position',[25 100 50 20],'TooltipString','Lower limit of latencies to include.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
                EPmain.handles.average.maxLatency = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.maxLatency,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.maxLatency=str2num(get(EPmain.handles.average.maxLatency,''String''));','ep(''start'');'],...
                    'Tag','maxLatency',...
                    'Position',[75 100 50 20],'TooltipString','Upper limit of latencies to include.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
                if ~isempty(EPmain.average.datasetList)
                    EPmain.handles.average.prestim = uicontrol('Style','text','HorizontalAlignment','left','String', EPdataset.dataset(EPmain.average.datasetList(EPmain.average.trialSpecDataset)).baseline,'FontSize',EPmain.fontsize,...
                        'Position',[125 100 50 20]);
                else
                    EPmain.handles.average.prestim = uicontrol('Style','text','HorizontalAlignment','left','String', '','FontSize',EPmain.fontsize,...
                        'Position',[125 100 50 20]);
                end
            case 3 %jitter
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Model','FontSize',EPmain.fontsize,...
                    'Position',[5 180 50 20]);
                
                EPmain.handles.average.channelDataset = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{EPdataset.dataset.dataName},...
                    'Value',EPmain.average.channelDataset,'Position',[45 180 150 20],'TooltipString','Model dataset for providing channel for jitter correction.',...
                    'Tag','channelDataset',...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.channelDataset,''Value'');','if tempVar ~=0,EPmain.average.channelDataset=tempVar;end;','if isempty(tempVar),EPmain.average.channelDataset=tempVar;end;','EPmain.average.channel=1;','ep(''start'');']);
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Channel','FontSize',EPmain.fontsize,...
                    'Position',[5 160 50 20]);
                
                EPmain.handles.average.channel = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',EPdataset.dataset(EPmain.average.channelDataset).chanNames,...
                    'Value',EPmain.average.channel,'Position',[45 160 150 20],'TooltipString','Channel used for jitter correction.',...
                    'Tag','channel',...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.channel,''Value'');','if tempVar ~=0,EPmain.average.channel=tempVar;end;','if isempty(tempVar),EPmain.average.channel=tempVar;end;','ep(''start'');']);
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Polarity','FontSize',EPmain.fontsize,...
                    'Position',[5 140 50 20]);
                
                EPmain.handles.average.peakPolarity = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{'Positive','Negative'},...
                    'Value',EPmain.average.peakPolarity,'Position',[45 140 150 20],'TooltipString','Use positive or negative peak for jitter correction.',...
                    'Tag','peakPolarity',...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.peakPolarity,''Value'');','if tempVar ~=0,EPmain.average.peakPolarity=tempVar;end;','if isempty(tempVar),EPmain.average.peakPolarity=tempVar;end;','ep(''start'');']);
                
                uicontrol('Style','text',...
                    'String','Range','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[5 120 50 20]);
                
                EPmain.handles.average.minLatency = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.minLatency,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.minLatency=str2num(get(EPmain.handles.average.minLatency,''String''));','ep(''start'');'],...
                    'Tag','minLatency',...
                    'Position',[45 120 50 20],'TooltipString','Lower limit of peaks to include.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
                EPmain.handles.average.maxLatency = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.maxLatency,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.maxLatency=str2num(get(EPmain.handles.average.maxLatency,''String''));','ep(''start'');'],...
                    'Tag','maxLatency',...
                    'Position',[100 120 50 20],'TooltipString','Upper limit of peaks to include.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
            case {4,5}
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Smoothing','FontSize',EPmain.fontsize,...
                    'Position',[25 180 65 20]);
                EPmain.handles.average.smoothing = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.smoothing,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.smoothing,''String'');','if tempVar ~=0,EPmain.average.smoothing=tempVar;end;','if isempty(tempVar),EPmain.average.smoothing=tempVar;end;','ep(''start'');'],...
                    'Tag','smoothing',...
                    'Position',[25 160 50 20]);
                
                if ~ismember(EPmain.average.method,[4 5]) || (EPmain.average.freqMethod ~= 1)
                    set(EPmain.handles.average.smoothing,'enable','off'); %smoothing only seems to be for multi-taper
                end
                
                set(EPmain.handles.average.dropNoise,'enable','off');
                set(EPmain.handles.average.dropSD,'enable','off');
                
            case 6 %RIDE
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'S latency','FontSize',EPmain.fontsize,...
                    'Position',[5 220 55 20]);
                EPmain.handles.average.sLatency = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.sLatency,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.sLatency,''String'');','if tempVar ~=0,EPmain.average.sLatency=tempVar;end;','if isempty(tempVar),EPmain.average.sLatency=tempVar;end;','ep(''start'');'],...
                    'Tag','sLatency',...
                    'Position',[75 220 50 20]);
                
                if isempty(strfind(menuItems{EPmain.average.freqMethod},'s'))
                    set(EPmain.handles.average.sLatency,'enable','off'); %s latency only relevant if s is being estimated.
                end
                
                uicontrol('Style','text',...
                    'String','S window','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[5 200 60 20]);
                
                EPmain.handles.average.sMin = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.sMin,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.sMin=str2num(get(EPmain.handles.average.sMin,''String''));','ep(''start'');'],...
                    'Tag','sMin',...
                    'Position',[75 200 50 20],'TooltipString','Start of window in ms for the stimulus-locked components.');
                
                EPmain.handles.average.sMax = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.sMax,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.sMax=str2num(get(EPmain.handles.average.sMax,''String''));','ep(''start'');'],...
                    'Tag','sMax',...
                    'Position',[130 200 50 20],'TooltipString','End of window in ms for the stimulus-locked components.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
                if isempty(strfind(menuItems{EPmain.average.freqMethod},'s'))
                    set(EPmain.handles.average.sMin,'enable','off');
                    set(EPmain.handles.average.sMax,'enable','off');
                end
                
                uicontrol('Style','text',...
                    'String','C1 window','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[5 180 65 20]);
                
                EPmain.handles.average.c1Min = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.c1Min,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.c1Min=str2num(get(EPmain.handles.average.c1Min,''String''));','ep(''start'');'],...
                    'Tag','c1Min',...
                    'Position',[75 180 50 20],'TooltipString','Start of window in ms for the first intermediate components.');
                
                EPmain.handles.average.c1Max = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.c1Max,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.c1Max=str2num(get(EPmain.handles.average.c1Max,''String''));','ep(''start'');'],...
                    'Tag','c1Max',...
                    'Position',[130 180 50 20],'TooltipString','End of window in ms for the first intermediate components.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
                if isempty(strfind(menuItems{EPmain.average.freqMethod},'c1'))
                    set(EPmain.handles.average.c1Min,'enable','off');
                    set(EPmain.handles.average.c1Max,'enable','off');
                end
                
                uicontrol('Style','text',...
                    'String','C2 window','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[5 160 65 20]);
                
                EPmain.handles.average.c2Min = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.c2Min,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.c2Min=str2num(get(EPmain.handles.average.c2Min,''String''));','ep(''start'');'],...
                    'Tag','c2Min',...
                    'Position',[75 160 50 20],'TooltipString','Start of window in ms for the second intermediate components.');
                
                EPmain.handles.average.c2Max = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.c2Max,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.c2Max=str2num(get(EPmain.handles.average.c2Max,''String''));','ep(''start'');'],...
                    'Tag','c2Max',...
                    'Position',[130 160 50 20],'TooltipString','End of window in ms for the second intermediate components.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
                if isempty(strfind(menuItems{EPmain.average.freqMethod},'c2'))
                    set(EPmain.handles.average.c2Min,'enable','off');
                    set(EPmain.handles.average.c2Max,'enable','off');
                end
                
                uicontrol('Style','text',...
                    'String','R window','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[5 140 60 20]);
                
                EPmain.handles.average.rMin = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.rMin,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.rMin=str2num(get(EPmain.handles.average.rMin,''String''));','ep(''start'');'],...
                    'Tag','rMin',...
                    'Position',[75 140 50 20],'TooltipString','Start of window in ms for the response-locked components.');
                
                EPmain.handles.average.rMax = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.rMax,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.rMax=str2num(get(EPmain.handles.average.rMax,''String''));','ep(''start'');'],...
                    'Tag','rMax',...
                    'Position',[130 140 50 20],'TooltipString','End of window in ms for the response-locked components.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
                if isempty(strfind(menuItems{EPmain.average.freqMethod},'r'))
                    set(EPmain.handles.average.rMin,'enable','off');
                    set(EPmain.handles.average.rMax,'enable','off');
                end
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'C cut-off','FontSize',EPmain.fontsize,...
                    'Position',[5 120 55 20]);
                EPmain.handles.average.cCutoff = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.cCutoff,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.cCutoff,''String'');','if tempVar ~=0,EPmain.average.cCutoff=tempVar;end;','if isempty(tempVar),EPmain.average.cCutoff=tempVar;end;','ep(''start'');'],...
                    'Tag','cCutoff',...
                    'Position',[65 120 40 20],'TooltipString','Low-pass filter setting for C latency estimation.');
                
                if isempty(strfind(menuItems{EPmain.average.freqMethod},'c'))
                    set(EPmain.handles.average.cCutoff,'enable','off'); %c cutoff only relevant if c is being estimated.
                end
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'C lag','FontSize',EPmain.fontsize,...
                    'Position',[105 120 40 20]);
                EPmain.handles.average.cLag = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.cLag,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.cLag,''String'');','if tempVar ~=0,EPmain.average.cLag=tempVar;end;','if isempty(tempVar),EPmain.average.cLag=tempVar;end;','ep(''start'');'],...
                    'Tag','cLag',...
                    'Position',[145 120 40 20],'TooltipString','+/- allowable variability for C latency estimation.');
                
                if isempty(strfind(menuItems{EPmain.average.freqMethod},'c'))
                    set(EPmain.handles.average.cLag,'enable','off'); %c lag only relevant if c is being estimated.
                end
                
            case 7 %ReSync
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Model','FontSize',EPmain.fontsize,...
                    'Position',[5 220 50 20]);
                
                EPmain.handles.average.rsChannelDataset = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{EPdataset.dataset.dataName},...
                    'Value',EPmain.average.rsChannelDataset,'Position',[50 220 150 20],'TooltipString','Model dataset for providing channel for ReSynch correction.',...
                    'Tag','rsChannelDataset',...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.rsChannelDataset,''Value'');','if tempVar ~=0,EPmain.average.rsChannelDataset=tempVar;end;','if isempty(tempVar),EPmain.average.rsChannelDataset=tempVar;end;','EPmain.average.channel=1;','ep(''start'');']);
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Wind1','FontSize',EPmain.fontsize,...
                    'Position',[5 200 35 20]);
                
                EPmain.handles.average.rsMin1 = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.rsMin{1},'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.rsMin{1}=str2num(get(EPmain.handles.average.rsMin1,''String''));','ep(''start'');'],...
                    'Tag','rsMin1',...
                    'Position',[40 200 40 20],'TooltipString','Start of window in ms for ReSynch.');
                
                EPmain.handles.average.rsMax1 = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.rsMax{1},'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.rsMax{1}=str2num(get(EPmain.handles.average.rsMax1,''String''));','ep(''start'');'],...
                    'Tag','rsMax1',...
                    'Position',[80 200 40 20],'TooltipString','End of window in ms for ReSynch.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
                EPmain.handles.average.rsChannel1 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',['global';EPdataset.dataset(EPmain.average.rsChannelDataset).chanNames],...
                    'Value',EPmain.average.rsChannel{1},'Position',[115 200 90 20],'TooltipString','Channel used for ReSynch correction.',...
                    'Tag','rsChannel1',...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.rsChannel1,''Value'');','if tempVar ~=0,EPmain.average.rsChannel{1}=tempVar;end;','if isempty(tempVar),EPmain.average.rsChannel{1}=tempVar;end;','ep(''start'');']);
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Wind2','FontSize',EPmain.fontsize,...
                    'Position',[5 180 35 20]);
                
                EPmain.handles.average.rsMin2 = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.rsMin{2},'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.rsMin{2}=str2num(get(EPmain.handles.average.rsMin2,''String''));','ep(''start'');'],...
                    'Tag','rsMin2',...
                    'Position',[40 180 40 20],'TooltipString','Start of window in ms for ReSynch.');
                
                EPmain.handles.average.rsMax2 = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.rsMax{2},'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.rsMax{2}=str2num(get(EPmain.handles.average.rsMax2,''String''));','ep(''start'');'],...
                    'Tag','rsMax2',...
                    'Position',[80 180 40 20],'TooltipString','End of window in ms for ReSynch.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
                EPmain.handles.average.rsChannel2 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',['global';EPdataset.dataset(EPmain.average.rsChannelDataset).chanNames],...
                    'Value',EPmain.average.rsChannel{2},'Position',[115 180 90 20],'TooltipString','Channel used for ReSynch correction.',...
                    'Tag','rsChannel2',...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.rsChannel2,''Value'');','if tempVar ~=0,EPmain.average.rsChannel{2}=tempVar;end;','if isempty(tempVar),EPmain.average.rsChannel{2}=tempVar;end;','ep(''start'');']);
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Wind3','FontSize',EPmain.fontsize,...
                    'Position',[5 160 35 20]);
                
                EPmain.handles.average.rsMin3 = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.rsMin{3},'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.rsMin{3}=str2num(get(EPmain.handles.average.rsMin3,''String''));','ep(''start'');'],...
                    'Tag','rsMin3',...
                    'Position',[40 160 40 20],'TooltipString','Start of window in ms for ReSynch.');
                
                EPmain.handles.average.rsMax3 = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.rsMax{3},'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.rsMax{3}=str2num(get(EPmain.handles.average.rsMax3,''String''));','ep(''start'');'],...
                    'Tag','rsMax3',...
                    'Position',[80 160 40 20],'TooltipString','End of window in ms for ReSynch.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
                EPmain.handles.average.rsChannel3 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',['global';EPdataset.dataset(EPmain.average.rsChannelDataset).chanNames],...
                    'Value',EPmain.average.rsChannel{3},'Position',[115 160 90 20],'TooltipString','Channel used for ReSynch correction.',...
                    'Tag','rsChannel3',...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.rsChannel3,''Value'');','if tempVar ~=0,EPmain.average.rsChannel{3}=tempVar;end;','if isempty(tempVar),EPmain.average.rsChannel{3}=tempVar;end;','ep(''start'');']);
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Wind4','FontSize',EPmain.fontsize,...
                    'Position',[5 140 35 20]);
                
                EPmain.handles.average.rsMin4 = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.rsMin{4},'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.rsMin{4}=str2num(get(EPmain.handles.average.rsMin4,''String''));','ep(''start'');'],...
                    'Tag','rsMin4',...
                    'Position',[40 140 40 20],'TooltipString','Start of window in ms for ReSynch.');
                
                EPmain.handles.average.rsMax4 = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.rsMax{4},'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.rsMax{4}=str2num(get(EPmain.handles.average.rsMax4,''String''));','ep(''start'');'],...
                    'Tag','rsMax4',...
                    'Position',[80 140 40 20],'TooltipString','End of window in ms for ReSynch.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
                EPmain.handles.average.rsChannel4 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',['global';EPdataset.dataset(EPmain.average.rsChannelDataset).chanNames],...
                    'Value',EPmain.average.rsChannel{4},'Position',[115 140 90 20],'TooltipString','Channel used for ReSynch correction.',...
                    'Tag','rsChannel4',...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.rsChannel4,''Value'');','if tempVar ~=0,EPmain.average.rsChannel{4}=tempVar;end;','if isempty(tempVar),EPmain.average.rsChannel{4}=tempVar;end;','ep(''start'');']);
                
                uicontrol('Style','text','HorizontalAlignment','left','String', 'Wind5','FontSize',EPmain.fontsize,...
                    'Position',[5 120 35 20]);
                
                EPmain.handles.average.rsMin5 = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.rsMin{5},'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.rsMin{5}=str2num(get(EPmain.handles.average.rsMin5,''String''));','ep(''start'');'],...
                    'Tag','rsMin5',...
                    'Position',[40 120 40 20],'TooltipString','Start of window in ms for ReSynch.');
                
                EPmain.handles.average.rsMax5 = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.average.rsMax{5},'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','EPmain.average.rsMax{5}=str2num(get(EPmain.handles.average.rsMax5,''String''));','ep(''start'');'],...
                    'Tag','rsMax5',...
                    'Position',[80 120 40 20],'TooltipString','End of window in ms for ReSynch.  Latency range specified as left to right-sided (e.g., 1 second of 250 samples could be -200 to 800ms).');
                
                EPmain.handles.average.rsChannel5 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',['global';EPdataset.dataset(EPmain.average.rsChannelDataset).chanNames],...
                    'Value',EPmain.average.rsChannel{5},'Position',[115 120 90 20],'TooltipString','Channel used for ReSynch correction.',...
                    'Tag','rsChannel5',...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.rsChannel5,''Value'');','if tempVar ~=0,EPmain.average.rsChannel{5}=tempVar;end;','if isempty(tempVar),EPmain.average.rsChannel{5}=tempVar;end;','ep(''start'');']);
                
        end
        
        if ~isempty(EPmain.average.datasetList)
            theString=[EPdataset.dataset(EPmain.average.datasetList(EPmain.average.trialSpecDataset)).trialSpecNames;'no ACC'];
            EPmain.handles.average.ACCspec = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',theString,...
                'Value',find(strcmp(EPmain.average.ACCspec,theString)),'Position',[5 80 85 20],'TooltipString','Trial spec to be used for accuracy numbers.',...
                'Tag','ACCspec',...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.ACCspec,''Value'');','tempVar2=get(EPmain.handles.average.ACCspec,''String'');','if tempVar ~=0,EPmain.average.ACCspec=tempVar2{tempVar};end;','if isempty(tempVar),EPmain.average.ACCspec=tempVar;end;','ep(''start'');']);
        else
            EPmain.handles.average.ACCspec = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String','no ACC',...
                'Tag','ACCspec',...
                'Value',1,'Position',[5 80 85 20],'TooltipString','Trial spec to be used for ACC information.');
        end
        
        EPmain.handles.average.dropBad= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Drop Bad Trials',...
            'CallBack',['global EPmain;','EPmain.average.dropBad=get(EPmain.handles.average.dropBad,''Value'');','ep(''start'');'],'FontSize',EPmain.fontsize,...
            'Tag','dropBad',...
            'Value',EPmain.average.dropBad,'Position',[5 100 120 20],'TooltipString','Drop bad trials.');           

        EPmain.handles.average.dropError= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Error',...
            'CallBack',['global EPmain;','EPmain.average.dropError=get(EPmain.handles.average.dropError,''Value'');','ep(''start'');'],'FontSize',EPmain.fontsize,...
            'Tag','dropError',...
            'Value',EPmain.average.dropError,'Position',[90 80 50 20],'TooltipString','Drop trials where ACC is recorded as an error trial with a 0 value.');                      
        
        EPmain.handles.average.dropTimeout= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Time',...
            'CallBack',['global EPmain;','EPmain.average.dropTimeout=get(EPmain.handles.average.dropTimeout,''Value'');','ep(''start'');'],'FontSize',EPmain.fontsize,...
            'Tag','dropTimeout',...
            'Value',EPmain.average.dropTimeout,'Position',[140 80 70 20],'TooltipString','Drop trials where ACC is recorded as a timeout trial with a 2 value.');

        if isempty(EPmain.average.datasetList) || strcmp(EPmain.average.ACCspec,'no ACC')
%             set(EPmain.handles.average.dropBad,'enable','off');
            set(EPmain.handles.average.dropError,'enable','off');
            set(EPmain.handles.average.dropTimeout,'enable','off');
        end
        
        if ~isempty(EPmain.average.datasetList)
            theString=[EPdataset.dataset(EPmain.average.datasetList(EPmain.average.trialSpecDataset)).trialSpecNames;'no RT'];
            EPmain.handles.average.RTspec = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',theString,...
                'Value',find(strcmp(EPmain.average.RTspec,theString)),'Position',[5 60 85 20],'TooltipString','Trial spec to be used for RT numbers.',...
                'Tag','RTspec',...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.RTspec,''Value'');','tempVar2=get(EPmain.handles.average.RTspec,''String'');','if tempVar ~=0,EPmain.average.RTspec=tempVar2{tempVar};end;','if isempty(tempVar),EPmain.average.RTspec=tempVar;end;','ep(''start'');']);
        else
            EPmain.handles.average.RTspec = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String','no RT',...
                'Tag','RTspec',...
                'Value',1,'Position',[5 60 85 20],'TooltipString','Trial spec to be used for RT numbers.');
        end
                
        EPmain.handles.average.RTmethod= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Median','Mean','Trimmed Mean'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.RTmethod,''Value'');','if tempVar ~=0,EPmain.average.RTmethod=tempVar;end;','if isempty(tempVar),EPmain.average.RTmethod=tempVar;end;','ep(''start'');'],...
            'Tag','RTmethod',...
            'Value',EPmain.average.RTmethod,'Position',[90 60 105 20],'TooltipString','Method for summarizing the reaction times.');
        
        uicontrol('Style','text',...
            'String','Min RT','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 40 40 20]);
        
        EPmain.handles.average.minRT= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.average.minRT,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.minRT,''String'');','if tempVar ~=0,EPmain.average.minRT=tempVar;end;','if isempty(tempVar),EPmain.average.minRT=tempVar;end;','ep(''start'');'],...
            'Tag','minRT',...
            'Position',[50 40 40 20],'TooltipString','RT less than this number is discarded from reaction time summary number.  Set to zero to disable this setting.');
        
        uicontrol('Style','text',...
            'String','Max RT SD','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[95 40 65 20]);
        
        EPmain.handles.average.maxRT= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.average.maxRT,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.average.maxRT,''String'');','if tempVar ~=0,EPmain.average.maxRT=tempVar;end;','if isempty(tempVar),EPmain.average.maxRT=tempVar;end;','ep(''start'');'],...
            'Tag','maxRT',...
            'Position',[160 40 40 20],'TooltipString','RT more than this number in standard deviation units is discarded from reaction time summary number.  Set to zero to disable this setting.');
        
        if isempty(EPmain.average.datasetList) || strcmp(EPmain.average.RTspec,'no RT')
            set(EPmain.handles.average.RTmethod,'enable','off');
            set(EPmain.handles.average.minRT,'enable','off');
            set(EPmain.handles.average.maxRT,'enable','off');
        end
        
        EPmain.handles.average.average = uicontrol('Style', 'pushbutton', 'String', 'Average','FontSize',EPmain.fontsize,...
            'Tag','average',...
            'Position', [70 0 60 35], 'Callback', 'ep(''averageData'')');
        
        if (EPmain.average.method==2) && isempty(EPmain.average.datasetList)
            set(EPmain.handles.average.average,'enable','off');
        end
        
        EPmain.handles.average.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Tag','done',...
            'Position', [130 0 60 35], 'Callback', ['global EPmain;','EPmain.average=rmfield(EPmain.average,''trialSpecDataset'');','EPmain.mode=''main'';','ep(''start'');']);
        
    case 'averageData'
        
        EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);
        typeNum = get(EPmain.handles.average.fileType,'value');
        switch typeNum
            case 1
                dataType='continuous';
                dataOutType='average';
            case 2
                dataType='single_trial';
                dataOutType='average';
            case 3
                dataType='average';
                dataOutType='grand_average';
            case 4
                dataType='grand_average';
                dataOutType='grand_average';
            case 5
                dataType='factors';
                dataOutType='grand_average';
        end
        
        importFormatNum = get(EPmain.handles.average.importFormat,'value');        
        [~,importFormatName,importFormat]=ep_fileFormats(dataType,EPmain.fileFormatReadList{importFormatNum});

        outputFormatNum = get(EPmain.handles.average.outputFormat,'value');
        [~,~,outputFormat]=ep_fileFormats(dataOutType,EPmain.fileFormatSaveList{outputFormatNum});
        
        if strcmp(importFormat,'ep_mat')
            dataType=[];
        end
        
        procedureNum = get(EPmain.handles.average.method,'value');
        methodNum = get(EPmain.handles.average.freqMethod,'value');
        cfg=[];
        switch procedureNum
            case 1
                averagingMethod='Average';
                switch methodNum
                    case 1
                        methodName='Mean';
                    case 2
                        methodName='Median';
                    case 3
                        methodName='Trimmed_Mean';
                end
                suffix='_avg';
            case 2
                averagingMethod='Latency-Lock';
                switch methodNum
                    case 1
                        methodName='Mean';
                    case 2
                        methodName='Median';
                    case 3
                        methodName='Trimmed_Mean';
                end
                suffix='_avg';
                
                cfg.minLatency=[];
                cfg.maxLatency=[];
                cfg.latencyName=[];
                cfg.minLatency=EPmain.average.minLatency;
                cfg.maxLatency=EPmain.average.maxLatency;
                if minLatency > maxLatency
                    msg{1}='Latency minimum cannot be larger than the latency maximum.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
                cfg.latencyName=EPdataset.dataset(EPmain.average.datasetList(EPmain.average.trialSpecDataset)).trialSpecNames{EPmain.average.latencySpec};
            case 3
                averagingMethod='Jitter-Correct';
                switch methodNum
                    case 1
                        methodName='Mean';
                    case 2
                        methodName='Median';
                    case 3
                        methodName='Trimmed_Mean';
                end
                suffix='_avg';
                
                cfg.jitterChan=EPdataset.dataset(EPmain.average.channelDataset).chanNames{EPmain.average.channel};
                cfg.jitterPolar=EPmain.average.peakPolarity;
                cfg.minLatency=[];
                cfg.maxLatency=[];
                cfg.minLatency=EPmain.average.minLatency;
                cfg.maxLatency=EPmain.average.maxLatency;
                if minLatency > maxLatency
                    msg{1}='Latency minimum cannot be larger than the latency maximum.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
            case 4
                averagingMethod='Frequency-Coherence';
                switch methodNum
                    case 1
                        methodName='multi-taper';
                    case 2
                        methodName='Hanning';
                end
                suffix='_coh';
            case 5
                averagingMethod='Frequency-Phase Lock';
                methodName='phase-lock wavelet';
                suffix='_plv';
            case 6
                averagingMethod='RIDE';
                menuItems=get(EPmain.handles.average.freqMethod,'string');
                methodName=menuItems{EPmain.average.freqMethod};
                suffix='_avg';
                sMin=EPmain.average.sMin;
                sMax=EPmain.average.sMax;
                c1Min=EPmain.average.c1Min;
                c1Max=EPmain.average.c1Max;
                c2Min=EPmain.average.c2Min;
                c2Max=EPmain.average.c2Max;
                rMin=EPmain.average.rMin;
                rMax=EPmain.average.rMax;
                sLatency=EPmain.average.sLatency;
                
                cfg.high_cutoff=EPmain.average.cCutoff;
                
                if sMin > sMax
                    msg{1}='S minimum cannot be larger than the S maximum.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
                if c1Min > c1Max
                    msg{1}='C1 minimum cannot be larger than the C1 maximum.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
                if c2Min > c2Max
                    msg{1}='C2 minimum cannot be larger than the C2 maximum.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
                if rMin > rMax
                    msg{1}='R minimum cannot be larger than the R maximum.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
                
                switch methodNum
                    case 1 %'s-c1-r'
                        cfg.comp.name = {'s','c','r'};
                        cfg.comp.twd = {[sMin sMax],[c1Min c1Max],[rMin rMax]};
                        cfg.comp.latency = {sLatency,'unknown',[]};
                        cfg.dur={0,EPmain.average.cLag,0};
                    case 2 %'s-r'
                        cfg.comp.name = {'s','r'};
                        cfg.comp.twd = {[sMin sMax],[rMin rMax]};
                        cfg.comp.latency = {sLatency,[]};
                        cfg.dur={0,0};
                    case 3 %'s-c1'
                        cfg.comp.name = {'s','c'};
                        cfg.comp.twd = {[sMin sMax],[c1Min c1Max]};
                        cfg.comp.latency = {sLatency,'unknown'};
                        cfg.dur={0,EPmain.average.cLag};
                    case 4 %'s-c1-c2-r'
                        cfg.comp.name = {'s','c1','c2','r'};
                        cfg.comp.twd = {[sMin sMax],[c1Min c1Max],[c2Min c2Max],[rMin rMax]};
                        cfg.comp.latency = {sLatency,'unknown','unknown',[]};
                        cfg.dur={0,EPmain.average.cLag,EPmain.average.cLag,0};
                    case 5 %'s-c1-c2'
                        cfg.comp.name = {'s','c1','c2'};
                        cfg.comp.twd = {[sMin sMax],[c1Min c1Max],[c2Min c2Max]};
                        cfg.comp.latency = {sLatency,'unknown','unknown'};
                        cfg.dur={0,EPmain.average.cLag,EPmain.average.cLag};
                    otherwise
                        disp('oops');
                end
            case 7
                averagingMethod='ReSynch';
                switch methodNum
                    case 1
                        methodName='Mean';
                    case 2
                        methodName='Median';
                    case 3
                        methodName='Trimmed_Mean';
                end
                cfg.rsChannel{1}=EPmain.average.rsChannel{1};
                cfg.rsMin{1}=EPmain.average.rsMin{1};
                cfg.rsMax{1}=EPmain.average.rsMax{1};
                cfg.rsChannel{2}=EPmain.average.rsChannel{2};
                cfg.rsMin{2}=EPmain.average.rsMin{2};
                cfg.rsMax{2}=EPmain.average.rsMax{2};
                cfg.rsChannel{3}=EPmain.average.rsChannel{3};
                cfg.rsMin{3}=EPmain.average.rsMin{3};
                cfg.rsMax{3}=EPmain.average.rsMax{3};
                cfg.rsChannel{4}=EPmain.average.rsChannel{4};
                cfg.rsMin{4}=EPmain.average.rsMin{4};
                cfg.rsMax{4}=EPmain.average.rsMax{4};
                cfg.rsChannel{5}=EPmain.average.rsChannel{5};
                cfg.rsMin{5}=EPmain.average.rsMin{5};
                cfg.rsMax{5}=EPmain.average.rsMax{5};
                suffix='_avg';
        end

        cfg.behav.ACC=EPmain.average.ACCspec;
        cfg.behav.RT=EPmain.average.RTspec;
        cfg.behav.codeCorrect=EPmain.average.codeCorrect;
        cfg.behav.codeError=EPmain.average.codeError;
        cfg.behav.codeTimeout=EPmain.average.codeTimeout;
        cfg.behav.dropBad=EPmain.average.dropBad;
        cfg.behav.dropError=EPmain.average.dropError;
        cfg.behav.dropTimeout=EPmain.average.dropTimeout;
        cfg.behav.minRT=EPmain.average.minRT;
        cfg.behav.maxRT=EPmain.average.maxRT;
        theString=get(EPmain.handles.average.RTmethod,'String');
        cfg.behav.RTmethod=theString{EPmain.average.RTmethod};
  
        switch EPmain.average.averageType
            case 1
                averageType='subject';
            case 2
                averageType='item';
        end
        
        smoothing=0;
        if methodNum > 3 %Frequency domain analysis
            smoothing=str2num(get(EPmain.handles.average.smoothing,'String'));
        end

        EPmain.average.importFormat=importFormatNum;
        EPmain.average.fileType=typeNum;
        EPmain.average.outputFormat=outputFormatNum;
        
        multiCellNumber=[];
        multiSubjectNumber=[];
        multiSessionNumber=[];
        if strcmp(averageType,'subject')
            if ~isempty(EPmain.average.cell)
                if ~isempty(findstr('-',EPmain.average.cell))
                    msg{1}='The cell field cannot have a dash.  If the intent was to specify a range of numbers, a colon is needed.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
                multiCellNumber=str2num(EPmain.average.cell);
                if ~isnumeric(multiCellNumber)
                    msg{1}='The cell field needs to be a number or numbers or empty.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
            end
            if ~isempty(EPmain.average.subject)
                if ~isempty(findstr('-',EPmain.average.subject))
                    msg{1}='The subject field cannot have a dash.  If the intent was to specify a range of numbers, a colon is needed.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
                multiSubjectNumber=str2num(EPmain.average.subject);
                if ~isnumeric(multiSubjectNumber)
                    msg{1}='The subject field needs to be a number or numbers or empty.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
            end
            if ~isempty(EPmain.average.session)
                if ~isempty(findstr('-',EPmain.average.session))
                    msg{1}='The session field cannot have a dash.  If the intent was to specify a range of numbers, a colon is needed.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
                multiSessionNumber=str2num(EPmain.average.session);
                if ~isnumeric(multiSessionNumber)
                    msg{1}='The session field needs to be a number or numbers or empty.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
            end
        end

        [sessionFiles, activeDirectory]=ep_getFilesUI(importFormat);
        if isempty(sessionFiles)
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
        if (sessionFiles{1}==0) & ~ischar(sessionFiles{1})
            msg{1}='No filenames selected. You have to click on a name.';
            [msg]=ep_errorMsg(msg);
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
        for iFile=1:size(sessionFiles,2)
            if ~isempty(multiCellNumber)
                [pathstr, name, ext] = fileparts(sessionFiles{iFile});
                if max(multiCellNumber) > length(name)
                    msg{1}=['The file name ' name ' is shorter than specified by the Cell field.'];
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
            end
            if ~isempty(multiSubjectNumber)
                [pathstr, name, ext] = fileparts(sessionFiles{iFile});
                if max(multiSubjectNumber) > length(name)
                    msg{1}=['The file name ' name ' is shorter than specified by the Subject field.'];
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
            end
            if ~isempty(multiSessionNumber)
                [pathstr, name, ext] = fileparts(sessionFiles{iFile});
                if max(multiSessionNumber) > length(name)
                    msg{1}=['The file name ' name ' is shorter than specified by the Session field.'];
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
            end
            sessionFiles{iFile}=[activeDirectory sessionFiles{iFile}];
        end
        
        sessionFiles=sort(sessionFiles);
        
        [pathstr, name, ext] = fileparts(sessionFiles{1});
        
        if strcmp(name(end-3:end),'_seg')
            name=name(1:end-4);
        end
        outFileName=[pathstr filesep name suffix ext];

        [outFileName, pathname] = uiputfile('*.*','Save:',outFileName);
        if outFileName == 0
            msg{1}='No output name selected.';
            [msg]=ep_errorMsg(msg);
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
        outFileNamePath=[pathname outFileName];
        [pathstr, name, ext] = fileparts(outFileNamePath);
        if exist(outFileNamePath,'file')
            delete(outFileNamePath); %user must have clicked "yes" to whether to replace existing file
        end
        
        if strcmp(averageType,'item')
            itemSpec=EPdataset.dataset(EPmain.average.datasetList(EPmain.average.trialSpecDataset)).trialSpecNames{EPmain.average.itemSpec};
        else
            itemSpec=[];
        end
        averagedData=ep_averageData(sessionFiles,importFormat,dataType,averagingMethod,EPmain.preferences.average.trimLevel,methodName,smoothing,multiSessionNumber,multiSubjectNumber,multiCellNumber,cfg,EPmain.preferences,averageType,itemSpec,EPmain.average.dropSD,EPmain.average.dropNoise);
        if EPtictoc.stop
            EPtictoc.stop=0;
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            ep('start');
        else
            if isfield(averagedData,'data') && ~isempty(averagedData.data)
                if EPmain.average.dropEvents
                    averagedData.events=cell(size(averagedData.events));
                end
                if isempty(EPmain.average.subject) && (length(averagedData.subNames)==1) && (size(sessionFiles,2) > 1)
                    disp('Note that you specified multiple files to be averaged but they were all assigned to the same subject.  Did you forget to use the Subject field?');
                end
                averagedData.dataName=name;
                
                if ~strcmp(importFormat,'ep_mat')
                    if length(sessionFiles)==1
                        theDescription=['Imported the file ' sessionFiles{1} '.'];
                    else
                        theDescription=['Imported the files.'];
                    end
                    averagedData.history=ep_addHistory(averagedData.history, EPmain.preferences.records, theDescription, EPmain.handles.hMainWindow, EPmain.preferences, EPmain.preferences.EPver,sessionFiles);
                end
                theDescription=['Averaged the data.'];
                averagedData.history=ep_addHistory(averagedData.history, EPmain.preferences.records, theDescription, EPmain.handles.hMainWindow, EPmain.preferences, EPmain.preferences.EPver,sessionFiles);
                [err]=ep_writeData(averagedData,outFileNamePath,EPmain.preferences.general.specSuffix,EPmain.preferences.general.subjectSpecSuffix,outputFormat);
                ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);ep('start');return;end
            end
            try
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            catch
                ep('start')
            end
            
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            ep('start');
        end
        
    case 'startRead'
        
        set(EPmain.handles.hMainWindow,'Name', 'Read Data');
        
        EPmain.handles.read.prefs = uicontrol('Style', 'pushbutton', 'String', '•','FontSize',EPmain.fontsize,...
            'Tag','prefs',...
            'Position', [5 485 10 10], 'Callback', ['global EPmain;','EPmain.prefReturn=''read'';','EPmain.mode=''preferenceGeneral'';','ep(''start'');']);
        
        uicontrol('Style','text',...
            'String','Format','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 430 50 20]);
        
        EPmain.handles.read.format = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatReadList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.read.format,''Value'');','if tempVar ~=0,EPmain.read.format=tempVar;end;','if isempty(tempVar),EPmain.read.format=tempVar;end;','ep(''start'');'],...
            'Tag','format',...
            'Value',EPmain.read.format,'Position',[50 430 150 20]);
        
        uicontrol('Style','text',...
            'String','Type','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 410 50 20]);
        
        EPmain.handles.read.fileType= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'continuous','single_trial','average','grand_average','factors'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.read.fileType,''Value'');','if tempVar ~=0,EPmain.read.fileType=tempVar;end;','if isempty(tempVar),EPmain.read.fileType=tempVar;end;','ep(''start'');'],...
            'Tag','fileType',...
            'Value',EPmain.read.fileType,'Position',[50 410 150 20]);
        
        uicontrol('Style','text',...
            'String','Mont','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[5 390 50 20]);
        
        EPmain.handles.read.importMontage = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.montageList,...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.read.importMontage,''Value'');','if tempVar ~=0,EPmain.read.importMontage=EPmain.montageList{tempVar};end;','if isempty(tempVar),EPmain.read.importMontage=EPmain.montageList{tempVar};end'],...
            'Tag','importMontage',...
            'Value',find(strcmp(EPmain.read.importMontage,EPmain.montageList)),'Position',[50 390 150 20]);
        
        [importSuffix,importFormatName,importFormat]=ep_fileFormats('average',EPmain.fileFormatReadList{EPmain.read.format});
        if strcmp(importFormat,'ep_mat')
            set(EPmain.handles.read.fileType,'enable','off');
            set(EPmain.handles.read.importMontage,'enable','off');
        end
        
        %         EPmain.handles.read.ced= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
%             'String','CED override',...
%             'CallBack',['global EPmain;','EPmain.read.ced=get(EPmain.handles.read.ced,''Value'');','ep(''start'');'],'FontSize',EPmain.fontsize,...
%             'Value',EPmain.read.ced,'Position',[10 370 150 20],'TooltipString','Override file''s internal electrode coordinates with CED file.');        
                
        uicontrol('Style','frame',...
            'Position',[5 265 195 90]);
        
        EPmain.handles.read.check= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Single File Mode',...
            'CallBack',['global EPmain;','EPmain.read.check=get(EPmain.handles.read.check,''Value'');','ep(''start'');'],'FontSize',EPmain.fontsize,...
            'Tag','check',...
            'Value',EPmain.read.check,'Position',[10 330 150 20]);
        
        if EPmain.read.check
            
            EPmain.handles.read.subject= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.read.subject,...
                'CallBack',['global EPmain;','EPmain.read.subject=get(EPmain.handles.read.subject,''String'');','ep(''start'');'],...
                'Tag','subject',...
                'Position',[10 310 45 20],'TooltipString','example 4:6');
            
            EPmain.handles.read.subjectLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Subject','HorizontalAlignment','left',...
                'Position',[65 310 45 20]);
            
            if isempty(EPmain.read.subject)
                set(EPmain.handles.read.subjectLabel,'enable','off');
            elseif isempty(str2num(EPmain.read.subject))
                set(EPmain.handles.read.subjectLabel,'enable','off');
            end
            
            EPmain.handles.read.cell= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.read.cell,...
                'CallBack',['global EPmain;','EPmain.read.cell=get(EPmain.handles.read.cell,''String'');','ep(''start'');'],...
                'Tag','cell',...
                'Position',[10 290 45 20],'TooltipString','example 7:9');
            
            EPmain.handles.read.cellLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Cell','HorizontalAlignment','left',...
                'Position',[65 290 45 20]);
            
            if isempty(EPmain.read.cell)
                set(EPmain.handles.read.cellLabel,'enable','off');
            elseif isempty(str2num(EPmain.read.cell))
                set(EPmain.handles.read.cellLabel,'enable','off');
            end
            
            EPmain.handles.read.session= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.read.session,...
                'CallBack',['global EPmain;','EPmain.read.session=get(EPmain.handles.read.session,''String'');','ep(''start'');'],...
                'Tag','session',...
                'Position',[110 290 50 20],'TooltipString','example 7:9');
            
            EPmain.handles.read.sessLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Sess','HorizontalAlignment','left',...
                'Position',[160 290 35 20]);
            
            if isempty(EPmain.read.session)
                set(EPmain.handles.read.sessLabel,'enable','off');
            elseif isempty(str2num(EPmain.read.session))
                set(EPmain.handles.read.sessLabel,'enable','off');
            end

            EPmain.handles.read.freq= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.read.freq,...
                'CallBack',['global EPmain;','EPmain.read.freq=get(EPmain.handles.read.freq,''String'');','ep(''start'');'],...
                'Tag','freq',...
                'Position',[10 270 45 20],'TooltipString','example 10:12');
            
            EPmain.handles.read.freqLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Freq','HorizontalAlignment','left',...
                'Position',[65 270 45 20]);
            
            if isempty(EPmain.read.freq)
                set(EPmain.handles.read.freqLabel,'enable','off');
            elseif isempty(str2num(EPmain.read.freq))
                set(EPmain.handles.read.freqLabel,'enable','off');
            end
            
            if (EPmain.read.fileType==2) || strcmp(importFormat,'ep_mat')
                EPmain.handles.read.trial= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                    'String',EPmain.read.trial,...
                    'CallBack',['global EPmain;','EPmain.read.trial=get(EPmain.handles.read.trial,''String'');','ep(''start'');'],...
                    'Tag','trial',...
                    'Position',[110 310 50 20],'TooltipString','example 10:12');
                
                EPmain.handles.read.trialLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Trial','HorizontalAlignment','left',...
                    'Position',[160 310 30 20]);
                
                if isempty(EPmain.read.trial)
                    set(EPmain.handles.read.trialLabel,'enable','off');
                elseif isempty(str2num(EPmain.read.trial))
                    set(EPmain.handles.read.trialLabel,'enable','off');
                end
            end
        end
        
        if ~isempty(EPdataset.dataset)
            for i=1:length(EPdataset.dataset)
                fileName=EPdataset.dataset(i).dataName;
                if strcmp(EPdataset.dataset(i).saved,'no')
                    fileName=['*' fileName];
                end
                tableData{i,1}=fileName;
            end
        else
            tableData=[];
        end
        
        tableNames{1}='data';
        columnEditable =  false;
        ColumnFormat{1}=[];
        
        EPmain.handles.read.readTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
            'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
            'CellSelectionCallback',@deleteReadData,'ForegroundColor','red',...
            'Tag','readTable',...
            'ColumnWidth',{300},'Position',[20 100 170 150]);
        
        EPmain.handles.read.read = uicontrol('Style', 'pushbutton', 'String', 'Read','FontSize',EPmain.fontsize,...
            'Tag','read',...
            'Position', [20 50 100 40], 'Callback', 'ep(''readData'')');
        
        EPmain.handles.read.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Tag','done',...
            'Position', [20 0 100 40], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        
    case 'readData'

        EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);
        
        typeNum = get(EPmain.handles.read.fileType,'value');
        switch typeNum
            case 1
                dataType='continuous';
            case 2
                dataType='single_trial';
            case 3
                dataType='average';
            case 4
                dataType='grand_average';
        end
        
        formatNum = get(EPmain.handles.read.format,'value');
        [importSuffix,importFormatName,importFormat]=ep_fileFormats(dataType,EPmain.fileFormatReadList{formatNum});
        
        if strcmp(importFormat,'ep_mat')
            dataType=''; %data type set by data file itself
        end
        
        EPmain.read.format=formatNum;
        EPmain.read.fileType=typeNum;
        
        EPmain.convertMode=0;
        theHandles=EPmain.handles.read;
        readFiles(theHandles,importFormat,dataType);

        try
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        catch
            ep('start')
        end
        
        ep('start');
        
    case 'startEdit'
        
        if isempty(EPoverview) || EPoverview.done
            doneFlag=1;
        else
            doneFlag=0;
        end
        EPoverview=[];
        tableData=[];
        
        if ~isempty(EPdataset.dataset)
            for i=1:length(EPdataset.dataset)
                fileName=EPdataset.dataset(i).dataName;
                if strcmp(EPdataset.dataset(i).saved,'no')
                    fileName=['*' fileName];
                end
                tableData{i,1}=fileName;
            end
        else
            tableData=[];
        end
        
        tableNames{1}='data';
        columnEditable =  false;
        ColumnFormat{1}=[];
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Click on dataset to edit.','FontSize',EPmain.fontsize,...
            'Position',[20 250 160 20]);
        
        EPmain.handles.edit.editTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
            'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
            'CellSelectionCallback',@pickEditData,...
            'Tag','editTable',...
            'ColumnWidth',{300},'Position',[20 100 170 150]);
        
        EPmain.handles.edit.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Tag','done',...
            'Position', [20 0 100 40], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        
        if isscalar(EPdataset.dataset) && ~doneFlag %if just one dataset, don't need to wait for it to be selected.
            EPoverview.dataset=1;
            EPoverview.mode=[];
            ep_editData
        end
        
    case 'startView'
        ep_viewPane;

    case 'viewWaves'
        EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);
        
        for iColor=1:EPmain.numColors 
            EPmain.view.dataset(iColor)=get(EPmain.handles.view.dataset(iColor),'value');
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                EPmain.view.cell(iColor)=get(EPmain.handles.view.cell(iColor),'value');
                EPmain.view.subject(iColor)=get(EPmain.handles.view.subject(iColor),'value');
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames)
                    EPmain.view.trial(iColor)=get(EPmain.handles.view.trial(iColor),'value');
                end
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)
                    EPmain.view.factor(iColor)=get(EPmain.handles.view.factor(iColor),'value');
                end
            end
        end
        
        EPwaves=[];
        EPwaves.mode=[];
        EPwaves.direction=EPmain.preferences.view.positive;
        
        eventLines=ep_collateEventLines('EPwaves');

        if EPmain.view.edited.bottomVolt
            theMin=str2num(get(EPmain.handles.view.bottomVolt,'string'));
        else
            theMin=[];
        end
        if EPmain.view.edited.topVolt
            theMax=str2num(get(EPmain.handles.view.topVolt,'string'));
        else
            theMax=[];
        end
        ep_showWaves(theMin,theMax,...
            str2num(get(EPmain.handles.view.startSamp,'string')),str2num(get(EPmain.handles.view.endSamp,'string')),...
            str2num(get(EPmain.handles.view.startHz,'string')),str2num(get(EPmain.handles.view.endHz,'string')),...
            str2num(get(EPmain.handles.view.marker1,'string')),str2num(get(EPmain.handles.view.marker2,'string')),EPmain.view.FFTunits,eventLines,1);
        if EPtictoc.stop
            if ishandle(EPwaves.handles.waves.hWaveWindow)
                close(EPwaves.handles.waves.hWaveWindow);
            end
            EPtictoc.stop=0;
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            ep('start');
        end

        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        ep('start')

    case 'viewTopos'
        EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);
        
        drawnow;
        
        for iColor=1:EPmain.numColors
            EPmain.view.dataset(iColor)=get(EPmain.handles.view.dataset(iColor),'value');
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                EPmain.view.cell(iColor)=get(EPmain.handles.view.cell(iColor),'value');
                EPmain.view.subject(iColor)=get(EPmain.handles.view.subject(iColor),'value');
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames)
                    EPmain.view.trial(iColor)=get(EPmain.handles.view.trial(iColor),'value');
                end
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)
                    EPmain.view.factor(iColor)=get(EPmain.handles.view.factor(iColor),'value');
                end
            end
        end
        
        eventLines=ep_collateEventLines('EPtopos');
        
        EPtopos=[];
        EPtopos.page=0;
        EPtopos.direction=EPmain.preferences.view.positive;
        if EPmain.view.edited.bottomVolt
            theMin=str2num(get(EPmain.handles.view.bottomVolt,'string'));
        else
            theMin=[];
        end
        if EPmain.view.edited.topVolt
            theMax=str2num(get(EPmain.handles.view.topVolt,'string'));
        else
            theMax=[];
        end
        
        ep_showTopos(str2num(get(EPmain.handles.view.startSamp,'string')),str2num(get(EPmain.handles.view.endSamp,'string')),...
            str2num(get(EPmain.handles.view.startHz,'string')),str2num(get(EPmain.handles.view.endHz,'string')),...
            str2num(get(EPmain.handles.view.marker1,'string')),str2num(get(EPmain.handles.view.marker2,'string')),EPmain.view.FFTunits,theMin,theMax,eventLines);
        if EPtictoc.stop
            EPtictoc.stop=0;
            close(EPtopos.handles.topos.topoWindow);
        end
        

    case 'viewEdit'

        if EPmain.view.dataset(1) > length(EPdataset.dataset)
            disp('Scan function requires the first color to not be set to none so that it can specify the dataset to be edited.')
            return
        end

        EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);
        
        for iColor=1:EPmain.numColors
            EPmain.view.dataset(iColor)=get(EPmain.handles.view.dataset(iColor),'value');
            if EPmain.view.dataset(iColor) <= length(EPdataset.dataset)
                EPmain.view.cell(iColor)=get(EPmain.handles.view.cell(iColor),'value');
                EPmain.view.subject(iColor)=get(EPmain.handles.view.subject(iColor),'value');
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).trialNames)
                    EPmain.view.trial(iColor)=get(EPmain.handles.view.trial(iColor),'value');
                end
                if ~isempty(EPdataset.dataset(EPmain.view.dataset(iColor)).facNames)
                    EPmain.view.factor(iColor)=get(EPmain.handles.view.factor(iColor),'value');
                end
            end
        end
        
        eventLines=ep_collateEventLines('EPwaves');
        %if an event type is selected, then there will at least be a flat line for the event.
        if length(EPmain.view.events) > 1
            for iColor=1:EPmain.numColors
                if isempty(eventLines{iColor})
                    eventLines{iColor}=0;
                end
            end
        end
        
        EPwaves.mode=[];
        EPwaves.direction=EPmain.preferences.view.positive;
        if strcmp(EPdataset.dataset(EPmain.view.dataset(1)).dataType,'continuous')
            ep_scanCont;
            if EPtictoc.stop
                EPtictoc.stop=0;
                if isfield(EPmain.handles,'scanContData') && ishandle(EPmain.handles.scanContData)
                    close(EPmain.handles.scanContData);
                end
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                ep('start');
            end
        else
            EPmain.view.midFigure=1;
            err=ep_showWaves(str2num(get(EPmain.handles.view.bottomVolt,'string')),str2num(get(EPmain.handles.view.topVolt,'string')),...
                str2num(get(EPmain.handles.view.startSamp,'string')),str2num(get(EPmain.handles.view.endSamp,'string')),...
                str2num(get(EPmain.handles.view.startHz,'string')),str2num(get(EPmain.handles.view.endHz,'string')),...
                str2num(get(EPmain.handles.view.marker1,'string')),str2num(get(EPmain.handles.view.marker2,'string')),EPmain.view.FFTunits,eventLines,0);
            
            if ~err %if the wave window did not error out
                EPmanualEdit=[];
                ep_manualEdit
            else
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            end
        end
        
    case 'startSampleTest'
        
        ep_tictoc('begin');
        
        set(EPmain.handles.hMainWindow,'Name', 'sampleTest');
        noSampTest=0;
        
        %initialize sampleTest parameters if entering in from the Main Menu either for the first time or list has changed.
        if ~isfield(EPmain.sampleTest,'dataset') ||  (EPmain.sampleTest.dataset > length(EPdataset.dataset)) || ~isfield(EPmain.sampleTest,'datasetName') || ~strcmp(EPmain.sampleTest.datasetName,EPdataset.dataset(EPmain.sampleTest.dataset).dataName) ||...
                (~isempty(EPmain.sampleTest.lastChange) && ~strcmp(EPmain.sampleTest.lastChange,EPdataset.dataset(EPmain.sampleTest.dataset).lastChange))
            EPmain.sampleTest.datasetList=[];
            EPmain.sampleTest.PCAlist=[];
            EPmain.sampleTest.AVElist=[];
            for iFile=1:length(EPdataset.dataset)
                if ~isempty(EPdataset.dataset(iFile).timeNames) && isempty(EPdataset.dataset(iFile).facVecT) && ~isempty(EPdataset.dataset(iFile).eloc) %cannot perform analysis on FFT data or temporal PCA data and need electrode coordinates.
                    switch EPdataset.dataset(iFile).dataType
                        case 'continuous'
                            EPmain.sampleTest.datasetList(end+1)=iFile;
                        case 'single_trial'
                            if length(unique(EPdataset.dataset(iFile).cellNames(find(strcmp(EPdataset.dataset(iFile).cellTypes,'SGL'))))) > 1
                                EPmain.sampleTest.datasetList(end+1)=iFile;
                            end
                        case 'average'
                            if length(EPdataset.dataset(iFile).subNames(find(strcmp(EPdataset.dataset(iFile).subTypes,'AVG')))) > 1
                                EPmain.sampleTest.datasetList(end+1)=iFile;
                            end
                    end
                end
            end
            EPmain.sampleTest.dataset=EPmain.sampleTest.datasetList(end); %the active dataset of EPdataset
            EPmain.sampleTest.datasetName=EPdataset.dataset(EPmain.sampleTest.dataset).dataName;
            theDataset=EPmain.sampleTest.dataset;
            if isfield(EPdataset.dataset(EPmain.sampleTest.dataset),'lastChange') && ~isempty(EPdataset.dataset(EPmain.sampleTest.dataset).lastChange)
                EPmain.sampleTest.lastChange=EPdataset.dataset(EPmain.sampleTest.dataset).lastChange;
            else
                EPmain.sampleTest.lastChange=cell(0);
            end
            EPmain.sampleTest.cell1=1;
            EPmain.sampleTest.cell2=2;
            EPmain.sampleTest.sub1=1;
            EPmain.sampleTest.sub2=1;
            if EPmain.sampleTest.dataset
                EPmain.sampleTest.cellNameList=unique(EPdataset.dataset(theDataset).cellNames(find(ismember(EPdataset.dataset(theDataset).cellTypes,{'SGL','CMB'}))));
                
                EPmain.sampleTest.subNameList=cell(0);
                EPmain.sampleTest.subList=cell(0);
                EPmain.sampleTest.subNameList{1}='-all-';
                EPmain.sampleTest.subList{1}=find(strcmp(EPdataset.dataset(theDataset).subTypes,'AVG'));
                for iSpec=1:length(EPdataset.dataset(theDataset).subjectSpecNames)
                    specLvLlist=cell(0);
                    for iSub=1:length(EPdataset.dataset(theDataset).subjectSpecs(:,iSpec))
                        theSpec=EPdataset.dataset(theDataset).subjectSpecs{iSub,iSpec};
                        if iscolumn(theSpec)
                            theSpec=theSpec';
                        end
                        if ~isempty(theSpec) && (isempty(specLvLlist) || ~any(strcmp(theSpec,specLvLlist)))
                            specLvLlist{end+1}=theSpec;
                            EPmain.sampleTest.subNameList{end+1}=[EPdataset.dataset(theDataset).subjectSpecNames{iSpec} '_' specLvLlist{end}];
                            EPmain.sampleTest.subList{end+1}=find(strcmp(specLvLlist{end},EPdataset.dataset(theDataset).subjectSpecs(:,iSpec)));
                        end
                    end
                end
                if ~isempty(EPdataset.dataset(theDataset).sessNums)
                    sessList=unique(EPdataset.dataset(theDataset).sessNums);
                    sessList=sessList(sessList>0);
                    EPmain.sampleTest.sessNameList=EPdataset.dataset(theDataset).sessNames(sessList);
                else
                    EPmain.sampleTest.sessNameList=cell(0);
                    EPmain.sampleTest.sess2=1;
                end
                EPmain.sampleTest.sess1=1;
                EPmain.sampleTest.sess2=min(2,length(EPmain.sampleTest.sessNameList));
                EPmain.sampleTest.datasetNameList={EPdataset.dataset(EPmain.sampleTest.datasetList).dataName}';
            else
                EPmain.sampleTest.cellNameList=cell(0);
                EPmain.sampleTest.subNameList=cell(0);
                EPmain.sampleTest.sessNameList=cell(0);
                EPmain.sampleTest.datasetNameList=cell(0);
            end
            EPmain.sampleTest.method=1;
            EPmain.sampleTest.contMethod=1;
            EPmain.sampleTest.test=1;
            EPmain.sampleTest.PCA=1;
            EPmain.sampleTest.AVE=1;
            EPmain.sampleTest.alpha=.05;
            EPmain.sampleTest.contiguous=4;
            EPmain.sampleTest.thresh=.5;
            EPmain.sampleTest.channel=1;
            EPmain.sampleTest.filterPass=1;
            EPmain.sampleTest.filterType=2;
            EPmain.sampleTest.filter1=[];
            EPmain.sampleTest.filter2=[];
            EPmain.sampleTest.filterOrder=6;
            EPmain.sampleTest.scale=300;
            EPmain.sampleTest.factor=1;
            EPmain.sampleTest.subject=1;
            EPmain.sampleTest.cell=1;
            EPmain.sampleTest.gridSize=67;
            EPmain.sampleTest.AVEdata=cell(0);
            EPmain.sampleTest.sampStart=[];
            EPmain.sampleTest.sampEnd=[];
            EPmain.sampleTest.changeFlag=1;
            EPmain.sampleTest.freqFlag=0;
            EPmain.sampleTest.freqBin=1;
            
            if EPmain.sampleTest.dataset
                EEGchans=find(strcmp('EEG',EPdataset.dataset(EPmain.sampleTest.dataset).chanTypes));
                EPmain.sampleTest.eloc=EPdataset.dataset(EPmain.sampleTest.dataset).eloc(EEGchans);
                EPmain.sampleTest.dataType=EPdataset.dataset(EPmain.sampleTest.dataset).dataType;
                EPmain.sampleTest.freqFlag=~isempty(EPdataset.dataset(EPmain.sampleTest.dataset).freqNames);
                
                maxRad=0.5;
                [y,x] = pol2cart(([EPmain.sampleTest.eloc.theta]/360)*2*pi,[EPmain.sampleTest.eloc.radius]);  % transform electrode locations from polar to cartesian coordinates
                y=-y; %flip y-coordinate so that nose is upwards.
                plotrad = min(1.0,max([EPmain.sampleTest.eloc.radius])*1.02);            % default: just outside the outermost electrode location
                plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
                x = x*(maxRad/plotrad);
                y = y*(maxRad/plotrad);
                
                xmin = min(-maxRad,min(x));
                xmax = max(maxRad,max(x));
                ymin = min(-maxRad,min(y));
                ymax = max(maxRad,max(y));
                
                EPmain.sampleTest.x=round(((x/(xmax-xmin))*EPmain.sampleTest.gridSize)+ceil(EPmain.sampleTest.gridSize/2));
                EPmain.sampleTest.y=round(((y/(ymax-ymin))*EPmain.sampleTest.gridSize)+ceil(EPmain.sampleTest.gridSize/2));
                
                for iFile=1:length(EPdataset.dataset)
                    if ~isempty(EPdataset.dataset(iFile).facNames) && ~isempty(EPdataset.dataset(iFile).eloc)
                        newEloc=EPdataset.dataset(iFile).eloc(find(strcmp('EEG',EPdataset.dataset(iFile).chanTypes)));
                        if length(EPmain.sampleTest.eloc) == length(newEloc)
                            if isequal(cellfun(@isempty,{EPmain.sampleTest.eloc.theta}),cellfun(@isempty,{newEloc.theta})) && ~any([EPmain.sampleTest.eloc.theta]-[newEloc.theta]) && ~any([EPmain.sampleTest.eloc.radius]-[newEloc.radius])
                                if (length(EPdataset.dataset(iFile).timeNames) <= length(EPdataset.dataset(EPmain.sampleTest.dataset).timeNames))
                                    if (length(EPdataset.dataset(iFile).freqNames)==length(EPdataset.dataset(EPmain.sampleTest.dataset).freqNames)) && all([EPdataset.dataset(iFile).freqNames] == [EPdataset.dataset(EPmain.sampleTest.dataset).freqNames])
                                        EPmain.sampleTest.PCAlist(end+1)=iFile;
                                    end
                                end
                            end
                        end
                    end
                end
                EPmain.sampleTest.PCAnameList={EPdataset.dataset(EPmain.sampleTest.PCAlist).dataName};
                
                disp('Loading in the compatible data in the working set.  The more data you have, the longer this will take.');
                for iFile=1:length(EPdataset.dataset)
                    if isempty(EPdataset.dataset(iFile).facNames) && ~isempty(EPdataset.dataset(iFile).eloc) && strcmp(EPdataset.dataset(iFile).dataType,'average')
                        newEloc=EPdataset.dataset(iFile).eloc(find(strcmp('EEG',EPdataset.dataset(iFile).chanTypes)));
                        if length(EPmain.sampleTest.eloc) == length(newEloc)
                            if isequal(cellfun(@isempty,{EPmain.sampleTest.eloc.theta}),cellfun(@isempty,{newEloc.theta})) && ~any([EPmain.sampleTest.eloc.theta]-[newEloc.theta]) && ~any([EPmain.sampleTest.eloc.radius]-[newEloc.radius])
                                if (length(EPdataset.dataset(iFile).timeNames) <= length(EPdataset.dataset(EPmain.sampleTest.dataset).timeNames))
                                    if (length(EPdataset.dataset(iFile).freqNames)==length(EPdataset.dataset(EPmain.sampleTest.dataset).freqNames)) && all([EPdataset.dataset(iFile).freqNames] == [EPdataset.dataset(EPmain.sampleTest.dataset).freqNames])
                                        EPmain.sampleTest.AVElist(end+1)=iFile;
                                        EPdata=ep_loadEPdataset(iFile);
                                        %convert virtual GAVEs to normal form so subject specs etc are available.
                                        EPdata=ep_combineData(EPdata,'convert',{[],[],[],[],[],[]},[],[],[]);
                                        if EPtictoc.stop;EPtictoc.stop=0;return;end
                                        if isempty(EPdata)
                                            disp('Error: sampleTest failed.')
                                            return
                                        end
                                        EPmain.sampleTest.AVEdata{end+1}=EPdata.data;
                                    end
                                end
                            end
                        end
                    end
                end
                EPmain.sampleTest.AVEnameList={EPdataset.dataset(EPmain.sampleTest.AVElist).dataName};
            end
        else
            EEGchans=[];
            EPmain.sampleTest.eloc=[];
        end
        
        if EPmain.sampleTest.dataset
            EPmain.handles.sampleTest.dataset = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',EPmain.sampleTest.datasetNameList,...
                'Value',find(EPmain.sampleTest.dataset==EPmain.sampleTest.datasetList),'Position',[5 480 160 20],...
                'Callback', @changeSampleTestDataset);
            
            if ~strcmp(EPdataset.dataset(EPmain.sampleTest.dataset).dataType,'continuous')
                if ~isempty(EPmain.sampleTest.cellNameList)
                    EPmain.handles.sampleTest.cell1 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.sampleTest.cellNameList,...
                        'Value',EPmain.sampleTest.cell1,'Position',[5 460 80 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.cell1,''Value'');','if tempVar ~=0,EPmain.sampleTest.cell1=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.cell1=tempVar;end;','ep(''start'');']);
                    
                    EPmain.handles.sampleTest.cell2 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.sampleTest.cellNameList,...
                        'Value',EPmain.sampleTest.cell2,'Position',[85 460 80 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.cell2,''Value'');','if tempVar ~=0,EPmain.sampleTest.cell2=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.cell2=tempVar;end;','ep(''start'');']);
                else
                    EPmain.handles.sampleTest.cell1 = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','Position',[5 460 80 20]);
                    
                    EPmain.handles.sampleTest.cell2 = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','Position',[85 460 80 20]);
                end
                
                if ~isempty(EPmain.sampleTest.subNameList)
                    EPmain.handles.sampleTest.sub1 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.sampleTest.subNameList,...
                        'Value',EPmain.sampleTest.sub1,'Position',[5 440 80 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.sub1,''Value'');','if tempVar ~=0,EPmain.sampleTest.sub1=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.sub1=tempVar;end;','ep(''start'');']);
                    
                    EPmain.handles.sampleTest.sub2 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.sampleTest.subNameList,...
                        'Value',EPmain.sampleTest.sub2,'Position',[85 440 80 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.sub2,''Value'');','if tempVar ~=0,EPmain.sampleTest.sub2=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.sub2=tempVar;end;','ep(''start'');']);
                else
                    EPmain.handles.sampleTest.sub1 = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','Position',[5 420 80 20]);
                    
                    EPmain.handles.sampleTest.sub2 = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','Position',[85 420 80 20]);
                end
                
                if ~isempty(EPmain.sampleTest.sessNameList)
                    EPmain.handles.sampleTest.sess1 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.sampleTest.sessNameList,...
                        'Value',EPmain.sampleTest.sess1,'Position',[5 420 80 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.sess1,''Value'');','if tempVar ~=0,EPmain.sampleTest.sess1=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.sess1=tempVar;end;','ep(''start'');']);
                    
                    EPmain.handles.sampleTest.sess2 = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.sampleTest.sessNameList,...
                        'Value',EPmain.sampleTest.sess2,'Position',[85 420 80 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.sess2,''Value'');','if tempVar ~=0,EPmain.sampleTest.sess2=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.sess2=tempVar;end;','ep(''start'');']);
                else
                    EPmain.handles.sampleTest.sess1 = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','Position',[5 420 80 20]);
                    
                    EPmain.handles.sampleTest.sess2 = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','Position',[85 420 80 20]);
                end
                
                uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Method','HorizontalAlignment','left',...
                    'Position',[10 400 150 20]);
                
                EPmain.handles.sampleTest.method = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{'sample','CWT','Template Woody','PCA Woody'},...
                    'Value',EPmain.sampleTest.method,'Position',[65 400 140 20],...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.method,''Value'');','if tempVar ~=0,EPmain.sampleTest.method=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.method=tempVar;end;','ep(''start'');']);
                
                uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Test','HorizontalAlignment','left',...
                    'Position',[10 380 150 20]);
                
                EPmain.handles.sampleTest.test = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{'t-test','jack-knife','t-map'},...
                    'Value',EPmain.sampleTest.test,'Position',[65 380 140 20],...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.test,''Value'');','if tempVar ~=0,EPmain.sampleTest.test=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.test=tempVar;end;','ep(''start'');']);
                
                switch EPmain.sampleTest.test
                    case {1,3} %t-test or t-map
                        
                        EPmain.handles.sampleTest.alphaLabel = uicontrol('Style','text','HorizontalAlignment','left','String', 'Alpha','FontSize',EPmain.fontsize,...
                            'Position',[10 340 70 20]);
                        
                        EPmain.handles.sampleTest.alpha = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.sampleTest.alpha),'FontSize',EPmain.fontsize,...
                            'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.sampleTest.alpha,''String''));','if tempVar ~=0,EPmain.sampleTest.alpha=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.alpha=tempVar;end;','ep(''start'');'],...
                            'Position',[10 320 70 20],'TooltipString','Threshold for sample by sample significance testing.');
                        
                        EPmain.handles.sampleTest.contiguousLabel=uicontrol('Style','text','HorizontalAlignment','left','String', 'Contiguous','FontSize',EPmain.fontsize,...
                            'Position',[80 340 70 20]);
                        
                        EPmain.handles.sampleTest.contiguous = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.sampleTest.contiguous,'FontSize',EPmain.fontsize,...
                            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.contiguous,''String'');','if tempVar ~=0,EPmain.sampleTest.contiguous=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.contiguous=tempVar;end;','ep(''start'');'],...
                            'Position',[80 320 70 20],'TooltipString','How many contiguous significant samples are required to consider to be significant.');
                        
                    case 2 %jack-knife
                        
                        EPmain.handles.sampleTest.threshLabel = uicontrol('Style','text','HorizontalAlignment','left','String', 'Threshold','FontSize',EPmain.fontsize,...
                            'Position',[10 340 70 20]);
                        
                        EPmain.handles.sampleTest.thresh = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.sampleTest.thresh),'FontSize',EPmain.fontsize,...
                            'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.sampleTest.thresh,''String''));','if tempVar ~=0,EPmain.sampleTest.thresh=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.thresh=tempVar;end;','ep(''start'');'],...
                            'Position',[10 320 70 20],'TooltipString','Threshold for jack-knife estimation of onset latency.');
                        
                        EPmain.handles.sampleTest.chanLabel=uicontrol('Style','text','HorizontalAlignment','left','String', 'Channel','FontSize',EPmain.fontsize,...
                            'Position',[80 340 70 20]);
                        
                        EPmain.handles.sampleTest.channel = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                            'String',EPdataset.dataset(EPmain.sampleTest.dataset).chanNames(find(strcmp(EPdataset.dataset(EPmain.sampleTest.dataset).chanTypes,'EEG'))),...
                            'Value',EPmain.sampleTest.channel,'Position',[80 320 100 20],...
                            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.channel,''Value'');','if tempVar ~=0,EPmain.sampleTest.channel=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.channel=tempVar;end;','ep(''start'');']);
                end
                
                EPmain.handles.sampleTest.filterPass= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{'Low Pass','High Pass','Band Pass','Band Stop','Notch'},...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.filterPass,''Value'');','if tempVar ~=0,EPmain.sampleTest.filterPass=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.filterPass=tempVar;end;','ep(''start'');'],...
                    'Value',EPmain.sampleTest.filterPass,'Position',[5 300 110 20]);
                
                EPmain.handles.sampleTest.filter1 = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.sampleTest.filter1),'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.sampleTest.filter1,''String''));','if tempVar ~=0,EPmain.sampleTest.filter1=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.filter1=tempVar;end;','ep(''start'');'],...
                    'Position',[110 300 40 20],'TooltipString','Lower frequency limit.');
                
                EPmain.handles.sampleTest.filter2 = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.sampleTest.filter2),'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.sampleTest.filter2,''String''));','if tempVar ~=0,EPmain.sampleTest.filter2=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.filter2=tempVar;end;','ep(''start'');'],...
                    'Position',[155 300 40 20],'TooltipString','Upper frequency limit.');
                
                EPmain.handles.sampleTest.filterType= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{'One-Pass Butterworth','Two-Pass Butterworth','One-Pass FIR','Two-Pass FIR','One-Pass FIRLS','Two-Pass FIRLS'},...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.filterType,''Value'');','if tempVar ~=0,EPmain.sampleTest.filterType=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.filterType=tempVar;end;','ep(''start'');'],...
                    'Value',EPmain.sampleTest.filterType,'Position',[5 280 160 20]);
                
                EPmain.handles.sampleTest.filterOrder = uicontrol('Style','edit','HorizontalAlignment','left','String', num2str(EPmain.sampleTest.filterOrder),'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','tempVar=str2num(get(EPmain.handles.sampleTest.filterOrder,''String''));','if tempVar ~=0,EPmain.sampleTest.filterOrder=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.filterOrder=tempVar;end;','ep(''start'');'],...
                    'Position',[160 280 30 20],'TooltipString','Order of the filter.');
                
                if ~isempty(EPdataset.dataset(EPmain.sampleTest.dataset).facVecT) %cannot filter temporal PCA data
                    set(EPmain.handles.sampleTest.filterPass,'enable','off');
                    set(EPmain.handles.sampleTest.filter1,'enable','off');
                    set(EPmain.handles.sampleTest.filter2,'enable','off');
                    set(EPmain.handles.sampleTest.filterType,'enable','off');
                    set(EPmain.handles.sampleTest.filterOrder,'enable','off');
                end
                
                EPmain.handles.sampleTest.scaleLabel = uicontrol('Style','text','HorizontalAlignment','left','String', 'Scale','FontSize',EPmain.fontsize,...
                    'Position',[10 260 50 20]);
                
                EPmain.handles.sampleTest.scale = uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.sampleTest.scale,'FontSize',EPmain.fontsize,...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.scale,''String'');','if tempVar ~=0,EPmain.sampleTest.scale=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.scale=tempVar;end;','ep(''start'');'],...
                    'Position',[70 260 50 20],'TooltipString','Width of the Mexican hat wavelet template.');
                
                if EPmain.sampleTest.method ~= 2
                    set(EPmain.handles.sampleTest.scaleLabel,'enable','off');
                    set(EPmain.handles.sampleTest.scale,'enable','off');
                end
            else
                EPmain.handles.sampleTest.contMethod = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{'Template Woody','PCA Woody'},...
                    'Value',EPmain.sampleTest.contMethod,'Position',[5 400 160 20],...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.contMethod,''Value'');','if tempVar ~=0,EPmain.sampleTest.contMethod=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.contMethod=tempVar;end;','ep(''start'');']);
            end
            
            %AVE template
            if (strcmp(EPdataset.dataset(EPmain.sampleTest.dataset).dataType,'continuous') && (EPmain.sampleTest.contMethod == 1)) || (~strcmp(EPdataset.dataset(EPmain.sampleTest.dataset).dataType,'continuous') && (EPmain.sampleTest.method == 3))
                EPmain.handles.sampleTest.AVElabel = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Average Template','HorizontalAlignment','left',...
                    'Position',[10 240 150 20]);
                
                if ~isempty(EPmain.sampleTest.AVElist)
                    EPmain.handles.sampleTest.AVE = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.sampleTest.AVEnameList,...
                        'Value',EPmain.sampleTest.AVE,'Position',[5 220 160 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.AVE,''Value'');','if tempVar ~=0,EPmain.sampleTest.AVE=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.AVE=tempVar;end;','EPmain.sampleTest.subject=1;','EPmain.sampleTest.cell=1;','ep(''start'');']);
                else
                    EPmain.handles.sampleTest.AVE = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','HorizontalAlignment','left',...
                        'Position',[5 220 160 20]);
                end
                
                EPmain.handles.sampleTest.subLabel = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Subject','HorizontalAlignment','left',...
                    'Position',[10 200 150 20]);
                
                EPmain.handles.sampleTest.cellLabel = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Cell','HorizontalAlignment','left',...
                    'Position',[10 180 150 20]);
                
                if ~isempty(EPmain.sampleTest.AVElist)
                    EPmain.handles.sampleTest.subject = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).subNames,...
                        'Value',EPmain.sampleTest.subject,'Position',[50 200 160 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.subject,''Value'');','if tempVar ~=0,EPmain.sampleTest.subject=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.subject=tempVar;end;','ep(''start'');']);
                    
                    EPmain.handles.sampleTest.cell = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).cellNames,...
                        'Value',EPmain.sampleTest.cell,'Position',[50 180 160 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.cell,''Value'');','if tempVar ~=0,EPmain.sampleTest.cell=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.cell=tempVar;end;','ep(''start'');']);
                else
                    noSampTest=1;
                    EPmain.handles.sampleTest.subject = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','HorizontalAlignment','left',...
                        'Position',[5 180 160 20]);
                    
                    EPmain.handles.sampleTest.cell = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','HorizontalAlignment','left',...
                        'Position',[5 140 160 20]);
                end
                %PCA Woody
            elseif (strcmp(EPdataset.dataset(EPmain.sampleTest.dataset).dataType,'continuous') && (EPmain.sampleTest.contMethod == 2)) || (~strcmp(EPdataset.dataset(EPmain.sampleTest.dataset).dataType,'continuous') && (EPmain.sampleTest.method == 4))
                EPmain.handles.sampleTest.PCAlabel = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','PCA','HorizontalAlignment','left',...
                    'Position',[10 240 150 20]);
                
                if ~isempty(EPmain.sampleTest.PCAlist)
                    EPmain.handles.sampleTest.PCA = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPmain.sampleTest.PCAnameList,...
                        'Value',EPmain.sampleTest.PCA,'Position',[5 220 160 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.PCA,''Value'');','if tempVar ~=0,EPmain.sampleTest.PCA=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.PCA=tempVar;end;','EPmain.sampleTest.factor=1;','ep(''start'');']);
                else
                    EPmain.handles.sampleTest.PCA = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','HorizontalAlignment','left',...
                        'Position',[5 220 160 20]);
                end
                
                EPmain.handles.sampleTest.factorLabel = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Factor','HorizontalAlignment','left',...
                    'Position',[10 200 150 20]);
                
                if ~isempty(EPmain.sampleTest.PCAlist)
                    EPmain.handles.sampleTest.factor = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facNames,...
                        'Value',EPmain.sampleTest.factor,'Position',[50 200 160 20],...
                        'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.factor,''Value'');','if tempVar ~=0,EPmain.sampleTest.factor=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.factor=tempVar;end;','ep(''start'');']);
                else
                    noSampTest=1;
                    EPmain.handles.sampleTest.factor = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none','HorizontalAlignment','left',...
                        'Position',[50 200 160 20]);
                end
            end
            
            if EPmain.sampleTest.freqFlag
                EPmain.handles.sampleTest.freqLabel = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Hz','HorizontalAlignment','left',...
                    'Position',[10 160 150 20]);
                EPmain.handles.sampleTest.freqScaleLabel = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Hz Scale','HorizontalAlignment','left',...
                    'Position',[10 140 150 20]);
                theFreqs=num2cell(EPdataset.dataset(EPmain.sampleTest.dataset).freqNames);
                theFreqs{end+1}='-all-';
                EPmain.handles.sampleTest.freqBin = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',theFreqs,...
                    'Value',EPmain.sampleTest.freqBin,'Position',[50 160 160 20],...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.freqBin,''Value'');','if tempVar ~=0,EPmain.sampleTest.freqBin=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.freqBin=tempVar;end;','ep(''start'');']);
                EPmain.handles.sampleTest.freqScale = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                    'String',{'cm-real','cm-imag','am','pw','dB'},...
                    'Value',EPmain.sampleTest.freqScale,'Position',[50 140 160 20],...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.sampleTest.freqScale,''Value'');','if tempVar ~=0,EPmain.sampleTest.freqScale=tempVar;end;','if isempty(tempVar),EPmain.sampleTest.freqScale=tempVar;end;','ep(''start'');']);
            end
            
            EEGchans=find(strcmp('EEG',EPdataset.dataset(EPmain.sampleTest.dataset).chanTypes));
            EEGchans=EEGchans(~cellfun(@isempty,{EPdataset.dataset(EPmain.sampleTest.dataset).eloc(EEGchans).radius}));
            
            if strcmp(EPmain.sampleTest.dataType,'continuous')
                switch EPmain.sampleTest.contMethod
                    case 1 %AVE template
                        if ~isempty(EPmain.sampleTest.AVElist)
                            templateData=EPmain.sampleTest.AVEdata{EPmain.sampleTest.AVE}(EEGchans,:,EPmain.sampleTest.cell,EPmain.sampleTest.subject,:,EPmain.sampleTest.freqBin);
                            if EPmain.sampleTest.freqFlag
                                if (EPmain.sampleTest.freqScale == 1)
                                    templateData=real(templateData);
                                end
                                if (EPmain.sampleTest.freqScale == 2)
                                    templateData=imag(templateData);
                                end
                                if (EPmain.sampleTest.freqScale > 2)
                                    templateData(EEGchans,:,:,:)=abs(templateData(EEGchans,:,:,:)); %convert complex number to real number
                                end
                                templateData(EEGchans,:,:,:)=templateData(EEGchans,:,:,:)/mean(diff(EPdataset.dataset(EPmain.sampleTest.dataset).freqNames)); %convert to spectral density
                                if EPmain.sampleTest.freqScale > 3
                                    templateData(EEGchans,:,:,:)=templateData(EEGchans,:,:,:).^2; %convert amplitude to power
                                end
                                if (EPmain.sampleTest.freqScale == 5)
                                    if ~all(templateData(EEGchans,:,:,:) >=0)
                                        disp('Some values negative so will take absolute value prior to log transform.  This can happen with PCA results due to minor imprecisions in PCA.');
                                    end
                                    templateData(EEGchans,:,:,:)=log10(abs(templateData(EEGchans,:,:,:)))*10; %convert to dB log scaling
                                    tempVar=templateData(EEGchans,:,:,:);
                                    tempVar(isinf(tempVar))=-flintmax;
                                    templateData(EEGchans,:,:,:)=tempVar; %log10 of zero is -inf.  Replace with maximum possible double-precision negative number.
                                end
                            end
                            if ~isempty(EPmain.sampleTest.sampEnd)
                                templateData2=templateData(:,EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd);
                            else
                                templateData2=templateData;
                            end
                            [A maxChan]=max(max(templateData2'));
                            [A maxPoint]=max(max(templateData2));
                            if ~isempty(EPmain.sampleTest.sampStart)
                                maxPoint=maxPoint+EPmain.sampleTest.sampStart-1;
                            end
                            EPmain.sampleTest.templateWaveform=templateData(maxChan,:)';
                            EPmain.sampleTest.templateTopo=templateData(:,maxPoint);
                            tMs=EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).timeNames;
                        else
                            templateData=zeros(length(EPdataset.dataset(EPmain.sampleTest.dataset).timeNames),length(EEGchans));
                            EPmain.sampleTest.templateWaveform=zeros(length(EPdataset.dataset(EPmain.sampleTest.dataset).timeNames),1);
                            EPmain.sampleTest.templateTopo=zeros(1,length(EEGchans));
                            tMs=EPdataset.dataset(EPmain.sampleTest.dataset).timeNames;
                        end
                    case 2 %PCA
                        if ~isempty(EPmain.sampleTest.PCAlist) && ~isempty(EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecT)
                            EPmain.sampleTest.templateWaveform=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecT(:,EPmain.sampleTest.factor);
                            tMs=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).timeNames;
                        else
                            EPmain.sampleTest.templateWaveform=zeros(length(EPdataset.dataset(EPmain.sampleTest.dataset).timeNames),1);
                            tMs=EPdataset.dataset(EPmain.sampleTest.dataset).timeNames;
                        end
                        if ~isempty(EPmain.sampleTest.PCAlist) && ~isempty(EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecS)
                            templateTopo=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecS(EEGchans,EPmain.sampleTest.factor);
                            if EPmain.sampleTest.test==2 %not just one channel
                                set(EPmain.handles.sampleTest.chanLabel,'enable','off');
                                set(EPmain.handles.sampleTest.channel,'enable','off');
                            end
                        else
                            templateTopo=zeros(length(EEGchans));
                        end
                    otherwise
                        disp('Oops: programmer error');
                        return;
                end
                if ~isempty(EPmain.sampleTest.AVElist)
                    EPmain.sampleTest.Fs=EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).Fs;
                    EPmain.sampleTest.baseline=EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).baseline;
                end
            else
                switch EPmain.sampleTest.method
                    case 1 %sample
                        EPmain.sampleTest.templateWaveform=zeros(length(EPdataset.dataset(EPmain.sampleTest.dataset).timeNames),1);
                        EPmain.sampleTest.templateTopo=zeros(length(EEGchans));
                        tMs=EPdataset.dataset(EPmain.sampleTest.dataset).timeNames;
                    case 2 %CWT
                        tMs=EPdataset.dataset(EPmain.sampleTest.dataset).timeNames;
                        theTau=tMs(floor(length(tMs)/2));
                        EPmain.sampleTest.templateWaveform=(1-16.*(((tMs-theTau)/(EPmain.sampleTest.scale/1000))/1000).^2).*exp(-8.*(((tMs-theTau)/(EPmain.sampleTest.scale/1000))/1000).^2); %The Mexican Hat wavelet
                        EPmain.sampleTest.templateTopo=zeros(length(EEGchans));
                    case 3 %AVE template
                        if ~isempty(EPmain.sampleTest.AVElist)
                            templateData=EPmain.sampleTest.AVEdata{EPmain.sampleTest.AVE}(EEGchans,EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd,EPmain.sampleTest.cell,EPmain.sampleTest.subject,:,EPmain.sampleTest.freqBin);
                            if EPmain.sampleTest.freqFlag
                                if (EPmain.sampleTest.freqScale == 1)
                                    templateData=real(templateData);
                                end
                                if (EPmain.sampleTest.freqScale == 2)
                                    templateData=imag(templateData);
                                end
                                if (EPmain.sampleTest.freqScale > 2)
                                    templateData(EEGchans,:,:,:)=abs(templateData(EEGchans,:,:,:)); %convert complex number to real number
                                end
                                templateData(EEGchans,:,:,:)=templateData(EEGchans,:,:,:)/mean(diff(EPdataset.dataset(EPmain.sampleTest.dataset).freqNames)); %convert to spectral density
                                if EPmain.sampleTest.freqScale > 3
                                    templateData(EEGchans,:,:,:)=templateData(EEGchans,:,:,:).^2; %convert amplitude to power
                                end
                                if (EPmain.sampleTest.freqScale == 5)
                                    if ~all(templateData(EEGchans,:,:,:) >=0)
                                        disp('Some values negative so will take absolute value prior to log transform.  This can happen with PCA results due to minor imprecisions in PCA.');
                                    end
                                    templateData(EEGchans,:,:,:)=log10(abs(templateData(EEGchans,:,:,:)))*10; %convert to dB log scaling
                                    tempVar=templateData(EEGchans,:,:,:);
                                    tempVar(isinf(tempVar))=-flintmax;
                                    templateData(EEGchans,:,:,:)=tempVar; %log10 of zero is -inf.  Replace with maximum possible double-precision negative number.
                                end
                            end
                            if ~isempty(EPmain.sampleTest.sampEnd)
                                templateData2=templateData(:,EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd);
                            else
                                templateData2=templateData;
                            end
                            [A maxChan]=max(max(templateData2'));
                            [A maxPoint]=max(max(templateData2));
                            if ~isempty(EPmain.sampleTest.sampStart)
                                maxPoint=maxPoint+EPmain.sampleTest.sampStart-1;
                            end
                            EPmain.sampleTest.templateWaveform=templateData(maxChan,:)';
                            EPmain.sampleTest.templateTopo=templateData(:,maxPoint);
                            tMs=EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).timeNames;
                        else
                            noSampTest=1;
                            templateData=zeros(length(EPdataset.dataset(EPmain.sampleTest.dataset).timeNames),length(EEGchans));
                            EPmain.sampleTest.templateWaveform=zeros(length(EPdataset.dataset(EPmain.sampleTest.dataset).timeNames),1);
                            EPmain.sampleTest.templateTopo=zeros(1,length(EEGchans));
                            tMs=EPdataset.dataset(EPmain.sampleTest.dataset).timeNames;
                        end
                    case 4 %PCA
                        if ~isempty(EPmain.sampleTest.PCAlist) && ~isempty(EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecT)
                            EPmain.sampleTest.templateWaveform=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecT(:,EPmain.sampleTest.factor);
                            tMs=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).timeNames;
                        else
                            EPmain.sampleTest.templateWaveform=zeros(length(EPdataset.dataset(EPmain.sampleTest.dataset).timeNames),1);
                            tMs=EPdataset.dataset(EPmain.sampleTest.dataset).timeNames;
                        end
                        if ~isempty(EPmain.sampleTest.PCAlist) && ~isempty(EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecS)
                            EPmain.sampleTest.templateTopo=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecS(EEGchans,EPmain.sampleTest.factor);
                            if EPmain.sampleTest.test==2 %not just one channel
                                set(EPmain.handles.sampleTest.chanLabel,'enable','off');
                                set(EPmain.handles.sampleTest.channel,'enable','off');
                            end
                        else
                            noSampTest=1;
                            EPmain.sampleTest.templateTopo=zeros(length(EEGchans));
                        end
                    otherwise
                        disp('Oops: programmer error');
                        return;
                end
            end
            if EPmain.sampleTest.changeFlag==1
                EPmain.sampleTest.changeFlag=0;
%                 if strcmp(EPmain.sampleTest.dataType,'continuous')
                    EPmain.sampleTest.sampStart=1;
                    EPmain.sampleTest.sampEnd=length(EPmain.sampleTest.templateWaveform);
%                 end
            end
            
            %Waveform plot of template
            EPmain.handles.sampleTest.templatePlot = axes('units','pixels','position',[5 88 100 50]);
            EPmain.handles.sampleTest.templateWaves = plot(tMs(EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd),EPmain.sampleTest.templateWaveform(EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd));
            axis([tMs(EPmain.sampleTest.sampStart) tMs(EPmain.sampleTest.sampEnd) min([-1; EPmain.sampleTest.templateWaveform]) max([1; EPmain.sampleTest.templateWaveform])]);
            
            %2D plot of template
            warning off
            [Xi,Yi,Zi] = griddata(EPmain.sampleTest.x,EPmain.sampleTest.y,EPmain.sampleTest.templateTopo,[1:EPmain.sampleTest.gridSize]',[1:EPmain.sampleTest.gridSize],'linear');
            warning on
            EPmain.handles.sampleTest.templateTopo = axes('units','pixels','position',[110 88 50 50]);
            EPmain.handles.sampleTest.templateTopoImage = imagesc(Zi);
            set(gca,'XTickLabel','','YTickLabel','');
            
            if strcmp(EPmain.sampleTest.dataType,'continuous')
                uicontrol('Style','text',...
                    'String','samples','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[5 55 50 20]);
                
                uicontrol('Style','text',...
                    'String','ms','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'Position',[75 55 50 20]);
                
                if ~isempty(EPmain.sampleTest.AVElist)
                    EPmain.handles.sampleTest.sampStart= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                        'String',EPmain.sampleTest.sampStart-EPmain.sampleTest.baseline-1,...
                        'TooltipString','First sample of template in relation to the event, where negative is before it.',...
                        'Position',[5 35 35 20],'Callback',@sampleTestSampStart);
                    
                    EPmain.handles.sampleTest.sampEnd= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                        'String',EPmain.sampleTest.sampEnd-EPmain.sampleTest.baseline,...
                        'TooltipString','Last sample of template in relation to the event, where negative is before it.',...
                        'Position',[40 35 35 20],'Callback',@sampleTestSampEnd);
                    
                    EPmain.handles.sampleTest.msStart= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                        'String',(EPmain.sampleTest.sampStart-EPmain.sampleTest.baseline-1)*(1000/EPmain.sampleTest.Fs),...
                        'TooltipString','Last ms of template in relation to the event, where negative is before it.',...
                        'Position',[75 35 35 20],'Callback',@sampleTestSampStart);
                    
                    EPmain.handles.sampleTest.msEnd= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                        'String',(EPmain.sampleTest.sampEnd-EPmain.sampleTest.baseline)*(1000/EPmain.sampleTest.Fs),...
                        'TooltipString','First ms of template in relation to the event, where negative is before it.',...
                        'Position',[110 35 35 20],'Callback',@sampleTestSampEnd);
                else
                    EPmain.handles.sampleTest.sampStart= uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none',...
                        'TooltipString','First sample of template in relation to the event, where negative is before it.',...
                        'Position',[5 35 35 20]);
                    
                    EPmain.handles.sampleTest.sampEnd= uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none',...
                        'TooltipString','Last sample of template in relation to the event, where negative is before it.',...
                        'Position',[40 35 35 20]);
                    
                    EPmain.handles.sampleTest.msStart= uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none',...
                        'TooltipString','Last ms of template in relation to the event, where negative is before it.',...
                        'Position',[75 35 35 20]);
                    
                    EPmain.handles.sampleTest.msEnd= uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','none',...
                        'TooltipString','First ms of template in relation to the event, where negative is before it.',...
                        'Position',[110 35 35 20]);
                end
            end
            
            EPmain.handles.sampleTest.sampleTest = uicontrol('Style', 'pushbutton', 'String', 'Run','FontSize',EPmain.fontsize,...
                'Position', [70 0 60 35], 'Callback', 'ep(''sampleTest'');');
            
            if noSampTest
                set(EPmain.handles.sampleTest.sampleTest,'enable','off');
            end
            
        else
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','No suitable data',...
                'Position',[5 480 160 20]);
            
        end
        
        ep_tictoc('end');
        
        EPmain.handles.sampleTest.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [130 0 60 35], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        
    case 'sampleTest'
        
        EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);
        ep_tictoc('begin');
        
        theDataset=EPmain.sampleTest.dataset;
        EPdata=ep_loadEPdataset(theDataset);
        
        if ~strcmp(EPmain.sampleTest.dataType,'continuous')
            templateTopo=[];
            templateWaveform=[];
            heightThresh=[];
            durThresh=[];
            waveletWidth=[];
            freq1=[];
            freq2=[];
            switch EPmain.sampleTest.method
                case 1
                    method='sample';
                    methodName='samp';
                    heightThresh=EPmain.sampleTest.alpha;
                    durThresh=EPmain.sampleTest.contiguous;
                case 2
                    method='CWT';
                    methodName='CWT';
                    heightThresh=EPmain.sampleTest.alpha;
                    durThresh=EPmain.sampleTest.contiguous;
                    waveletWidth=EPmain.sampleTest.scale;
                case 3
                    method='Template Woody';
                    methodName='TW';
                    templateTopo=EPmain.sampleTest.templateTopo;
                    templateWaveform=EPmain.sampleTest.templateWaveform;
                    heightThresh=EPmain.sampleTest.alpha;
                    durThresh=EPmain.sampleTest.contiguous;
                case 4
                    method='PCA Woody';
                    methodName='PW';
                    facVecT=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecT;
                    if ~isempty(facVecT)
                        templateWaveform=facVecT(EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd,EPmain.sampleTest.factor);
                    else
                        templateWaveform=[];
                    end
                    facVecS=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecS;
                    if ~isempty(facVecS)
                        templateTopo=facVecS(:,EPmain.sampleTest.factor);
                    else
                        templateTopo=[];
                    end
                    heightThresh=EPmain.sampleTest.alpha;
                    durThresh=EPmain.sampleTest.contiguous;
                otherwise
                    disp('Programmer error');
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return;
            end
            
            switch EPmain.sampleTest.test
                case 1
                    test='t-test';
                case 2
                    test='jack-knife';
                case 3
                    test='t-map';
                otherwise
                    disp('Programmer error');
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return;
            end
            
            cellList1=find(ismember(EPdataset.dataset(theDataset).cellNames,EPmain.sampleTest.cellNameList(EPmain.sampleTest.cell1)));
            cellList2=find(ismember(EPdataset.dataset(theDataset).cellNames,EPmain.sampleTest.cellNameList(EPmain.sampleTest.cell2)));
            bothCellList=union(cellList1,cellList2);
            cellList1=find(ismember(bothCellList,cellList1));
            cellList2=find(ismember(bothCellList,cellList2));            
            
            subList1=EPmain.sampleTest.subList{EPmain.sampleTest.sub1};
            subList2=EPmain.sampleTest.subList{EPmain.sampleTest.sub2};
            if ~isempty(EPdata.sessNums)
                subList1=intersect(find(EPdata.sessNums==(EPmain.sampleTest.sess1)),subList1);
                subList2=intersect(find(EPdata.sessNums==(EPmain.sampleTest.sess2)),subList2);
            end
            bothSubList=union(subList1,subList2);
            subList1=find(ismember(bothSubList,subList1));
            subList2=find(ismember(bothSubList,subList2));            
            
            if ~isempty(EPdataset.dataset(theDataset).freqNames)
                if EPmain.sampleTest.freqBin > length(EPdataset.dataset(theDataset).freqNames)
                    freqList=[1:length(EPdataset.dataset(theDataset).freqNames)];
                else
                    freqList=EPmain.sampleTest.freqBin;
                end
            else
                freqList=[];
            end

            EPdataST=ep_selectData(EPdata,{[],[],bothCellList,[],[],freqList}); %no need to select subjects as sampTest handles that.
            if isempty(EPdataST)
                msg{1}='No regular data to analyze.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            EEGchans=find(strcmp('EEG',EPdataST.chanTypes));
            if isempty(EEGchans)
                msg{1}='No EEG channels to analyze.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
                        
            if isempty(EPdataset.dataset(EPmain.sampleTest.dataset).facVecT) && ~isempty(EPmain.sampleTest.filter1) %cannot filter temporal PCA data
                if (EPmain.sampleTest.filterPass==3) && (EPmain.sampleTest.filter1 >= EPmain.sampleTest.filter2)
                    msg{1}='The start of the bandpass window must come before the end of the bandpass window.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
                
                if (EPmain.sampleTest.filterPass==4) && (EPmain.sampleTest.filter1 >= EPmain.sampleTest.filter2)
                    msg{1}='The start of the bandstop window must come before the end of the bandstop window.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
                
                if (length(EPmain.sampleTest.filter1) > 1) || (length(EPmain.sampleTest.filter2) > 1) || (length(EPmain.sampleTest.filterOrder) > 1)
                    msg{1}='The filter setting in each field needs to be a single number.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
                
                if (EPmain.sampleTest.filter1 < 0) || (EPmain.sampleTest.filterOrder < 0)
                    msg{1}='Filter settings cannot be negative numbers.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
                
                if strcmp(test,'jack-knife') && ~isempty(freqList)
                    msg{1}='Jack-knife can only be performed on a single frequency band.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
                
                cfg=[];
                switch EPmain.sampleTest.filterType
                    case 1
                        filterDirection='onepass';
                        filterType='but';
                    case 2
                        filterDirection='twopass';
                        filterType='but';
                    case 3
                        filterDirection='onepass';
                        filterType='fir';
                    case 4
                        filterDirection='twopass';
                        filterType='fir';
                    case 5
                        filterDirection='onepass';
                        filterType='firls';
                    case 6
                        filterDirection='twopass';
                        filterType='firls';
                end
                switch EPmain.sampleTest.filterPass
                    case 1 %low pass
                        cfg.lpfilter='yes';
                        cfg.lpfreq=EPmain.sampleTest.filter1;
                        cfg.lpfiltdir=filterDirection;
                        cfg.lpfilttype=filterType;
                        cfg.lpfiltord=EPmain.sampleTest.filterOrder;
                    case 2 %high pass
                        cfg.hpfilter='yes';
                        cfg.hpfreq=EPmain.sampleTest.filter1;
                        cfg.hpfiltdir=filterDirection;
                        cfg.hpfilttype=filterType;
                        cfg.hpfiltord=EPmain.sampleTest.filterOrder;
                    case 3 %band pass
                        cfg.bpfilter='yes';
                        cfg.bpfreq=[EPmain.sampleTest.filter1 EPmain.sampleTest.filter2];
                        cfg.bpfiltdir=filterDirection;
                        cfg.bpfilttype=filterType;
                        cfg.bpfiltord=EPmain.sampleTest.filterOrder;
                    case 4 %band stop
                        cfg.bsfilter='yes';
                        cfg.bsfreq=[EPmain.sampleTest.filter1 EPmain.sampleTest.filter2];
                        cfg.bsfiltdir=filterDirection;
                        cfg.bsfilttype=filterType;
                        cfg.bsfiltord=EPmain.sampleTest.filterOrder;
                    case 5 %notch
                        cfg.dftfilter='yes';
                        cfg.dftfreq=EPmain.sampleTest.filter1;
                        cfg.dftfiltdir=filterDirection;
                        cfg.dftfilttype=filterType;
                        cfg.dftfiltord=EPmain.sampleTest.filterOrder;
                end
                EPdataST=ep_filterData(EPdataST,cfg,EEGchans);
            end
            
            if EPmain.sampleTest.cell1~=EPmain.sampleTest.cell2
                theMode='cells';
            elseif EPmain.sampleTest.sub1~=EPmain.sampleTest.sub2
                theMode='subjects';
            elseif EPmain.sampleTest.sess1~=EPmain.sampleTest.sess2
                theMode='sessions';
            else
                msg{1}='Two different aspects of the data must be chosen to be contrasted.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            
            [outputData, sampLat1, sampLat2, sampAmp1, sampAmp2]=ep_sampleTest(EPdataST,theMode,cellList1,cellList2,subList1,subList2,method,heightThresh,durThresh,waveletWidth,templateWaveform,templateTopo,test,EPmain.sampleTest.thresh,EPmain.sampleTest.channel);
            if EPtictoc.stop
                EPtictoc.stop=0;
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                ep('start');
            end
            if isempty(outputData)
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            
            %add the results to the dataset if t-test
            if any(strcmp(test,{'t-test','t-map'}))
                if strcmp(test,'t-map')
                    methodName=[method '-map'];
                end
                if ~isempty(EPmain.sampleTest.sessNameList)
                    sessName1=['_' EPmain.sampleTest.sessNameList{EPmain.sampleTest.sess1}];
                    sessName2=['_' EPmain.sampleTest.sessNameList{EPmain.sampleTest.sess2}];
                else
                    sessName1='';
                    sessName2='';
                end
                if EPmain.sampleTest.sub1==1
                    subName1='';
                else
                    subName1=['_' EPmain.sampleTest.subNameList{EPmain.sampleTest.sub1}];
                end
                if EPmain.sampleTest.sub2==1
                    subName2='';
                else
                    subName2=['_' EPmain.sampleTest.subNameList{EPmain.sampleTest.sub2}];
                end
                
                sampleTestCellName=[methodName ':' EPmain.sampleTest.cellNameList{EPmain.sampleTest.cell1} subName1 sessName1 '-' EPmain.sampleTest.cellNameList{EPmain.sampleTest.cell2} subName2 sessName2];
                
                if ~any(strcmp(sampleTestCellName,EPdata.cellNames))
                    EPadd=[];
                    EPadd.cellTypes{1}='STS';
                    if strcmp(EPdata.dataType,'single_trial')
                        EPadd.trialNames=1;
                    else
                        EPadd.trialNames=[];
                    end
                    EPadd.cellNames{1}=sampleTestCellName;
                    [EPdata]=ep_addData(EPdata,EPadd,'cells');
                    if isempty(EPdata)
                        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                        return
                    end
                else
                    disp('Overwriting existing sample test results');
                end
                if isempty(freqList)
                    theFreq=1;
                else
                    theFreq=freqList;
                end
                if strcmp(EPdata.dataType,'average')
                    if ~any(strcmp('sampleTest',EPdata.subNames))
                        EPadd=[];
                        EPadd.subTypes{1}='GAV';
                        EPadd.subNames{1}='sampleTest';
                        [EPdata]=ep_addData(EPdata,EPadd,'subjects');
                        if isempty(EPdata)
                            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                            return
                        end
                    end
                    EPdata.data(:,:,strcmp(sampleTestCellName,EPdata.cellNames),strcmp('sampleTest',EPdata.subNames),:,theFreq,:)=outputData;
                else %single_trial
                    EPdata.data(:,:,strcmp(sampleTestCellName,EPdata.cellNames),:,:,theFreq,:)=outputData;
                end
                
                %add latency information if Woody PCA and if filter was spatial PCA and thus a single latency for all the channels
                if strcmp('PCA Woody',method) && ~isempty(templateTopo) && strcmp('single_trial',EPdata.dataType) && (isscalar(freqList))
                    cellNames=unique(EPdataST.cellNames);
                    goodCellTrials1=find(EPdata.analysis.badTrials(1,strcmp(cellNames{1},EPdata.cellNames))==0);
                    goodCellTrials2=find(EPdata.analysis.badTrials(1,strcmp(cellNames{2},EPdata.cellNames))==0);
                    ce11Trials1=find(strcmp(cellNames{1},EPdata.cellNames));
                    ce11Trials2=find(strcmp(cellNames{2},EPdata.cellNames));
                    goodTrials1=ce11Trials1(goodCellTrials1);
                    goodTrials2=ce11Trials2(goodCellTrials2);
                    
                    trialSpecName=['Woody PCA latency: ' EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facNames{EPmain.sampleTest.factor}];
                    
                    if ~any(strcmp(trialSpecName,EPdata.trialSpecNames))
                        EPdata.trialSpecNames{end+1}=trialSpecName;
                        EPdata.trialSpecs(:,end+1)=cell(length(EPdata.cellNames),1);
                    end
                    
                    EPdata.trialSpecs(goodTrials1,find(strcmp(trialSpecName,EPdata.trialSpecNames)))=num2cell(sampLat1);
                    EPdata.trialSpecs(goodTrials2,find(strcmp(trialSpecName,EPdata.trialSpecNames)))=num2cell(sampLat2);
                    
                    trialSpecName=['Woody PCA amplitude: ' EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facNames{EPmain.sampleTest.factor}];
                    if ~any(strcmp(trialSpecName,EPdata.trialSpecNames))
                        EPdata.trialSpecNames{end+1}=trialSpecName;
                        EPdata.trialSpecs(:,end+1)=cell(length(EPdata.cellNames),1);
                    end
                    EPdata.trialSpecs(goodTrials1,find(strcmp(trialSpecName,EPdata.trialSpecNames)))=num2cell(sampAmp1);
                    EPdata.trialSpecs(goodTrials2,find(strcmp(trialSpecName,EPdata.trialSpecNames)))=num2cell(sampAmp2);
                end
            end
        else %continuous
            
            switch EPmain.sampleTest.method
                case 1
                    method='Template Woody';
                    methodName='TW';
                    templateTopo=EPmain.sampleTest.templateTopo;
                    templateWaveform=EPmain.sampleTest.templateWaveform(EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd);
                    sampleTestChanName=[methodName ':' EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).subNames{EPmain.sampleTest.subject} '-' EPdataset.dataset(EPmain.sampleTest.AVElist(EPmain.sampleTest.AVE)).cellNames{EPmain.sampleTest.cell},[],[]];
                case 2
                    method='PCA Woody';
                    methodName='PW';
                    sampleTestChanName=[methodName ':' EPdata.facNames{EPmain.sampleTest.factor}];
                    facVecT=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecT;
                    if ~isempty(facVecT)
                        templateWaveform=facVecT(EPmain.sampleTest.sampStart:EPmain.sampleTest.sampEnd,EPmain.sampleTest.factor);
                    else
                        templateWaveform=[];
                    end
                    facVecS=EPdataset.dataset(EPmain.sampleTest.PCAlist(EPmain.sampleTest.PCA)).facVecS;
                    if ~isempty(facVecS)
                        templateTopo=facVecS(:,EPmain.sampleTest.factor);
                    else
                        templateTopo=[];
                    end
                otherwise
                    disp('Programmer error');
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return;
            end
            
            test='';
            
            [outputData, sampLat1, sampLat2, sampAmp1, sampAmp2]=ep_sampleTest(EPdata,'',1,1,1,1,method,EPmain.sampleTest.alpha,EPmain.sampleTest.contiguous,EPmain.sampleTest.scale,templateWaveform,templateTopo,test,EPmain.sampleTest.thresh,EPmain.sampleTest.channel);

            if EPtictoc.stop
                EPtictoc.stop=0;
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                ep('start');
            end
            if isempty(outputData)
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            
            if ~any(strcmp(sampleTestChanName,EPdata.chanNames))
                EPadd=[];
                EPadd.chanTypes{1}='REG';
                EPadd.chanNames{1}=sampleTestChanName;
                [EPdata]=ep_addData(EPdata,EPadd,'channels');
                if isempty(EPdata)
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
            else
                disp('Overwriting existing Woody waveform');
            end
            
            EPdata.data(strcmp(sampleTestChanName,EPdata.chanNames),:,:,:,:,:,:)=outputData;
        end
        
        [err]=ep_checkEPfile(EPdata);
        if ~err
            theDescription=['Added ' method ' sample test.'];
            EPdata.history=ep_addHistory(EPdata.history, EPmain.preferences.records, theDescription, EPmain.handles.hMainWindow, EPmain.preferences, EPmain.preferences.EPver);
            
            %update the copy of the dataset in the working set
            ep_tictoc('ioStart');
            disp('Saving results to the working set.');
            whichData=EPmain.sampleTest.dataset;
            delete([EPdataset.EPwork filesep 'EPwork' filesep EPdataset.dataset(whichData).dataName '.mat']);
            EPdataset.dataset=EPdataset.dataset(setdiff([1:length(EPdataset.dataset)],whichData));
            ep_saveEPdataset(EPdata,length(EPdataset.dataset)+1,'no');
            EPdataset.dataset=[EPdataset.dataset(1:whichData-1) EPdataset.dataset(end) EPdataset.dataset(whichData:end-1)];
            if isfield(EPdataset.dataset(EPmain.sampleTest.dataset),'lastChange') && ~isempty(EPdataset.dataset(EPmain.sampleTest.dataset).lastChange)
                EPmain.sampleTest.lastChange=EPdataset.dataset(EPmain.sampleTest.dataset).lastChange;
            else
                EPmain.sampleTest.lastChange=cell(0);
            end
            ep_tictoc('ioFinish');
        end
        ep_tictoc('end');
        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        
    case 'startPCA'
        
        set(EPmain.handles.hMainWindow,'Name', 'PCA');

        figure(EPmain.handles.hMainWindow) %as of Matlab2024a, sometimes focus going to scree window and then putting the controls there instead.
        
        EPmain.handles.pca.prefs = uicontrol('Style', 'pushbutton', 'String', '•','FontSize',EPmain.fontsize,...
            'Position', [5 485 10 10], 'Callback', ['global EPmain;','EPmain.prefReturn=''PCA'';','EPmain.mode=''preferencePCA'';','ep(''start'');']);
        
        h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Mode','HorizontalAlignment','left',...
            'Position',[25 470 50 20]);
        
        EPmain.handles.pca.mode = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'spatial','temporal','frequency','cross-verification'},...
            'CallBack',['global EPmain;','tempVar=get(EPmain.handles.pca.mode,''Value'');','if tempVar ~=0,EPmain.pca.mode=tempVar;end;','if isempty(tempVar),EPmain.pca.mode=tempVar;end;EPmain.pca.crossVerifyPCA=[];','ep(''start'');'],...
            'Value',EPmain.pca.mode,'Position',[20 450 160 20]);
        
        if EPmain.pca.mode == 4
            %cross-verification mode
            
            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','PCA Datasets','HorizontalAlignment','left',...
                'Position',[25 430 150 20]);
            
            crossVerifyTableData=cell(0);
            EPmain.pca.PCAdatasets=[];
            if ~isempty(EPdataset.dataset)
                for i=1:length(EPdataset.dataset)
                    if ~isempty(EPdataset.dataset(i).facNames)
                        fileName=EPdataset.dataset(i).dataName;
                        if strcmp(EPdataset.dataset(i).saved,'no')
                            fileName=['*' fileName];
                        end
                        EPmain.pca.PCAdatasets(end+1)=i;
                        crossVerifyTableData{length(EPmain.pca.PCAdatasets),1}=fileName;
                    end
                end
            else
                crossVerifyTableData=cell(0);
            end
            
            tableNames{1}='data';
            columnEditable =  false;
            ColumnFormat{1}=[];
            
            EPmain.handles.crossVerify.hTable = uitable('Data',crossVerifyTableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
                'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
                'CellSelectionCallback',@crossVerifyData,...
                'ColumnWidth',{300},'Position',[20 280 170 150]);
            
            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Chosen:','HorizontalAlignment','left',...
                'Position',[25 250 45 20]);            

            if isempty(EPmain.pca.crossVerifyPCA)
                theText='none';
            else
                theText=EPdataset.dataset(EPmain.pca.crossVerifyPCA).dataName;
            end
            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String',theText,'HorizontalAlignment','left',...
                'Position',[75 250 150 20]);            

            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Title of PCA','HorizontalAlignment','left',...
                'Position',[25 210 150 20]);
            
            EPmain.handles.pca.name= uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'String',EPmain.pca.name,...
                'Callback', ['global EPmain;','EPmain.pca.name=get(EPmain.handles.pca.name,''String'');'],...
                'Position',[25 190 150 20]);
            
            if ~isempty(EPdataset.dataset)
                for i=1:length(EPdataset.dataset)
                    fileName=EPdataset.dataset(i).dataName;
                    if strcmp(EPdataset.dataset(i).saved,'no')
                        fileName=['*' fileName];
                    end
                    tableData{i,1}=fileName;
                end
            else
                tableData=[];
            end
            
            tableData=cell(0);
            EPmain.pca.targetDatasets=[];
            if ~isempty(EPdataset.dataset)
                for i=1:length(EPdataset.dataset)
                    dropFlag=0;
                    if ~isempty(EPmain.pca.crossVerifyPCA)
                        if ~isempty(EPdataset.dataset(EPmain.pca.crossVerifyPCA).facVecT)
                            if ~isempty(EPdataset.dataset(i).facVecT)
                                dropFlag=1;
                            end
                            if size(EPdataset.dataset(EPmain.pca.crossVerifyPCA).facVecT,1) ~= length(EPdataset.dataset(i).timeNames)
                                dropFlag=1;
                            end
                        end
                        if ~isempty(EPdataset.dataset(EPmain.pca.crossVerifyPCA).facVecS)
                            if ~isempty(EPdataset.dataset(i).facVecS)
                                dropFlag=1;
                            end
                            if size(EPdataset.dataset(EPmain.pca.crossVerifyPCA).facVecS,1) ~= length(EPdataset.dataset(i).chanNames(strcmp('EEG',EPdataset.dataset(i).chanTypes)))
                                dropFlag=1;
                            end
                        end
                        if ~isempty(EPdataset.dataset(EPmain.pca.crossVerifyPCA).facVecF)
                            if ~isempty(EPdataset.dataset(i).facVecF)
                                dropFlag=1;
                            end
                            if size(EPdataset.dataset(EPmain.pca.crossVerifyPCA).facVecF,1) ~= length(EPdataset.dataset(i).freqNames)
                                dropFlag=1;
                            end
                        end
                    end
                    if ~dropFlag
                        fileName=EPdataset.dataset(i).dataName;
                        if strcmp(EPdataset.dataset(i).saved,'no')
                            fileName=['*' fileName];
                        end
                        EPmain.pca.targetDatasets(end+1)=i;
                        tableData{length(EPmain.pca.targetDatasets),1}=fileName;
                    end
                end
            else
                tableData=cell(0);
            end

            tableNames{1}='data';
            columnEditable =  false;
            ColumnFormat{1}=[];
            
            EPmain.handles.pca.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
                'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
                'CellSelectionCallback',@pickPCAdata,...
                'ColumnWidth',{300},'Position',[20 40 170 150]);
        else
            
            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Rotation','HorizontalAlignment','left',...
                'Position',[25 430 50 20]);
            
            rotationList=ep_doPCA;
            EPmain.handles.pca.rotation= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',rotationList,...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.pca.rotation,''Value'');','if tempVar ~=0,EPmain.pca.rotation=tempVar;end;','if isempty(tempVar),EPmain.pca.rotation=tempVar;end;','EPmain.pca.rotFlag=1;','ep(''start'');'],...
                'Value',EPmain.pca.rotation,'Position',[20 410 160 20]);
            
            rotationNum = get(EPmain.handles.pca.rotation,'value');
            rotationName = get(EPmain.handles.pca.rotation,'string');
            PCArotation=rotationName{rotationNum};
            
            if strcmp(PCArotation,'SOBI') && (EPmain.pca.mode~=1)
                EPmain.pca.mode=1;
                set(EPmain.handles.pca.mode, 'Value', EPmain.pca.mode)
                disp('SOBI can only be used in spatial mode.');
            end

            rotoptLabel = uicontrol('Style','text',...
                'String','Parameter','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[25 390 150 20]);
            
            EPmain.handles.pca.rotopt= uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'String',EPmain.pca.rotopt,...
                'CallBack',['global EPmain;','EPmain.pca.rotopt=get(EPmain.handles.pca.rotopt,''String'');','ep(''start'');'],...
                'TooltipString','Parameter controlling the rotation method.',...
                'Position',[25 370 50 20]);
            
            if strcmp(PCArotation,{'Promax'})
            h = uicontrol('Style','text',...
                'String','Algorithm','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[110 390 150 20]);

                theList={'SAS';'SPSS'};
                EPmain.handles.pca.theAlgorithm= uicontrol('Style','popupmenu','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                    'String',theList,...
                    'CallBack',['global EPmain;','tempVar=get(EPmain.handles.pca.theAlgorithm,''String'');','EPmain.pca.theAlgorithm=tempVar{get(EPmain.handles.pca.theAlgorithm,''Value'')};','if strcmp(EPmain.pca.theAlgorithm,''SAS''),EPmain.pca.rotopt=3;else EPmain.pca.rotopt=4;end;','ep(''start'');'],...
                    'TooltipString','The version of the rotation algorithm rotation to use.',...
                    'Value',find(strcmp(EPmain.pca.theAlgorithm,theList)),'Position',[110 370 80 20]);

%                 verList=ver;
%                 if isempty(find(strcmp('Statistics and Machine Learning Toolbox',{verList.Name})))
%                     set(EPmain.handles.pca.theAlgorithm,'enable','off');
%                     EPmain.pca.theAlgorithm='SAS';
%                 end
            end

            h = uicontrol('Style','text',...
                'String','Decomposition','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[25 350 150 20]);

            EPmain.handles.pca.decomp= uicontrol('Style','popupmenu','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'String',{'SVD','NIPALS'},...
                'CallBack',['global EPmain;','EPmain.pca.decomp=get(EPmain.handles.pca.decomp,''Value'');','ep(''start'');'],...
                'TooltipString','The method by which the relationship matrix is to be decomposed.',...
                'Value',EPmain.pca.decomp,'Position',[25 330 80 20]);
            
            verList=ver;
            if isempty(find(strcmp('Deep Learning Toolbox',{verList.Name})))
                set(EPmain.handles.pca.decomp,'enable','off');
                EPmain.pca.decomp=1;
            end
            
            h = uicontrol('Style','text',...
                'String','Relationships','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[25 310 80 20]);
            
            EPmain.handles.pca.rel = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',{'correlation','covariance'},...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.pca.rel,''Value'');','if tempVar ~=0,EPmain.pca.rel=tempVar;end;','if isempty(tempVar),EPmain.pca.rel=tempVar;end;','ep(''start'');'],...
                'TooltipString','The relationship matrix to be decomposed.',...
                'Value',EPmain.pca.rel,'Position',[20 290 80 20]);
            
            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Loadings','HorizontalAlignment','left',...
                'Position',[110 310 80 20]);
            
            EPmain.handles.pca.loadings = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',{'none','Kaiser','covariance','C-M'},...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.pca.loadings,''Value'');','if tempVar ~=0,EPmain.pca.loadings=tempVar;end;','if isempty(tempVar),EPmain.pca.loadings=tempVar;end;','ep(''start'');'],...
                'TooltipString','Weighting for the factor loadings.',...
                'Value',EPmain.pca.loadings,'Position',[110 290 80 20]);
            
            if any(strcmp(PCArotation,{'Infomax','JADE','SOBI','fastICA'}))
                set(EPmain.handles.pca.loadings,'enable','off');
            end
            switch PCArotation
                case 'Promax'
                    set(rotoptLabel,'String','kappa');
                case 'Oblimin'
                    set(rotoptLabel,'String','delta');
                case 'SOBI'
                    set(rotoptLabel,'String','number of lags');
                case 'Extended-Infomax'
                    set(rotoptLabel,'String','training blocks');
                case 'Geomin'
                    set(rotoptLabel,'String','epsilon');
                otherwise
                    set(EPmain.handles.pca.rotopt,'enable','off');
            end
            if EPmain.pca.rotFlag
                EPmain.pca.rotFlag=0;
                switch PCArotation
                    case 'Promax'
                        EPmain.pca.rotopt=3;
                    case 'Oblimin'
                        EPmain.pca.rotopt=0;
                    case 'SOBI'
                        EPmain.pca.rotopt=100;
                    case 'Extended-Infomax'
                        EPmain.pca.rotopt=1;
                    case 'Geomin'
                        EPmain.pca.rotopt=.01;
                end
                set(EPmain.handles.pca.rotopt,'Value',EPmain.pca.rotopt);
                ep('start');
            end
            
            EPmain.handles.pca.parametric= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
                'String','Parametric Analysis',...
                'Value',EPmain.pca.parametric,'Position',[20 270 150 20]);
            
            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','# Factors (0=scree)','HorizontalAlignment','left',...
                'Position',[25 250 150 20]);
            
            EPmain.handles.pca.facNum= uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'String',EPmain.pca.facNum,...
                'Callback', ['global EPmain;','EPmain.pca.facNum=str2num(get(EPmain.handles.pca.facNum,''String''));'],...
                'Position',[25 230 150 20]);
            
            h = uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Title of PCA','HorizontalAlignment','left',...
                'Position',[25 210 150 20]);
            
            EPmain.handles.pca.name= uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'String',EPmain.pca.name,...
                'Callback', ['global EPmain;','EPmain.pca.name=get(EPmain.handles.pca.name,''String'');'],...
                'Position',[25 190 150 20]);
            
            if ~isempty(EPdataset.dataset)
                for i=1:length(EPdataset.dataset)
                    fileName=EPdataset.dataset(i).dataName;
                    if strcmp(EPdataset.dataset(i).saved,'no')
                        fileName=['*' fileName];
                    end
                    tableData{i,1}=fileName;
                end
            else
                tableData=[];
            end
            
            tableNames{1}='data';
            columnEditable =  false;
            ColumnFormat{1}=[];
            
            EPmain.handles.pca.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
                'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
                'CellSelectionCallback',@pickPCAdata,...
                'ColumnWidth',{300},'Position',[20 40 170 150]);
            
        end
        
        EPmain.handles.pca.hQuitRead = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', ['global EPmain;','EPmain.mode=''main'';EPmain.pca.crossVerifyPCA=[];','ep(''start'');']);
        
    case 'startWindow'
        
        set(EPmain.handles.hMainWindow,'Name', 'Window Data');
        
        EPmain.handles.window.prefs = uicontrol('Style', 'pushbutton', 'String', '•','FontSize',EPmain.fontsize,...
            'Position', [5 485 10 10], 'Callback', ['global EPmain;','EPmain.prefReturn=''window'';','EPmain.mode=''preferenceWindow'';','ep(''start'');']);
        
        %initialize Window Data parameters if entering in from the Main Menu either for the first time or list has changed.
        if ~isfield(EPmain.window,'measure') || (EPmain.window.dataset > length(EPdataset.dataset)) || ~strcmp(EPmain.window.datasetName,EPdataset.dataset(EPmain.window.dataset).dataName) ||...
                (length(EPmain.window.dataNames) ~= length(EPdataset.dataset)) || (~isempty(EPmain.window.lastChange) && ~strcmp(EPmain.window.lastChange,EPdataset.dataset(EPmain.window.dataset).lastChange))
            
            if isfield(EPmain.window,'dataset')
                if (EPmain.window.dataset > length(EPdataset.dataset)) || (~isempty(EPmain.window.datasetName) && ~strcmp(EPmain.window.datasetName,EPdataset.dataset(EPmain.window.dataset).dataName))
                    disp('The working set has been changed so reinitializing the Window Data pane.')
                end
            end
            
            EPmain.window.minFacVar=EPmain.preferences.window.minFacVar;
            EPmain.window.measure=1;
            EPmain.window.chanGrp=1;
            EPmain.window.specSelect(1)=false;
            EPmain.window.factor=1;
            EPmain.window.FFTunits=4;
            EPmain.window.sampAdapt=0;
            if ~isfield(EPmain.window,'dataset')
                EPmain.window.dataset=1;
            end

            EPmain.window.dataNames=cell(length(EPdataset.dataset),1);
            EPmain.window.factorData=[];
            EPmain.window.aveData=[];
            factorDataNum=0;
            if ~isempty(EPdataset.dataset)
                for i=1:length(EPdataset.dataset)
                    if ~strcmp(EPdataset.dataset(i).dataType,'continuous') && any(ismember({'EEG','BSC'},EPdataset.dataset(i).chanTypes)) &&...
                            ~isempty(find(ismember({'RAW','AVG'},EPdataset.dataset(i).subTypes))) &&...
                            ~isempty(find(strcmp('SGL',EPdataset.dataset(i).cellTypes)))
                        EPmain.window.aveData(end+1)=i; %keep track of which datasets are suitable for windowing
                        EPmain.window.dataNames{i,1}=EPdataset.dataset(i).dataName;
                    end
                    
                    if isfield(EPdataset.dataset(i).pca,'PCAmode')
                        if isfield(EPdataset.dataset(i).pca,'PCAmode2')
                            if strcmp(EPdataset.dataset(i).pca.PCAmode2,'spat')
                                factorDataNum=factorDataNum+1;
                                EPmain.window.factorData(factorDataNum).name=EPdataset.dataset(i).dataName;
                                EPmain.window.factorData(factorDataNum).FacPat=EPdataset.dataset(i).pca.FacPatST;
                                facNames=EPdataset.dataset(i).facNames(find(strcmp('SGL',EPdataset.dataset(i).facTypes)));
                                if length(facNames)==size(EPmain.window.factorData(factorDataNum).FacPat,2)
                                    EPmain.window.factorData(factorDataNum).facNames=facNames;
                                else
                                    for fac =1:size(EPmain.window.factorData(factorDataNum).FacPat,2)
                                        EPmain.window.factorData(factorDataNum).facNames{fac}=['fac' num2str(fac)];
                                    end
                                end
                            end
                        end
                        if strcmp(EPdataset.dataset(i).pca.PCAmode,'spat')
                            factorDataNum=factorDataNum+1;
                            EPmain.window.factorData(factorDataNum).name=EPdataset.dataset(i).dataName;
                            EPmain.window.factorData(factorDataNum).FacPat=EPdataset.dataset(i).pca.FacPat;
                            facNames=EPdataset.dataset(i).facNames(find(strcmp('SGL',EPdataset.dataset(i).facTypes)));
                            if length(facNames)==size(EPmain.window.factorData(factorDataNum).FacPat,2)
                                EPmain.window.factorData(factorDataNum).facNames=facNames;
                            else
                                for fac =1:size(EPmain.window.factorData(factorDataNum).FacPat,2)
                                    EPmain.window.factorData(factorDataNum).facNames{fac}=['fac' num2str(fac)];
                                end
                            end
                        end
                    end
                end
                EPmain.window.dataset=EPmain.window.aveData(end); %the active dataset
                EPmain.window.datasetName=EPdataset.dataset(EPmain.window.dataset).dataName;
                if isfield(EPdataset.dataset(EPmain.window.dataset),'lastChange') && ~isempty(EPdataset.dataset(EPmain.window.dataset).lastChange)
                    EPmain.window.lastChange=EPdataset.dataset(EPmain.window.dataset).lastChange;
                else
                    EPmain.window.lastChange=cell(0);
                end
            else
                EPmain.window.dataNames='none';
                EPmain.window.lastChange=cell(0);
            end
            EPmain.window.sampStart=1;
            EPmain.window.sampEnd=length(EPdataset.dataset(EPmain.window.dataset).timeNames);
            if isempty(EPmain.window.sampEnd)
                EPmain.window.sampEnd=1;
            end
            EPmain.window.HzStart=1;
            EPmain.window.HzEnd=length(EPdataset.dataset(EPmain.window.dataset).freqNames);
            if isempty(EPmain.window.HzEnd)
                EPmain.window.HzEnd=1;
            end
            cellNames=EPdataset.dataset(EPmain.window.dataset).cellNames(find(strcmp('SGL',EPdataset.dataset(EPmain.window.dataset).cellTypes)));
            [u i]=unique(cellNames,'first');
            EPmain.window.inCells=cellNames(sort(i));
            EPmain.window.outCells=EPmain.window.inCells;
        end
        
%         uicontrol('Style','text','FontSize',EPmain.fontsize,...
%             'String','Dataset','HorizontalAlignment','left',...
%             'Position',[25 470 140 20]);
        
        EPmain.handles.window.dataset = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.window.dataNames(EPmain.window.aveData),...
            'Value',find(EPmain.window.dataset==EPmain.window.aveData),'Position',[5 450 190 20],...
            'TooltipString','Dataset to be windowed.',...
            'Callback', @changeWindowDataset);
        
        if isempty(EPdataset.dataset)
            set(EPmain.handles.window.dataset,'enable','off');
        end
        
        uicontrol('Style','text','FontSize',EPmain.fontsize,...
            'String','Measure','HorizontalAlignment','left',...
            'Position',[25 430 140 20]);
        
        if isempty(EPdataset.dataset(EPmain.window.dataset).timeNames) && ~isempty(EPdataset.dataset(EPmain.window.dataset).freqNames)
            measures={'mean','minHzPeak','maxHzPeak','minHzLatency','maxHzLatency','minHzCentroid','maxHzCentroid','meanStdVar','meanSME','meanNoise','meanERA'};
        elseif ~isempty(EPdataset.dataset(EPmain.window.dataset).timeNames) && isempty(EPdataset.dataset(EPmain.window.dataset).freqNames)
            measures={'mean','minpeak','maxpeak','minlatency','maxlatency','mincentroid','maxcentroid','meanStdVar','meanSME','meanNoise','meanERA'};
        else
            measures={'mean','minpeak','maxpeak','minlatency','maxlatency','mincentroid','maxcentroid','minHzPeak','maxHzPeak','minHzLatency','maxHzLatency','minHzCentroid','maxHzCentroid','meanStdVar','meanSME','meanNoise','meanERA'};
        end
        if ~isempty(EPdataset.dataset(EPmain.window.dataset).trialSpecNames) && any(ismember({'RT','ACC'},EPdataset.dataset(EPmain.window.dataset).trialSpecNames))
            measures{end+1}='behavioral';
        end
        if any(strcmp('BSC',EPdataset.dataset(EPmain.window.dataset).chanTypes))
            measures{end+1}='BOSCepisodes';
        end
        
        EPmain.handles.window.measure= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',measures,...
            'Value',EPmain.window.measure,'Position',[5 410 110 20],...
            'Callback', ['global EPmain;','tempVar=get(EPmain.handles.window.measure,''Value'');','if tempVar ~=0,EPmain.window.measure=tempVar;end;','if isempty(tempVar),EPmain.window.measure=tempVar;end;','ep(''start'');']);
        
        if strcmp(measures(EPmain.window.measure),'behavioral')
            
            
            
            
        else
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','samples','HorizontalAlignment','left',...
                'Position',[110 430 55 20]);
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','ms','HorizontalAlignment','left',...
                'Position',[165 430 40 20]);
            
            EPmain.handles.window.sampAdapt= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.window.sampAdapt,...
                'Position',[110 410 40 20],'Callback',@windowSampAdapt,...
                'TooltipString','Samples surrounding peak that are averaged together for peak measure (0 for peak only).');
            
            EPmain.handles.window.msAdapt= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',round(EPmain.window.sampAdapt*(1000/EPdataset.dataset(EPmain.window.dataset).Fs)),...
                'Position',[160 410 40 20],'Callback',@windowSampAdapt,...
                'TooltipString','Ms surrounding peak that are averaged together for peak measure (0 for peak only).');
            
            if ~any(strcmp(measures{EPmain.window.measure},{'minpeak','maxpeak'}))
                set(EPmain.handles.window.sampAdapt,'enable','off');
                set(EPmain.handles.window.msAdapt,'enable','off');
            end
            
            uicontrol('Style','text',...
                'String','samples','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[25 390 140 20]);
            
            if strcmp(EPdataset.dataset(EPmain.window.dataset).timeUnits,'per')
                theUnit='%';
            else
                theUnit='ms';
            end
            uicontrol('Style','text',...
                'String',theUnit,'HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[120 390 140 20]);
            
            EPmain.handles.window.sampStart= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.window.sampStart,...
                'Position',[25 370 40 20],'Callback',@windowSampStart);
            
            EPmain.handles.window.sampEnd= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.window.sampEnd,...
                'Position',[60 370 40 20],'Callback',@windowSampEnd);
            
            EPmain.handles.window.msStart= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',round((EPmain.window.sampStart-EPdataset.dataset(EPmain.window.dataset).baseline-1)*(1000/EPdataset.dataset(EPmain.window.dataset).Fs)),...
                'Position',[100 370 40 20],'Callback',@windowSampStart);
            
            EPmain.handles.window.msEnd= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',round((EPmain.window.sampEnd-EPdataset.dataset(EPmain.window.dataset).baseline)*(1000/EPdataset.dataset(EPmain.window.dataset).Fs)),...
                'Position',[140 370 40 20],'Callback',@windowSampEnd);
            
            if isempty(EPdataset.dataset(EPmain.window.dataset).timeNames)
                set(EPmain.handles.window.sampStart,'enable','off');
                set(EPmain.handles.window.sampEnd,'enable','off');
                set(EPmain.handles.window.msStart,'enable','off');
                set(EPmain.handles.window.msEnd,'enable','off');
            end
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Bins','HorizontalAlignment','left',...
                'Position',[25 350 140 20]);
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Hz','HorizontalAlignment','left',...
                'Position',[120 350 140 20]);
            
            EPmain.handles.window.binStart= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.window.HzStart,...
                'Position',[25 330 40 20],'Callback',@windowHzStart);
            
            EPmain.handles.window.binEnd= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                'String',EPmain.window.HzEnd,...
                'Position',[60 330 40 20],'Callback',@windowHzEnd);
            
            HzStart=0;
            HzEnd=0;
            if ~isempty(EPdataset.dataset(EPmain.window.dataset).freqNames)
                HzStart=round(EPdataset.dataset(EPmain.window.dataset).freqNames(EPmain.window.HzStart)*10)/10;
                HzEnd=round(EPdataset.dataset(EPmain.window.dataset).freqNames(EPmain.window.HzEnd)*10)/10;
            end
            EPmain.handles.window.HzStart= uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String',sprintf('%5.1f',HzStart),...
                'Position',[100 330 40 20]);
            
            EPmain.handles.window.HzEnd= uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String',sprintf('%5.1f',HzEnd),...
                'Position',[140 330 40 20]);
            
            h = uicontrol('Style','text','HorizontalAlignment','left','String', 'units','FontSize',EPmain.fontsize,...
                'ForegroundColor','black','Position',[160 310 30 20]);
            EPmain.handles.window.FFTunits = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'Value',EPmain.window.FFTunits,'Position',[155 290 65 20],...
                'String',{'cm','am','pw','dB'},...
                'Callback', ['global EPmain;','EPmain.window.FFTunits=get(EPmain.handles.window.FFTunits,''value'');','ep(''start'')'],...
                'TooltipString','Units for spectral data.');
            
            if isempty(EPdataset.dataset(EPmain.window.dataset).freqNames)
                set(EPmain.handles.window.HzStart,'enable','off');
                set(EPmain.handles.window.HzEnd,'enable','off');
                set(EPmain.handles.window.binStart,'enable','off');
                set(EPmain.handles.window.binEnd,'enable','off');
                set(EPmain.handles.window.FFTunits,'enable','off');
            elseif ~isempty(EPdataset.dataset(EPmain.window.dataset).relNames)
                set(EPmain.handles.window.FFTunits,'enable','off');
            end
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Baseline:','HorizontalAlignment','left',...
                'Position',[10 280 70 20]);
            
            EPmain.handles.window.baseline= uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String',(EPdataset.dataset(EPmain.window.dataset).baseline)*(1000/EPdataset.dataset(EPmain.window.dataset).Fs),'HorizontalAlignment','left',...
                'Position',[60 280 30 20]);
            
            if isfield(EPchanGrp,'group')
                if isempty(EPchanGrp.group)
                    chanGrpNames='-all-';
                else
                    if isempty({EPchanGrp.group.name})
                        chanGrpNames='-all-';
                    else
                        chanGrpNames=[{EPchanGrp.group.name} '-all-'];
                    end
                end
            else
                chanGrpNames='-all-';
            end
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Channels','HorizontalAlignment','left',...
                'Position',[5 310 55 20]);
            
            EPmain.handles.window.chanGrp = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',chanGrpNames,...
                'Value',EPmain.window.chanGrp,'Position',[60 310 100 20],...
                'Callback', ['global EPmain;','tempVar=get(EPmain.handles.window.chanGrp,''Value'');','if tempVar ~=0,EPmain.window.chanGrp=tempVar;end;','if isempty(tempVar),EPmain.window.chanGrp=tempVar;end;','ep(''start'');']);
            
            EPmain.handles.window.channels = uicontrol('Style', 'pushbutton', 'String', 'Channels','FontSize',EPmain.fontsize,...
                'Position', [90 275 60 35], 'Callback', ['global EPmain EPdataset;,ep_chanGrp(ep_loadEPdataset(EPmain.window.dataset),EPmain.window.factorData);']);
            
            if isempty(EPdataset.dataset)
                set(EPmain.handles.window.chanGrp,'enable','off');
            end
            
            uicontrol('Style','text','FontSize',EPmain.fontsize,...
                'String','Factor','HorizontalAlignment','left',...
                'Position',[5 255 40 20]);
            
            if ~isempty(EPdataset.dataset(EPmain.window.dataset).facNames(find(strcmp('SGL',EPdataset.dataset(EPmain.window.dataset).facTypes))))
                EPmain.handles.window.factor = uicontrol('Style', 'popupmenu', 'String', EPdataset.dataset(EPmain.window.dataset).facNames(find(strcmp('SGL',EPdataset.dataset(EPmain.window.dataset).facTypes))),'FontSize',EPmain.fontsize,...
                    'Value',EPmain.window.factor,'Position', [50 240 100 35],...
                    'Callback', ['global EPmain;','EPmain.window.factor=get(EPmain.handles.window.factor,''Value'');','ep(''start'');']);
            else
                EPmain.handles.window.factors = uicontrol('Style', 'popupmenu', 'String', 'none','FontSize',EPmain.fontsize,...
                    'Position', [70 240 100 35],'enable','off');
            end
            
        end

%         EPmain.handles.window.chanGrp = uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
%             'String','SME',...
%             'Value',EPmain.window.SME,'Position',[150 255 50 20],...
%             'TooltipString','Implement SME data quality measure.');
% 
        %cell selection
        if any(ismember({'EEG','BSC'},EPdataset.dataset(EPmain.window.dataset).chanTypes)) && ~isempty(find(ismember({'RAW','AVG'},EPdataset.dataset(EPmain.window.dataset).subTypes)))...
                && ~isempty(find(strcmp('SGL',EPdataset.dataset(EPmain.window.dataset).cellTypes))) && ~(~isempty(EPdataset.dataset(EPmain.window.dataset).facNames) && isempty(find(strcmp('SGL',EPdataset.dataset(EPmain.window.dataset).cellTypes))))
            cellNames=unique(EPdataset.dataset(EPmain.window.dataset).cellNames(find(strcmp('SGL',EPdataset.dataset(EPmain.window.dataset).cellTypes))),'first');
            cellNames{end+1}='none';
            for i=1:length(cellNames)-1
                tableData(i,1)=EPmain.window.outCells(i);
                tableData(i,2)=EPmain.window.inCells(i);
            end
            
            tableNames{1}='outCells';
            tableNames{2}='inCells';
            columnEditable=[true true];
            ColumnFormat{1}='char';
            ColumnFormat{2}='char';
            
            EPmain.handles.window.cellTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
                'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
                'ColumnWidth',{50 50},...
                'Position',[10 125 150 130],...
                'CellEditCallback', @windowCellTable);

            EPmain.handles.window.loadWindowtable = uicontrol('Style', 'pushbutton', 'String', 'Load','FontSize',EPmain.fontsize,...
                'Position', [165 190 40 20], 'Callback', @loadWindowtable);

            EPmain.handles.window.saveWindowtable = uicontrol('Style', 'pushbutton', 'String', 'Save','FontSize',EPmain.fontsize,...
                'Position', [165 170 40 20], 'Callback', @saveWindowtable);

            EPmain.handles.window.clearWindowtable = uicontrol('Style', 'pushbutton', 'String', 'Reset','FontSize',EPmain.fontsize,...
                'Position', [165 150 40 20], 'Callback', @resetWindowtable);

            EPmain.handles.window.undoWindowtable = uicontrol('Style', 'pushbutton', 'String', 'Undo','FontSize',EPmain.fontsize,...
                'Position', [165 130 40 20], 'Callback', @undoWindowtable);

            if isempty(EPmain.window.undo)
                set(EPmain.handles.window.undoWindowtable,'enable','off');
            end

        else
            uicontrol('Style','text','HorizontalAlignment','left','String', 'This average file cannot be windowed.',...
                'Position',[10 125 190 130]);
        end
        
        %subject spec selection
        subjectSpecNames=EPdataset.dataset(EPmain.window.dataset).subjectSpecNames;
        if length(EPmain.window.specSelect) ~= length(subjectSpecNames)
            EPmain.window.specSelect=repmat(false,length(subjectSpecNames),1);
        end
        if ~isempty(subjectSpecNames)
            tableData=[];
            for i=1:length(subjectSpecNames)
                tableData{i,1}=EPmain.window.specSelect(i);
                tableData{i,2}=subjectSpecNames{i};
            end
        else
            tableData=[];
        end
        
        tableNames{1}='Select';
        tableNames{2}='Specs';
        columnEditable=[true false];
        ColumnFormat{1}='logical';
        ColumnFormat{2}='char';
        
        EPmain.handles.window.specsTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
            'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
            'ColumnWidth',{40 100},...
            'Position',[10 40 190 80],...
            'CellEditCallback', @windowSpecsTable);
        
        EPmain.handles.window.window = uicontrol('Style', 'pushbutton', 'String', 'Window','FontSize',EPmain.fontsize,...
            'Position', [10 0 60 35], 'Callback', 'ep(''WindowData'')');
        
        EPmain.handles.window.automatic = uicontrol('Style', 'pushbutton', 'String', 'AutoPCA','FontSize',EPmain.fontsize,...
            'Position', [70 0 60 35], 'Callback', 'ep(''AutoPCA'')');
        
        if isempty(EPdataset.dataset(EPmain.window.dataset).facNames) || strcmp(measures(EPmain.window.measure),'behavioral')
            set(EPmain.handles.window.automatic,'enable','off');
        end
        
        EPmain.handles.window.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [130 0 60 35], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        
    case 'WindowData'
        %output the windowed measure to a file for later ANOVA
        
        measureList=get(EPmain.handles.window.measure,'String');
        measureNum = EPmain.window.measure;
        measure=measureList{measureNum};
        
        EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);
        if  ~strcmp(measure,'behavioral')
            set(EPmain.handles.window.channels,'enable','off');
            set(EPmain.handles.window.automatic,'enable','off');
            set(EPmain.handles.window.FFTunits,'enable','off');
        end
        drawnow

        [FileName,PathName,FilterIndex] = uiputfile('*.txt','Output Measures File',EPdataset.dataset(EPmain.window.dataset).dataName);
        
        if isnumeric(FileName)
            if FileName == 0
                msg{1}='No file name specified.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
        end
        
        [pathstr, name, ext] = fileparts(FileName);
        
        if ~strcmp(ext,'.txt')
            ext=[ext '.txt'];
        end
        
        ep_tictoc('begin');
                
        EPdata=ep_stripAdds(ep_loadEPdataset(EPmain.window.dataset));
        sampStart = EPmain.window.sampStart;
        sampEnd = EPmain.window.sampEnd;
        HzStart = EPmain.window.HzStart;
        HzEnd = EPmain.window.HzEnd;
        factor = EPmain.window.factor;
        theChanGrp = EPmain.window.chanGrp;
        
        if isempty(EPchanGrp) || theChanGrp > length(EPchanGrp.group)
            chanGrp.channel=[1:length(EPdata.chanNames)];
            chanGrp.name='-all-';
            chanGrp.areaName=[EPdata.chanNames; 'none'];
        else
            chanGrp=EPchanGrp.group(theChanGrp);
        end
        
        areaList=setdiff(unique(chanGrp.channel),length(chanGrp.areaName));
        if length(areaList) > 1 && ~isempty(EPdata.facVecS)
            warndlg('Each spatial PCA factor is a virtual electrode covering the entire head.  The electrode regions would be 100% correlated so the ANOVA would fail.');
        end
        
        if ~isempty(EPdataset.dataset(EPmain.window.dataset).relNames)
            disp('For coherence data, if there is one channel in an area then the correlations with all other channels will be used.');
            disp('The absolute value will be taken since otherwise they would tend to sum to zero (at least if average referenced).');
            disp('If a single common reference is used, it will be excluded since coherence with a flat channel produces a not-a-number result.');
            disp('If there is more than one channel in an area, then only the correlations between the channels in that area will be used.');
            disp('The absolute value will not be taken since the sign is meaningful and will not necessarily sum to zero.');
        end
        
        inputcells=[];
        outCellNames=[];
        for iCell=1:length(EPmain.window.inCells)
            if ~isempty(EPmain.window.inCells{iCell}) && ~isempty(EPmain.window.outCells{iCell})
                if ~any(strcmp(EPmain.window.outCells(iCell),outCellNames))
                    outCellNames{end+1}=EPmain.window.outCells{iCell};
                    inputcells{end+1}=[];
                end
                theCell=find(strcmp(EPmain.window.inCells(iCell),EPdata.cellNames));
                if isempty(theCell)
                    disp(['warning: The cell ' EPmain.window.inCells{iCell} ' is missing from the dataset and will be ignored.']);
                else
                    inputcells{find(strcmp(EPmain.window.outCells(iCell),outCellNames))}(end+1:end+length(theCell))=theCell;
                end
            end
        end
        
        emptyCells=find(cellfun(@isempty,inputcells));
        if ~isempty(emptyCells)
            disp('warning: some outCells were not assigned inCells and are being dropped.')
            inputcells(emptyCells)=[];
            outCellNames(emptyCells)=[];
        end
        
        subjectSpecs =find(EPmain.window.specSelect);
        
        if  strcmp(measure,'behavioral')
            if any(strcmp('RT',EPdataset.dataset(EPmain.window.dataset).trialSpecNames))
                ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'RT');
                if EPtictoc.stop
                    EPtictoc.stop=0;
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    ep('start');
                end
            end
            if any(strcmp('ACC',EPdataset.dataset(EPmain.window.dataset).trialSpecNames))
                ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'ACC');
                if EPtictoc.stop
                    EPtictoc.stop=0;
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    ep('start');
                end
            end
            if ~isempty(EPdata.taskSpecs) && strcmp('average',EPdataset.dataset(EPmain.window.dataset).dataType)
                for iMeasure=1:length(EPdata.taskMeasNames)
                    ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, 'task', [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,EPdata.taskMeasNames{iMeasure});
                    if EPtictoc.stop
                        EPtictoc.stop=0;
                        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                        ep('start');
                    end
                end
            end
        elseif any(strcmp(measure,{'meanSME','meanERA'}))
            if ~strcmp(EPdata.dataType,'average')
                msg{1}=['The ' measure ' option is only intended for average files.'];
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            if isempty(EPdata.history)
                msg{1}='The average file does not have a history record with which to recreate the averaging process.  An average file generated by version 2.96 or later is required.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            averageRecord=0;
            for iRecord=size(EPdata.history,1):-1:1
                if strcmp(EPdata.history{iRecord,2},'Averaged the data.')
                    averageRecord=iRecord;
                    break
                end
            end
            if averageRecord==0
                msg{1}='The average file does not have the averaging step in its history record with which to recreate the averaging process.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end

            [EPdata]=ep_stripAdds(EPdata);
            if isempty(EPdata.data)
                msg{1}='Error: The file had no data left after additions were removed.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end

            newDir = uigetdir(pwd,'Select directory with the single-trial subject files used to generate the average file.');
            if newDir == 0
                msg{1}='No directory selected.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end

            inputFormat='';
            formatType='';
            fileType='';
            averagingMethod='';
            methodName='';
            multiSubjectNumber=[];
            multiSessionNumber=[];
            multiCellNumber=[];
            cfg=[];
            if ~isfield(EPdata.history{averageRecord,4},'average')
                cfg.behav.codeCorrect=EPdata.history{averageRecord,4}.codeCorrect;
                cfg.behav.codeError=EPdata.history{averageRecord,4}.codeError;
                cfg.behav.codeTimeout=EPdata.history{averageRecord,4}.codeTimeout;
            else
                cfg.behav.codeCorrect=EPdata.history{averageRecord,4}.average.codeCorrect;
                cfg.behav.codeError=EPdata.history{averageRecord,4}.average.codeError;
                cfg.behav.codeTimeout=EPdata.history{averageRecord,4}.average.codeTimeout;
            end
            cfg.behav.ACC='no ACC';
            cfg.behav.RT='no RT';

            for iChildren=1:length(EPdata.history{averageRecord,3}.children)
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'importFormat')
                    formatType=EPdata.history{averageRecord,3}.children(iChildren).String{EPdata.history{averageRecord,3}.children(iChildren).Value};
                end
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'fileType')
                    fileType=EPdata.history{averageRecord,3}.children(iChildren).String{EPdata.history{averageRecord,3}.children(iChildren).Value};
                end
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'method')
                    averagingMethod=EPdata.history{averageRecord,3}.children(iChildren).String{EPdata.history{averageRecord,3}.children(iChildren).Value};
                end
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'freqMethod')
                    methodName=EPdata.history{averageRecord,3}.children(iChildren).String{EPdata.history{averageRecord,3}.children(iChildren).Value};
                end
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'subject')
                    multiSubjectNumber=str2num(EPdata.history{averageRecord,3}.children(iChildren).String);
                    if isnan(multiSubjectNumber)
                        multiSubjectNumber=[];
                    end
                end
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'session')
                    multiSessionNumber=str2num(EPdata.history{averageRecord,3}.children(iChildren).String);
                    if isnan(multiSessionNumber)
                        multiSessionNumber=[];
                    end
                end
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'cell')
                    multiCellNumber=str2num(EPdata.history{averageRecord,3}.children(iChildren).String);
                    if isnan(multiCellNumber)
                        multiCellNumber=[];
                    end
                end
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'averageType')
                    averageType=EPdata.history{averageRecord,3}.children(iChildren).String{EPdata.history{averageRecord,3}.children(iChildren).Value};
                end
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'maxRT')
                    cfg.behav.maxRT=str2double(EPdata.history{averageRecord,3}.children(iChildren).String);
                end
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'minRT')
                    cfg.behav.minRT=str2double(EPdata.history{averageRecord,3}.children(iChildren).String);
                end
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'RTmethod')
                    cfg.behav.RTmethod=EPdata.history{averageRecord,3}.children(iChildren).String{EPdata.history{averageRecord,3}.children(iChildren).Value};
                end
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'dropBad')
                    cfg.behav.dropBad=EPdata.history{averageRecord,3}.children(iChildren).Value;
                end
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'dropError')
                    cfg.behav.dropError=EPdata.history{averageRecord,3}.children(iChildren).Value;
                end
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'dropTimeout')
                    cfg.behav.dropTimeout=EPdata.history{averageRecord,3}.children(iChildren).Value;
                end
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'ACC')
                    cfg.behav.ACC=EPdata.history{averageRecord,3}.children(iChildren).String;
                end
                if strcmp(EPdata.history{averageRecord,3}.children(iChildren).Tag,'RT')
                    cfg.behav.RT=EPdata.history{averageRecord,3}.children(iChildren).String;
                end
            end
            if ~isempty(formatType)
                [~,~,inputFormat]=ep_fileFormats(fileType,formatType);
            end

            if ~strcmp(averagingMethod,'Average') || ~strcmp(methodName,'Mean') || ~strcmp(averageType,'subject')
                msg{1}=[measure ' is currently only available for average subject files that used the Average procedure and the Mean method.'];
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end

            %Determine which original files are in the designated directory.
            fileList=cell(0);
            for iFile=1:length(EPdata.history{averageRecord,5})
                [pathstr, name, ext] = fileparts(EPdata.history{averageRecord,5}(iFile));
                if ~exist([newDir filesep name ext],'file')
                    disp(['The file ' name ext ' is not present.'])
                else
                    fileList{end+1,1}=[newDir filesep name ext];
                end
            end

            if isempty(fileList)
                msg{1}='None of the original files are in the designated directory.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end

            %output the trial-level data in temp files, including just the retained trials
            SMElist=ep_averageData(fileList,inputFormat,fileType,'Average',[],measure,[],multiSessionNumber,multiSubjectNumber,multiCellNumber,cfg,EPdata.history{averageRecord,4},'subject',[],3,1);
            if EPtictoc.stop
                EPtictoc.stop=0;
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end

            if isempty(SMElist)
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end

            %further transforms and edits
            if size(EPdata.history,1) > averageRecord
                for iRecord = averageRecord+1:size(EPdata.history,1)
                    if contains(EPdata.history{iRecord,2},'Transformed the file')
                        transformSettings.importFormat=1;
                        transformSettings.outputFormat=1;
                        transformSettings.refChan1=[];
                        transformSettings.refChan2=[];
                        transformSettings.baselineStart=[];
                        transformSettings.baselineEnd=[];
                        transformSettings.preStim=[];
                        transformSettings.smoothing=[];
                        transformSettings.detrend=[];
                        transformSettings.mainsFix=[];
                        transformSettings.cfgFilter=[];
                        for iChildren=1:length(EPdata.history{iRecord,3}.children)
                            if strcmp(EPdata.history{iRecord,3}.children(iChildren).Tag,'reference')
                                referenceMethod=EPdata.history{iRecord,3}.children(iChildren).String{EPdata.history{iRecord,3}.children(iChildren).Value};
                            end
                            if strcmp(EPdata.history{iRecord,3}.children(iChildren).Tag,'refChan1')
                                transformSettings.refChan1=str2double(EPdata.history{iRecord,3}.children(iChildren).String);
                            end
                            if strcmp(EPdata.history{iRecord,3}.children(iChildren).Tag,'refChan2')
                                transformSettings.refChan2=str2double(EPdata.history{iRecord,3}.children(iChildren).String);
                            end
                            if strcmp(EPdata.history{iRecord,3}.children(iChildren).Tag,'preStim')
                                transformSettings.preStim=str2double(EPdata.history{iRecord,3}.children(iChildren).String);
                            end
                            if strcmp(EPdata.history{iRecord,3}.children(iChildren).Tag,'baselineStart')
                                transformSettings.baselineStart=str2double(EPdata.history{iRecord,3}.children(iChildren).String);
                            end
                            if strcmp(EPdata.history{iRecord,3}.children(iChildren).Tag,'baselineEnd')
                                transformSettings.baselineEnd=str2double(EPdata.history{iRecord,3}.children(iChildren).String);
                            end
                            if strcmp(EPdata.history{iRecord,3}.children(iChildren).Tag,'mainsFix')
                                transformSettings.mainsFix=EPdata.history{iRecord,3}.children(iChildren).String{EPdata.history{iRecord,3}.children(iChildren).Value};
                            end
                            if strcmp(EPdata.history{iRecord,3}.children(iChildren).Tag,'detrend')
                                transformSettings.detrend=EPdata.history{iRecord,3}.children(iChildren).Value;
                            end
                            if strcmp(EPdata.history{iRecord,3}.children(iChildren).Tag,'domain')
                                domainName=EPdata.history{iRecord,3}.children(iChildren).String{EPdata.history{iRecord,3}.children(iChildren).Value};
                                if strcmp(domainName,'no change')
                                    domainName='Time';
                                end
                            end
                            if strcmp(EPdata.history{iRecord,3}.children(iChildren).Tag,'method')
                                methodName=EPdata.history{iRecord,3}.children(iChildren).String{EPdata.history{iRecord,3}.children(iChildren).Value};
                            end
                            if strcmp(EPdata.history{iRecord,3}.children(iChildren).Tag,'dataMode')
                                dataMode=EPdata.history{iRecord,3}.children(iChildren).String{EPdata.history{iRecord,3}.children(iChildren).Value};
                            end
                        end
                        SMElist2=ep_transformData(SMElist,'ep_mat','single_trial','ep_mat',referenceMethod,transformSettings,domainName,methodName,dataMode);
                        if EPtictoc.stop
                            EPtictoc.stop=0;
                            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                            return
                        end
                    end

                    EPdataAVG=EPdata;

                    if strcmp(EPdata.history{iRecord,2},'Sampling rate halved using Samples Pane of the Edit Function.')
                        for iFile=1:length(SMElist2)
                            theFile=SMElist2{iFile};
                            load('-mat', theFile, 'EPdata');
                            newtimeNames=EPdata.timeNames(1:2:end);
                            EPdata=ep_interpTime(EPdata,newtimeNames);
                            save('-mat', theFile, 'EPdata');
                        end
                    end
                    if strcmp(EPdata.history{iRecord,2},'Sampling rate doubled using Samples Pane of the Edit Function.')
                        for iFile=1:length(SMElist2)
                            theFile=SMElist2{iFile};
                            load('-mat', theFile, 'EPdata');
                            sampleSize=median(diff(EPdata.timeNames))/2;
                            newtimeNames=[EPdata.timeNames(1):sampleSize:EPdata.timeNames(end)+sampleSize]';
                            EPdata=ep_interpTime(EPdata,newtimeNames);
                            save('-mat', theFile, 'EPdata');
                        end
                    end
                    if strcmp(EPdata.history{iRecord,2},'Sampling rate changed using Samples Pane of the Edit Function.')
                        for iChildren=1:length(EPdata.history{iRecord,3}.children)
                            if strcmp(EPdata.history{iRecord,3}.children(iChildren).Tag,'FsResample')
                                newRate=EPdata.history{iRecord,3}.children(iChildren).String{EPdata.history{iRecord,3}.children(iChildren).Value};
                            end
                        end
                        for iFile=1:length(SMElist2)
                            theFile=SMElist2{iFile};
                            load('-mat', theFile, 'EPdata');
                            sampleSize=median(diff(EPdata.timeNames))*(EPdata.Fs/newRate);
                            newtimeNames=[EPdata.timeNames(1):sampleSize:EPdata.timeNames(end)+sampleSize]';
                            newNum=floor((newRate/EPdata.Fs)*length(EPdata.timeNames));
                            newtimeNames=newtimeNames(1:newNum);
                            EPdata=ep_interpTime(EPdata,newtimeNames);
                            save('-mat', theFile, 'EPdata');
                        end
                    end
                end
            else
                EPdataAVG=EPdata;
                SMElist2=SMElist;
            end

            if strcmp(measure,'meanSME')
                %generate the SME output
                EPdataSME=ep_newFile(EPdataAVG,length(EPdataAVG.cellNames),length(SMElist2));
                EPdataSME.dataType='average';
                EPdataSME.cellNames=EPdataAVG.cellNames;
                EPdataSME.cellTypes(:) = {'SGL'};
                EPdataSME.cellTypes=EPdataSME.cellTypes(:);
                for iFile=1:length(SMElist2)
                    theFile=SMElist2{iFile};
                    load('-mat', theFile, 'EPdata');
                    [ALLEEG]=ep_ep2alleeg(EPdata);
                    %epoch_list.good_bep_indx(numCells,1)=col vectors with the integers of the trials going into each average
                    epoch_list=[];
                    for iCell=1:length(ALLEEG)
                        theCell=find(strcmp(ALLEEG(iCell).condition,EPdataAVG.cellNames));
                        epoch_list.good_bep_indx=[1:size(ALLEEG(iCell).data,3)];
                        smeData = sme_analytic(ALLEEG(iCell),epoch_list,(([sampStart sampEnd]-EPdata.baseline)*round(1000/EPdata.Fs)));
                        if length(EPdataSME.chanNames)==length(EPdata.chanNames)
                            EPdataSME.data(:,:,theCell,iFile)=repmat(smeData,1,length(EPdata.timeNames));
                        else
                            for iChan=1:length(EPdataSME.chanNames)
                                theChan=find(strcmp(EPdataSME.chanNames{iChan},EPdata.chanNames));
                                if ~isempty(theChan)
                                    EPdataSME.data(iChan,:,theCell,iFile)=repmat(smeData(theChan),1,length(EPdata.timeNames));
                                end
                            end
                        end
                    end
                    EPdataSME.subNames{iFile,1}=EPdata.subNames{1};
                    EPdataSME.subTypes{iFile,1}='AVG';
                end

                if ~isempty(EPdataSME.freqNames) && (EPmain.window.FFTunits ==1)
                    outputdata1=ep_windowData(EPdataSME, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'real');
                    outputdata2=ep_windowData(EPdataSME, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'imag');
                    outputdata=complex(outputdata1,outputdata2);
                else
                    outputdata=ep_windowData(EPdataSME, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'normal');
                end
            else
                %meanERA

                if ~isempty(EPdataAVG.freqNames) && (EPmain.window.FFTunits ==1)
                    outputdata1=ep_windowData(EPdataAVG, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'real',SMElist2);
                    outputdata2=ep_windowData(EPdataAVG, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'imag',SMElist2);
                    outputdata=complex(outputdata1,outputdata2);
                else
                    outputdata=ep_windowData(EPdataAVG, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'normal',SMElist2);
                end
            end

            for iFile=1:length(SMElist)
                delete(SMElist{iFile});
            end
            for iFile=1:length(SMElist2)
                delete(SMElist2{iFile});
            end

        else
            if ~isempty(EPdata.freqNames) && (EPmain.window.FFTunits ==1)
                outputdata1=ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'real');
                outputdata2=ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'imag');
                outputdata=complex(outputdata1,outputdata2);
            else
                outputdata=ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, sampStart, sampEnd, measure, [PathName FileName], factor, HzStart, HzEnd, EPmain.window.FFTunits, EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, EPmain.window.sampAdapt,'normal');
            end
            if EPtictoc.stop
                EPtictoc.stop=0;
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                ep('start');
            end            
            if EPmain.preferences.window.adds && ~isempty(outputdata) && ~any(strcmp(measure,{'meanStdVar','meanSME','meanNoise','meanERA'}))
                windowAdds(inputcells, outCellNames, false);
            end
        end
        
        disp('done');
        ep_tictoc('end');
        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        ep('start')
        
    case 'AutoPCA'
        %output the windowed measure to a file for later ANOVA for an entire PCA dataset, setting channels and windows to
        %the peaks.
        
        ep_tictoc('begin');
        
        EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);
        
        EPdata=ep_stripAdds(ep_loadEPdataset(EPmain.window.dataset));
        
        numChans=length(EPdata.chanNames);
        numPoints=length(EPdata.timeNames);
        numFacs=length(EPdata.facNames);
        numFreqs=length(EPdata.freqNames);
        
        inputcells=[];
        outCellNames=[];
        for iCell=1:length(EPmain.window.inCells)
            if ~isempty(EPmain.window.inCells{iCell}) && ~isempty(EPmain.window.outCells{iCell})
                if ~any(strcmp(EPmain.window.outCells(iCell),outCellNames))
                    outCellNames{end+1}=EPmain.window.outCells{iCell};
                    inputcells{end+1}=[];
                end
                theCell=find(strcmp(EPmain.window.inCells(iCell),EPdata.cellNames));
                if isempty(theCell)
                    disp(['warning: The cell ' EPmain.window.inCells{iCell} ' is missing from the dataset and will be ignored.']);
                else
                    inputcells{find(strcmp(EPmain.window.outCells(iCell),outCellNames))}(end+1:end+length(theCell))=theCell;
                end
            end
        end
        
        emptyCells=find(cellfun(@isempty,inputcells));
        if ~isempty(emptyCells)
            disp('warning: some outCells were not assigned inCells and are being dropped.')
            inputcells(emptyCells)=[];
            outCellNames(emptyCells)=[];
        end
        
        subjectSpecs =find(EPmain.window.specSelect);
                
        disp('Starting AutoPCA.  Window size and adjoining samples option for peak measures ignored.');
        if isfield(EPdata,'facVar')
            factorList=find(EPdata.facVar >= EPmain.window.minFacVar);
            disp(['There are ' num2str(length(factorList)) ' factors that meet the minimum variance criterion of: ' num2str(EPmain.window.minFacVar) '.']);
        else
            factorList=[1:length(EPdata.facNames)];
        end
        
        if isempty(factorList)
            disp('Procedure aborted as there are no factors to process.');
        else
            [FileName,PathName,FilterIndex] = uiputfile('*.*','PCA ANOVA file root name');
            
            if isnumeric(FileName)
                if FileName == 0
                    msg{1}='No file name specified.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
            end
            
            [pathstr, name, ext] = fileparts(FileName);
            
            if ~strcmp(ext,'.txt')
                ext=[ext '.txt'];
            end
            
            for theFactor=factorList
                gave=ep_expandFacs(EPdata,find(strcmp('EEG',EPdata.chanTypes)),[],[],[],theFactor,[]);
                ep_tictoc;if EPtictoc.stop;return;end
                if ~isempty(EPdata.facVecT)
                    [C peakLatency]=max(abs(EPdata.facVecT(:,theFactor)));
                elseif numPoints
                    [C peakLatency]=max(max(reshape(shiftdim(abs(mean(gave,4)),1),numPoints,[])'));
                else
                    peakLatency=1;
                end
                
                startLatency=peakLatency;
                endLatency=peakLatency;
                
                if ~isempty(EPdata.facVecS)
                    [C peakChan]=max(abs(EPdata.facVecS(:,theFactor)));
                else
                    [C peakChan]=max(max(reshape(abs(mean(gave,4)),numChans,[])'));
                end
                
                if ~isempty(EPdata.facVecF)
                    [C peakFreq]=max(abs(EPdata.facVecF(:,theFactor)));
                else
                    if numFreqs
                        [C peakFreq]=max(max(reshape(shiftdim(abs(mean(gave,4)),5),numFreqs,[])'));
                    else
                        peakFreq=1;
                    end
                end
                
                measureNum = get(EPmain.handles.window.measure,'Value');
                switch measureNum
                    case 1
                        measure='mean';
                    case 2
                        measure='minpeak';
                    case 3
                        measure='maxpeak';
                    case 4
                        measure='minlatency';
                    case 5
                        measure='maxlatency';
                    case 6
                        measure='mincentroid';
                        startLatency=1;
                        endLatency=numPoints;
                        disp('Note: AutoPCA will use the entire epoch for the centroid measure, which will not work well if there are multiple ERP components.  For better results, manually window using a smaller window customized for each factor.');
                    case 7
                        measure='maxcentroid';
                        startLatency=1;
                        endLatency=numPoints;
                        disp('Note: AutoPCA will use the entire epoch for the centroid measure, which will not work well if there are multiple ERP components.  For better results, manually window using a smaller window customized for each factor.');
                end
                
                chanGrp=[];
                chanGrp.areaName(1)=EPdata.chanNames(peakChan);
                chanGrp.areaName{2}='none';
                chanGrp.name='autoPCA';
                chanGrp.channel=zeros(length(EPdata.chanNames),1);
                chanGrp.channel(peakChan)=1;
                
                if ~isreal(EPdata.data)
                    outputdata=ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, startLatency, endLatency, measure, [PathName FileName '-' EPdata.facNames{theFactor} ext], theFactor, peakFreq, peakFreq, EPmain.window.FFTunits,EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, 0,'real');
                    outputdata=ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, startLatency, endLatency, measure, [PathName FileName '-' EPdata.facNames{theFactor} ext], theFactor, peakFreq, peakFreq, EPmain.window.FFTunits,EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, 0,'imag');
                else
                    outputdata=ep_windowData(EPdata, chanGrp, inputcells, outCellNames, subjectSpecs, startLatency, endLatency, measure, [PathName FileName '-' EPdata.facNames{theFactor} ext], theFactor, peakFreq, peakFreq, EPmain.window.FFTunits,EPmain.preferences.window.chanGrp, EPmain.preferences.anova.missing, 0,'normal');
                end
            end
            
            if EPmain.preferences.window.adds
                windowAdds(inputcells, outCellNames, true);
            end
            
            disp('done');
        end
        
        ep_tictoc('end');
        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        ep('start')
        
    case 'startANOVA'
        %Set up ANOVA pane of main window.
        
        ep_ANOVApane;
        
    case 'loadANOVA'
        %Load in an initial ANOVA file in order to set up the structure of the analysis.
        
        [ANOVAfile, pathname] = uigetfile({'*.txt;*.csv'},'ANOVA File');
        
        if ANOVAfile == 0
            return;
        end

        EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);
        [ANOVAheader, ANOVAdata] = ep_loadANOVA([pathname filesep ANOVAfile]);
        if isempty(ANOVAdata)
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
        
        numCols=size(ANOVAdata.data,2);
        numSpecs=length(find(strcmp('spec',ANOVAdata.areaNames)));
        
        EPmain.anova.data.columnNames=cell(0);
        if ~strcmp(ANOVAdata.name,'behavioral')
            if ANOVAdata.ANOVAver==1
                for iCol=1:numCols+numSpecs
                    EPmain.anova.data.columnNames{iCol}=[ANOVAdata.cellNames{iCol} '-' ANOVAdata.areaNames{iCol}];
                end
            else
                for iCol=1:numCols+numSpecs
                    if iCol<=length(ANOVAdata.sessNames)
                        sessName=[ANOVAdata.sessNames{iCol} '-'];
                    else
                        sessName='';
                    end
                    EPmain.anova.data.columnNames{iCol}=[ sessName ANOVAdata.cellNames{iCol} '-' ANOVAdata.areaNames{iCol}];
                end
            end
        else
            for iCol=1:numCols+numSpecs
                EPmain.anova.data.columnNames{iCol}=[ANOVAdata.cellNames{iCol}];
            end
        end
        
        for iFac=1:6
            EPmain.anova.data.between(iFac)=numCols+numSpecs+1;
            EPmain.anova.data.factor{iFac}='';
            EPmain.anova.data.levels{iFac}='';
            EPmain.anova.data.betweenName{iFac}='';
            EPmain.anova.data.covariate(iFac)=1;
        end
        
        EPmain.anova.data.leftColumn=1;
        
        EPmain.anova.data.rightColumn=numCols;
        
        EPmain.anova.data.data=ANOVAdata.data;

        EPmain.anova.data.origData=ANOVAdata.origData;
        
        EPmain.anova.data.betweenLvl=ANOVAdata.betweenLvl;
        
        EPmain.anova.data.sessNames=ANOVAdata.sessNames;

        EPmain.anova.data.cellNames=ANOVAdata.cellNames;
        
        EPmain.anova.data.areaNames=ANOVAdata.areaNames;
        
        EPmain.anova.data.name=ANOVAdata.name;
        
        EPmain.anova.data.totalLevels=0;
        
        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        ep('start');
        
    case 'runANOVA'
        %run the ANOVAs on the data
        
        EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);
        
        %initial error checking of settings
        
        verList=ver;
        if (EPmain.anova.method==2) && isempty(find(strcmp('Statistics and Machine Learning Toolbox',{verList.Name})))
                msg{1}='Statistics and Machine Learning Toolbox not installed so conventional ANOVAs not available.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
        end
        switch EPmain.anova.method
            case 1 %robust
                ANOVAmethod=1;
                trimOption=1;
            case 2 %conventional
                ANOVAmethod=0;
                trimOption=0;
            otherwise
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                error('oops programmer error.')
        end
        
        for iGroup=1:6
            if xor(isempty(EPmain.anova.data.betweenName{iGroup}),(EPmain.anova.data.between(iGroup)>length(EPmain.anova.data.columnNames)))
                msg{1}='All between-group independent variables need to have a three letter label specified.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            
            if ~isempty(EPmain.anova.data.factor{iGroup}) && isempty(EPmain.anova.data.levels{iGroup})
                msg{1}='There is a within-group factor with a name but no levels.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
        end
        
        %determine the within group factors and which are electrode factors
        factorNames=cell(0);
        levelNames=cell(0);
        for i=1:6
            if ~isempty(EPmain.anova.data.factor{i})
                factorNames{end+1}=EPmain.anova.data.factor{i};
                levelNames{end+1}=EPmain.anova.data.levels{i};
            end
        end
        
        repFactor=1;        
        if (length(EPmain.anova.data.name)<10) || ~strcmp(EPmain.anova.data.name(1:10),'behavioral')
            elecFactors=ones(length(factorNames),1);
            sessFactors=ones(length(factorNames),1);
            for theFactor =length(factorNames):-1:1
                factorArea=cell(length(levelNames{theFactor}),1);
                factorSess=cell(length(levelNames{theFactor}),1);
                lvlCount=1;
                repCount=0;
                for theCol=EPmain.anova.data.leftColumn:EPmain.anova.data.rightColumn
                    repCount=repCount+1;
                    if repCount > repFactor
                        repCount=1;
                        lvlCount=lvlCount+1;
                        if lvlCount > length(levelNames{theFactor})
                            lvlCount=1;
                        end
                    end
                    factorArea{lvlCount}{end+1}=EPmain.anova.data.areaNames{theCol};
                    if theCol<=length(EPmain.anova.data.sessNames)
                        factorSess{lvlCount}{end+1}=EPmain.anova.data.sessNames{theCol};
                    end
                end
                repFactor=repFactor*length(levelNames{theFactor});
                for theLevel=1:length(levelNames{theFactor})
                    if length(unique(factorArea{theLevel})) > 1
                        elecFactors(theFactor)=0; %not an electrode factor
                    end
                    if isempty(factorSess{theLevel}) || (length(unique(factorSess{theLevel})) > 1)
                        sessFactors(theFactor)=0; %not a sess factor
                    end
                end
            end
        else
            elecFactors=zeros(length(factorNames),1);
            sessFactors=zeros(length(factorNames),1);
        end
        
        if isempty(factorNames)
            factorNames{1}='   ';
        end
        if isempty(levelNames)
            levelNames{1}=' ';
        end
        
        %Select input and output files
        
        [outFileName, pathname] = uiputfile('*.html*','ANOVA output:');
        if outFileName == 0
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            ep('start');
            return %user hit cancel on file requestor
        end
        outFileName=[pathname outFileName];
        if exist(outFileName,'file')
            delete(outFileName); %user must have clicked "yes" to whether to replace existing file
        end
        [outPathstr, outFileName, outExt] = fileparts(outFileName);
        if ~strcmp(outExt,'.html')
            outExt='.html';
        end
        
        outfid=fopen([outPathstr filesep outFileName outExt],'w');
        
        [ANOVAfiles, pathname] = uigetfile({'*.txt;*.csv'},'Open:','MultiSelect','on');
        activeDirectory=pathname;
        if ~iscell(ANOVAfiles)
            tempVar=ANOVAfiles;
            ANOVAfiles=[];
            ANOVAfiles{1}=tempVar;
        end
        if ANOVAfiles{1}==0
            msg{1}='No filenames selected. You have to click on a name';
            [msg]=ep_errorMsg(msg);
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
        if ~iscell(ANOVAfiles)
            tempVar=ANOVAfiles;
            ANOVAfiles=[];
            ANOVAfiles{1}=tempVar;
        end
        for theFile=1:size(ANOVAfiles,2)
            ANOVAfiles{theFile}=[activeDirectory ANOVAfiles{theFile}];
        end
        
        ANOVAfiles=sort(ANOVAfiles);
        
        %start processing each ANOVA file
        
        fprintf(outfid,'<font face="Courier">');
        fprintf(outfid,'<font size="2">');
        disp('Commencing ANOVA run.  This may take some time.  You may monitor its progress by opening the html output file with a browser and periodically reloading it.');
        ep_tictoc('begin');

        lastName='';
        cellTables=cell(0,4); %the cell array containing all the cell tables (SD, SE, CI, N)
        plotTables=cell(0,3);
        tableCounter=0;
        for theFile=1:length(ANOVAfiles)
            ep_tictoc;if EPtictoc.stop;ep('start');return;end
            if theFile==length(ANOVAfiles)
                lastFlag=1;
            else
                lastFlag=0;
            end
            [pathstr, fileName, ext] = fileparts(ANOVAfiles{theFile});
            [ANOVAheader, ANOVAdata2] = ep_loadANOVA(ANOVAfiles{theFile},EPmain.anova.data);
            if isempty(ANOVAdata2)
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            if (theFile > 1) && isfield(ANOVAdata,'EPdata')
                ANOVAdata2.EPdata=ANOVAdata.EPdata;
            end
            ANOVAdata=ANOVAdata2;
            
            fprintf(outfid,'%s </BR>',fileName);
            disp(fileName);
            
            numCols=size(ANOVAdata.data,2);
            numSpecs=length(find(strcmp('spec',ANOVAdata.areaNames)));
            
            if EPmain.anova.data.rightColumn > numCols
                disp(['Mismatch between number of columns of data specified for analysis ( ' num2str(EPmain.anova.data.rightColumn) ' ) and number in the data file ( ' num2str(numCols) ' ).']);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            
%             ANOVAareas=cell(numCols,1);
%             oldAreas=cell(numCols,1);
%             ANOVAareaList=unique(ANOVAdata.areaNames,'stable');
%             oldList=unique(EPmain.anova.data.areaNames,'stable');
%             for iCol=1:numCols
%                 ANOVAareas{iCol}=find(strcmp(ANOVAdata.areaNames{iCol},ANOVAareaList));
%                 if iCol<=length(EPmain.anova.data.areaNames)
%                     oldAreas{iCol}=find(strcmp(EPmain.anova.data.areaNames{iCol},oldList));
%                 end
%             end
            
%             if ~isequal(EPmain.anova.data.cellNames,ANOVAdata.cellNames) || ~isequal(ANOVAareas,oldAreas)
%             %if (any(~ismember(EPmain.anova.data.cellNames,ANOVAcellNames)) && isempty(findstr('autoPCA',ANOVAdata.changrp))) || (length(EPmain.anova.data.columnNames) ~= numCols)
%                 msg{1}=['The column names for ' fileName ' were different from the ANOVA structure that was set up.'];
%                 [msg]=ep_errorMsg(msg);
%                 set(EPmain.handles.anova.contrast,'enable','on');
%                 set(EPmain.handles.anova.load,'enable','on');
%                 set(EPmain.handles.anova.view,'enable','on');
%                 set(EPmain.handles.anova.run,'enable','on');
%                 set(EPmain.handles.anova.done,'enable','on');
%                 return;
%             end


            numDF=length(levelNames{1})-1;
            for i=2:length(levelNames)
                numDF=numDF*(length(levelNames{i})-1);
            end
            denDF=ANOVAdata.subjects-2*floor(EPmain.preferences.anova.trimming*ANOVAdata.subjects);
            
            if numDF > denDF
                msg{1}=['Product of (number of levels-1) of all within factors (' num2str(numDF) ') cannot exceed the number of participants minus trimming (' num2str(denDF) ').'];
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return;
            end
                
            factorGroupNames=cell(0);
            levelGroupNames=cell(0);
            ANOVAinfo.covariate=[];
            ANOVAinfo.betweenLvl=cell(0);
            for iFactor=1:6
                if ~isempty(EPmain.anova.data.betweenName{iFactor})
                    factorGroupNames{end+1}=EPmain.anova.data.betweenName{iFactor};
                    theNames=EPmain.anova.data.betweenLvl(:,EPmain.anova.data.between(iFactor)-numCols);
                    firstLetter=cell(length(theNames),1);
                    for iName=1:length(firstLetter)
                        firstLetter{iName}=theNames{iName}(1); %first letter only
                    end
                    levelGroupNames{end+1}=cell2mat(unique(firstLetter)');
                    ANOVAinfo.covariate(end+1,1)=EPmain.anova.data.covariate(iFactor);
                    ANOVAinfo.betweenLvl(:,end+1)=theNames;
                end
            end
            ANOVAinfo.leftColumn=EPmain.anova.data.leftColumn;
            ANOVAinfo.rightColumn=EPmain.anova.data.rightColumn;
            
            numComps=EPmain.anova.numComps;
            if isempty(numComps) || (numComps==0)
                numComps=1;
            end
            alpha.uncorrected=.05; %threshold for declaring statistical significance.
%             alpha.corrected=alpha.uncorrected/numComps; %with Bonferroni correction
            alpha.corrected=1-(1-alpha.uncorrected)^(1/numComps); %with Dunn-Šidák correction
            ANOVAheader{end+1}=numComps;
            Y=ANOVAdata.data(:,EPmain.anova.data.leftColumn:EPmain.anova.data.rightColumn);
            ep_tictoc;if EPtictoc.stop;ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);ep('start');return;end
            [cellTables, plotTables, tableCounter, GAVsubsTable] = ep_ADF(Y, ANOVAdata, trimOption, EPmain.preferences.anova, ANOVAmethod, factorNames, levelNames, factorGroupNames, levelGroupNames, elecFactors, sessFactors, alpha, 0, 1, 1, 0, 0, outfid,[],[], ANOVAheader, EPmain.anova.allPosthoc, ANOVAinfo, cellTables, plotTables, tableCounter, lastFlag);
            if EPtictoc.stop
                EPtictoc.stop=0;
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                break
            end
            %if there are between-group factors, then add grand averages corresponding to their levels to the original data.
            if EPmain.preferences.anova.adds && ~strcmp(ANOVAdata.name,'behavioral') && ~any(strcmp(ANOVAdata.measure,{'meanStdVar','meanSME','meanNoise'})) && ~isempty(levelGroupNames)
                
                %search through active datasets to see if one has a name matching that in the ANOVA file
                if ~isempty(EPdataset.dataset)
                    ANOVAname=ANOVAdata.name;
                    if strfind(ANOVAdata.name,'Dataset Name: ')
                        ANOVAname=ANOVAname(15:end);
                    end
                    if ~isempty(ANOVAname)
                        whichData=find(strcmp(ANOVAname,{EPdataset.dataset.dataName}));
                        if isempty(whichData)
                            disp('The ANOVA dataset is not present in the working set so corresponding grand averages will not be added.');
                        else
                            if ~strcmp(lastName,ANOVAname)
                                lastName=ANOVAname;
                                ANOVAdata.EPdata=ep_loadEPdataset(whichData);
                            end
                            ANOVAAdds(ANOVAdata, levelGroupNames, whichData, GAVsubsTable);
                        end
                    end
                end
            end
        end
        fclose(outfid);
        disp('Finished ANOVA run.');
        ep_tictoc('end');
        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        ep('start');
        
    case 'startSave'
        %Set up save pane of main window.
        
        set(EPmain.handles.hMainWindow,'Name', 'Save Data');
        
        EPmain.handles.save.prefs = uicontrol('Style', 'pushbutton', 'String', '•','FontSize',EPmain.fontsize,...
            'Position', [5 485 10 10], 'Callback', ['global EPmain;','EPmain.prefReturn=''save'';','EPmain.mode=''preferenceGeneral'';','ep(''start'');']);
        
        uicontrol('Style','text',...
            'String','Save File Format','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[25 450 100 20]);
        
        EPmain.handles.save.format = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatSaveList,...
            'Value',EPmain.save.format,'Position',[20 420 150 20],...
            'Callback', ['global EPmain;','tempVar=get(EPmain.handles.save.format,''Value'');','if tempVar ~=0,EPmain.save.format=tempVar;end;','if isempty(tempVar),EPmain.save.format=tempVar;end;','ep(''start'');']);
        
        uicontrol('Style','text',...
            'String','','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[20 390 70 20]);
        
        uicontrol('Style','text',...
            'String','Single','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[80 390 50 20]);

        uicontrol('Style','text',...
            'String','Combined','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[120 390 60 20]);
        
        uicontrol('Style','text',...
            'String','Channels','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[20 370 70 20]);
        
        EPmain.handles.save.SGLchan= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.SGLchan,'Position',[80 370 50 20],...
            'Callback', ['global EPmain;','EPmain.save.SGLchan=get(EPmain.handles.save.SGLchan,''Value'');'],'TooltipString','Single channels.');
        
        EPmain.handles.save.REGchan= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.REGchan,'Position',[120 370 50 20],...
            'Callback', ['global EPmain;','EPmain.save.REGchan=get(EPmain.handles.save.REGchan,''Value'');'],'TooltipString','Regional channels.');
                
        uicontrol('Style','text',...
            'String','Cells','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[20 350 70 20]);
        
        EPmain.handles.save.SGLcell= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.SGLcell,'Position',[80 350 50 20],...
            'Callback', ['global EPmain;','EPmain.save.SGLcell=get(EPmain.handles.save.SGLcell,''Value'');'],'TooltipString','Single cells.');
        
        EPmain.handles.save.CMBcell= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.CMBcell,'Position',[120 350 50 20],...
            'Callback', ['global EPmain;','EPmain.save.CMBcell=get(EPmain.handles.save.CMBcell,''Value'');'],'TooltipString','Combination cells.');
        
        uicontrol('Style','text',...
            'String','Trials','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[20 330 70 20]);
        
        EPmain.handles.save.RAW= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.RAW,'Position',[80 330 90 20],'FontSize',EPmain.fontsize,...
            'Callback', ['global EPmain;','EPmain.save.RAW=get(EPmain.handles.save.RAW,''Value'');'],'TooltipString','Single trial data.');
        
        uicontrol('Style','text',...
            'String','Averages','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[20 310 70 20]);
        
        EPmain.handles.save.AVG= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.AVG,'Position',[80 310 50 20],...
            'Callback', ['global EPmain;','EPmain.save.AVG=get(EPmain.handles.save.AVG,''Value'');'],'TooltipString','Subject averages.');
        
        EPmain.handles.save.GAV= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.GAV,'Position',[120 310 50 20],...
            'Callback', ['global EPmain;','EPmain.save.GAV=get(EPmain.handles.save.GAV,''Value'');'],'TooltipString','Grand averages.');
        
        uicontrol('Style','text',...
            'String','Factors','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'Position',[20 290 70 20]);
        
        EPmain.handles.save.SGLfac= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.SGLfac,'Position',[80 290 50 20],...
            'Callback', ['global EPmain;','EPmain.save.SGLfac=get(EPmain.handles.save.SGLfac,''Value'');'],'TooltipString','Single Factors (including two-step factors).');
        
        EPmain.handles.save.CMBfac= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.CMBfac,'Position',[120 290 50 20],...
            'Callback', ['global EPmain;','EPmain.save.CMBfac=get(EPmain.handles.save.CMBfac,''Value'');'],'TooltipString','Combined Factors (grand factors, not two-step factors).');
        
        EPmain.handles.save.CMBfac= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.save.CMBfac,'Position',[120 290 50 20],...
            'Callback', ['global EPmain;','EPmain.save.CMBfac=get(EPmain.handles.save.CMBfac,''Value'');'],'TooltipString','Combined Factors (grand factors, not two-step factors).');
         
        
        EPmain.handles.save.batch = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Single','Convert'},...
            'Value',EPmain.save.batch,'Position',[20 250 150 20],...
            'Callback', ['global EPmain;','tempVar=get(EPmain.handles.save.batch,''Value'');','if tempVar ~=0,EPmain.save.batch=tempVar;end;','if isempty(tempVar),EPmain.save.batch=tempVar;end;','ep(''start'');']);
       
        if EPmain.save.batch == 1
            
            if ~isempty(EPdataset.dataset)
                for iFile=1:length(EPdataset.dataset)
                    fileName=EPdataset.dataset(iFile).dataName;
                    if strcmp(EPdataset.dataset(iFile).saved,'no')
                        fileName=['*' fileName];
                    end
                    tableData{iFile,1}=false;
                    tableData{iFile,2}=fileName;
                end
            else
                tableData=[];
            end
            
            tableNames={'';''};
            columnEditable =  [true,false];
            ColumnFormat={[],[]};
            
            EPmain.handles.save.hTable = uitable('Data',tableData,'ColumnName',tableNames,'FontSize',EPmain.fontsize,...
                'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,...
                'CellSelectionCallback',{@saveData,[],0},...
                'ColumnWidth',{20,280},'Position',[20 60 170 150]);
            
        else
            uicontrol('Style','text',...
                'String','Read File Format','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[25 230 100 20]);
            
            EPmain.handles.save.readFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',EPmain.fileFormatReadList,...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.save.readFormat,''Value'');','if tempVar ~=0,EPmain.save.readFormat=tempVar;end;','if isempty(tempVar),EPmain.save.readFormat=tempVar;end;','ep(''start'');'],...
                'Value',EPmain.save.readFormat,'Position',[20 210 150 20]);
            
            uicontrol('Style','text',...
                'String','File Type','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
                'Position',[25 190 100 20]);
            
            EPmain.handles.save.fileType= uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                'String',{'continuous','single_trial','average','grand_average','factors'},...
                'CallBack',['global EPmain;','tempVar=get(EPmain.handles.save.fileType,''Value'');','if tempVar ~=0,EPmain.save.fileType=tempVar;end;','if isempty(tempVar),EPmain.save.fileType=tempVar;end;','ep(''start'');'],...
                'Value',EPmain.save.fileType,'Position',[20 170 150 20]);
            
            if EPmain.preferences.general.BVheader && (strcmp(EPmain.fileFormatReadList{EPmain.save.readFormat},'BrainVision (.eeg/.dat/.seg)')) && (EPmain.save.fileType==1)
                EPmain.handles.save.EPheaderFix= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
                    'Value',EPmain.save.EPheaderFix,'Position',[20 270 160 20],'String','Fix damaged EP header',...
                    'Callback', ['global EPmain;','EPmain.save.EPheaderFix=get(EPmain.handles.save.EPheaderFix,''Value'');'],'TooltipString','Fix partly missing EP header if needed.');
            end
            
            if (strcmp(EPmain.fileFormatReadList{EPmain.save.readFormat},'EGI Simple Binary (.sbin)')) && (EPmain.save.fileType==1)
                EPmain.handles.save.NS4fix= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
                    'Value',EPmain.save.NS4fix,'Position',[20 270 150 20],'String','Salvage NS4 files',...
                    'Callback', ['global EPmain;','EPmain.save.NS4fix=get(EPmain.handles.save.NS4fix,''Value'');'],'TooltipString','Reconstitute NS4 files.  See tutorial.');
            else
                EPmain.save.NS4fix=0;
            end
            
            [importSuffix,importFormatName,importFormat]=ep_fileFormats('average',EPmain.fileFormatReadList{EPmain.save.readFormat});
            if strcmp(importFormat,'ep_mat')
                set(EPmain.handles.save.fileType,'enable','off');
            end
            
            uicontrol('Style','frame',...
                'Position',[5 80 195 90]);
            
            EPmain.handles.save.check= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
                'String','Single File Mode',...
                'CallBack',['global EPmain;','EPmain.save.check=get(EPmain.handles.save.check,''Value'');','ep(''start'');'],'FontSize',EPmain.fontsize,...
                'Value',EPmain.save.check,'Position',[10 145 150 20]);
            
            if EPmain.save.check
                
                EPmain.handles.save.subject= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                    'String',EPmain.save.subject,...
                    'CallBack',['global EPmain;','EPmain.save.subject=get(EPmain.handles.save.subject,''String'');','ep(''start'');'],...
                    'Position',[10 125 50 20],'TooltipString','example 4:6');
                
                EPmain.handles.save.subjectLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Sub','HorizontalAlignment','left',...
                    'Position',[70 125 50 20]);
                
                if isempty(EPmain.save.subject)
                    set(EPmain.handles.save.subjectLabel,'enable','off');
                elseif isempty(str2num(EPmain.save.subject))
                    set(EPmain.handles.save.subjectLabel,'enable','off');
                end
                
                EPmain.handles.save.cell= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                    'String',EPmain.save.cell,...
                    'CallBack',['global EPmain;','EPmain.save.cell=get(EPmain.handles.save.cell,''String'');','ep(''start'');'],...
                    'Position',[10 105 50 20],'TooltipString','example 7:9');
                
                EPmain.handles.save.cellLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Cell','HorizontalAlignment','left',...
                    'Position',[70 105 50 20]);
                
                if isempty(EPmain.save.cell)
                    set(EPmain.handles.save.cellLabel,'enable','off');
                elseif isempty(str2num(EPmain.save.cell))
                    set(EPmain.handles.save.cellLabel,'enable','off');
                end
                
                EPmain.handles.save.freq= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                    'String',EPmain.save.freq,...
                    'CallBack',['global EPmain;','EPmain.save.freq=get(EPmain.handles.save.freq,''String'');','ep(''start'');'],...
                    'Position',[10 85 50 20],'TooltipString','example 10:12');
                
                EPmain.handles.save.freqLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Freq','HorizontalAlignment','left',...
                    'Position',[70 85 50 20]);
                
                if isempty(EPmain.save.freq)
                    set(EPmain.handles.save.freqLabel,'enable','off');
                elseif isempty(str2num(EPmain.save.freq))
                    set(EPmain.handles.save.freqLabel,'enable','off');
                end
                
                if (EPmain.save.fileType==2) || strcmp(importFormat,'ep_mat')
                    EPmain.handles.save.trial= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                        'String',EPmain.save.trial,...
                        'CallBack',['global EPmain;','EPmain.save.trial=get(EPmain.handles.save.trial,''String'');','ep(''start'');'],...
                        'Position',[110 125 50 20],'TooltipString','example 10:12');
                    
                    EPmain.handles.save.trialLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                        'String','Trial','HorizontalAlignment','left',...
                        'Position',[160 125 30 20]);
                    
                    if isempty(EPmain.save.trial)
                        set(EPmain.handles.save.trialLabel,'enable','off');
                    elseif isempty(str2double(EPmain.save.trial))
                        set(EPmain.handles.save.trialLabel,'enable','off');
                    end
                end
                
                EPmain.handles.save.session= uicontrol('Style','edit','FontSize',EPmain.fontsize,...
                    'String',EPmain.save.session,...
                    'CallBack',['global EPmain;','EPmain.save.session=get(EPmain.handles.save.session,''String'');','ep(''start'');'],...
                    'Position',[110 105 50 20],'TooltipString','example 7:9');
                
                EPmain.handles.save.sessionLabel=uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','Sess','HorizontalAlignment','left',...
                    'Position',[160 105 35 20]);
                
                if isempty(EPmain.save.session)
                    set(EPmain.handles.save.sessionLabel,'enable','off');
                elseif isempty(str2double(EPmain.save.session))
                    set(EPmain.handles.save.sessionLabel,'enable','off');
                end
            end
            
            EPmain.handles.save.convert = uicontrol('Style', 'pushbutton', 'String', 'Convert','FontSize',EPmain.fontsize,...
                'Position', [20 40 80 40], 'Callback', @convertFiles);
        end
        
        EPmain.handles.save.done = uicontrol('Style', 'pushbutton', 'String', 'Main','FontSize',EPmain.fontsize,...
            'Position', [20 0 80 40], 'Callback', ['global EPmain;','EPmain.mode=''main'';','EPmain.save.EPheaderFix=0;','EPmain.save.NS4fix=0;','ep(''start'');']);
        
        EPmain.handles.save.save = uicontrol('Style', 'pushbutton', 'String', 'Save','FontSize',EPmain.fontsize,...
            'Position', [120 0 80 40], 'Callback', @saveSelectData);

        if EPmain.save.batch || ~any([tableData{:,1}])
            set(EPmain.handles.save.save,'enable','off');
        end
        
    case 'startPreferenceMain'
        
        set(EPmain.handles.hMainWindow,'Name', 'Preferences');
        
        EPmain.handles.preferences.hGeneral = uicontrol('Style', 'pushbutton', 'String', 'Files',...
            'Position', [20 450 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferenceGeneral'';','ep(''start'');']);
        
        EPmain.handles.preferences.hPreprocess = uicontrol('Style', 'pushbutton', 'String', 'Preprocess',...
            'Position', [20 420 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferencePreprocess'';','ep(''start'');']);
        
        EPmain.handles.preferences.hAverage = uicontrol('Style', 'pushbutton', 'String', 'Average',...
            'Position', [20 390 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferenceAverage'';','ep(''start'');']);
        
        EPmain.handles.preferences.hTransform = uicontrol('Style', 'pushbutton', 'String', 'Transform',...
            'Position', [20 360 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferenceTransform'';','ep(''start'');']);
        
        EPmain.handles.preferences.hView = uicontrol('Style', 'pushbutton', 'String', 'View',...
            'Position', [20 330 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferenceView'';','ep(''start'');']);
        
        EPmain.handles.preferences.hPCA = uicontrol('Style', 'pushbutton', 'String', 'PCA',...
            'Position', [20 300 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferencePCA'';','ep(''start'');']);
        
        EPmain.handles.preferences.hWindow = uicontrol('Style', 'pushbutton', 'String', 'Window',...
            'Position', [20 270 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferenceWindow'';','ep(''start'');']);
        
        EPmain.handles.preferences.hANOVA = uicontrol('Style', 'pushbutton', 'String', 'ANOVA',...
            'Position', [20 240 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferenceANOVA'';','ep(''start'');']);
        
        EPmain.handles.preferences.advanced = uicontrol('Style', 'pushbutton', 'String', 'Addvanced',...
            'Position', [20 210 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferenceAdvanced'';','ep(''start'');']);
        
        EPmain.handles.preferences.records = uicontrol('Style', 'pushbutton', 'String', 'Records',...
            'Position', [20 180 100 30], 'Callback', ['global EPmain;','EPmain.mode=''preferenceRecords'';','ep(''start'');']);
        
        EPmain.handles.preferences.done = uicontrol('Style', 'pushbutton', 'String', 'Main',...
            'Position', [20 100 100 40], 'Callback', ['global EPmain;','EPmain.mode=''main'';','ep(''start'');']);
        
        EPmain.handles.preferences.reset = uicontrol('Style', 'pushbutton', 'String', 'Reset',...
            'Position', [20 50 100 40], 'Callback', @resetPrefs);
        
        EPmain.handles.preferences.save = uicontrol('Style', 'pushbutton', 'String', 'Save',...
            'Position', [20 0 100 40], 'Callback', 'ep(''savePrefs'');');
        
        EPmain.tempColor=EPmain.preferences.view.color;
        EPmain.prefReturn='preferenceMain';
        
    case 'startPreferenceGeneral'
        
        set(EPmain.handles.hMainWindow,'Name', 'Files Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Session Import Format','FontSize',EPmain.fontsize,...
            'Position',[20 470 150 20]);
        
        EPmain.handles.preferences.general.sessionImportFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatReadList,...
            'Value',EPmain.preferences.general.sessionImportFormat,'Position',[20 450 150 20],'TooltipString','Default file format for reading session files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Session Ouput Format','FontSize',EPmain.fontsize,...
            'Position',[20 430 150 20]);
        
        EPmain.handles.preferences.general.sessionOutputFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatSaveList,...
            'Value',EPmain.preferences.general.sessionOutputFormat,'Position',[20 410 150 20],'TooltipString','Default file format for saving session files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'File Import Format','FontSize',EPmain.fontsize,...
            'Position',[20 390 150 20]);
        
        EPmain.handles.preferences.general.importFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatReadList,...
            'Value',EPmain.preferences.general.importFormat,'Position',[20 370 150 20],'TooltipString','Default file format for reading average files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'File Ouput Format','FontSize',EPmain.fontsize,...
            'Position',[20 350 150 20]);
        
        EPmain.handles.preferences.general.outputFormat = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',EPmain.fileFormatSaveList,...
            'Value',EPmain.preferences.general.outputFormat,'Position',[20 330 150 20],'TooltipString','Default file format for saving average files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'First Text Row','FontSize',EPmain.fontsize,...
            'Position',[20 310 90 20]);
        
        EPmain.handles.preferences.general.firstRow= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.general.firstRow),'FontSize',EPmain.fontsize,...
            'Position',[120 310 70 20],'TooltipString','First row of data (as opposed to header rows) when reading in text files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Last Text Row','FontSize',EPmain.fontsize,...
            'Position',[20 290 90 20]);
        
        EPmain.handles.preferences.general.lastRow= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.general.lastRow),'FontSize',EPmain.fontsize,...
            'Position',[120 290 70 20],'TooltipString','Last row of data when reading in text files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'First Text Col','FontSize',EPmain.fontsize,...
            'Position',[20 270 90 20]);
        
        EPmain.handles.preferences.general.firstCol= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.general.firstCol),'FontSize',EPmain.fontsize,...
            'Position',[120 270 70 20],'TooltipString','First column of data when reading in text files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Last Text Col','FontSize',EPmain.fontsize,...
            'Position',[20 250 150 20]);
        
        EPmain.handles.preferences.general.lastCol= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.general.lastCol),'FontSize',EPmain.fontsize,...
            'Position',[120 250 70 20],'TooltipString','Last column of data when reading in text files (0=last).');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Sample rate','FontSize',EPmain.fontsize,...
            'Position',[20 230 150 20]);
        
        EPmain.handles.preferences.general.sampleRate= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.general.sampleRate),'FontSize',EPmain.fontsize,...
            'Position',[120 230 70 20],'TooltipString','Sample rate when reading in text files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'segment suffix','FontSize',EPmain.fontsize,...
            'Position',[20 210 90 20]);
        
        EPmain.handles.preferences.general.segSuffix= uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.preferences.general.segSuffix,'FontSize',EPmain.fontsize,...
            'Position',[120 210 70 20],'TooltipString','Suffix for output segment files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'segment suffix','FontSize',EPmain.fontsize,...
            'Position',[20 190 90 20]);
        
        EPmain.handles.preferences.general.specSuffix= uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.preferences.general.specSuffix,'FontSize',EPmain.fontsize,...
            'Position',[120 190 70 20],'TooltipString','Suffix for event data file to merge into file being read.  Should have both an underscore and a dot suffix (e.g., _evt.txt)');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Sub spec suffix','FontSize',EPmain.fontsize,...
            'Position',[20 170 90 20]);
        
        EPmain.handles.preferences.general.subjectSpecSuffix= uicontrol('Style','edit','HorizontalAlignment','left','String', EPmain.preferences.general.subjectSpecSuffix,'FontSize',EPmain.fontsize,...
            'Position',[120 170 70 20],'TooltipString','Suffix for subject spec text file to merge into file being read.  Should have both an underscore and a dot suffix (e.g., _sub.txt)');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Text File Orientation','FontSize',EPmain.fontsize,...
            'Position',[20 150 150 20]);
        
        EPmain.handles.preferences.general.orientation = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Chan Cols','Chan Rows'},...
            'Value',EPmain.preferences.general.orientation,'Position',[20 130 150 20],'TooltipString','When reading in text files, are channels the rows or the columns?');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Montage','FontSize',EPmain.fontsize,...
            'Position',[5 90 50 20]);
        
        EPmain.handles.preferences.general.defaultMontage = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,'String',EPmain.montageList,...
            'Value',find(strcmp(EPmain.preferences.general.defaultMontage,EPmain.montageList)),'Position',[60 90 150 20],'TooltipString','Default montage setting when importing data files.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'BV header','FontSize',EPmain.fontsize,...
            'Position',[20 70 150 20]);
        
        EPmain.handles.preferences.general.BVheader = uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'Value',EPmain.preferences.general.BVheader,'Position',[180 70 150 20],'TooltipString','When loading BrainVision files, interpret event codes as being encoded EP Toolkit header and TRSP information?');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'In Memory','FontSize',EPmain.fontsize,...
            'Position',[20 50 150 20]);
        
        EPmain.handles.preferences.general.numEEG = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.general.numEEG,'Position',[120 50 70 20],'TooltipString','Number of most recently accessed working set to keep in RAM to save on loading time.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'No internal CED','FontSize',EPmain.fontsize,...
            'Position',[20 30 100 20]);
        
        EPmain.handles.preferences.general.noInternal=uicontrol('Style','checkbox','HorizontalAlignment','left','value', EPmain.preferences.general.noInternal,'FontSize',EPmain.fontsize,...
            'Position',[130 30 20 20],'TooltipString','No Internal CED means ignoring the electrode coordinates stored in mff or set files. Instead, the EP Toolkit will ask for a CED file.');
        
        EPmain.handles.preferences.general.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [10 0 100 30], 'Callback', ['global EPmain;','EPmain.mode=EPmain.prefReturn;','ep(''start'');']);
        
        EPmain.handles.preferences.general.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [110 0 100 30], 'Callback', 'ep(''getGeneralPrefs'');');
                
    case 'getGeneralPrefs' %retrieve preference settings from the general preference input fields
        
        tempVar=EPmain.preferences.general;
        
        EPmain.preferences.general.sessionImportFormat=get(EPmain.handles.preferences.general.sessionImportFormat,'Value');
        
        EPmain.preferences.general.sessionOutputFormat=get(EPmain.handles.preferences.general.sessionOutputFormat,'Value');
        
        EPmain.preferences.general.importFormat=get(EPmain.handles.preferences.general.importFormat,'Value');
        
        EPmain.preferences.general.outputFormat=get(EPmain.handles.preferences.general.outputFormat,'Value');
        
        EPmain.preferences.general.firstRow=str2num(get(EPmain.handles.preferences.general.firstRow,'String'));
        
        EPmain.preferences.general.lastRow=str2num(get(EPmain.handles.preferences.general.lastRow,'String'));

        EPmain.preferences.general.firstCol=str2num(get(EPmain.handles.preferences.general.firstCol,'String'));
        
        EPmain.preferences.general.lastCol=str2num(get(EPmain.handles.preferences.general.lastCol,'String'));
        
        EPmain.preferences.general.sampleRate=str2num(get(EPmain.handles.preferences.general.sampleRate,'String'));

        EPmain.preferences.general.segSuffix=get(EPmain.handles.preferences.general.segSuffix,'String');
        
        EPmain.preferences.general.specSuffix=get(EPmain.handles.preferences.general.specSuffix,'String');
        
        EPmain.preferences.general.subjectSpecSuffix=get(EPmain.handles.preferences.general.subjectSpecSuffix,'String');

        EPmain.preferences.general.orientation=get(EPmain.handles.preferences.general.orientation,'Value');
        
        EPmain.preferences.general.defaultMontage=EPmain.montageList{get(EPmain.handles.preferences.general.defaultMontage,'Value')};
        
        EPmain.preferences.general.BVheader=get(EPmain.handles.preferences.general.BVheader,'Value');
        
        EPmain.preferences.general.numEEG=str2num(get(EPmain.handles.preferences.general.numEEG,'String'));

        EPmain.preferences.general.noInternal=get(EPmain.handles.preferences.general.noInternal,'Value');
        
        [prefErr prefMissing]=checkPrefs;
        
        if prefErr
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.general=tempVar;
        else
            EPmain.mode=EPmain.prefReturn;
            
            EPmain.average.importFormat=EPmain.preferences.general.sessionImportFormat;
            EPmain.average.outputFormat=EPmain.preferences.general.outputFormat;
            EPmain.transform.importFormat=EPmain.preferences.general.importFormat;
            EPmain.transform.outputFormat=EPmain.preferences.general.outputFormat;
            
            ep('start');
        end
        
    case 'startPreferencePreprocess'
        
        set(EPmain.handles.hMainWindow,'Name', 'Preprocess Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'No Figure','FontSize',EPmain.fontsize,...
            'Position',[20 470 100 20]);
        
        EPmain.handles.preferences.preprocess.noFigure= uicontrol('Style','checkbox','HorizontalAlignment','left','value', EPmain.preferences.preprocess.noFigure,'FontSize',EPmain.fontsize,...
            'Position',[130 470 20 20],'TooltipString','No Figure is an option to not provide a summary figure for the artifact correction process.  While a very useful figure, it requires substantial memory and so dropping it can be helpful when encountering recalcitrant memory problems. ');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Moving Window','FontSize',EPmain.fontsize,...
            'Position',[20 450 100 20]);
        
        EPmain.handles.preferences.preprocess.window= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.window),'FontSize',EPmain.fontsize,...
            'Position',[130 450 70 20],'TooltipString','Moving Window is the number of milliseconds over which the artifact correction routines average the data in a form of low pass filtering.  The larger the number, the less sensitive it is to high frequency spikes.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Chan Min-Max','FontSize',EPmain.fontsize,...
            'Position',[20 430 100 20]);
        
        EPmain.handles.preferences.preprocess.minmax= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.minmax),'FontSize',EPmain.fontsize,...
            'Position',[130 430 70 20],'TooltipString',['Chan Min-Max ' char(181) 'v is the maximum allowed change in voltage levels for a channel during a trial before it is deemed to be a bad channel for that trial.']);
        
        uicontrol('Style','text','HorizontalAlignment','left','String', '% Bad Channel','FontSize',EPmain.fontsize,...
            'Position',[20 410 100 20]);
        
        EPmain.handles.preferences.preprocess.badnum= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.badnum),'FontSize',EPmain.fontsize,...
            'Position',[130 410 70 20],'TooltipString','% Bad Channel is the maximum number of channels allowed to be bad in a trial before it is deemed to be a bad trial.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', '# Neigh Chans','FontSize',EPmain.fontsize,...
            'Position',[20 390 100 20]);
        
        EPmain.handles.preferences.preprocess.neighbors= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.neighbors),'FontSize',EPmain.fontsize,...
            'Position',[130 390 70 20],'TooltipString','# Neigh Chans is the number of channels considered to be a neighbor for purposes of the artifact correction algorithms.  The electrode coordinates are then used to figure out which channels are to be used.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', ['Neigh Diff ' char(181) 'v'],'FontSize',EPmain.fontsize,...
            'Position',[20 370 100 20]);
        
        EPmain.handles.preferences.preprocess.maxneighbor= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.maxneighbor),'FontSize',EPmain.fontsize,...
            'Position',[130 370 70 20],'TooltipString',['Neigh Diff ' char(181) 'v is the maximum voltage difference allowed between a channel and its neighbors before it is deemed to be a bad channel.']);
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Bad Chan Corr','FontSize',EPmain.fontsize,...
            'Position',[20 350 100 20]);
        
        EPmain.handles.preferences.preprocess.badchan= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%.3f', EPmain.preferences.preprocess.badchan),'FontSize',EPmain.fontsize,...
            'Position',[130 350 70 20],'TooltipString','Bad Chan Corr is the correlation criterion for determining whether a channel is a globally bad channel over the entire session.  If its best correlation with any neighbor is lower than this criterion then it is judged to be globally bad.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Blink Corr','FontSize',EPmain.fontsize,...
            'Position',[20 330 100 20]);
        
        EPmain.handles.preferences.preprocess.blink= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%.3f', EPmain.preferences.preprocess.blink),'FontSize',EPmain.fontsize,...
            'Position',[130 330 70 20],'TooltipString','Blink Corr is the correlation criterion for determining whether an ICA factor matches the blink template and should therefore be subtracted from the data.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', '% Bad Trial','FontSize',EPmain.fontsize,...
            'Position',[20 310 100 20]);
        
        EPmain.handles.preferences.preprocess.badtrials= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.badtrials),'FontSize',EPmain.fontsize,...
            'Position',[130 310 70 20],'TooltipString','% Bad Trial Chan is the maximum number of trials a channel is allowed to be judged bad before it is deemed to be globally bad.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Size of Chunks','FontSize',EPmain.fontsize,...
            'Position',[20 290 100 20]);
        
        EPmain.handles.preferences.preprocess.chunkSize= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.chunkSize),'FontSize',EPmain.fontsize,...
            'Position',[130 290 70 20],'TooltipString','Size of Chunks is the number of time points that are read into each chunk (about 100,000 per GB of available RAM seems to generally work).');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Warn Trials/Cell','FontSize',EPmain.fontsize,...
            'Position',[20 270 100 20]);
        
        EPmain.handles.preferences.preprocess.minTrialsPerCell= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.minTrialsPerCell),'FontSize',EPmain.fontsize,...
            'Position',[130 270 70 20],'TooltipString','Warn Trials/Cell is the minimum number of good trials that is considered to be sufficient for a cell.  Any cells dropping below this number will trigger a warning in the artifact correction log.  There is no other effect of this setting.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Bad Neighbors','FontSize',EPmain.fontsize,...
            'Position',[20 250 100 20]);
        
        EPmain.handles.preferences.preprocess.noadjacent= uicontrol('Style','checkbox','HorizontalAlignment','left','value', EPmain.preferences.preprocess.noadjacent,'FontSize',EPmain.fontsize,...
            'Position',[130 250 70 20],'TooltipString','Bad Neighbors is an option where if two neighboring channels are marked as being locally bad then the trial is also marked bad (because this typically means that a movement artifact of some sort is present, as opposed to isolated bad channels).');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', [char(181) 'v Move Fac'],'FontSize',EPmain.fontsize,...
            'Position',[20 230 100 20]);
        
        EPmain.handles.preferences.preprocess.trialminmax= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.trialminmax),'FontSize',EPmain.fontsize,...
            'Position',[130 230 70 20],'TooltipString',[char(181) 'v Move Fac is the maximum voltage difference (maximum-minimum) allowed by a factor by the movement artifact correction step.  Factors exceeding this limit are deeemed to reflect artifacts and are subtracted from the data.']);
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Move Corr Facs','FontSize',EPmain.fontsize,...
            'Position',[20 210 100 20]);
        
        EPmain.handles.preferences.preprocess.movefacs= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.movefacs),'FontSize',EPmain.fontsize,...
            'Position',[130 210 70 20],'TooltipString','Move Corr Facs is the number of factors to be retained by the movement artifact correction routine.  A larger number results in a more accurate but slower process.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'EMG Ratio','FontSize',EPmain.fontsize,...
            'Position',[20 190 100 20]);
        
        EPmain.handles.preferences.preprocess.EMGratio= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.EMGratio),'FontSize',EPmain.fontsize,...
            'Position',[130 190 70 20],'TooltipString','The minimum ratio of signal power to EMG noise to retain during EMG correction.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'EMG Threshold','FontSize',EPmain.fontsize,...
            'Position',[20 170 100 20]);
        
        EPmain.handles.preferences.preprocess.EMGthresh= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.EMGthresh),'FontSize',EPmain.fontsize,...
            'Position',[130 170 70 20],'TooltipString','The Hz threshold considered to be the lower bound of possible EMG frquencies during EMG correction.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'EOG channels','FontSize',EPmain.fontsize,...
            'Position',[20 150 100 20]);
        
        EPmain.handles.preferences.preprocess.EOGchans= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d %d %d %d %d %d', EPmain.preferences.preprocess.EOGchans),'FontSize',EPmain.fontsize,...
            'Position',[130 150 70 20],'TooltipString','EOG chans if automatic identification is not working [LUV RUV LLV RLV LH RH].  Enter in -1 for a missing electrode.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Spike Pot','FontSize',EPmain.fontsize,...
            'Position',[20 130 100 20]);
        
        EPmain.handles.preferences.preprocess.sacPot= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%.3f', EPmain.preferences.preprocess.sacPot),'FontSize',EPmain.fontsize,...
            'Position',[130 130 70 20],'TooltipString','Threshold for saccade potential detection.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Saturation','FontSize',EPmain.fontsize,...
            'Position',[20 110 100 20]);
        
        EPmain.handles.preferences.preprocess.saturation= uicontrol('Style','edit','HorizontalAlignment','left','String', sprintf('%d', EPmain.preferences.preprocess.saturation),'FontSize',EPmain.fontsize,...
            'Position',[130 110 70 20],'TooltipString','Saturation is the voltage level at which the amplifier reaches the maximum number that it is capable of reporting.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'fMRI correct','FontSize',EPmain.fontsize,...
            'Position',[20 90 100 20]);
        
        EPmain.handles.preferences.preprocess.fMRI = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'fMRIb OBS','AMRI ICA'},...
            'Value',EPmain.preferences.preprocess.fMRI,'Position',[125 90 100 20],'TooltipString','Algorithm for correcting for fMRI artifacts.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Blink Rotation','FontSize',EPmain.fontsize,...
            'Position',[20 70 100 20]);
        
        rotationList=ep_doPCA;
        EPmain.handles.preferences.preprocess.blinkRotation = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',rotationList,...
            'Value',find(strcmp(EPmain.preferences.preprocess.blinkRotation,rotationList)),'Position',[125 70 100 20],'TooltipString','Rotation for blink correction.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Saccade Method','FontSize',EPmain.fontsize,...
            'Position',[20 50 100 20]);
        
        rotationList(2:end+1)=rotationList;
        rotationList{1}='regression';
        EPmain.handles.preferences.preprocess.saccRotation = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',rotationList,...
            'Value',find(strcmp(EPmain.preferences.preprocess.saccRotation,rotationList)),'Position',[125 50 100 20],'TooltipString','Method for saccade correction.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'SP Rotation','FontSize',EPmain.fontsize,...
            'Position',[20 30 100 20]);
        
        rotationList=ep_doPCA;
        rotationList(2:end+1)=rotationList;
        rotationList{1}='vector';
        EPmain.handles.preferences.preprocess.SProtation = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',rotationList,...
            'Value',find(strcmp(EPmain.preferences.preprocess.SProtation,rotationList)),'Position',[125 30 100 20],'TooltipString','Rotation for SP correction.');
        
        EPmain.handles.preferences.preprocess.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 0 80 30], 'Callback', ['global EPmain;','EPmain.mode=EPmain.prefReturn;','ep(''start'');']);
        
        EPmain.handles.preferences.preprocess.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [100 0 80 30], 'Callback', 'ep(''getPreprocessPrefs'');');
        
    case 'getPreprocessPrefs' %retrieve preference settings from the preprocessing preference input fields
        
        tempVar=EPmain.preferences.preprocess;
        
        EPmain.preferences.preprocess.noFigure=get(EPmain.handles.preferences.preprocess.noFigure,'Value');
        
        EPmain.preferences.preprocess.window=str2num(get(EPmain.handles.preferences.preprocess.window,'String'));
        
        EPmain.preferences.preprocess.minmax=str2num(get(EPmain.handles.preferences.preprocess.minmax,'String'));
        
        EPmain.preferences.preprocess.badnum=str2num(get(EPmain.handles.preferences.preprocess.badnum,'String'));
        
        EPmain.preferences.preprocess.neighbors=str2num(get(EPmain.handles.preferences.preprocess.neighbors,'String'));
        
        EPmain.preferences.preprocess.maxneighbor=str2num(get(EPmain.handles.preferences.preprocess.maxneighbor,'String'));
        
        EPmain.preferences.preprocess.badchan=str2num(get(EPmain.handles.preferences.preprocess.badchan,'String'));
        
        EPmain.preferences.preprocess.blink=str2num(get(EPmain.handles.preferences.preprocess.blink,'String'));
        
        EPmain.preferences.preprocess.badtrials=str2num(get(EPmain.handles.preferences.preprocess.badtrials,'String'));
        
        EPmain.preferences.preprocess.chunkSize=str2num(get(EPmain.handles.preferences.preprocess.chunkSize,'String'));
        
        EPmain.preferences.preprocess.minTrialsPerCell=str2num(get(EPmain.handles.preferences.preprocess.minTrialsPerCell,'String'));
        
        EPmain.preferences.preprocess.noadjacent=get(EPmain.handles.preferences.preprocess.noadjacent,'Value');
        
        EPmain.preferences.preprocess.trialminmax=str2num(get(EPmain.handles.preferences.preprocess.trialminmax,'String'));
        
        EPmain.preferences.preprocess.movefacs=str2num(get(EPmain.handles.preferences.preprocess.movefacs,'String'));
        
        EPmain.preferences.preprocess.EMGratio=str2num(get(EPmain.handles.preferences.preprocess.EMGratio,'String'));
        
        EPmain.preferences.preprocess.EMGthresh=str2num(get(EPmain.handles.preferences.preprocess.EMGthresh,'String'));
        
        EPmain.preferences.preprocess.EOGchans=str2num(get(EPmain.handles.preferences.preprocess.EOGchans,'String'));
        
        EPmain.preferences.preprocess.sacPot=str2num(get(EPmain.handles.preferences.preprocess.sacPot,'String'));

        EPmain.preferences.preprocess.saturation=str2num(get(EPmain.handles.preferences.preprocess.saturation,'String'));

        EPmain.preferences.preprocess.fMRI=get(EPmain.handles.preferences.preprocess.fMRI,'Value');

        rotationList=get(EPmain.handles.preferences.preprocess.blinkRotation,'String');
        EPmain.preferences.preprocess.blinkRotation=rotationList{get(EPmain.handles.preferences.preprocess.blinkRotation,'Value')};

        rotationList=get(EPmain.handles.preferences.preprocess.saccRotation,'String');
        EPmain.preferences.preprocess.saccRotation=rotationList{get(EPmain.handles.preferences.preprocess.saccRotation,'Value')};

        rotationList=get(EPmain.handles.preferences.preprocess.SProtation,'String');
        EPmain.preferences.preprocess.SProtation=rotationList{get(EPmain.handles.preferences.preprocess.SProtation,'Value')};

        [prefErr prefMissing]=checkPrefs;
        
        if prefErr
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.preprocess=tempVar;
        else
            EPmain.mode=EPmain.prefReturn;
            ep('start');
        end
        
    case 'startPreferenceAverage'
        
        set(EPmain.handles.hMainWindow,'Name', 'Average Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Averaging Method','FontSize',EPmain.fontsize,...
            'Position',[20 470 150 20]);
        
        EPmain.handles.preferences.average.method = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Mean','Median','Trimmed Mean'},...
            'Value',EPmain.preferences.average.method,'Position',[20 450 150 20],'TooltipString','Central tendency estimator used for averaging procedure.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Trimming Level','FontSize',EPmain.fontsize,...
            'Position',[20 430 150 20]);
        
        EPmain.handles.preferences.average.trimLevel = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.average.trimLevel,'Position',[20 410 150 20],'TooltipString','If trimmed mean chosen, proportion of trials trimmed from each tail of the distribution.');
         
        uicontrol('Style','text','HorizontalAlignment','left','String', 'C','FontSize',EPmain.fontsize,...
            'Position',[20 390 10 20]);
        
        EPmain.handles.preferences.average.codeCorrect = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.average.codeCorrect,'Position',[30 390 20 20],'TooltipString','Value for Correct Response.');
          
        uicontrol('Style','text','HorizontalAlignment','left','String', 'E','FontSize',EPmain.fontsize,...
            'Position',[50 390 10 20]);
        
        EPmain.handles.preferences.average.codeError = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.average.codeError,'Position',[60 390 20 20],'TooltipString','Value for Error Response.');
          
        uicontrol('Style','text','HorizontalAlignment','left','String', 'T','FontSize',EPmain.fontsize,...
            'Position',[80 390 10 20]);
        
        EPmain.handles.preferences.average.codeTimeout = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.average.codeTimeout,'Position',[90 390 20 20],'TooltipString','Value for Response Timeout.');
      
        EPmain.handles.preferences.average.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', ['global EPmain;','EPmain.mode=EPmain.prefReturn;','ep(''start'');']);
        
        EPmain.handles.preferences.average.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', 'ep(''getAveragePrefs'');');
        
    case 'getAveragePrefs' %retrieve preference settings from the general preference input fields
        
        tempVar=EPmain.preferences.average;
        
        EPmain.preferences.average.method=get(EPmain.handles.preferences.average.method,'Value');
        
        EPmain.preferences.average.trimLevel=str2double(get(EPmain.handles.preferences.average.trimLevel,'String'));
        
        EPmain.preferences.average.codeCorrect=get(EPmain.handles.preferences.average.codeCorrect,'String');
        
        EPmain.preferences.average.codeError=get(EPmain.handles.preferences.average.codeError,'String');
        
        EPmain.preferences.average.codeTimeout=get(EPmain.handles.preferences.average.codeTimeout,'String');
        
        [prefErr prefMissing]=checkPrefs;
        
        if prefErr
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.average=tempVar;
        else
            EPmain.mode=EPmain.prefReturn;
            ep('start');
        end
        
    case 'startPreferenceTransform'
        
        set(EPmain.handles.hMainWindow,'Name', 'Transform Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Reference','FontSize',EPmain.fontsize,...
            'Position',[20 470 150 20]);
        
        EPmain.handles.preferences.transform.reference = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Average','Traditional','none'},...
            'Value',EPmain.preferences.transform.reference,'Position',[20 450 150 20],'TooltipString','Default type of referencing scheme to apply when transforming data.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Reference Channel 1','FontSize',EPmain.fontsize,...
            'Position',[20 430 150 20]);
        
        EPmain.handles.preferences.transform.refChan1 = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.transform.refChan1,'Position',[20 410 150 20],'TooltipString','If traditional reference chosen, default channel to use.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Reference Channel 2','FontSize',EPmain.fontsize,...
            'Position',[20 390 150 20]);
        
        EPmain.handles.preferences.transform.refChan2 = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.transform.refChan2,'Position',[20 370 150 20],'TooltipString','If traditional reference chosen, default second channel to use, if any.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Baseline Start','FontSize',EPmain.fontsize,...
            'Position',[20 350 150 20]);
        
        EPmain.handles.preferences.transform.baselineStart = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.transform.baselineStart,'Position',[20 330 150 20],'TooltipString','Default msec start (left side of sample) of period to use for baseline correction when transforming data.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Baseline End','FontSize',EPmain.fontsize,...
            'Position',[20 310 150 20]);
        
        EPmain.handles.preferences.transform.baselineEnd = uicontrol('Style','edit','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.transform.baselineEnd,'Position',[20 290 150 20],'TooltipString','Default msec end (right side of sample) of period to use for baseline correction when transforming data.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Fix Mains Noise','FontSize',EPmain.fontsize,...
            'Position',[20 270 150 20]);
        
        theString={'None','50','60'};
        EPmain.handles.preferences.transform.mainsFix = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,'Value',find(strcmp(EPmain.preferences.transform.mainsFix,theString)),...
            'String',theString,'Position',[20 250 150 20],'TooltipString','Fix mains noise at this frequency.');
        
        EPmain.handles.preferences.transform.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', ['global EPmain;','EPmain.mode=EPmain.prefReturn;','ep(''start'');']);
        
        EPmain.handles.preferences.transform.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', 'ep(''getTransformPrefs'');');
        
    case 'getTransformPrefs' %retrieve preference settings from the transform preference input fields
        
        tempVar=EPmain.preferences.transform;
        
        EPmain.preferences.transform.reference=get(EPmain.handles.preferences.transform.reference,'Value');
        
        EPmain.preferences.transform.refChan1=str2num(get(EPmain.handles.preferences.transform.refChan1,'String'));
        
        EPmain.preferences.transform.refChan2=str2num(get(EPmain.handles.preferences.transform.refChan2,'String'));
        
        EPmain.preferences.transform.baselineStart=str2num(get(EPmain.handles.preferences.transform.baselineStart,'String'));
        
        EPmain.preferences.transform.baselineEnd=str2num(get(EPmain.handles.preferences.transform.baselineEnd,'String'));
        
        theString=get(EPmain.handles.preferences.transform.mainsFix,'string');
        EPmain.preferences.transform.mainsFix=theString{get(EPmain.handles.preferences.transform.mainsFix,'Value')};

        [prefErr prefMissing]=checkPrefs;
        
        if prefErr
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.transform=tempVar;
        else
            EPmain.mode=EPmain.prefReturn;
            ep('start');
        end
        
    case 'startPreferencePCA'
        
        set(EPmain.handles.hMainWindow,'Name', 'PCA Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Mode','FontSize',EPmain.fontsize,...
            'Position',[20 470 150 20]);
        
        EPmain.handles.preferences.pca.mode = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'spatial','temporal'},...
            'Value',EPmain.preferences.pca.mode,'Position',[20 450 150 20],'TooltipString','When conducting PCA, whether default mode is temporal or spatial.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Rotation','FontSize',EPmain.fontsize,...
            'Position',[20 430 150 20]);
        
        EPmain.handles.preferences.pca.rotation = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Varimax','Promax','Infomax','Quartimax','Quartimin','Oblimin','Crawford-Ferguson','minimum entropy','invariant pattern simplicity','tandem II criterion','Geomin','McCammon minimum entropy','Variable-Oblimin'},...
            'Value',EPmain.preferences.pca.rotation,'Position',[20 410 150 20],'TooltipString','When conducting PCA, default rotation.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Parameter','FontSize',EPmain.fontsize,...
            'Position',[20 390 150 20]);
        
        EPmain.handles.preferences.pca.rotopt = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.pca.rotopt,'Position',[20 370 150 20],'TooltipString','When conducting PCA, default parameter for rotations that have a parameter.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Rel matrix','FontSize',EPmain.fontsize,...
            'Position',[20 350 150 20]);
        
        EPmain.handles.preferences.pca.rel = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'correlation','covariance'},...
            'Value',EPmain.preferences.pca.rel,'Position',[20 330 150 20],'TooltipString','When conducting PCA, the default relationship matrix.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Loadings','FontSize',EPmain.fontsize,...
            'Position',[20 310 150 20]);
        
        EPmain.handles.preferences.pca.loadings = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'none','Kaiser','covariance','C-M'},...
            'Value',EPmain.preferences.pca.loadings,'Position',[20 290 150 20],'TooltipString','When conducting PCA, the default factor loading weighting scheme.');
        
        EPmain.handles.preferences.pca.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', ['global EPmain;','EPmain.mode=EPmain.prefReturn;','ep(''start'');']);
        
        EPmain.handles.preferences.pca.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', 'ep(''getPCAPrefs'');');
        
    case 'getPCAPrefs' %retrieve preference settings from the PCA preference input fields
        
        tempVar=EPmain.preferences.pca;
        
        EPmain.preferences.pca.mode=get(EPmain.handles.preferences.pca.mode,'Value');
        
        EPmain.preferences.pca.rotation=get(EPmain.handles.preferences.pca.rotation,'Value');
        
        EPmain.preferences.pca.rotopt=str2num(get(EPmain.handles.preferences.pca.rotopt,'String'));
        
        EPmain.preferences.pca.rel=get(EPmain.handles.preferences.pca.rel,'Value');
        
        EPmain.preferences.pca.loadings=get(EPmain.handles.preferences.pca.loadings,'Value');
        
        [prefErr prefMissing]=checkPrefs;
        
        if prefErr
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.pca=tempVar;
        else
            EPmain.mode=EPmain.prefReturn;
            ep('start');
        end
        
    case 'startPreferenceView'
        
        set(EPmain.handles.hMainWindow,'Name', 'View Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Plot positive...','FontSize',EPmain.fontsize,...
            'Position',[20 470 150 20]);
        
        EPmain.handles.preferences.view.positive = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'up','down'},...
            'Value',EPmain.preferences.view.positive,'Position',[20 450 50 20],'TooltipString','When plotting waveforms, whether positive is up or down.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Label Font Size','FontSize',EPmain.fontsize,...
            'Position',[20 430 150 20]);
        
        EPmain.handles.preferences.view.labelSize = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.view.labelSize,'Position',[20 410 50 20],'TooltipString','Font size of the labels in the waveform plots.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Line Sizes','FontSize',EPmain.fontsize,...
            'Position',[20 390 150 20]);
        
        EPmain.handles.preferences.view.lineSize = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.view.lineSize,'Position',[20 370 50 20],'TooltipString','Width of the marking lines in the waveform plots.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Topo Map Electrodes','FontSize',EPmain.fontsize,...
            'Position',[20 350 150 20]);
        
        theList={'on','off','labels','numbers'};
        EPmain.handles.preferences.view.topoElectrodes = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',theList,...
            'Value',find(strcmp(EPmain.preferences.view.topoElectrodes,theList)),'Position',[20 330 100 20],'TooltipString','When plotting topo maps, how to mark electrodes.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String','R','FontSize',EPmain.fontsize,...
            'Position',[25 310 30 20]);
        
        uicontrol('Style','text','HorizontalAlignment','left','String','G','FontSize',EPmain.fontsize,...
            'Position',[60 310 30 20]);
        
        uicontrol('Style','text','HorizontalAlignment','left','String','B','FontSize',EPmain.fontsize,...
            'Position',[95 310 30 20]);
        
        uicontrol('Style','text','HorizontalAlignment','left','String','width','FontSize',EPmain.fontsize,...
            'Position',[130 310 35 20]);
        
        uicontrol('Style','text','HorizontalAlignment','left','String','style','FontSize',EPmain.fontsize,...
            'Position',[170 310 35 20]);
        
        for iColor=1:EPmain.maxColors
            
            EPmain.handles.preferences.view.color(iColor).label=uicontrol('Style','pushbutton','HorizontalAlignment','left','String',num2str(iColor),'FontSize',EPmain.fontsize,...
                'ForegroundColor',EPmain.preferences.view.color(iColor).RGB,'Position',[10 310-(20*iColor) 15 20],...
                'Callback', ['global EPmain;','EPmain.preferences.view.color(' num2str(iColor) ').RGB=uisetcolor(EPmain.preferences.view.color(' num2str(iColor) ').RGB);','ep(''start'');']);            
            EPmain.handles.preferences.view.color(iColor).R = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,'Callback',{@changeColor,iColor,1},...
                'String',EPmain.preferences.view.color(iColor).RGB(1),'Position',[25 310-(20*iColor) 30 20],'TooltipString',['Default: ' num2str(EPmain.defaultColor(iColor).RGB(1))]);
            
            EPmain.handles.preferences.view.color(iColor).G = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,'Callback',{@changeColor,iColor,2},...
                'String',EPmain.preferences.view.color(iColor).RGB(2),'Position',[60 310-(20*iColor) 30 20],'TooltipString',['Default: ' num2str(EPmain.defaultColor(iColor).RGB(2))]);
            
            EPmain.handles.preferences.view.color(iColor).B = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,'Callback',{@changeColor,iColor,3},...
                'String',EPmain.preferences.view.color(iColor).RGB(3),'Position',[95 310-(20*iColor) 30 20],'TooltipString',['Default: ' num2str(EPmain.defaultColor(iColor).RGB(3))]);
            
            EPmain.handles.preferences.view.color(iColor).lineSize = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,'Callback',{@changeColor,iColor,4},...
                'String',EPmain.preferences.view.color(iColor).lineSize,'Position',[130 310-(20*iColor) 30 20],'TooltipString','Thickness of waveform lines.');
            
            EPmain.handles.preferences.view.color(iColor).lineStyle = uicontrol('Style','popupmenu','HorizontalAlignment','left','FontSize',EPmain.fontsize,'Callback',{@changeColor,iColor,5},...
                'String',EPmain.lineStyleList,'Value',find(strcmp(EPmain.preferences.view.color(iColor).lineStyle,EPmain.lineStyleList)),...
                'Position',[170 310-(20*iColor) 30 20],'TooltipString','Style of waveform lines.');
        end
        
        EPmain.handles.preferences.view.reset = uicontrol('Style', 'pushbutton', 'String', 'Reset Colors','FontSize',EPmain.fontsize,...
            'Position', [20 100 100 40], 'Callback', @resetColors);
        
        EPmain.handles.preferences.view.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', ['global EPmain;','EPmain.mode=EPmain.prefReturn;','EPmain.preferences.view.color=EPmain.tempColor;','ep(''start'');']);
        
        EPmain.handles.preferences.view.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', 'ep(''getViewPrefs'');');
        
    case 'getViewPrefs' %retrieve preference settings from the view preference input fields
        
        tempVar=EPmain.preferences.view;
        
        EPmain.preferences.view.positive=get(EPmain.handles.preferences.view.positive,'Value');
        
        EPmain.preferences.view.labelSize=str2double(get(EPmain.handles.preferences.view.labelSize,'String'));
        
        EPmain.preferences.view.lineSize=str2double(get(EPmain.handles.preferences.view.lineSize,'String'));

        theList=get(EPmain.handles.preferences.view.topoElectrodes,'String');
        EPmain.preferences.view.topoElectrodes=theList{get(EPmain.handles.preferences.view.topoElectrodes,'Value')};

        %color settings already registered and could have been rolled back if cancel button clicked.
        
        [prefErr prefMissing]=checkPrefs;
        
        if prefErr
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.view=tempVar;
        else
            EPmain.mode=EPmain.prefReturn;
            ep('start');
        end
        
    case 'startPreferenceWindow'
        
        set(EPmain.handles.hMainWindow,'Name', 'Window Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Min Factor Variance','FontSize',EPmain.fontsize,...
            'Position',[20 470 150 20]);
        
        EPmain.handles.preferences.window.minFacVar = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.window.minFacVar,'Position',[20 450 150 20],'TooltipString','When using autoPCA option, the minimum percent of variance for a factor to be included.');
        
        EPmain.handles.preferences.window.adds= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Adds',...
            'Value',EPmain.preferences.window.adds,'Position',[20 430 150 20],'TooltipString','When windowing data, whether to add summary cell and channel waveforms corresponding to the analysis to the dataset.');
        
        EPmain.handles.preferences.window.chanGrp = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',{'Collapse First','Measure First'},...
            'Value',EPmain.preferences.window.chanGrp,'Position',[20 410 150 20],'TooltipString','Collapse channels then measure or vice versa.');
        
        EPmain.handles.preferences.window.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', ['global EPmain;','EPmain.mode=EPmain.prefReturn;','ep(''start'');']);
        
        EPmain.handles.preferences.window.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', 'ep(''getWindowPrefs'');');
        
    case 'getWindowPrefs' %retrieve preference settings from the window preference input fields
        
        tempVar=EPmain.preferences.window;
        
        EPmain.preferences.window.minFacVar=str2num(get(EPmain.handles.preferences.window.minFacVar,'String'));
        
        EPmain.preferences.window.adds=get(EPmain.handles.preferences.window.adds,'Value');
        
        EPmain.preferences.window.chanGrp=get(EPmain.handles.preferences.window.chanGrp,'Value');
        
        [prefErr prefMissing]=checkPrefs;
        
        if prefErr
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.window=tempVar;
        else
            EPmain.mode=EPmain.prefReturn;
            ep('start');
        end
        
    case 'startPreferenceANOVA'
        
        set(EPmain.handles.hMainWindow,'Name', 'ANOVA Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Trimming','FontSize',EPmain.fontsize,...
            'Position',[20 470 150 20]);
        
        EPmain.handles.preferences.anova.trimming = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.anova.trimming,'Position',[20 450 150 20],'TooltipString','Proportion of distribution of each cell''s tail to trim during ANOVA (0 to .5).');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Bootstrap Samples','FontSize',EPmain.fontsize,...
            'Position',[20 430 150 20]);
        
        EPmain.handles.preferences.anova.bootstrap = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.anova.bootstrap,'Position',[20 410 150 20],'TooltipString','Number of times to run bootstrap to generate ANOVA distributions.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Bootstrap Reps','FontSize',EPmain.fontsize,...
            'Position',[20 390 150 20]);
        
        EPmain.handles.preferences.anova.reps = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.anova.reps,'Position',[20 370 150 20],'TooltipString','Number of times to run ANOVA reps to determine p-value variability.  Needs to be an odd number.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Seed','FontSize',EPmain.fontsize,...
            'Position',[20 350 150 20]);
        
        EPmain.handles.preferences.anova.seed = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.anova.seed,'Position',[20 330 150 20],'TooltipString','The starting seed number to use for random number generator for bootstrapping.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Missing','FontSize',EPmain.fontsize,...
            'Position',[20 310 150 20]);
        
        EPmain.handles.preferences.anova.missing = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.anova.missing,'Position',[20 290 150 20],'TooltipString','Value used to indicate a missing number.');
        
        EPmain.handles.preferences.anova.adds= uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Adds',...
            'Value',EPmain.preferences.anova.adds,'Position',[20 270 150 20],'TooltipString','When conducting ANOVA, whether to add summary subject waveforms corresponding to the analysis to the dataset.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Epsilon','FontSize',EPmain.fontsize,...
            'Position',[20 250 150 20]);
        
        theString=ep_enumVar('epsilon');
        EPmain.handles.preferences.anova.epsilon = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',theString,...
            'Value',find(strcmp(EPmain.preferences.anova.epsilon,theString)),'Position',[20 230 150 20],'TooltipString','Epsilon correction factor when using conventional ANOVAs.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Post-Hoc','FontSize',EPmain.fontsize,...
            'Position',[20 210 150 20]);
        
        theString=ep_enumVar('postHoc');
        EPmain.handles.preferences.anova.posthoc = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
            'String',theString,...
            'Value',find(strcmp(EPmain.preferences.anova.posthoc,theString)),'Position',[20 190 150 20],'TooltipString','Between-group post-hoc when using conventional ANOVAs.');

%         uicontrol('Style','text','HorizontalAlignment','left','String', 'Degrees of Freedom','FontSize',EPmain.fontsize,...
%             'Position',[20 170 150 20]);
%         
%         theString={'ADF','epsilon'};
%         EPmain.handles.preferences.anova.adf = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
%             'String',theString,...
%             'Value',find(strcmp(EPmain.preferences.anova.adf,theString)),'Position',[20 150 150 20],'TooltipString','ANOVA output presents DF adjustment directly or separately as epsilon.');

        EPmain.handles.preferences.anova.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', ['global EPmain;','EPmain.mode=EPmain.prefReturn;','ep(''start'');']);
        
        EPmain.handles.preferences.anova.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', 'ep(''getANOVAPrefs'');');
        
    case 'getANOVAPrefs' %retrieve preference settings from the ANOVA preference input fields
        
        tempVar=EPmain.preferences.anova;
        
        EPmain.preferences.anova.trimming=str2num(get(EPmain.handles.preferences.anova.trimming,'String'));
        
        EPmain.preferences.anova.bootstrap=str2num(get(EPmain.handles.preferences.anova.bootstrap,'String'));
        
        EPmain.preferences.anova.reps=str2num(get(EPmain.handles.preferences.anova.reps,'String'));
        
        EPmain.preferences.anova.seed=str2num(get(EPmain.handles.preferences.anova.seed,'String'));
        
        EPmain.preferences.anova.missing=str2num(get(EPmain.handles.preferences.anova.missing,'String'));
        
        EPmain.preferences.anova.adds=get(EPmain.handles.preferences.anova.adds,'Value');
        
        theString=get(EPmain.handles.preferences.anova.epsilon,'string');
        EPmain.preferences.anova.epsilon=theString{get(EPmain.handles.preferences.anova.epsilon,'Value')};
        
        theString=get(EPmain.handles.preferences.anova.posthoc,'string');
        EPmain.preferences.anova.posthoc=theString{get(EPmain.handles.preferences.anova.posthoc,'Value')};
%         
%         theString=get(EPmain.handles.preferences.anova.adf,'string');
%         EPmain.preferences.anova.adf=theString{get(EPmain.handles.preferences.anova.adf,'Value')};
        
        [prefErr prefMissing]=checkPrefs;
        
        if prefErr
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.anova=tempVar;
        else
            EPmain.mode=EPmain.prefReturn;
            ep('start');
        end
        
    case 'startPreferenceAdvanced' %using handles2 to avoid weird bug in Matlab 2020a
        
        set(EPmain.handles.hMainWindow,'Name', 'Advanced Pref');
        
        EPmain.handles2.preferences.advanced.parallel = uicontrol('Style','checkbox','FontSize',EPmain.fontsize,...
            'String','Parallel','Value',EPmain.preferences.advanced.parallel,'Position',[20 470 150 20],'TooltipString','Enable use of Matlab Distributed Computing Toolbox.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Monitor Position (X & Y)','FontSize',EPmain.fontsize,...
            'Position',[20 450 200 20]);
        
        EPmain.handles2.preferences.advanced.monitor1 = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.advanced.monitor(1),'Position',[20 430 50 20],'TooltipString','Horizontal monitor position.  Zero to disable.');
        
        EPmain.handles2.preferences.advanced.monitor2 = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.advanced.monitor(2),'Position',[80 430 50 20],'TooltipString','Vertical monitor position.  Zero to disable.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Monitor Resolution (X & Y)','FontSize',EPmain.fontsize,...
            'Position',[20 410 200 20]);
        
        EPmain.handles2.preferences.advanced.monitor3 = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.advanced.monitor(3),'Position',[20 390 50 20],'TooltipString','Horizontal monitor resolution.  Zero to disable.');
        
        EPmain.handles2.preferences.advanced.monitor4 = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.advanced.monitor(4),'Position',[80 390 50 20],'TooltipString','Vertical monitor resolution.  Zero to disable.');

        uicontrol('Style','text','HorizontalAlignment','left','String', 'RAM','FontSize',EPmain.fontsize,...
            'Position',[20 370 150 20]);
        
        EPmain.handles2.preferences.advanced.RAM = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.advanced.RAM,'Position',[20 350 150 20],'TooltipString','How much RAM you have installed on your computer, in GB.  Zero to disable.');
        
        EPmain.handles2.preferences.advanced.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', ['global EPmain;','EPmain.mode=EPmain.prefReturn;','ep(''start'');']);
        
        EPmain.handles2.preferences.advanced.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', 'ep(''getAdvancedPrefs'');');
        
    case 'getAdvancedPrefs' %retrieve preference settings from the advanced preference input fields
        
        tempVar=EPmain.preferences.advanced;
        
        EPmain.preferences.advanced.parallel=get(EPmain.handles2.preferences.advanced.parallel,'Value');
        
        EPmain.preferences.advanced.monitor(1)=str2num(get(EPmain.handles2.preferences.advanced.monitor1,'String'));
        
        EPmain.preferences.advanced.monitor(2)=str2num(get(EPmain.handles2.preferences.advanced.monitor2,'String'));
        
        EPmain.preferences.advanced.monitor(3)=str2num(get(EPmain.handles2.preferences.advanced.monitor3,'String'));
        
        EPmain.preferences.advanced.monitor(4)=str2num(get(EPmain.handles2.preferences.advanced.monitor4,'String'));
        
        EPmain.preferences.advanced.RAM=str2num(get(EPmain.handles2.preferences.advanced.RAM,'String'));
        
        [prefErr prefMissing]=checkPrefs;
        
        if prefErr
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.advanced=tempVar;
        else
            if all(EPmain.preferences.advanced.monitor)
                EPmain.scrsz=EPmain.preferences.advanced.monitor;
            end
            EPmain.mode=EPmain.prefReturn;
            ep('start');
        end
        
    case 'startPreferenceRecords'
        
       set(EPmain.handles.hMainWindow,'Name', 'Records Pref');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'User','FontSize',EPmain.fontsize,...
            'Position',[20 450 200 20]);
        
        EPmain.handles2.preferences.records.user = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.records.user, 'Position',[20 430 150 20],'TooltipString','Name of the user.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Lab','FontSize',EPmain.fontsize,...
            'Position',[20 410 200 20]);
        
        EPmain.handles2.preferences.records.lab = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.records.lab,'Position',[20 390 150 20],'TooltipString','Name of the lab.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Institution','FontSize',EPmain.fontsize,...
            'Position',[20 370 200 20]);
        
        EPmain.handles2.preferences.records.institution = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.records.institution,'Position',[20 350 150 20],'TooltipString','Name of the institution.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Project','FontSize',EPmain.fontsize,...
            'Position',[20 330 200 20]);
        
        EPmain.handles2.preferences.records.project = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.records.project, 'Position',[20 310 150 20],'TooltipString','Name of the project.');
        
        uicontrol('Style','text','HorizontalAlignment','left','String', 'Experiment','FontSize',EPmain.fontsize,...
            'Position',[20 290 200 20]);
        
        EPmain.handles2.preferences.records.experiment = uicontrol('Style','edit','HorizontalAlignment','left','FontSize',EPmain.fontsize,...
            'String',EPmain.preferences.records.experiment, 'Position',[20 270 150 20],'TooltipString','Name of the experiment.');
        
        EPmain.handles2.preferences.records.cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel','FontSize',EPmain.fontsize,...
            'Position', [20 50 100 40], 'Callback', ['global EPmain;','EPmain.mode=EPmain.prefReturn;','ep(''start'');']);
        
        EPmain.handles2.preferences.records.done = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
            'Position', [20 0 100 40], 'Callback', 'ep(''getRecordsPrefs'');');
        
    case 'getRecordsPrefs' %retrieve preference settings from the records preference input fields
        
        tempVar=EPmain.preferences.records;
        
        EPmain.preferences.records.user=get(EPmain.handles2.preferences.records.user,'String');
        
        EPmain.preferences.records.lab=get(EPmain.handles2.preferences.records.lab,'String');
        
        EPmain.preferences.records.institution=get(EPmain.handles2.preferences.records.institution,'String');
        
        EPmain.preferences.records.project=get(EPmain.handles2.preferences.records.project,'String');
         
        EPmain.preferences.records.experiment=get(EPmain.handles2.preferences.records.experiment,'String');
       
        [prefErr prefMissing]=checkPrefs;
        
        if prefErr
            msg{1}=['New preferences are incorrect.  Either correct or cancel.'];
            [msg]=ep_errorMsg(msg);
            EPmain.preferences.records=tempVar;
        else
            EPmain.mode=EPmain.prefReturn;
            ep('start');
        end
        
    case 'savePrefs' %save current preference settings
        
        %         PathName=userpath;
%         if any(strcmp(PathName(end),{';',':'}))
%             PathName=PathName(1:end-1);
%         end
        
        prefs=EPmain.preferences;
%         if exist('~/Library/Preferences/EPprefs.mat','file')
%             eval('save ~/Library/Preferences/EPprefs.mat prefs');
%         elseif exist('EPprefs.mat','file')
%             location = which('EPprefs.mat');
%             eval(['save ''' location ''' prefs']);
%         elseif exist('~/Library/Preferences','dir')
%             resetPrefs;
%             prefs=EPmain.preferences;
%             eval('save ~/Library/Preferences/EPprefs.mat prefs');
         if exist([EPdataset.EPwork filesep 'EPwork' filesep 'EPprefs.mat'],'file')
             eval(['save ''' EPdataset.EPwork filesep 'EPwork' filesep 'EPprefs.mat'' prefs']);
        else
%            [FileName,PathName,FilterIndex] = uiputfile('','Save Preferences File','EPprefs');
            eval(['save ''' EPdataset.EPwork 'EPprefs'' prefs']);
%             if ~strcmp(PathName,path)
%                 addpath(PathName,'-end'); %adds directory with preferences file to Matlab's path
%             end
        end
        
        EPmain.mode='main';
        ep('start');
        
    otherwise
        disp('oops - programmer mistake.  Invalid EP mode.')
end

figure(EPmain.handles.hMainWindow)
drawnow


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function deleteReadData(src,eventdata)
global EPdataset EPmain EPeeg

if isempty(eventdata.Indices) %if just deselecting
    return;
end
theDataset=eventdata.Indices(1);

set(EPmain.handles.read.readTable,'enable','off');
set(EPmain.handles.read.read,'enable','off');
set(EPmain.handles.read.done,'enable','off');
drawnow

ep_tictoc('ioStart');
if exist('EPeeg','var') && ~isempty(EPeeg)
    theEEG=EPdataset.dataset(theDataset).inMemory;
    if theEEG>0
        EPeeg(theEEG)=[];
        largerList=find([EPdataset.dataset.inMemory]>theEEG);
        for iEEG=1:length(largerList)
            theLarger=largerList(iEEG);
            EPdataset.dataset(theLarger).inMemory=EPdataset.dataset(theLarger).inMemory-1;
        end
    end
end
delete([EPdataset.EPwork filesep 'EPwork' filesep EPdataset.dataset(theDataset).dataName '.mat']);
EPdataset.dataset=EPdataset.dataset(setdiff([1:length(EPdataset.dataset)],theDataset));
eval(['save ''' EPdataset.EPwork filesep 'EPwork' filesep 'EPdataset'' EPdataset'], '-v7.3');
ep_tictoc('ioFinish');

set(EPmain.handles.read.readTable,'enable','on');
set(EPmain.handles.read.read,'enable','on');
set(EPmain.handles.read.done,'enable','on');
drawnow

ep('start');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pickEditData(src,eventdata)
global EPoverview EPmain

if isempty(eventdata.Indices) %if just deselecting
    return;
end

if isempty(EPoverview) || EPoverview.done
    EPoverview=[];
end

EPoverview.dataset=eventdata.Indices(1);

EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);

EPoverview.mode=[];
ep_editData


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pickPCAdata(src,eventdata)
global EPmain EPdataset EPscree EPtictoc

if isempty(eventdata.Indices) %if just deselecting
    return;
end

ep_tictoc('begin');

EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);

EPscree=[];
screeFigure=findobj('Name', 'ScreeWindow');
if ~isempty(screeFigure)
    close(screeFigure)
end

theDataset=eventdata.Indices(1);

if EPmain.pca.mode == 4
    if isempty(EPmain.pca.crossVerifyPCA)
        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        return
    end
    theDataset=EPmain.pca.targetDatasets(theDataset);
end

theDataFull=ep_loadEPdataset(theDataset);

if ~isempty(theDataFull.facNames) && isempty(theDataFull.facVecS) && isempty(theDataFull.facVecT) && isempty(theDataFull.facVecF)
    warndlg('Error: This factor file does not have the original PCA information so it cannot be further factored.');
    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
    return;    
end

theData=ep_stripAdds(theDataFull);
if isempty(theData)
    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
    return;    
end

if length(theData.subNames) ~= length(theDataFull.subNames)
    disp(['Stripping off ' num2str(length(theDataFull.subNames)-length(theData.subNames)) ' subject adds.']);
end

if length(theData.cellNames) ~= length(theDataFull.cellNames)
    disp(['Stripping off ' num2str(length(theDataFull.cellNames)-length(theData.cellNames)) ' cell adds.']);
end

EEGchans=find(strcmp('EEG',theData.chanTypes));
numChans=length(theData.chanNames);

if length(EEGchans) ~= numChans
    disp(['Stripping off ' num2str(numChans-length(EEGchans)) ' non-EEG channels and regional EEG channels.']);
end

theData=ep_selectData(theData,{EEGchans,[],[],[],[],[]});

if isempty(theData) || isempty(theData.data)
    warndlg(['Error: The file had no data left after additions were removed.']);
    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
    return;
end

if strcmp(theData.dataType,'continuous') && (length(theData.timeNames) > 1000) && (EPmain.pca.mode == 2)
    warndlg(['Error: Attempting a temporal PCA on continuous data would likely crash the computer due to excessive memory usage.']);
    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
    return;
end

PCAname=get(EPmain.handles.pca.name,'String');
if isempty(PCAname)
    PCAname='PCA';
end
sameName=1;
suffix=0;
PCAnameSuffix=PCAname;
while sameName
    sameName=0;
    for i=1:length(EPdataset.dataset)
        if strcmp(EPdataset.dataset(i).dataName,PCAnameSuffix)
            sameName=1;
        end
    end
    if sameName
        suffix=suffix+1;
        PCAnameSuffix=[PCAname '-' num2str(suffix)];
    end
end

ep_tictoc;if EPtictoc.stop;ep('start');return;end
if EPmain.pca.mode == 4
    crossVerifyPCA=ep_loadEPdataset(EPmain.pca.crossVerifyPCA);
    crossVerifyPCA=ep_stripAdds(crossVerifyPCA);
    
    %check to see if the data has been edited since the PCA
    if isfield(crossVerifyPCA.pca,'PCAmode')
        switch crossVerifyPCA.pca.PCAmode
            case 'spat'
                if size(crossVerifyPCA.pca.FacPat,1) ~= length(crossVerifyPCA.chanNames)
                    warndlg(['Error: The channels of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return;
                end
            case 'temp'
                if size(crossVerifyPCA.pca.FacPat,1) ~= length(crossVerifyPCA.timeNames)
                    warndlg(['Error: The time points of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return;
                end
            case 'freq'
                if size(crossVerifyPCA.pca.FacPat,1) ~= length(crossVerifyPCA.freqNames)
                    warndlg(['Error: The frequencies of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return;
                end
        end
    end
    if isfield(crossVerifyPCA.pca,'PCAmode2')
        switch crossVerifyPCA.pca.PCAmode2
            case 'spat'
                if size(crossVerifyPCA.pca.FacPatST,1) ~= length(crossVerifyPCA.chanNames)
                    warndlg(['Error: The channels of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return;
                end
            case 'temp'
                if size(crossVerifyPCA.pca.FacPatST,1) ~= length(crossVerifyPCA.timeNames)
                    warndlg(['Error: The time points of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return;
                end
            case 'freq'
                if size(crossVerifyPCA.pca.FacPatST,1) ~= length(crossVerifyPCA.freqNames)
                    warndlg(['Error: The frequencies of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return;
                end
        end
    end
    if isfield(crossVerifyPCA.pca,'PCAmode3')
        switch crossVerifyPCA.pca.PCAmode
            case 'spat'
                if size(crossVerifyPCA.pca.FacPat3,1) ~= length(crossVerifyPCA.chanNames)
                    warndlg(['Error: The channels of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return;
                end
            case 'temp'
                if size(crossVerifyPCA.pca.FacPat3,1) ~= length(crossVerifyPCA.timeNames)
                    warndlg(['Error: The time points of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return;
                end
            case 'freq'
                if size(crossVerifyPCA.pca.FacPat3,1) ~= length(crossVerifyPCA.freqNames)
                    warndlg(['Error: The frequencies of the chosen PCA cross-verification file have been changed since the PCA was conducted.']);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return;
                end
        end
    end
    
    if isfield(crossVerifyPCA.pca,'PCAmode')
        CVpca=crossVerifyPCA.pca.FacCof;
        FactorResults = ep_doPCA(crossVerifyPCA.pca.PCAmode, crossVerifyPCA.pca.ROTATION, crossVerifyPCA.pca.RotOpt, crossVerifyPCA.pca.decomp, crossVerifyPCA.pca.MAT_TYPE, crossVerifyPCA.pca.numFacs, theData, crossVerifyPCA.pca.LOADING, 'N', [], CVpca, crossVerifyPCA.pca.theAlgorithm);
        if EPtictoc.stop
            EPtictoc.stop=0;
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            ep('start');
        end
        if isempty(FactorResults)
            msg{1}='Error: PCA was not successful.';
            [msg]=ep_errorMsg(msg);
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
        [theData] = ep_PCAoutput(FactorResults, [], []);
        if isempty(theData)
            msg{1}='PCA data reconstruction aborted.';
            [msg]=ep_errorMsg(msg);
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
        theData.dataName=PCAnameSuffix;
    end
    if isfield(crossVerifyPCA.pca,'PCAmode2')
        CVpca=crossVerifyPCA.pca.FacCofST;
        [FactorResults] = ep_doPCAst(theData.pca, crossVerifyPCA.pca.ROTATION2, crossVerifyPCA.pca.RotOpt2, crossVerifyPCA.pca.decomp2, crossVerifyPCA.pca.MAT_TYPE2, crossVerifyPCA.pca.numFacs2, crossVerifyPCA.pca.LOADING2,crossVerifyPCA.pca.PCAmode2,CVpca, crossVerifyPCA.pca.theAlgorithm2);
        if EPtictoc.stop
            EPtictoc.stop=0;
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            ep('start');
        end
        if isempty(FactorResults)
            msg{1}='Error: PCA was not successful.';
            [msg]=ep_errorMsg(msg);
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
        [theData] = ep_PCAoutput(FactorResults, [], []);
        if isempty(theData)
            msg{1}='PCA data reconstruction aborted.';
            [msg]=ep_errorMsg(msg);
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
        theData.dataName=PCAnameSuffix;
    end
    if isfield(crossVerifyPCA.pca,'PCAmode3')
        CVpca=crossVerifyPCA.pca.FacCof3;
        [FactorResults] = ep_doPCAst(theData.pca, crossVerifyPCA.pca.ROTATION3, crossVerifyPCA.pca.RotOpt3, crossVerifyPCA.pca.decomp3, crossVerifyPCA.pca.MAT_TYPE3, crossVerifyPCA.pca.numFacs3, crossVerifyPCA.pca.LOADING3,crossVerifyPCA.pca.PCAmode3,CVpca, crossVerifyPCA.pca.theAlgorithm3);
        if EPtictoc.stop
            EPtictoc.stop=0;
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
        if isempty(FactorResults)
            msg{1}='Error: PCA was not successful.';
            [msg]=ep_errorMsg(msg);
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
        [theData] = ep_PCAoutput(FactorResults, [], []);
        if isempty(theData)
            msg{1}='PCA data reconstruction aborted.';
            [msg]=ep_errorMsg(msg);
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
        theData.dataName=PCAnameSuffix;
    end
    
    ep_saveEPdataset(theData,length(EPdataset.dataset)+1,'no');
    
else
    modeNum = get(EPmain.handles.pca.mode,'value');
    switch modeNum
        case 1
            PCAmode='spat';
        case 2
            PCAmode='temp';
        case 3
            PCAmode='freq';
    end
    EPmain.pca.mode=modeNum;
    
    rotationNum = get(EPmain.handles.pca.rotation,'value');
    rotationName = get(EPmain.handles.pca.rotation,'string');
    PCArotation=rotationName{rotationNum};
    
    switch PCArotation
        case 'Varimax'
            disp('Kaiser, H. F. (1958). The varimax criterion for analytic rotation in factor analysis. Psychometrika, 23, 187-200.');
        case 'Promax'
            disp('Hendrickson, A. E. & White, P. O. (1964). Promax: A quick method for rotation to oblique simple structure. The British Journal of Statistical Psychology, 17, 65-70.');
        case 'Infomax'
            disp('Bell, A. J. & Sejnowski, T. J. (1995). An information-maximisation approach to blind separation and blind deconvolution. Neural Computation, 7(6), 1129-1159.');
        case 'Extended-Infomax'
            disp('Lee, T.-W., Girolami, M., & Sejnowski, T. (1999). Independent component analysis using an extended infomax algorithm for mixed sub-Gaussian and super-Gaussian sources. Neural Computation, 11(2), 609–633.');
        case 'Quartimax'
            disp('Carroll, J. B. (1953). Analytical solution for approximating simple structure in factor analysis. Psychometrika, 18, 23-38.');
        case 'Quartimin'
            disp('Carroll, J. B. (1953). Analytical solution for approximating simple structure in factor analysis. Psychometrika, 18, 23-38.');
        case 'Oblimin'
            disp('Carroll, J. B. (1957). Biquartimin Criterion for Rotation to Oblique Simple Structure in Factor Analysis. Science, 126(3283), 1114-1115.');
        case 'Crawford-Ferguson'
            disp('Crawford, C. B. & Ferguson, G. A. (1970). A general rotation criterion and its use in orthogonal rotation. Psychometrika, 35, 321-332.');
        case 'minimum entropy'
            disp('Jennrich, R. I. (2004). Rotation to simple loadings using component loss functions: The orthogonal case. Psychometrika, 69, 257-273.');
        case 'invariant pattern simplicity'
            disp('Bentler, P. M. (1977). Factor simplicity index and transformations. Psychometrika, 42, 277-295.');
        case 'tandem II criterion'
            disp('Comrey, A. L. (1967). Tandem criteria for analytic rotation in factor analysis. Psychometrika, 32, 277-295.');
        case 'Geomin'
            disp('Yates, A. (1987). Multivariate exploratory data analysis: A perspective on exploratory factor analysis. Albany, New York: State University of New York Press.');
        case 'McCammon minimum entropy'
            disp('McCammon, R. B. (1966). Principal components analysis and its application in large-scale correlation studies. Journal of Geology, 74, 721-733.');
        case 'Variable-Oblimin'
            disp('Carroll, J. B. (1957). Biquartimin Criterion for Rotation to Oblique Simple Structure in Factor Analysis. Science, 126(3283), 1114-1115.');
            disp('Gorsuch, R. L.  (1983).  Factor analysis, 2nd edition.  Hillsdale, NJ:Lawrence Erlbaum Associates.');
        case 'JADE'
            disp('Cardoso, J. F., & Souloumiac, A. (1993). Blind beamforming for non-gaussian signals. IEE Proceedings F (Radar and Signal Processing), 140(6), 362–370.');
        case 'SOBI'
            disp('Belouchrani, A., Abed-Meraim, K, Cardoso, J. F., & Moulines, E. (1993). Second-order blind separation of temporally correlated sources. Proc. Int. Conf. Digital Signal Processing, 346–351.');
        case 'fastICA'
            disp('Hyvarinen, A. (1999). Fast and robust fixed-point algorithms for independent component analysis. IEEE Transactions on Neural Networks, 10(3), 626–634.');
    end
            
    if strcmp(PCArotation,'Infomax')
        disp('Using EEGlab function runica to perform Infomax rotation.');
    end
    
    EPmain.pca.rotation=rotationNum;
    
    EPmain.pca.rotopt=str2num(get(EPmain.handles.pca.rotopt,'String'));
    
    decompList=get(EPmain.handles.pca.decomp,'String');
    PCAdecomp=decompList{get(EPmain.handles.pca.decomp,'value')};
    
    relNum = get(EPmain.handles.pca.rel,'value');
    switch relNum
        case 1
            PCArel='COR';
        case 2
            PCArel='COV';
    end
    EPmain.pca.rel=relNum;
    
    loadingNum = get(EPmain.handles.pca.loadings,'value');
    switch loadingNum
        case 1
            PCAloading='N';
        case 2
            PCAloading='K';
        case 3
            PCAloading='C';
        case 4
            PCAloading='W';
    end
    EPmain.pca.loading=loadingNum;
    
    parametric = get(EPmain.handles.pca.parametric,'value');
    EPmain.pca.parametric=parametric;
    
    EPmain.pca.facNum=str2num(get(EPmain.handles.pca.facNum,'String'));
    
    if parametric
        [parametricFile, pathname] = uigetfile('*.*','Parametric File');
        if parametricFile ~= 0
            parametricData=importdata([pathname parametricFile],sprintf('\t'),1);
            if isempty(parametricData)
                msg{1}=['Was unable to retrieve the parametric data.  It needs to be a tab-delimited text file in which the first row are the labels.'];
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
        else
            msg{1}=['Parametric PCA cancelled.  If using this option was not your intent, uncheck the box.'];
            [msg]=ep_errorMsg(msg);
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
    else
        parametricData.data=[];
        parametricData.colheaders=[];
    end
    
    disp('Starting the PCA.');
    ep_tictoc;if EPtictoc.stop;ep('start');return;end
    
    if strcmp(PCAmode,'temp')  && isempty(theData.timeNames)
        warndlg(['Error: The file has no time information with which to perform a temporal PCA.']);
        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        return;
    end
    
    if strcmp(PCAmode,'spat')  && isempty(theData.chanNames)
        warndlg(['Error: The file has no spatial information with which to perform a spatial PCA.']);
        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        return;
    end
    
    if strcmp(PCAmode,'freq')  && isempty(theData.freqNames)
        warndlg(['Error: The file has no frequency information with which to perform a frequency PCA.']);
        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        return;
    end
    
    if any(any(theData.analysis.badTrials)) && strcmp(theData.dataType,'single_trial')
        disp('Conducting PCA despite presence of bad data.  Observations with bad data will be left out from PCA calculations')
        disp('so they will not affect results but doing so will underweight their aspects of the dataset.');
    end
    
    if any(any(any(theData.analysis.badChans < 0))) && strcmp(theData.dataType,'single_trial')
        disp('Conducting PCA despite presence of bad data.  Observations with bad data will be left out from PCA calculations')
        disp('so they will not affect results but doing so will underweight their aspects of the dataset.');
    end
    
    if any(any(any(isnan(theData.analysis.badChans)))) && strcmp(theData.dataType,'average')
        disp('Conducting PCA despite presence of bad data.  Observations with bad data will be left out from PCA calculations')
        disp('so they will not affect results but doing so will underweight their aspects of the dataset.');
    end
    
    if strcmp('PCAdecomp','NIPALS')
        disp('For NIPALS: Filippo Amato (2020). NIPALS (https://www.github.com/fulepo/NIPALS), GitHub. Retrieved March 14, 2020.')
    end
    
    if ~isfield(theData.pca,'PCAmode')
        if EPmain.pca.facNum==0 %scree test only needs unrotated solution
            disp('Scree: The PCA of the data.');
            FactorResults = ep_doPCA(PCAmode, 'unrotated', EPmain.pca.rotopt, PCAdecomp, PCArel, 1, theData, PCAloading);
            if isempty(FactorResults)
                msg{1}='Error: PCA was not successful.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            if ~isempty(theData.noise) && any(any(any(any(any(theData.noise)))))
                disp('Scree: The PCA of the noise data.');
                noiseResults = ep_doPCA(PCAmode, 'unrotated', EPmain.pca.rotopt, PCAdecomp, PCArel, 1, theData.noise, PCAloading);
                if isempty(noiseResults)
                    msg{1}='Error: PCA was not successful.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
            else
                noiseResults = [];
            end
            disp('Scree: The PCA of the random data.');
            randResults = ep_doPCA(PCAmode, 'unrotated', EPmain.pca.rotopt, PCAdecomp, PCArel, 1, randn(size(theData.data)), PCAloading);
            if isempty(randResults)
                msg{1}='Error: PCA was not successful.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
        else
            FactorResults = ep_doPCA(PCAmode, PCArotation, EPmain.pca.rotopt, PCAdecomp, PCArel, EPmain.pca.facNum, theData, PCAloading, 'N', [], [], EPmain.pca.theAlgorithm);
            if EPtictoc.stop
                EPtictoc.stop=0;
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            if isempty(FactorResults)
                msg{1}='Error: PCA was not successful.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            
            %jack-knife analysis
            numSubs=length(theData.subNames);
            if numSubs > 1
                disp('Computing Jack-Knife PCA loadings.');
                fprintf('%60s\n',' ' );
                jackLoadings=zeros(size(FactorResults.FacPat,1),size(FactorResults.FacPat,2),numSubs);
                jackSD=zeros(length(FactorResults.varSD),numSubs);
                for theSubject=1:numSubs
                    ep_tictoc;if EPtictoc.stop;ep('start');return;end
                    fprintf('%s%-60s',repmat(sprintf('\b'),1,60),sprintf('%s%4d%s%4d.','Working on subject# ', theSubject, ' of ', numSubs))
                    theDataJK=theData;
                    subjectList=setdiff([1:numSubs],theSubject);
                    theDataJK.data=theDataJK.data(:,:,:,subjectList,:,:,:);
                    theDataJK.subNames=theDataJK.subNames(subjectList);
                    theDataJK.subTypes=theDataJK.subTypes(subjectList);
                    if ~isempty(theDataJK.subjectSpecs)
                        theDataJK.subjectSpecs=theDataJK.subjectSpecs(subjectList,:);
                    end
                    theDataJK.avgNum=theDataJK.avgNum(subjectList,:);
                    theDataJK.covNum=theDataJK.covNum(subjectList,:);
                    theDataJK.subNum=theDataJK.subNum(subjectList,:);
                    theDataJK.events=theDataJK.events(subjectList,:);
                    theDataJK.analysis.badChans=theDataJK.analysis.badChans(subjectList,:,:);
                    theDataJK.analysis.moveTrial=theDataJK.analysis.moveTrial(subjectList,:);
                    theDataJK.analysis.blinkTrial=theDataJK.analysis.blinkTrial(subjectList,:);
                    theDataJK.analysis.saccadeTrial=theDataJK.analysis.saccadeTrial(subjectList,:);
                    theDataJK.analysis.saccadeOnset=theDataJK.analysis.saccadeOnset(subjectList,:);
                    theDataJK.analysis.badTrials=theDataJK.analysis.badTrials(subjectList,:);
                    FactorResultsJK = ep_doPCA(PCAmode, PCArotation, EPmain.pca.rotopt, PCAdecomp, PCArel, EPmain.pca.facNum, theDataJK, PCAloading, 'N', [], [], EPmain.pca.theAlgorithm);
                    if ~isempty(FactorResultsJK)
                        jackLoadings(:,:,theSubject)=FactorResultsJK.FacPat;
                        jackSD(:,theSubject)=FactorResultsJK.varSD;
                    end
                end
                FactorResults.jack.FacPat=jackLoadings;
                FactorResults.jack.varSD=jackSD;
                fprintf('%60s\n',' ' );
            end
        end
    elseif isfield(theData.pca,'PCAmode3')
        msg{1}='These data have already undergone a three-step PCA process.';
        [msg]=ep_errorMsg(msg);
        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        return
    elseif isfield(theData.pca,'PCAmode') || isfield(theData.pca,'PCAmode2')
        if isfield(theData.pca,'PCAmode')
            if strcmp(PCAmode,theData.pca.PCAmode)
                if strcmp(PCAmode,'temp')
                    msg{1}='These data have already undergone the temporal PCA step.';
                elseif strcmp(PCAmode,'spat')
                    msg{1}='These data have already undergone the spatial PCA step.';
                elseif strcmp(PCAmode,'freq')
                    msg{1}='These data have already undergone the frequency PCA step.';
                else
                    msg{1}='Error: PCA mode not recognized.';
                end
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
        end
        if isfield(theData.pca,'PCAmode2')
            if strcmp(PCAmode,theData.pca.PCAmode2)
                if strcmp(PCAmode,'temp')
                    msg{1}='These data have already undergone the temporal PCA step.';
                elseif strcmp(PCAmode,'spat')
                    msg{1}='These data have already undergone the spatial PCA step.';
                elseif strcmp(PCAmode,'freq')
                    msg{1}='These data have already undergone the frequency PCA step.';
                else
                    msg{1}='Error: PCA mode not recognized.';
                end
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
        end
        if EPmain.pca.facNum==0
            disp('Scree: The PCA of the data.');
            FactorResults = ep_doPCAst(theData.pca, 'unrotated', EPmain.pca.rotopt, PCAdecomp, PCArel, 1, PCAloading,PCAmode); %scree test only needs unrotated solution
            if EPtictoc.stop
                EPtictoc.stop=0;
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                ep('start');
            end
            if isempty(FactorResults)
                msg{1}='Error: PCA was not successful.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            if isfield(theData.pca,'noise') && ~isempty(theData.pca.noise) && any(any(any(any(any(theData.pca.noise)))))
                %perform two-step PCA of noise data.  Since it is meant to show results when the signal is removed,
                %use the first step scoring coefficients so that they are within the same "window" as the regular data.
                
                disp('Scree: The PCA of the noise data.');
                noiseResults=theData.pca;
                
                numchan=size(noiseResults.noise,1);
                timePoints=size(noiseResults.noise,2);
                numCells=size(noiseResults.noise,3);
                numSubs=size(noiseResults.noise,4);
                
                if strcmp(noiseResults.PCAmode,'temp')
                    data2=zeros(numchan*numCells*numSubs,timePoints);
                    for cell = 1:numCells
                        for iSub = 1:numSubs
                            data2((iSub-1)*numchan+(cell-1)*numSubs*numchan+1:iSub*numchan+(cell-1)*numSubs*numchan,:)=squeeze(noiseResults.noise(:,:,cell,iSub));
                        end
                    end
                elseif strcmp(noiseResults.PCAmode,'spat')
                    data2=zeros(timePoints*numCells*numSubs,numchan);
                    for cell = 1:numCells
                        for iSub = 1:numSubs
                            data2((iSub-1)*timePoints+(cell-1)*numSubs*timePoints+1:iSub*timePoints+(cell-1)*numSubs*timePoints,:)=squeeze(noiseResults.noise(:,:,cell,iSub))';
                        end
                    end
                elseif strcmp(noiseResults.PCAmode,'asis')
                    data2=noiseResults.noise;
                else
                    error('PCAmode must be set to either temp or spat or asis');
                end
                
                noiseResults.FacScr=data2*noiseResults.FacCof;
                noiseResults.FacScr=(noiseResults.FacScr)*inv(diag(std(noiseResults.FacScr))); %Standardize factor scores, not mean corrected.
                if isfield(theData.pca,'PCAmode2')
                    noiseResults = ep_doPCAst(noiseResults, theData.pca.ROTATION2, theData.pca.RotOpt2, theData.pca.decomp2, theData.pca.MAT_TYPE2, theData.pca.numFacs2, theData.pca.LOADING2, theData.pca.PCAmode2, [], theData.pca.theAlgorithm2);
                    if isempty(noiseResults)
                        msg{1}='Error: PCA was not successful.';
                        [msg]=ep_errorMsg(msg);
                        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                        return
                    end
                end
                noiseResults = ep_doPCAst(noiseResults, 'unrotated', EPmain.pca.rotopt, PCAdecomp, PCArel, 1, PCAloading,PCAmode);
                if isempty(noiseResults)
                    msg{1}='Error: PCA was not successful.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
            else
                noiseResults = [];
            end
            disp('Scree: The PCA of the random data.');
            randResults = ep_doPCA(theData.pca.PCAmode, theData.pca.ROTATION, theData.pca.RotOpt, theData.pca.decomp, theData.pca.MAT_TYPE, theData.pca.numFacs, randn(theData.pca.numchan,theData.pca.timepoints,theData.pca.numCells,theData.pca.numSubs,1,theData.pca.numFreqs), theData.pca.LOADING, 'N', [], [], theData.pca.theAlgorithm);
            if isempty(randResults)
                msg{1}='Error: PCA was not successful.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            if isfield(theData.pca,'PCAmode2')
                randResults = ep_doPCAst(randResults, theData.pca.ROTATION2, theData.pca.RotOpt2, theData.pca.decomp2, theData.pca.MAT_TYPE2, theData.pca.numFacs2, theData.pca.LOADING2, theData.pca.PCAmode2, [], theData.pca.theAlgorithm2);
                if isempty(randResults)
                    msg{1}='Error: PCA was not successful.';
                    [msg]=ep_errorMsg(msg);
                    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                    return
                end
            end
            randResults = ep_doPCAst(randResults, 'unrotated', EPmain.pca.rotopt, PCAdecomp, PCArel, 1, PCAloading, PCAmode);
            if isempty(randResults)
                msg{1}='Error: PCA was not successful.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
        else
            [FactorResults] = ep_doPCAst(theData.pca, PCArotation, EPmain.pca.rotopt, PCAdecomp, PCArel, EPmain.pca.facNum, PCAloading,PCAmode, [], EPmain.pca.theAlgorithm);
            if EPtictoc.stop
                EPtictoc.stop=0;
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                ep('start');
            end
            if isempty(FactorResults)
                msg{1}='Error: PCA was not successful.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
            
            %jack-knife analysis
            numSubs=length(theData.subNames);
            numCells=length(theData.cellNames);
            if isfield(theData.pca,'PCAmode2')
                numVars=size(FactorResults.FacScrST,1)/(numSubs*numCells);
            else
                numVars=size(FactorResults.FacScr,1)/(numSubs*numCells);
            end
            if numSubs > 1
                disp('Computing Jack-Knife PCA loadings.');
                fprintf('%60s\n',' ' );
                jackLoadingsST=zeros(size(FactorResults.FacPatST,1),size(FactorResults.FacPatST,2),numSubs);
                jackSDST=zeros(size(FactorResults.varSDST,1),size(FactorResults.varSDST,2),numSubs);
                for theSubject=1:numSubs
                    ep_tictoc;if EPtictoc.stop;ep('start');return;end
                    fprintf('%s%-60s',repmat(sprintf('\b'),1,60),sprintf('%s%4d%s%4d.','Working on subject# ', theSubject, ' of ', numSubs))
                    theDataJK=theData.pca;
                    subjectVec=ones(numSubs,1);
                    subjectVec(theSubject)=0;
                    subjectList=kron(ones(numCells,1),subjectVec);
                    subjectList=logical(kron(subjectList,ones(numVars,1)));
                    if isfield(theData.pca,'PCAmode2')
                        theDataJK.FacScrST=theDataJK.FacScrST(subjectList,:);
                        theDataJK.badObsST=theDataJK.badObsST(subjectList,:);
                    else
                        theDataJK.FacScr=theDataJK.FacScr(subjectList,:);
                        theDataJK.badObs=theDataJK.badObs(subjectList);
                    end
                    theDataJK.numSubs=numSubs-1;
                    [FactorResultsJK] = ep_doPCAst(theDataJK, PCArotation, EPmain.pca.rotopt, PCAdecomp, PCArel, EPmain.pca.facNum, PCAloading,PCAmode, [], EPmain.pca.theAlgorithm);
                    if ~isempty(FactorResultsJK)
                        jackLoadingsST(:,:,theSubject)=FactorResultsJK.FacPatST;
                        jackSDST(:,:,theSubject)=FactorResultsJK.varSDST;
                    end
                end
                FactorResults.jack.FacPatST=jackLoadingsST;
                FactorResults.jack.varSDST=jackSDST;
                fprintf('%60s\n',' ' );
            end
        end
    else
        msg{1}='These data have already undergone a PCA process but is corrupted in some manner.';
        [msg]=ep_errorMsg(msg);
        ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
        return
    end
    
    ep_tictoc;if EPtictoc.stop;ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);ep('start');return;end
    if EPmain.pca.facNum==0
        ep_scree('start',FactorResults,noiseResults,randResults);
    else
        if parametric
            if size(parametricData.data,1) ~= (size(FactorResults.FacScr,1)/FactorResults.numchan)
                msg{1}='The number of observations in the parameters file does not match the data.  It needs to be a tab-delimited text file in which the first row are the labels.';
                [msg]=ep_errorMsg(msg);
                ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
                return
            end
        end
        
        [PCAoutput] = ep_PCAoutput(FactorResults, parametricData.data, parametricData.colheaders,1,EPmain.preferences.anova.missing);
        if isempty(PCAoutput)
            msg{1}='PCA data reconstruction aborted.';
            [msg]=ep_errorMsg(msg);
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            return
        end
        
        PCAoutput.dataName=PCAnameSuffix;
        ep_saveEPdataset(PCAoutput,length(EPdataset.dataset)+1,'no');
    end
    
end

ep_tictoc('end');
ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resetPrefs(src,eventdata)
%reset the permanent preference settings

global EPmain EP_VER

if ~isfield(EPmain,'preferences')
    EPmain.preferences=[];
end
if ~isfield(EPmain.preferences,'general')
    EPmain.preferences.general=[];
end
if ~isfield(EPmain.preferences.general,'sessionImportFormat')
    EPmain.preferences.general.sessionImportFormat=1; %EP format
end
if ~isfield(EPmain.preferences.general,'sessionOutputFormat')
    EPmain.preferences.general.sessionOutputFormat=1; %EP format
end
if ~isfield(EPmain.preferences.general,'importFormat')
    EPmain.preferences.general.importFormat=1; %EP format
end
if ~isfield(EPmain.preferences.general,'outputFormat')
    EPmain.preferences.general.outputFormat=1; %EP format
end
if ~isfield(EPmain.preferences.general,'firstRow')
    EPmain.preferences.general.firstRow=1;
end
if ~isfield(EPmain.preferences.general,'lastRow')
    EPmain.preferences.general.lastRow=0;
end
if ~isfield(EPmain.preferences.general,'firstCol')
    EPmain.preferences.general.firstCol=1;
end
if ~isfield(EPmain.preferences.general,'lastCol')
    EPmain.preferences.general.lastCol=0;
end
if ~isfield(EPmain.preferences.general,'sampleRate')
    EPmain.preferences.general.sampleRate=250;
end
if ~isfield(EPmain.preferences.general,'SMIsuffix')
    EPmain.preferences.general.SMIsuffix='_smi.txt';
end
if ~isfield(EPmain.preferences.general,'segSuffix')
    EPmain.preferences.general.segSuffix='_seg';
end
if ~isfield(EPmain.preferences.general,'specSuffix')
    EPmain.preferences.general.specSuffix='_evt.txt';
end
if ~isfield(EPmain.preferences.general,'subjectSpecSuffix')
    EPmain.preferences.general.subjectSpecSuffix='_sub.txt';
end
if ~isfield(EPmain.preferences.general,'orientation')
    EPmain.preferences.general.orientation=1;
end
if ~isfield(EPmain.preferences.general,'defaultMontage')
    EPmain.preferences.general.defaultMontage='Generic 10-05';
end
if ~isfield(EPmain.preferences.general,'BVheader')
    EPmain.preferences.general.BVheader=1;
end
if ~isfield(EPmain.preferences.general,'numEEG')
    EPmain.preferences.general.numEEG=10;
end
if ~isfield(EPmain.preferences.general,'noInternal')
    EPmain.preferences.general.noInternal=0;
end

if ~isfield(EPmain.preferences,'preprocess')
    EPmain.preferences.preprocess=[];
end
if ~isfield(EPmain.preferences.preprocess,'noFigure')
    EPmain.preferences.preprocess.noFigure=0;
end
if ~isfield(EPmain.preferences.preprocess,'saturation')
    EPmain.preferences.preprocess.saturation=1000;
end
if ~isfield(EPmain.preferences.preprocess,'window')
    EPmain.preferences.preprocess.window=80;
end
if ~isfield(EPmain.preferences.preprocess,'minmax')
    EPmain.preferences.preprocess.minmax=100;
end
if ~isfield(EPmain.preferences.preprocess,'badnum')
    EPmain.preferences.preprocess.badnum=10;
end
if ~isfield(EPmain.preferences.preprocess,'neighbors')
    EPmain.preferences.preprocess.neighbors=6;
end
if ~isfield(EPmain.preferences.preprocess,'maxneighbor')
    EPmain.preferences.preprocess.maxneighbor=30;
end
if ~isfield(EPmain.preferences.preprocess,'badchan')
    EPmain.preferences.preprocess.badchan=.4;
end
if ~isfield(EPmain.preferences.preprocess,'blink')
    EPmain.preferences.preprocess.blink=.9;
end
if ~isfield(EPmain.preferences.preprocess,'badtrials')
    EPmain.preferences.preprocess.badtrials=20;
end
if ~isfield(EPmain.preferences.preprocess,'chunkSize')
    EPmain.preferences.preprocess.chunkSize=20000000;
end
if ~isfield(EPmain.preferences.preprocess,'minTrialsPerCell')
    EPmain.preferences.preprocess.minTrialsPerCell=15;
end
if ~isfield(EPmain.preferences.preprocess,'noadjacent')
    EPmain.preferences.preprocess.noadjacent=0;
end
if ~isfield(EPmain.preferences.preprocess,'trialminmax')
    EPmain.preferences.preprocess.trialminmax=200;
end
if ~isfield(EPmain.preferences.preprocess,'movefacs')
    EPmain.preferences.preprocess.movefacs=20;
end
if ~isfield(EPmain.preferences.preprocess,'EOGchans')
    EPmain.preferences.preprocess.EOGchans=[];
end
if ~isfield(EPmain.preferences.preprocess,'sacPot')
    EPmain.preferences.preprocess.sacPot=2;
end
if ~isfield(EPmain.preferences.preprocess,'fMRI')
    EPmain.preferences.preprocess.fMRI=1;
end
if ~isfield(EPmain.preferences.preprocess,'blinkRotation')
    EPmain.preferences.preprocess.blinkRotation='Infomax';
end
if ~isfield(EPmain.preferences.preprocess,'saccRotation')
    EPmain.preferences.preprocess.saccRotation='regression';
end
if ~isfield(EPmain.preferences.preprocess,'SProtation')
    EPmain.preferences.preprocess.SProtation='vector';
end
if ~isfield(EPmain.preferences.preprocess,'EMGratio')
    EPmain.preferences.preprocess.EMGratio=9;
end
if ~isfield(EPmain.preferences.preprocess,'EMGthresh')
    EPmain.preferences.preprocess.EMGthresh=15;   
end

if ~isfield(EPmain.preferences,'average')
    EPmain.preferences.average=[];
end
if ~isfield(EPmain.preferences.average,'method')
    EPmain.preferences.average.method=1;  
end
if ~isfield(EPmain.preferences.average,'trimLevel')
    EPmain.preferences.average.trimLevel=.25;  
end
if ~isfield(EPmain.preferences.average,'codeCorrect')
    EPmain.preferences.average.codeCorrect='1';  
end
if ~isfield(EPmain.preferences.average,'codeError')
    EPmain.preferences.average.codeError='0';  
end
if ~isfield(EPmain.preferences.average,'codeTimeout')
    EPmain.preferences.average.codeTimeout='2';  
end

if ~isfield(EPmain.preferences,'transform')
    EPmain.preferences.transform=[];
end
if ~isfield(EPmain.preferences.transform,'reference')
    EPmain.preferences.transform.reference=1; %average reference
end
if ~isfield(EPmain.preferences.transform,'refChan1')
    EPmain.preferences.transform.refChan1=57; 
end
if ~isfield(EPmain.preferences.transform,'refChan2')
    EPmain.preferences.transform.refChan2=100; 
end
if ~isfield(EPmain.preferences.transform,'baselineStart')
    EPmain.preferences.transform.baselineStart=-200;
end
if ~isfield(EPmain.preferences.transform,'baselineEnd')
    EPmain.preferences.transform.baselineEnd=0; 
end
if ~isfield(EPmain.preferences.transform,'mainsFix')
    EPmain.preferences.transform.mainsFix='None'; 
end

if ~isfield(EPmain.preferences,'view')
    EPmain.preferences.view=[];
end
if ~isfield(EPmain.preferences.view,'positive')
    EPmain.preferences.view.positive=1;
end
if ~isfield(EPmain.preferences.view,'lineSize')
    EPmain.preferences.view.lineSize=1;
end
if ~isfield(EPmain.preferences.view,'labelSize')
    EPmain.preferences.view.labelSize=12;
end
if ~isfield(EPmain.preferences.view,'topoElectrodes')
    EPmain.preferences.view.topoElectrodes='on';
end
if ~isfield(EPmain.preferences.view,'color')
    EPmain.preferences.view.color=[];
end
if ~isfield(EPmain.preferences.view.color,'RGB')
    for iColor=1:EPmain.maxColors
        EPmain.preferences.view.color(iColor).RGB=EPmain.defaultColor(iColor).RGB;
    end
end
if ~isfield(EPmain.preferences.view.color,'lineSize')
    for iColor=1:EPmain.maxColors
        EPmain.preferences.view.color(iColor).lineSize=2;
    end
end
if ~isfield(EPmain.preferences.view.color,'lineStyle')
    for iColor=1:EPmain.maxColors
        EPmain.preferences.view.color(iColor).lineStyle=EPmain.lineStyleList{1};
    end
end

if ~isfield(EPmain.preferences,'pca')
    EPmain.preferences.pca=[];
end
if ~isfield(EPmain.preferences.pca,'mode')
    EPmain.preferences.pca.mode=2;
end
if ~isfield(EPmain.preferences.pca,'rotation')
    EPmain.preferences.pca.rotation=2;
end
if ~isfield(EPmain.preferences.pca,'rotopt')
    EPmain.preferences.pca.rotopt=3;
end
if ~isfield(EPmain.preferences.pca,'rel')
    EPmain.preferences.pca.rel=2;
end
if ~isfield(EPmain.preferences.pca,'loadings')
    EPmain.preferences.pca.loadings=2;
end

if ~isfield(EPmain.preferences,'window')
    EPmain.preferences.window=[];
end
if ~isfield(EPmain.preferences.window,'minFacVar')
    EPmain.preferences.window.minFacVar=.005;
end
if ~isfield(EPmain.preferences.window,'adds')
    EPmain.preferences.window.adds=1;
end
if ~isfield(EPmain.preferences.window,'chanGrp')
    EPmain.preferences.window.chanGrp=1;
end

if ~isfield(EPmain.preferences,'anova')
    EPmain.preferences.anova=[];
end
if ~isfield(EPmain.preferences.anova,'trimming')
    EPmain.preferences.anova.trimming=.05;
end
if ~isfield(EPmain.preferences.anova,'missing')
    EPmain.preferences.anova.missing=-999;
end
if ~isfield(EPmain.preferences.anova,'seed')
    EPmain.preferences.anova.seed=1000;
end
if ~isfield(EPmain.preferences.anova,'bootstrap')
    EPmain.preferences.anova.bootstrap=499999;
end
if ~isfield(EPmain.preferences.anova,'reps')
    EPmain.preferences.anova.reps=11;
end
if ~isfield(EPmain.preferences.anova,'adds')
    EPmain.preferences.anova.adds=1;
end
if ~isfield(EPmain.preferences.anova,'epsilon')
    theString=ep_enumVar('epsilon');
    EPmain.preferences.anova.epsilon=theString{1};
end
if ~isfield(EPmain.preferences.anova,'posthoc')
    theString=ep_enumVar('postHoc');
    EPmain.preferences.anova.posthoc=theString{1};
end
% if ~isfield(EPmain.preferences.anova,'adf')
%     EPmain.preferences.anova.adf='ADF';
% end

if ~isfield(EPmain.preferences,'advanced')
    EPmain.preferences.advanced=[];
end
if ~isfield(EPmain.preferences.advanced,'parallel')
    EPmain.preferences.advanced.parallel=1;
end
if ~isfield(EPmain.preferences.advanced,'monitor')
    EPmain.preferences.advanced.monitor=[0 0 0 0];
end
if ~isfield(EPmain.preferences.advanced,'RAM')
    EPmain.preferences.advanced.RAM=0;
end

if ~isfield(EPmain.preferences,'records')
    EPmain.preferences.records=[];
end
if ~isfield(EPmain.preferences.records,'user')
    EPmain.preferences.records.user='';
end
if ~isfield(EPmain.preferences.records,'lab')
    EPmain.preferences.records.lab='';
end
if ~isfield(EPmain.preferences.records,'institution')
    EPmain.preferences.records.institution='';
end
if ~isfield(EPmain.preferences.records,'project')
    EPmain.preferences.records.project='';
end
if ~isfield(EPmain.preferences.records,'experiment')
    EPmain.preferences.records.experiment='';
end

EPmain.preferences.EPver=EP_VER;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [err prefMissing]=checkPrefs(src,eventdata) %check the preference settings
global EPmain EP_VER

err=[];
prefMissing=[];

if ~isfield(EPmain.preferences,'EPver')
    disp(['Updating preferences to version: ' EP_VER]);
    err = [err ;100];
end    

if isfield(EPmain.preferences,'general')
    if ~isfield(EPmain.preferences.general,'sessionImportFormat')
        prefMissing = [prefMissing ;14];
        disp('The EPmain.preferences.general.sessionImportFormat field is missing.');
    elseif isempty(EPmain.preferences.general.sessionImportFormat)
        err = [err ;14];
        disp('The EPmain.preferences.general.sessionImportFormat field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.sessionImportFormat)
        err = [err ;14];
        disp('The EPmain.preferences.general.sessionImportFormat field is not a number.');
    elseif EPmain.preferences.general.sessionImportFormat < 1 || EPmain.preferences.general.sessionImportFormat > length(EPmain.fileFormatReadList)
        err = [err ;14];
        disp(['The EPmain.preferences.general.sessionImportFormat field is not within the range of 1 to ' num2str(length(EPmain.fileFormatSaveList)) '.']);
    end

    if ~isfield(EPmain.preferences.general,'sessionOutputFormat')
        prefMissing = [prefMissing ;15];
        disp('The EPmain.preferences.general.sessionOutputFormat field is missing.');
    elseif isempty(EPmain.preferences.general.sessionOutputFormat)
        err = [err ;15];
        disp('The EPmain.preferences.general.sessionOutputFormat field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.sessionOutputFormat)
        err = [err ;15];
        disp('The EPmain.preferences.general.sessionOutputFormat field is not a number.');
    elseif EPmain.preferences.general.sessionOutputFormat < 1 || EPmain.preferences.general.sessionOutputFormat > length(EPmain.fileFormatSaveList)
        err = [err ;15];
        disp(['The EPmain.preferences.general.sessionOutputFormat field is not within the range of 1 to ' num2str(length(EPmain.fileFormatSaveList)) '.']);
    end
    if ~isfield(EPmain.preferences.general,'importFormat')
        prefMissing = [prefMissing ;19];
        disp('The EPmain.preferences.general.importFormat field is missing.');
    elseif isempty(EPmain.preferences.general.importFormat)
        err = [err ;19];
        disp('The EPmain.preferences.general.importFormat field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.importFormat)
        err = [err ;19];
        disp('The EPmain.preferences.general.importFormat field is not a number.');
    elseif EPmain.preferences.general.importFormat < 1 || EPmain.preferences.general.importFormat > length(EPmain.fileFormatReadList)
        err = [err ;19];
        disp(['The EPmain.preferences.general.importFormat field is not within the range of 1 to ' num2str(length(EPmain.fileFormatSaveList)) '.']);
    end
    
    if ~isfield(EPmain.preferences.general,'outputFormat')
        prefMissing = [prefMissing ;20];
        disp('The EPmain.preferences.general.outputFormat field is missing.');
    elseif isempty(EPmain.preferences.general.outputFormat)
        err = [err ;20];
        disp('The EPmain.preferences.general.outputFormat field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.outputFormat)
        err = [err ;20];
        disp('The EPmain.preferences.general.outputFormat field is not a number.');
    elseif EPmain.preferences.general.outputFormat < 1 || EPmain.preferences.general.outputFormat > length(EPmain.fileFormatSaveList)
        err = [err ;20];
        disp(['The EPmain.preferences.general.outputFormat field is not within the range of ' num2str(length(EPmain.fileFormatSaveList)) '.']);
    end
    
    if ~isfield(EPmain.preferences.general,'firstRow')
        prefMissing = [prefMissing ;47];
        disp('The EPmain.preferences.general.firstRow field is missing.');
    elseif isempty(EPmain.preferences.general.firstRow)
        err = [err ;47];
        disp('The EPmain.preferences.general.firstRow field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.firstRow)
        err = [err ;47];
        disp('The EPmain.preferences.general.firstRow field is not a number.');
    elseif EPmain.preferences.general.firstRow < 1
        err = [err ;47];
        disp('The EPmain.preferences.general.firstRow field is less than 1.');
    end
    
    if ~isfield(EPmain.preferences.general,'lastRow')
        prefMissing = [prefMissing ;58];
        disp('The EPmain.preferences.general.lastRow field is missing.');
    elseif isempty(EPmain.preferences.general.lastRow)
        err = [err ;58];
        disp('The EPmain.preferences.general.lastRow field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.lastRow)
        err = [err ;58];
        disp('The EPmain.preferences.general.lastRow field is not a number.');
    elseif EPmain.preferences.general.lastRow < 0
        err = [err ;58];
        disp('The EPmain.preferences.general.lastRow field is less than zero.');
    end
    
    if ~isfield(EPmain.preferences.general,'firstCol')
        prefMissing = [prefMissing ;48];
        disp('The EPmain.preferences.general.firstCol field is missing.');
    elseif isempty(EPmain.preferences.general.firstCol)
        err = [err ;48];
        disp('The EPmain.preferences.general.firstCol field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.firstCol)
        err = [err ;48];
        disp('The EPmain.preferences.general.firstCol field is not a number.');
    elseif EPmain.preferences.general.firstCol < 1
        err = [err ;48];
        disp('The EPmain.preferences.general.firstCol field is less than 1.');
    end
    
    if ~isfield(EPmain.preferences.general,'lastCol')
        prefMissing = [prefMissing ;49];
        disp('The EPmain.preferences.general.lastCol field is missing.');
    elseif isempty(EPmain.preferences.general.lastCol)
        err = [err ;49];
        disp('The EPmain.preferences.general.lastCol field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.lastCol)
        err = [err ;49];
        disp('The EPmain.preferences.general.lastCol field is not a number.');
    elseif EPmain.preferences.general.lastCol < 0
        err = [err ;49];
        disp('The EPmain.preferences.general.lastCol field is less than zero.');
    end
    
    if ~isfield(EPmain.preferences.general,'sampleRate')
        prefMissing = [prefMissing ;63];
        disp('The EPmain.preferences.general.sampleRate field is missing.');
    elseif isempty(EPmain.preferences.general.sampleRate)
        err = [err ;63];
        disp('The EPmain.preferences.general.sampleRate field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.sampleRate)
        err = [err ;63];
        disp('The EPmain.preferences.general.sampleRate field is not a number.');
    elseif EPmain.preferences.general.sampleRate < 0
        err = [err ;63];
        disp('The EPmain.preferences.general.sampleRate field is less than zero.');
    end

    if ~isfield(EPmain.preferences.general,'segSuffix')
        prefMissing = [prefMissing ;60];
        disp('The EPmain.preferences.general.segSuffix field is missing.');
    elseif isempty(EPmain.preferences.general.segSuffix)
        err = [err ;60];
        disp('The EPmain.preferences.general.segSuffix field is empty.');
    elseif isnumeric(EPmain.preferences.general.segSuffix)
        err = [err ;60];
        disp('The EPmain.preferences.general.segSuffix field is a number.');
    elseif ~strcmp(EPmain.preferences.general.segSuffix(1),'_')
        err = [err ;60];
        disp('The first character of the EPmain.preferences.general.segSuffix field is not an underscore.');
    end
    
    if ~isfield(EPmain.preferences.general,'specSuffix')
        prefMissing = [prefMissing ;61];
        disp('The EPmain.preferences.general.specSuffix field is missing.');
    elseif isempty(EPmain.preferences.general.specSuffix)
        err = [err ;61];
        disp('The EPmain.preferences.general.specSuffix field is empty.');
    elseif isnumeric(EPmain.preferences.general.specSuffix)
        err = [err ;61];
        disp('The EPmain.preferences.general.specSuffix field is a number.');
    elseif ~strcmp(EPmain.preferences.general.specSuffix(1),'_')
        err = [err ;61];
        disp('The first character of the EPmain.preferences.general.specSuffix field is not an underscore.');
    end
    
    if ~isfield(EPmain.preferences.general,'subjectSpecSuffix')
        prefMissing = [prefMissing ;62];
        disp('The EPmain.preferences.general.subjectSpecSuffix field is missing.');
    elseif isempty(EPmain.preferences.general.subjectSpecSuffix)
        err = [err ;62];
        disp('The EPmain.preferences.general.subjectSpecSuffix field is empty.');
    elseif isnumeric(EPmain.preferences.general.subjectSpecSuffix)
        err = [err ;62];
        disp('The EPmain.preferences.general.subjectSpecSuffix field is a number.');
    elseif ~strcmp(EPmain.preferences.general.subjectSpecSuffix(1),'_')
        err = [err ;62];
        disp('The first character of the EPmain.preferences.general.subjectSpecSuffix field is not an underscore.');
    end

    if ~isfield(EPmain.preferences.general,'orientation')
        prefMissing = [prefMissing ;50];
        disp('The EPmain.preferences.general.orientation field is missing.');
    elseif isempty(EPmain.preferences.general.orientation)
        err = [err ;50];
        disp('The EPmain.preferences.general.orientation field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.orientation)
        err = [err ;50];
        disp('The EPmain.preferences.general.orientation field is not a number.');
    elseif EPmain.preferences.general.orientation < 1 || EPmain.preferences.general.orientation > 2
        err = [err ;50];
        disp('The EPmain.preferences.general.orientation field is not within the range of 1 to 2.');
    end
    
    if ~isfield(EPmain.preferences.general,'defaultMontage')
        prefMissing = [prefMissing ;83];
        disp('The EPmain.preferences.general.defaultMontage field is missing.');
    elseif isempty(EPmain.preferences.general.defaultMontage)
        err = [err ;83];
        disp('The EPmain.preferences.general.defaultMontage field is empty.');
    elseif isnumeric(EPmain.preferences.general.defaultMontage)
        err = [err ;83];
        disp('The EPmain.preferences.general.defaultMontage field is a number.');
    end
    
    if ~isfield(EPmain.preferences.general,'BVheader')
        prefMissing = [prefMissing ;54];
        disp('The EPmain.preferences.general.BVheader field is missing.');
    elseif isempty(EPmain.preferences.general.BVheader)
        err = [err ;54];
        disp('The EPmain.preferences.general.BVheader field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.BVheader)
        err = [err ;54];
        disp('The EPmain.preferences.general.BVheader field is not a number.');
    elseif EPmain.preferences.general.BVheader ~= 0 && EPmain.preferences.general.BVheader ~= 1
        err = [err ;54];
        disp('The EPmain.preferences.general.BVheader field does not equal zero or one.');
    end

    if ~isfield(EPmain.preferences.general,'numEEG')
        prefMissing = [prefMissing ;64];
        disp('The EPmain.preferences.general.numEEG field is missing.');
    elseif isempty(EPmain.preferences.general.numEEG)
        err = [err ;64];
        disp('The EPmain.preferences.general.numEEG field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.numEEG)
        err = [err ;64];
        disp('The EPmain.preferences.general.numEEG field is not a number.');
    elseif EPmain.preferences.general.numEEG < 0
        err = [err ;64];
        disp('The EPmain.preferences.general.numEEG field is a negative number.');
    end
    
    if ~isfield(EPmain.preferences.general,'noInternal')
        prefMissing = [prefMissing ;74];
        disp('The EPmain.preferences.general.noInternal field is missing.');
    elseif isempty(EPmain.preferences.general.noInternal)
        err = [err ;74];
        disp('The EPmain.preferences.general.noInternal field is empty.');
    elseif ~isnumeric(EPmain.preferences.general.noInternal)
        err = [err ;74];
        disp('The EPmain.preferences.general.noInternal field is not a number.');
    elseif EPmain.preferences.general.noInternal ~= 0 && EPmain.preferences.general.noInternal ~= 1
        err = [err ;74];
        disp('The EPmain.preferences.general.noInternal field does not equal zero or one.');
    end
    
else
    prefMissing = [prefMissing ;17];
    disp('The EPmain.preferences.general field is missing.');
end

if isfield(EPmain.preferences,'preprocess')
    if ~isfield(EPmain.preferences.preprocess,'saturation')
        prefMissing = [prefMissing ;1];
        disp('The EPmain.preferences.preprocess.saturation field is missing.');
    elseif isempty(EPmain.preferences.preprocess.saturation)
        err = [err ;1];
        disp('The EPmain.preferences.preprocess.saturation field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.saturation)
        err = [err ;1];
        disp('The EPmain.preferences.preprocess.saturation field is not a number.');
    elseif EPmain.preferences.preprocess.saturation == 0
        err = [err ;1];
        disp('The EPmain.preferences.preprocess.saturation field is zero.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'window')
        prefMissing = [prefMissing ;2];
        disp('The EPmain.preferences.preprocess.window field is missing.');
    elseif isempty(EPmain.preferences.preprocess.window)
        err = [err ;2];
        disp('The EPmain.preferences.preprocess.window field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.window)
        err = [err ;2];
        disp('The EPmain.preferences.preprocess.window field is not a number.');
    elseif EPmain.preferences.preprocess.window < 1
        err = [err ;2];
        disp('The EPmain.preferences.preprocess.window field is less than one.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'minmax')
        prefMissing = [prefMissing ;3];
        disp('The EPmain.preferences.preprocess.minmax field is missing.');
    elseif isempty(EPmain.preferences.preprocess.minmax)
        err = [err ;3];
        disp('The EPmain.preferences.preprocess.minmax field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.minmax)
        err = [err ;3];
        disp('The EPmain.preferences.preprocess.minmax field is not a number.');
    elseif EPmain.preferences.preprocess.saturation < 1
        err = [err ;3];
        disp('The EPmain.preferences.preprocess.minmax field is less than one.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'badnum')
        prefMissing = [prefMissing ;4];
        disp('The EPmain.preferences.preprocess.badnum field is missing.');
    elseif isempty(EPmain.preferences.preprocess.badnum)
        err = [err ;4];
        disp('The EPmain.preferences.preprocess.badnum field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.badnum)
        err = [err ;4];
        disp('The EPmain.preferences.preprocess.badnum field is not a number.');
    elseif EPmain.preferences.preprocess.badnum < 0 || EPmain.preferences.preprocess.badnum > 100
        err = [err ;4];
        disp('The EPmain.preferences.preprocess.badnum field is not between 0 and 100.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'neighbors')
        prefMissing = [prefMissing ;6];
        disp('The EPmain.preferences.preprocess.neighbors field is missing.');
    elseif isempty(EPmain.preferences.preprocess.neighbors)
        err = [err ;6];
        disp('The EPmain.preferences.preprocess.neighbors field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.neighbors)
        err = [err ;6];
        disp('The EPmain.preferences.preprocess.neighbors field is not a number.');
    elseif EPmain.preferences.preprocess.neighbors < 0
        err = [err ;6];
        disp('The EPmain.preferences.preprocess.neighbors field is less than zero.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'maxneighbor')
        prefMissing = [prefMissing ;7];
        disp('The EPmain.preferences.preprocess.maxneighbor field is missing.');
    elseif isempty(EPmain.preferences.preprocess.maxneighbor)
        err = [err ;7];
        disp('The EPmain.preferences.preprocess.maxneighbor field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.maxneighbor)
        err = [err ;7];
        disp('The EPmain.preferences.preprocess.maxneighbor field is not a number.');
    elseif EPmain.preferences.preprocess.maxneighbor < 0
        err = [err ;7];
        disp('The EPmain.preferences.preprocess.maxneighbor field is less than zero.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'badchan')
        prefMissing = [prefMissing ;8];
        disp('The EPmain.preferences.preprocess.badchan field is missing.');
    elseif isempty(EPmain.preferences.preprocess.badchan)
        err = [err ;8];
        disp('The EPmain.preferences.preprocess.badchan field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.badchan)
        err = [err ;8];
        disp('The EPmain.preferences.preprocess.badchan field is not a number.');
    elseif EPmain.preferences.preprocess.badchan < 0 || EPmain.preferences.preprocess.badchan > 1
        err = [err ;8];
        disp('The EPmain.preferences.preprocess.badchan field is not between 0 and 1.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'blink')
        prefMissing = [prefMissing ;9];
        disp('The EPmain.preferences.preprocess.blink field is missing.');
    elseif isempty(EPmain.preferences.preprocess.blink)
        err = [err ;9];
        disp('The EPmain.preferences.preprocess.blink field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.blink)
        err = [err ;9];
        disp('The EPmain.preferences.preprocess.blink field is not a number.');
    elseif EPmain.preferences.preprocess.blink < 0 || EPmain.preferences.preprocess.blink > 1
        err = [err ;9];
        disp('The EPmain.preferences.preprocess.blink field is not between 0 and 1.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'badtrials')
        prefMissing = [prefMissing ;10];
        disp('The EPmain.preferences.preprocess.badtrials field is missing.');
    elseif isempty(EPmain.preferences.preprocess.badtrials)
        err = [err ;10];
        disp('The EPmain.preferences.preprocess.badtrials field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.badtrials)
        err = [err ;10];
        disp('The EPmain.preferences.preprocess.badtrials field is not a number.');
    elseif EPmain.preferences.preprocess.badtrials < 0 || EPmain.preferences.preprocess.badtrials > 100
        err = [err ;10];
        disp('The EPmain.preferences.preprocess.badtrials field is not between 0 and 100.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'chunkSize')
        prefMissing = [prefMissing ;11];
        disp('The EPmain.preferences.preprocess.chunkSize field is missing.');
    elseif isempty(EPmain.preferences.preprocess.chunkSize)
        err = [err ;11];
        disp('The EPmain.preferences.preprocess.chunkSize field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.chunkSize)
        err = [err ;11];
        disp('The EPmain.preferences.preprocess.chunkSize field is not a number.');
    elseif EPmain.preferences.preprocess.chunkSize < 1
        err = [err ;11];
        disp('The EPmain.preferences.preprocess.chunkSize field is less than one.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'minTrialsPerCell')
        prefMissing = [prefMissing ;12];
        disp('The EPmain.preferences.preprocess.minTrialsPerCell field is missing.');
    elseif isempty(EPmain.preferences.preprocess.minTrialsPerCell)
        err = [err ;12];
        disp('The EPmain.preferences.preprocess.minTrialsPerCell field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.minTrialsPerCell)
        err = [err ;12];
        disp('The EPmain.preferences.preprocess.minTrialsPerCell field is not a number.');
    elseif EPmain.preferences.preprocess.minTrialsPerCell < 0
        err = [err ;12];
        disp('The EPmain.preferences.preprocess.minTrialsPerCell field is less than zero.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'noadjacent')
        prefMissing = [prefMissing ;13];
        disp('The EPmain.preferences.preprocess.noadjacent field is missing.');
    elseif isempty(EPmain.preferences.preprocess.noadjacent)
        err = [err ;13];
        disp('The EPmain.preferences.preprocess.noadjacent field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.noadjacent)
        err = [err ;13];
        disp('The EPmain.preferences.preprocess.noadjacent field is not a number.');
    elseif EPmain.preferences.preprocess.noadjacent ~= 0 && EPmain.preferences.preprocess.noadjacent ~= 1
        err = [err ;13];
        disp('The EPmain.preferences.preprocess.noadjacent field does not equal zero or one.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'trialminmax')
        prefMissing = [prefMissing ;18];
        disp('The EPmain.preferences.preprocess.trialminmax field is missing.');
    elseif isempty(EPmain.preferences.preprocess.trialminmax)
        err = [err ;18];
        disp('The EPmain.preferences.preprocess.trialminmax field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.trialminmax)
        err = [err ;18];
        disp('The EPmain.preferences.preprocess.trialminmax field is not a number.');
    elseif EPmain.preferences.preprocess.trialminmax < 0
        err = [err ;18];
        disp('The EPmain.preferences.preprocess.trialminmax field is less than zero.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'movefacs')
        prefMissing = [prefMissing ;21];
        disp('The EPmain.preferences.preprocess.movefacs field is missing.');
    elseif isempty(EPmain.preferences.preprocess.movefacs)
        err = [err ;21];
        disp('The EPmain.preferences.preprocess.movefacs field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.movefacs)
        err = [err ;21];
        disp('The EPmain.preferences.preprocess.movefacs field is not a number.');
    elseif EPmain.preferences.preprocess.movefacs < 0
        err = [err ;21];
        disp('The EPmain.preferences.preprocess.movefacs field is less than zero.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'noFigure')
        prefMissing = [prefMissing ;51];
        disp('The EPmain.preferences.preprocess.noFigure field is missing.');
    elseif isempty(EPmain.preferences.preprocess.noFigure)
        err = [err ;51];
        disp('The EPmain.preferences.preprocess.noFigure field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.noFigure)
        err = [err ;51];
        disp('The EPmain.preferences.preprocess.noFigure field is not a number.');
    elseif EPmain.preferences.preprocess.noFigure ~= 0 && EPmain.preferences.preprocess.noFigure ~= 1
        err = [err ;51];
        disp('The EPmain.preferences.preprocess.noFigure field does not equal zero or one.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'EMGratio')
        prefMissing = [prefMissing ;53];
        disp('The EPmain.preferences.preprocess.EMGratio field is missing.');
    elseif isempty(EPmain.preferences.preprocess.EMGratio)
        err = [err ;53];
        disp('The EPmain.preferences.preprocess.EMGratio field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.EMGratio)
        err = [err ;53];
        disp('The EPmain.preferences.preprocess.EMGratio field is not a number.');
    elseif EPmain.preferences.preprocess.EMGratio < 0
        err = [err ;53];
        disp('The EPmain.preferences.preprocess.EMGratio field is negative.');
    end
    
    if ~isfield(EPmain.preferences.preprocess,'EMGthresh')
        prefMissing = [prefMissing ;54];
        disp('The EPmain.preferences.preprocess.EMGthresh field is missing.');
    elseif isempty(EPmain.preferences.preprocess.EMGthresh)
        err = [err ;54];
        disp('The EPmain.preferences.preprocess.EMGthresh field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.EMGthresh)
        err = [err ;54];
        disp('The EPmain.preferences.preprocess.EMGthresh field is not a number.');
    elseif EPmain.preferences.preprocess.EMGthresh < 0
        err = [err ;54];
        disp('The EPmain.preferences.preprocess.EMGthresh field is less than zero.');
    end

    if ~isfield(EPmain.preferences.preprocess,'EOGchans')
        prefMissing = [prefMissing ;55];
        disp('The EPmain.preferences.preprocess.EOGchans field is missing.');
    elseif ~isnumeric(EPmain.preferences.preprocess.EOGchans)
        err = [err ;55];
        disp('The EPmain.preferences.preprocess.EOGchans field is not a number.');
    elseif length(EPmain.preferences.preprocess.EOGchans) ~= 6 && ~isempty(EPmain.preferences.preprocess.EOGchans)
        err = [err ;55];
        disp('The EPmain.preferences.preprocess.EOGchans field does not have six numbers in it.');
    end

    if ~isfield(EPmain.preferences.preprocess,'sacPot')
        prefMissing = [prefMissing ;63];
        disp('The EPmain.preferences.preprocess.sacPot field is missing.');
    elseif isempty(EPmain.preferences.preprocess.sacPot)
        err = [err ;63];
        disp('The EPmain.preferences.preprocess.sacPot field is empty.');
    elseif ~isnumeric(EPmain.preferences.preprocess.sacPot)
        err = [err ;63];
        disp('The EPmain.preferences.preprocess.sacPot field is not a number.');
    elseif EPmain.preferences.preprocess.sacPot < 0
        err = [err ;63];
        disp('The EPmain.preferences.preprocess.sacPot field is less than zero.');
    end

    if ~isfield(EPmain.preferences.preprocess,'fMRI')
        prefMissing = [prefMissing ;58];
        disp('The EPmain.preferences.preprocess.fMRI field is missing.');
    elseif ~isnumeric(EPmain.preferences.preprocess.fMRI)
        err = [err ;58];
        disp('The EPmain.preferences.preprocess.fMRI field is not a number.');
    elseif EPmain.preferences.preprocess.fMRI < 1 || EPmain.preferences.preprocess.fMRI > 2
        err = [err ;58];
        disp('The EPmain.preferences.preprocess.fMRI field is not between 1 and 2.');
    end

    rotationList=ep_doPCA;
    if ~isfield(EPmain.preferences.preprocess,'blinkRotation')
        prefMissing = [prefMissing ;67];
        disp('The EPmain.preferences.preprocess.blinkRotation field is missing.');
    elseif ~ischar(EPmain.preferences.preprocess.blinkRotation)
        err = [err ;67];
        disp('The EPmain.preferences.preprocess.blinkRotation field is not a string.');
    elseif ~any(strcmp(EPmain.preferences.preprocess.blinkRotation,rotationList))
        err = [err ;67];
        disp('The EPmain.preferences.preprocess.blinkRotation field is not a recognized rotation name.');
    end
    rotationList(2:end+1)=rotationList;
    rotationList{1}='regression';
    if ~isfield(EPmain.preferences.preprocess,'saccRotation')
        prefMissing = [prefMissing ;70];
        disp('The EPmain.preferences.preprocess.saccRotation field is missing.');
    elseif ~ischar(EPmain.preferences.preprocess.saccRotation)
        err = [err ;70];
        disp('The EPmain.preferences.preprocess.saccRotation field is not a string.');
    elseif ~any(strcmp(EPmain.preferences.preprocess.saccRotation,rotationList))
        err = [err ;70];
        disp('The EPmain.preferences.preprocess.saccRotation field is not a recognized method name.');
    end
    rotationList=ep_doPCA;
    rotationList(2:end+1)=rotationList;
    rotationList{1}='vector';
    if ~isfield(EPmain.preferences.preprocess,'SProtation')
        prefMissing = [prefMissing ;73];
        disp('The EPmain.preferences.preprocess.SProtation field is missing.');
    elseif ~ischar(EPmain.preferences.preprocess.SProtation)
        err = [err ;73];
        disp('The EPmain.preferences.preprocess.SProtation field is not a string.');
    elseif ~any(strcmp(EPmain.preferences.preprocess.SProtation,rotationList))
        err = [err ;73];
        disp('The EPmain.preferences.preprocess.SProtation field is not a recognized method name.');
    end
else
    prefMissing = [prefMissing ;16];
    disp('The EPmain.preferences.preprocess field is missing.');
end

if isfield(EPmain.preferences,'average')
    if ~isfield(EPmain.preferences.average,'method')
        prefMissing = [prefMissing ;23];
        disp('The EPmain.preferences.average.method field is missing.');
    elseif isempty(EPmain.preferences.average.method)
        err = [err ;23];
        disp('The EPmain.preferences.average.method field is empty.');
    elseif ~isnumeric(EPmain.preferences.average.method)
        err = [err ;23];
        disp('The EPmain.preferences.average.method field is not a number.');
    elseif EPmain.preferences.average.method < 1 && EPmain.preferences.average.method > 3
        err = [err ;23];
        disp('The EPmain.preferences.average.method field is not within the range of 1 to 3.');
    end
    if ~isfield(EPmain.preferences.average,'trimLevel')
        prefMissing = [prefMissing ;24];
        disp('The EPmain.preferences.average.trimLevel field is missing.');
    elseif isempty(EPmain.preferences.average.trimLevel)
        err = [err ;24];
        disp('The EPmain.preferences.average.trimLevel field is empty.');
    elseif ~isnumeric(EPmain.preferences.average.trimLevel)
        err = [err ;24];
        disp('The EPmain.preferences.average.trimLevel field is not a number.');
    elseif ~(EPmain.preferences.average.trimLevel < .5 && EPmain.preferences.average.trimLevel > 0)
        err = [err ;24];
        disp('The EPmain.preferences.average.trimLevel field is not between 0 and .5.');
    end
    if ~isfield(EPmain.preferences.average,'codeCorrect')
        prefMissing = [prefMissing ;80];
        disp('The EPmain.preferences.average.codeCorrect field is missing.');
    elseif isempty(EPmain.preferences.average.codeCorrect)
        err = [err ;80];
        disp('The EPmain.preferences.average.codeCorrect field is empty.');
    elseif isnumeric(EPmain.preferences.average.codeCorrect)
        err = [err ;80];
        disp('The EPmain.preferences.average.codeCorrect field is a number.');
    end
    if ~isfield(EPmain.preferences.average,'codeError')
        prefMissing = [prefMissing ;81];
        disp('The EPmain.preferences.average.codeError field is missing.');
    elseif isempty(EPmain.preferences.average.codeError)
        err = [err ;81];
        disp('The EPmain.preferences.average.codeError field is empty.');
    elseif isnumeric(EPmain.preferences.average.codeError)
        err = [err ;81];
        disp('The EPmain.preferences.average.codeError field is a number.');
    end
    if ~isfield(EPmain.preferences.average,'codeTimeout')
        prefMissing = [prefMissing ;82];
        disp('The EPmain.preferences.average.codeTimeout field is missing.');
    elseif isempty(EPmain.preferences.average.codeTimeout)
        err = [err ;82];
        disp('The EPmain.preferences.average.codeTimeout field is empty.');
    elseif isnumeric(EPmain.preferences.average.codeTimeout)
        err = [err ;82];
        disp('The EPmain.preferences.average.codeTimeout field is a number.');
    end
else
    prefMissing = [prefMissing ;22];
    disp('The EPmain.preferences.average field is missing.');
end

if isfield(EPmain.preferences,'transform')
    if ~isfield(EPmain.preferences.transform,'reference')
        prefMissing = [prefMissing ;26];
        disp('The EPmain.preferences.transform.reference field is missing.');
    elseif isempty(EPmain.preferences.transform.reference)
        err = [err ;26];
        disp('The EPmain.preferences.transform.reference field is empty.');
    elseif ~isnumeric(EPmain.preferences.transform.reference)
        err = [err ;26];
        disp('The EPmain.preferences.transform.reference field is not a number.');
    elseif EPmain.preferences.transform.reference < 1 && EPmain.transform.transform.reference > 3
        err = [err ;26];
        disp('The EPmain.preferences.transform.reference field is not within the range of 1 to 3.');
    end
    if ~isfield(EPmain.preferences.transform,'refChan1')
        prefMissing = [prefMissing ;27];
        disp('The EPmain.preferences.transform.refChan1 field is missing.');
    elseif isempty(EPmain.preferences.transform.refChan1)
        err = [err ;27];
        disp('The EPmain.preferences.transform.refChan1 field is empty.');
    elseif ~isnumeric(EPmain.preferences.transform.refChan1)
        err = [err ;27];
        disp('The EPmain.preferences.transform.refChan1 field is not a number.');
    end
    if ~isfield(EPmain.preferences.transform,'refChan2')
        prefMissing = [prefMissing ;28];
        disp('The EPmain.preferences.transform.refChan2 field is missing.');
    elseif isempty(EPmain.preferences.transform.refChan2)
        err = [err ;28];
        disp('The EPmain.preferences.transform.refChan2 field is empty.');
    elseif ~isnumeric(EPmain.preferences.transform.refChan2)
        err = [err ;28];
        disp('The EPmain.preferences.transform.refChan2 field is not a number.');
    end
    if ~isfield(EPmain.preferences.transform,'baselineStart')
        prefMissing = [prefMissing ;29];
        disp('The EPmain.preferences.transform.baselineStart field is missing.');
    elseif isempty(EPmain.preferences.transform.baselineStart)
        err = [err ;29];
        disp('The EPmain.preferences.transform.baselineStart field is empty.');
    elseif ~isnumeric(EPmain.preferences.transform.baselineStart)
        err = [err ;29];
        disp('The EPmain.preferences.transform.baselineStart field is not a number.');
    end
    if ~isfield(EPmain.preferences.transform,'baselineEnd')
        prefMissing = [prefMissing ;56];
        disp('The EPmain.preferences.transform.baselineEnd field is missing.');
    elseif isempty(EPmain.preferences.transform.baselineEnd)
        err = [err ;56];
        disp('The EPmain.preferences.transform.baselineEnd field is empty.');
    elseif ~isnumeric(EPmain.preferences.transform.baselineEnd)
        err = [err ;56];
        disp('The EPmain.preferences.transform.baselineEnd field is not a number.');
    end
    if ~isfield(EPmain.preferences.transform,'mainsFix')
        prefMissing = [prefMissing ;80];
        disp('The EPmain.preferences.transform.mainsFix field is missing.');
    elseif isempty(EPmain.preferences.transform.mainsFix)
        err = [err ;80];
        disp('The EPmain.preferences.transform.mainsFix field is empty.');
    elseif ~ischar(EPmain.preferences.transform.mainsFix)
        err = [err ;80];
        disp('The EPmain.preferences.transform.mainsFix field is not a character string.');
    end
else
    prefMissing = [prefMissing ;25];
    disp('The EPmain.preferences.transform field is missing.');
end

if isfield(EPmain.preferences,'view')
    if ~isfield(EPmain.preferences.view,'positive')
        prefMissing = [prefMissing ;37];
        disp('The EPmain.preferences.view.positive field is missing.');
    elseif isempty(EPmain.preferences.view.positive)
        err = [err ;37];
        disp('The EPmain.preferences.view.positive field is empty.');
    elseif ~isnumeric(EPmain.preferences.view.positive)
        err = [err ;37];
        disp('The EPmain.preferences.view.positive field is not a number.');
    elseif EPmain.preferences.view.positive < 1 && EPmain.preferences.view.positive > 2
        err = [err ;37];
        disp('The EPmain.preferences.view.positive field is not within the range of 1 to 2.');
    end
    
    if ~isfield(EPmain.preferences.view,'topoElectrodes')
        prefMissing = [prefMissing ;91];
        disp('The EPmain.preferences.view.topoElectrodes field is missing.');
    elseif isempty(EPmain.preferences.view.topoElectrodes)
        err = [err ;91];
        disp('The EPmain.preferences.view.topoElectrodes field is empty.');
    elseif ~ischar(EPmain.preferences.view.topoElectrodes)
        err = [err ;91];
        disp('The EPmain.preferences.view.topoElectrodes field is not a character string.');
    end

    if isfield(EPmain.preferences.view,'color')
        if length(EPmain.preferences.view.color)==EPmain.maxColors
            if ~isfield(EPmain.preferences.view.color,'lineSize')
                prefMissing = [prefMissing ;65];
                disp('The EPmain.preferences.view.color.lineSize field is missing.');
            else
                for iColor=1:EPmain.maxColors
                    if isempty([EPmain.preferences.view.color(iColor).lineSize])
                        err = [err ;65];
                        disp(['The EPmain.preferences.view.color(' num2str(iColor) ').lineSize field is empty.']);
                    elseif ~isnumeric(EPmain.preferences.view.color(iColor).lineSize)
                        err = [err ;65];
                        disp(['The EPmain.preferences.view.color(' num2str(iColor) ').lineSize field is not a number.']);
                    elseif EPmain.preferences.view.color(iColor).lineSize < 0
                        err = [err ;65];
                        disp(['The EPmain.preferences.view.color(' num2str(iColor) ').lineSize field is not positive.']);
                    end
                end
            end
            if ~isfield(EPmain.preferences.view.color,'lineStyle')
                prefMissing = [prefMissing ;90];
                disp('The EPmain.preferences.view.color.lineStyle field is missing.');
            else
                for iColor=1:EPmain.maxColors
                    if isempty([EPmain.preferences.view.color(iColor).lineStyle])
                        err = [err ;90];
                        disp(['The EPmain.preferences.view.color(' num2str(iColor) ').lineStyle field is empty.']);
                    elseif isnumeric(EPmain.preferences.view.color(iColor).lineStyle)
                        err = [err ;90];
                        disp(['The EPmain.preferences.view.color(' num2str(iColor) ').lineStyle field is a number.']);
                    elseif EPmain.preferences.view.color(iColor).lineStyle < 0
                        err = [err ;90];
                        disp(['The EPmain.preferences.view.color(' num2str(iColor) ').lineStyle field is not positive.']);
                    end
                end
            end
            if ~isfield(EPmain.preferences.view.color,'RGB')
                prefMissing = [prefMissing ;69];
                disp('The EPmain.preferences.view.color.RGB field is missing.');
            else
                for iColor=1:EPmain.maxColors
                    if isempty([EPmain.preferences.view.color(iColor).RGB])
                        err = [err ;69];
                        disp(['The EPmain.preferences.view.color(' num2str(iColor) ').RGB field is empty.']);
                    elseif length(EPmain.preferences.view.color(iColor).RGB)~=3
                        err = [err ;69];
                        disp(['The EPmain.preferences.view.color(' num2str(iColor) ').RGB field is not length three.']);
                    elseif any(~isnumeric(EPmain.preferences.view.color(iColor).RGB))
                        err = [err ;69];
                        disp(['An EPmain.preferences.view.color(' num2str(iColor) ').RGB field is not a number.']);
                    elseif any(EPmain.preferences.view.color(iColor).RGB < 0)
                        err = [err ;69];
                        disp(['An EPmain.preferences.view.color(' num2str(iColor) ').RGB field is not positive.']);
                    end
                end
            end
        else
            err = [err ;68];
            disp('The EPmain.preferences.view.color field is the wrong length.');
        end
    else
        prefMissing = [prefMissing ;68];
        disp('The EPmain.preferences.view.color field is missing.');
    end
    
    if ~isfield(EPmain.preferences.view,'lineSize')
        prefMissing = [prefMissing ;79];
        disp('The EPmain.preferences.view.lineSize field is missing.');
    elseif isempty(EPmain.preferences.view.lineSize)
        err = [err ;79];
        disp('The EPmain.preferences.view.lineSize field is empty.');
    elseif ~isnumeric(EPmain.preferences.view.lineSize)
        err = [err ;79];
        disp('The EPmain.preferences.view.lineSize field is not a number.');
    elseif EPmain.preferences.view.lineSize < 0
        err = [err ;79];
        disp('The EPmain.preferences.view.lineSize field is not equal to or larger than 0.');
    end
    
    if ~isfield(EPmain.preferences.view,'labelSize')
        prefMissing = [prefMissing ;66];
        disp('The EPmain.preferences.view.labelSize field is missing.');
    elseif isempty(EPmain.preferences.view.labelSize)
        err = [err ;66];
        disp('The EPmain.preferences.view.labelSize field is empty.');
    elseif ~isnumeric(EPmain.preferences.view.labelSize)
        err = [err ;66];
        disp('The EPmain.preferences.view.labelSize field is not a number.');
    elseif EPmain.preferences.view.labelSize < 0
        err = [err ;66];
        disp('The EPmain.preferences.view.labelSize field is not equal to or larger than 0.');
    end
else
    prefMissing = [prefMissing ;36];
    disp('The EPmain.preferences.view field is missing.');
end

if isfield(EPmain.preferences,'pca')
    if ~isfield(EPmain.preferences.pca,'mode')
        prefMissing = [prefMissing ;31];
        disp('The EPmain.preferences.pca.mode field is missing.');
    elseif isempty(EPmain.preferences.pca.mode)
        err = [err ;31];
        disp('The EPmain.preferences.pca.mode field is empty.');
    elseif ~isnumeric(EPmain.preferences.pca.mode)
        err = [err ;31];
        disp('The EPmain.preferences.pca.mode field is not a number.');
    elseif EPmain.preferences.pca.mode < 1 && EPmain.preferences.pca.mode > 2
        err = [err ;31];
        disp('The EPmain.preferences.pca.mode field is not within the range of 1 to 2.');
    end
    if ~isfield(EPmain.preferences.pca,'rotation')
        prefMissing = [prefMissing ;32];
        disp('The EPmain.preferences.pca.rotation field is missing.');
    elseif isempty(EPmain.preferences.pca.rotation)
        err = [err ;32];
        disp('The EPmain.preferences.pca.rotation field is empty.');
    elseif ~isnumeric(EPmain.preferences.pca.rotation)
        err = [err ;32];
        disp('The EPmain.preferences.pca.rotation field is not a number.');
    elseif EPmain.preferences.pca.rotation < 1 && EPmain.preferences.pca.rotation > 13
        err = [err ;32];
        disp('The EPmain.preferences.pca.rotation field is not within the range of 1 to 13.');
    end
    if ~isfield(EPmain.preferences.pca,'rotopt')
        prefMissing = [prefMissing ;33];
        disp('The EPmain.preferences.pca.rotopt field is missing.');
    elseif isempty(EPmain.preferences.pca.rotopt)
        err = [err ;33];
        disp('The EPmain.preferences.pca.rotopt field is empty.');
    elseif ~isnumeric(EPmain.preferences.pca.rotopt)
        err = [err ;33];
        disp('The EPmain.preferences.pca.rotopt field is not a number.');
    end
    if ~isfield(EPmain.preferences.pca,'rel')
        prefMissing = [prefMissing ;34];
        disp('The EPmain.preferences.pca.rel field is missing.');
    elseif isempty(EPmain.preferences.pca.rel)
        err = [err ;34];
        disp('The EPmain.preferences.pca.rel field is empty.');
    elseif ~isnumeric(EPmain.preferences.pca.rel)
        err = [err ;34];
        disp('The EPmain.preferences.pca.rel field is not a number.');
    elseif EPmain.preferences.pca.rel < 1 && EPmain.preferences.pca.rel > 2
        err = [err ;34];
        disp('The EPmain.preferences.pca.rel field is not within the range of 1 to 2.');
    end
    if ~isfield(EPmain.preferences.pca,'loadings')
        prefMissing = [prefMissing ;35];
        disp('The EPmain.preferences.pca.loadings field is missing.');
    elseif isempty(EPmain.preferences.pca.loadings)
        err = [err ;35];
        disp('The EPmain.preferences.pca.loadings field is empty.');
    elseif ~isnumeric(EPmain.preferences.pca.loadings)
        err = [err ;35];
        disp('The EPmain.preferences.pca.loadings field is not a number.');
    elseif EPmain.preferences.pca.loadings < 1 && EPmain.preferences.pca.loadings > 4
        err = [err ;35];
        disp('The EPmain.preferences.pca.loadings field is not within the range of 1 to 4.');
    end
else
    prefMissing = [prefMissing ;30];
    disp('The EPmain.preferences.pca field is missing.');
end

if isfield(EPmain.preferences,'window')
    if ~isfield(EPmain.preferences.window,'minFacVar')
        prefMissing = [prefMissing ;39];
        disp('The EPmain.preferences.window.minFacVar field is missing.');
    elseif isempty(EPmain.preferences.window.minFacVar)
        err = [err ;39];
        disp('The EPmain.preferences.window.minFacVar field is empty.');
    elseif ~isnumeric(EPmain.preferences.window.minFacVar)
        err = [err ;39];
        disp('The EPmain.preferences.window.minFacVar field is not a number.');
    elseif EPmain.preferences.window.minFacVar < 0 && EPmain.preferences.window.minFacVar > 1
        err = [err ;39];
        disp('The EPmain.preferences.window.minFacVar field is not within the range of 0 to 1.');
    end
    if ~isfield(EPmain.preferences.window,'adds')
        prefMissing = [prefMissing ;46];
        disp('The EPmain.preferences.window.adds field is missing.');
    elseif isempty(EPmain.preferences.window.adds)
        err = [err ;46];
        disp('The EPmain.preferences.window.adds field is empty.');
    elseif ~isnumeric(EPmain.preferences.window.adds)
        err = [err ;46];
        disp('The EPmain.preferences.window.adds field is not a number.');
    elseif EPmain.preferences.window.adds ~= 0 && EPmain.preferences.window.adds ~=1
        err = [err ;46];
        disp('The EPmain.preferences.window.adds field must be either 0 or 1.');
    end
    if ~isfield(EPmain.preferences.window,'chanGrp')
        prefMissing = [prefMissing ;57];
        disp('The EPmain.preferences.window.chanGrp field is missing.');
    elseif isempty(EPmain.preferences.window.chanGrp)
        err = [err ;57];
        disp('The EPmain.preferences.window.chanGrp field is empty.');
    elseif ~isnumeric(EPmain.preferences.window.chanGrp)
        err = [err ;57];
        disp('The EPmain.preferences.window.chanGrp field is not a number.');
    elseif EPmain.preferences.window.chanGrp ~= 1 && EPmain.preferences.window.chanGrp ~=2
        err = [err ;57];
        disp('The EPmain.preferences.window.adds field must be either 1 or 2.');
    end
else
    prefMissing = [prefMissing ;38];
    disp('The EPmain.preferences.window field is missing.');
end

if isfield(EPmain.preferences,'anova')
    if ~isfield(EPmain.preferences.anova,'trimming')
        prefMissing = [prefMissing ;41];
        disp('The EPmain.preferences.anova.trimming field is missing.');
    elseif isempty(EPmain.preferences.anova.trimming)
        err = [err ;41];
        disp('The EPmain.preferences.anova.trimming field is empty.');
    elseif ~isnumeric(EPmain.preferences.anova.trimming)
        err = [err ;41];
        disp('The EPmain.preferences.anova.trimming field is not a number.');
    elseif EPmain.preferences.anova.trimming < 0 || EPmain.preferences.anova.trimming > .5
        err = [err ;41];
        disp('The EPmain.preferences.anova.trimming field is not within the range of 0 to .5.');
    end
    if ~isfield(EPmain.preferences.anova,'bootstrap')
        prefMissing = [prefMissing ;42];
        disp('The EPmain.preferences.anova.bootstrap field is missing.');
    elseif isempty(EPmain.preferences.anova.bootstrap)
        err = [err ;42];
        disp('The EPmain.preferences.anova.bootstrap field is empty.');
    elseif ~isnumeric(EPmain.preferences.anova.bootstrap)
        err = [err ;42];
        disp('The EPmain.preferences.anova.bootstrap field is not a number.');
    end
    if ~isfield(EPmain.preferences.anova,'reps')
        prefMissing = [prefMissing ;62];
        disp('The EPmain.preferences.anova.reps field is missing.');
    elseif isempty(EPmain.preferences.anova.reps)
        err = [err ;62];
        disp('The EPmain.preferences.anova.reps field is empty.');
    elseif ~isnumeric(EPmain.preferences.anova.reps)
        err = [err ;62];
        disp('The EPmain.preferences.anova.reps field is not a number.');
    elseif ceil(EPmain.preferences.anova.reps) ~= EPmain.preferences.anova.reps
        err = [err ;62];
        disp('The EPmain.preferences.anova.reps field is not an integer.');
    elseif ceil(EPmain.preferences.anova.reps/2) == EPmain.preferences.anova.reps/2
        err = [err ;62];
        disp('The EPmain.preferences.anova.reps field is not an odd number.');
    end
    if ~isfield(EPmain.preferences.anova,'seed')
        prefMissing = [prefMissing ;43];
        disp('The EPmain.preferences.anova.seed field is missing.');
    elseif isempty(EPmain.preferences.anova.seed)
        err = [err ;43];
        disp('The EPmain.preferences.anova.seed field is empty.');
    elseif ~isnumeric(EPmain.preferences.anova.seed)
        err = [err ;43];
        disp('The EPmain.preferences.anova.seed field is not a number.');
    end
    if ~isfield(EPmain.preferences.anova,'missing')
        prefMissing = [prefMissing ;44];
        disp('The EPmain.preferences.anova.missing field is missing.');
    elseif isempty(EPmain.preferences.anova.missing)
        err = [err ;44];
        disp('The EPmain.preferences.anova.missing field is empty.');
    elseif ~isnumeric(EPmain.preferences.anova.missing)
        err = [err ;44];
        disp('The EPmain.preferences.anova.missing field is not a number.');
    end
    if ~isfield(EPmain.preferences.anova,'adds')
        prefMissing = [prefMissing ;45];
        disp('The EPmain.preferences.anova.adds field is missing.');
    elseif isempty(EPmain.preferences.anova.adds)
        err = [err ;45];
        disp('The EPmain.preferences.anova.adds field is empty.');
    elseif ~isnumeric(EPmain.preferences.anova.adds)
        err = [err ;45];
        disp('The EPmain.preferences.anova.adds field is not a number.');
    elseif EPmain.preferences.anova.adds ~= 0 && EPmain.preferences.anova.adds ~=1
        err = [err ;45];
        disp('The EPmain.preferences.anova.adds field must be either 0 or 1.');
    end
    if ~isfield(EPmain.preferences.anova,'epsilon')
        prefMissing = [prefMissing ;77];
        disp('The EPmain.preferences.anova.epsilon field is missing.');
    elseif isempty(EPmain.preferences.anova.epsilon)
        err = [err ;77];
        disp('The EPmain.preferences.anova.epsilon field is empty.');
    elseif isnumeric(EPmain.preferences.anova.epsilon)
        err = [err ;77];
        disp('The EPmain.preferences.anova.epsilon field is a number.');
    elseif ~ep_enumVar('epsilon',EPmain.preferences.anova.epsilon)
        err = [err ;77];
        disp('The EPmain.preferences.anova.epsilon field is an invalid value.');
        theString=ep_enumVar('epsilon');
        EPmain.preferences.anova.epsilon=theString{1};
    end
    if ~isfield(EPmain.preferences.anova,'posthoc')
        prefMissing = [prefMissing ;78];
        disp('The EPmain.preferences.anova.posthoc field is missing.');
    elseif isempty(EPmain.preferences.anova.posthoc)
        err = [err ;78];
        disp('The EPmain.preferences.anova.posthoc field is empty.');
    elseif isnumeric(EPmain.preferences.anova.posthoc)
        err = [err ;78];
        disp('The EPmain.preferences.anova.posthoc field is a number.');
    elseif strcmp(EPmain.preferences.anova.posthoc,'ADF')
        err = [err ;78];
        disp('The EPmain.preferences.anova.posthoc field is an invalid value.');
        EPmain.preferences.anova.posthoc='';
    elseif ~ep_enumVar('postHoc',EPmain.preferences.anova.posthoc)
        err = [err ;78];
        disp('The EPmain.preferences.anova.posthoc field is an invalid value.');
        theString=ep_enumVar('postHoc');
        EPmain.preferences.anova.posthoc=theString{1};
    end
    %     if ~isfield(EPmain.preferences.anova,'adf')
    %         prefMissing = [prefMissing ;79];
%         disp('The EPmain.preferences.anova.adf field is missing.');
%     elseif isempty(EPmain.preferences.anova.adf)
%         err = [err ;79];
%         disp('The EPmain.preferences.anova.adf field is empty.');
%     elseif isnumeric(EPmain.preferences.anova.adf)
%         err = [err ;79];
%         disp('The EPmain.preferences.anova.adf field is a number.');
%     end
else
    prefMissing = [prefMissing ;40];
    disp('The EPmain.preferences.anova field is missing.');
end

if isfield(EPmain.preferences,'advanced')
    if ~isfield(EPmain.preferences.advanced,'parallel')
        prefMissing = [prefMissing ;71];
        disp('The EPmain.preferences.advanced.parallel field is missing.');
    elseif isempty(EPmain.preferences.advanced.parallel)
        err = [err ;71];
        disp('The EPmain.preferences.advanced.parallel field is empty.');
    elseif ~isnumeric(EPmain.preferences.advanced.parallel)
        err = [err ;71];
        disp('The EPmain.preferences.advanced.parallel field is not a number.');
    elseif EPmain.preferences.advanced.parallel ~= 0 && EPmain.preferences.advanced.parallel ~=1
        err = [err ;71];
        disp('The EPmain.preferences.advanced.parallel field must be either 0 or 1.');
    end
    if ~isfield(EPmain.preferences.advanced,'monitor')
        prefMissing = [prefMissing ;75];
        disp('The EPmain.preferences.advanced.monitor field is missing.');
    elseif isempty(EPmain.preferences.advanced.monitor)
        err = [err ;75];
        disp('The EPmain.preferences.advanced.monitor field is empty.');
    elseif length(EPmain.preferences.advanced.monitor)~=4
        err = [err ;75];
        disp('The EPmain.preferences.advanced.monitor field must be four numbers (horizontal position, vertical position, horizontal resolution, vertical resolution.');
    elseif any(~isnumeric(EPmain.preferences.advanced.monitor))
        err = [err ;75];
        disp('A EPmain.preferences.advanced.monitor field is not a number.');
    elseif any(EPmain.preferences.advanced.monitor<0)
        err = [err ;75];
        disp('EPmain.preferences.advanced.monitor fields may not be negative.');
    end
    if ~isfield(EPmain.preferences.advanced,'RAM')
        prefMissing = [prefMissing ;76];
        disp('The EPmain.preferences.advanced.RAM field is missing.');
    elseif isempty(EPmain.preferences.advanced.RAM)
        err = [err ;76];
        disp('The EPmain.preferences.advanced.RAM field is empty.');
    elseif ~isnumeric(EPmain.preferences.advanced.RAM)
        err = [err ;76];
        disp('The EPmain.preferences.advanced.RAM field is not a number.');
    elseif length(EPmain.preferences.advanced.RAM)~=1
        err = [err ;76];
        disp('The EPmain.preferences.advanced.RAM field must be one number.');
    elseif EPmain.preferences.advanced.RAM<0
        err = [err ;76];
        disp('The EPmain.preferences.advanced.RAM field may not be negative.');
    end
else
    prefMissing = [prefMissing ;72];
    disp('The EPmain.preferences.advanced field is missing.');
end

if isfield(EPmain.preferences,'records')
    if ~isfield(EPmain.preferences.records,'user')
        prefMissing = [prefMissing ;84];
        disp('The EPmain.preferences.records.user field is missing.');
    elseif isnumeric(EPmain.preferences.records.user)
        err = [err ;84];
        disp('The EPmain.preferences.records.user field is a number.');
    end
    if ~isfield(EPmain.preferences.records,'lab')
        prefMissing = [prefMissing ;85];
        disp('The EPmain.preferences.records.lab field is missing.');
    elseif isnumeric(EPmain.preferences.records.lab)
        err = [err ;85];
        disp('The EPmain.preferences.records.lab field is a number.');
    end
    if ~isfield(EPmain.preferences.records,'institution')
        prefMissing = [prefMissing ;86];
        disp('The EPmain.preferences.records.institution field is missing.');
    elseif isnumeric(EPmain.preferences.records.institution)
        err = [err ;86];
        disp('The EPmain.preferences.records.institution field is a number.');
    end
    if ~isfield(EPmain.preferences.records,'project')
        prefMissing = [prefMissing ;87];
        disp('The EPmain.preferences.records.project field is missing.');
    elseif isnumeric(EPmain.preferences.records.project)
        err = [err ;87];
        disp('The EPmain.preferences.records.project field is a number.');
    end
    if ~isfield(EPmain.preferences.records,'experiment')
        prefMissing = [prefMissing ;88];
        disp('The EPmain.preferences.records.experiment field is missing.');
    elseif isnumeric(EPmain.preferences.records.experiment)
        err = [err ;88];
        disp('The EPmain.preferences.records.experiment field is a number.');
    end
else
    prefMissing = [prefMissing ;89];
    disp('The EPmain.preferences.records field is missing.');
end

if err
    beep
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function splash(src,eventdata) %present splash page

global EPmain

EPmain.handles.hSPlash=figure('Name', 'About EP Toolkit', 'NumberTitle', 'off', 'MenuBar', 'none', 'Position', [EPmain.scrsz(1) EPmain.scrsz(4)-550 570 413]);


try
    EPver=ver('EP_Toolkit');
catch
    EPver.Version='unavailable'; %workaround for bug in earlier version of Matlab
end

copyright=['ERP PCA Toolkit version: ' EPver.Version ' '...
    'Copyright (C) 1999-2025  Joseph Dien. ' ...
    sprintf('\n')...
    'This program is free software: you can redistribute it and/or modify '...
    'it under the terms of the GNU General Public License as published by '...
    'the Free Software Foundation, either version 3 of the License, or '...
    '(at your option) any later version. '...
    'This program is distributed in the hope that it will be useful, '...
    'but WITHOUT ANY WARRANTY; without even the implied warranty of '...
    'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE (especially medical).  See the '...
    'GNU General Public License for more details. '...
    'You should have received a copy of the GNU General Public License '...
    'along with this program.  If not, see <http://www.gnu.org/licenses/>. '...
    sprintf('\n')...
    'Note also that Dr. Dien does not have grant funding for this software '...
    'project and cannot therefore guarantee the same level of support '...
    'available for software projects that have full-time staff allocated '...
    'to them.  Nonetheless, Dr. Dien will endeavor to help users with '...
    'bug reports and questions to the best of his ability, given his time constraints. '...
    'He will be happy to help with implementing support for new file types. '...
    'Others are welcome to contribute code to this software project. '...
    'See files in the documentation directory for more information.'...
    sprintf('\n')...
    '***Again, this software is intended solely for scientific research and should never be used for medical purposes.***'];

uicontrol('Style','text','FontSize',EPmain.fontsize,...
    'String',copyright,...
    'HorizontalAlignment','left','Position',[20 50 500 350]);

if ~isfield(EPmain,'mode')
    EPmain.handles.splash.agree = uicontrol('Style', 'pushbutton', 'String', 'Agree','FontSize',EPmain.fontsize,...
        'Position', [20 0 100 40], 'Callback', ['global EPmain;','EPmain.mode=''agree'';','ep']);
    
    EPmain.handles.splash.disagree = uicontrol('Style', 'pushbutton', 'String', 'Disagree','FontSize',EPmain.fontsize,...
        'Position', [120 0 100 40], 'Callback', ['global EPmain;','EPmain.mode=''disagree'';','ep']);
else
    uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',EPmain.fontsize,...
        'Position', [120 0 100 40], 'Callback', ['close(''About EP Toolkit'');','ep(''start'');']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowSampAdapt(src,eventdata) %change period around peak that is averaged together for peak measures

global EPmain EPdataset

windowLength=EPmain.window.sampEnd-EPmain.window.sampStart+1;
Fs=EPdataset.dataset(EPmain.window.dataset).Fs;

if src==EPmain.handles.window.sampAdapt
    sampAdapt=str2num(get(EPmain.handles.window.sampAdapt,'String'));
else
    sampAdapt=round((str2num(get(EPmain.handles.window.msAdapt,'String'))/(1000/Fs)));
end

if (sampAdapt*2+1) > windowLength
    sampAdapt=floor((windowLength-1)/2);
    disp(['For current window size, max peak width is: ' num2str(sampAdapt) '.']);
end

if sampAdapt < 0
    sampAdapt=0;
    disp('Peak width cannot be negative.');
end

set(EPmain.handles.window.sampAdapt,'String',sampAdapt);
set(EPmain.handles.window.msAdapt,'String',round(sampAdapt*(1000/Fs)));

EPmain.window.sampAdapt=sampAdapt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowSampStart(src,eventdata) %change start of window

global EPmain EPdataset

numPoints=length(EPdataset.dataset(EPmain.window.dataset).timeNames);
Fs=EPdataset.dataset(EPmain.window.dataset).Fs;
baseline=EPdataset.dataset(EPmain.window.dataset).baseline;

if src==EPmain.handles.window.sampStart
    sampStart=str2num(get(EPmain.handles.window.sampStart,'String'));
    sampEnd=str2num(get(EPmain.handles.window.sampEnd,'String'));
else
    sampStart=round((str2num(get(EPmain.handles.window.msStart,'String'))/(1000/Fs))+baseline+1);
    sampEnd=round((str2num(get(EPmain.handles.window.msEnd,'String'))/(1000/Fs))+baseline);
end

if sampStart > numPoints
    sampStart=numPoints;
end

if sampStart < 1
    sampStart=1;
end

if sampStart > sampEnd
    sampEnd = sampStart;
end

set(EPmain.handles.window.sampStart,'String',sampStart);
set(EPmain.handles.window.sampEnd,'String',sampEnd);
set(EPmain.handles.window.msStart,'String',(sampStart-baseline-1)*(1000/Fs));
set(EPmain.handles.window.msEnd,'String',(sampEnd-baseline)*(1000/Fs));

EPmain.window.sampStart=sampStart;
EPmain.window.sampEnd=sampEnd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowSampEnd(src,eventdata) %change end of window

global EPmain EPdataset

numPoints=length(EPdataset.dataset(EPmain.window.dataset).timeNames);
Fs=EPdataset.dataset(EPmain.window.dataset).Fs;
baseline=EPdataset.dataset(EPmain.window.dataset).baseline;

if src==EPmain.handles.window.sampEnd
    sampStart=str2num(get(EPmain.handles.window.sampStart,'String'));
    sampEnd=str2num(get(EPmain.handles.window.sampEnd,'String'));
else
    sampStart=round((str2num(get(EPmain.handles.window.msStart,'String'))/(1000/Fs))+baseline+1);
    sampEnd=round((str2num(get(EPmain.handles.window.msEnd,'String'))/(1000/Fs))+baseline);
end

if sampEnd > numPoints
    sampEnd=numPoints;
end

if sampEnd < 1
    sampEnd=1;
end

if sampStart > sampEnd
    sampStart = sampEnd;
end

set(EPmain.handles.window.sampStart,'String',sampStart);
set(EPmain.handles.window.sampEnd,'String',sampEnd);
set(EPmain.handles.window.msStart,'String',(sampStart-baseline-1)*(1000/Fs));
set(EPmain.handles.window.msEnd,'String',(sampEnd-baseline)*(1000/Fs));

EPmain.window.sampStart=sampStart;
EPmain.window.sampEnd=sampEnd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowHzStart(src,eventdata) %change start of Hz window

global EPmain EPdataset

HzStart=str2num(get(EPmain.handles.window.binStart,'String'));
HzEnd=str2num(get(EPmain.handles.window.binEnd,'String'));
numFreqs=length(EPdataset.dataset(EPmain.window.dataset).freqNames);

if HzStart > numFreqs
    HzStart=numFreqs;
end

if HzStart < 1
    HzStart=1;
end

if HzStart > HzEnd
    HzEnd = HzStart;
end

set(EPmain.handles.window.binStart,'String',HzStart);
set(EPmain.handles.window.binEnd,'String',HzEnd);

EPmain.window.HzStart=HzStart;
EPmain.window.HzEnd=HzEnd;

set(EPmain.handles.window.HzStart,'String',sprintf('%5.1f',round(EPdataset.dataset(EPmain.window.dataset).freqNames(HzStart)*10)/10));
set(EPmain.handles.window.HzEnd,'String',sprintf('%5.1f',round(EPdataset.dataset(EPmain.window.dataset).freqNames(HzEnd)*10)/10));
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowHzEnd(src,eventdata) %change end of Hz window

global EPmain EPdataset

HzStart=str2num(get(EPmain.handles.window.binStart,'String'));
HzEnd=str2num(get(EPmain.handles.window.binEnd,'String'));
numFreqs=length(EPdataset.dataset(EPmain.window.dataset).freqNames);

if HzEnd > numFreqs
    HzEnd=numFreqs;
end

if HzEnd < 1
    HzEnd=1;
end

if HzStart > HzEnd
    HzStart = HzEnd;
end

set(EPmain.handles.window.binStart,'String',HzStart);
set(EPmain.handles.window.binEnd,'String',HzEnd);

EPmain.window.HzStart=HzStart;
EPmain.window.HzEnd=HzEnd;

set(EPmain.handles.window.HzStart,'String',sprintf('%5.1f',round(EPdataset.dataset(EPmain.window.dataset).freqNames(HzStart)*10)/10));
set(EPmain.handles.window.HzEnd,'String',sprintf('%5.1f',round(EPdataset.dataset(EPmain.window.dataset).freqNames(HzEnd)*10)/10));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowSpecsTable(src,eventdata) %change specs table data in Window pane

global EPmain

specTable=get(EPmain.handles.window.specsTable,'Data');
EPmain.window.specSelect=cell2mat(specTable(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowCellTable(src,eventdata) %change cell table data in Window pane

global EPmain

cellTable=get(EPmain.handles.window.cellTable,'Data');

EPmain.window.undo.outCells=EPmain.window.outCells;
EPmain.window.undo.inCells=EPmain.window.inCells;

EPmain.window.outCells=cellTable(:,1);
EPmain.window.inCells=cellTable(:,2);

set(EPmain.handles.window.undoWindowtable,'enable','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveData(src,eventdata,caller,batchFlag) 
%respond to click on save table

global EPmain EPdataset

if isempty(caller) %callback from table of file names
    if isempty(eventdata.Indices) %if just deselecting
        return;
    end
    if eventdata.Indices(2) == 1 %if just checking a checkbox
        tableData=get(EPmain.handles.save.hTable,'Data');
        if any([tableData{:,1}])
            set(EPmain.handles.save.save,'enable','on');
        else
            set(EPmain.handles.save.save,'enable','off');
        end
        return;
    end
    whichData=eventdata.Indices(1);
    EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);
else
    whichData=caller;
end

ep_tictoc('begin');

ep_tictoc('ioStart');
EPdata=ep_loadEPdataset(whichData);
ep_tictoc('ioFinish');

[~,~,formatCode]=ep_fileFormats(EPdata.dataType,EPmain.fileFormatSaveList{get(EPmain.handles.save.format,'value')});
[fileSuffix,~]=ep_fileExtensions(formatCode);

if ~batchFlag
    [outFileName, pathname] = uiputfile('*.*','Save:',[EPdataset.dataset(whichData).dataName fileSuffix{1}]);
    if outFileName == 0
        ep_tictoc('end');
        if isempty(caller) %callback from table of file names
            ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
            ep_tictoc('end');
            ep('start');
        end
        return %user hit cancel on file requestor
    end
else
    [outFileName, pathname, ext] = fileparts([EPdataset.dataset(whichData).dataName fileSuffix{1}]);
    pathname=[pathname ext];
end

if exist([pathname outFileName],'file')
    delete([pathname outFileName]); %user must have clicked "yes" to whether to replace existing file
end


outFileStem=outFileName;
periodList = strfind(outFileStem,'.');
if ~isempty(periodList)
    outFileStem=outFileStem(1:periodList(end)-1);
end

EPmain.convertMode=0;

saveFlag=saveFile(EPdata,[pathname outFileName]);

if saveFlag && ~batchFlag
    if ~strcmp(EPdataset.dataset(whichData).dataName,outFileStem) %if changed the dataset name
        button = questdlg('Did you want to change the name of the dataset in the working set as well?','Rename dataset?','Yes','No','Yes');
        if strcmp(button,'Yes')
            EPdata.dataName=outFileStem;
        end
    end
    ep_saveEPdataset(EPdata,whichData,'yes');
end

if isempty(caller) %callback from table of file names
    ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);
    ep_tictoc('end');
    ep('start');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveSelectData(src,eventdata) 
%saving all the checked files

global EPmain

EPmain.handleList=ep_disableGUI(EPmain.handles.hMainWindow);
tableData=get(EPmain.handles.save.hTable,'Data');

if ~any(cell2mat(tableData(:,1)))
    disp('Error: No files are selected.')
else
    batchFlag=1;
    for iFile=1:size(tableData,1)
        if tableData{iFile,1}
            tableData{iFile,1}=false;
            saveData(src,eventdata,iFile,batchFlag);
        end
    end
end

set(EPmain.handles.save.hTable,'Data',tableData);
ep_disableGUI(EPmain.handles.hMainWindow,EPmain.handleList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function saveFlag=saveFile(EPdata, fileName) 
%save a file

global EPmain EPtictoc

saveFlag=0;

set(EPmain.handles.save.done,'enable','off');
drawnow

[~,~,formatCode]=ep_fileFormats(EPdata.dataType,EPmain.fileFormatSaveList{get(EPmain.handles.save.format,'value')});

if strcmp(formatCode,'egi_egis') && ~strcmp(EPdata.dataType,'single_trial')
    formatCode='egi_egia';
end

if strcmp(formatCode,'egi_egia') && EPdata.baseline ==0
    choice = questdlg('EGIS average file with a zero baseline will cause problems for NetStation.  Save anyway?', ...
        'Options', ...
        'Yes','Cancel','Yes');
    % Handle response
    switch choice
        case 'Yes'
        case 'Cancel'
            set(EPmain.handles.save.done,'enable','on');
            drawnow
            return
    end
end

EPmain.save.format=get(EPmain.handles.save.format,'value');

if ~EPmain.save.SGLchan && ~EPmain.save.REGchan
    EPmain.save.SGLchan=1; %at least some type of channel needs to be kept
    disp('at least some type of channel needs to be kept so keeping single channels.');
end

if ~EPmain.save.SGLcell && ~EPmain.save.CMBcell
    EPmain.save.SGLcell=1; %at least some type of cell needs to be kept
    disp('at least some type of cell needs to be kept so keeping single channels.');
end

if ~EPmain.save.RAW && ~EPmain.save.AVG && ~EPmain.save.GAV
    EPmain.save.AVG=1; %at least some type of subject needs to be kept
    disp('at least some type of subject needs to be kept so keeping single channels.');
end

if ~EPmain.save.SGLfac && ~EPmain.save.CMBfac
    EPmain.save.SGLfac=1; %at least some type of factor needs to be kept
    disp('at least some type of factor needs to be kept (if any) so keeping single channels.');
end

adds=[];
if EPmain.save.SGLchan
    adds{end+1}='SGLchan';
end
if EPmain.save.REGchan
    adds{end+1}='REGchan';
end
if EPmain.save.SGLcell
    adds{end+1}='SGLcell';
end
if EPmain.save.CMBcell
    adds{end+1}='CMBcell';
    adds{end+1}='STScell';
end
if EPmain.save.RAW
    adds{end+1}='RAW';
end
if EPmain.save.AVG
    if any(strcmp(formatCode,{'egi_egis','egi_egia','egi_sbin','eeglab_set','text','neuromag_fif'})) && ~isempty(EPdata.facNames) &&  (length(EPdata.subNames)>1)
        choice = questdlg('Averaged Data option will result in an additional file for each subject.  Do anyway?', ...
            'Options', ...
            'Yes','No','Cancel','Yes');
        % Handle response
        switch choice
            case 'Yes'
                adds{end+1}='AVG';
            case 'No'
            case 'Cancel'
                set(EPmain.handles.save.done,'enable','on');
                drawnow
                return
        end
    else
        adds{end+1}='AVG';
    end
end
if EPmain.save.GAV
    adds{end+1}='GAV';
end
if EPmain.save.SGLfac
    adds{end+1}='SGLfac';
end
if EPmain.save.CMBfac
    adds{end+1}='CMBfac';
end

err=ep_writeData(EPdata,fileName,EPmain.preferences.general.specSuffix,EPmain.preferences.general.subjectSpecSuffix,formatCode,adds);
ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end

if ~err
    saveFlag=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function windowAdds(inputcells, outCellNames, autoPCA)
%add regional channels and combined cells to the active dataset that correspond to the windowing procedure.

global EPmain EPdataset EPchanGrp EPtictoc

EPdata=ep_loadEPdataset(EPmain.window.dataset);
tempGAVsubs=EPdata.GAVsubs; %save the GAVsubs info so it can be restored and to keep track of which ones are virtual.
numSubs=length(EPdata.subNames);
numVsubs=max(0,size(tempGAVsubs,1)-1);
numRsubs=numSubs-numVsubs;
numCells=length(EPdata.cellNames);
numVcells=max(0,size(EPdata.GAVsubs,2)-1);
numRcells=numCells-numVcells;

%convert virtual GAVEs to normal form so subject specs etc are available.
EPdata=ep_combineData(EPdata,'convert',{[],[],[],[],[],[]},[],[],[]);
if EPtictoc.stop;EPtictoc.stop=0;return;end
if isempty(EPdata)
    disp('Error: windowAdds failed.')
    return
end
numChans=length(EPdata.chanNames);
numPoints=length(EPdata.timeNames);
numCells=length(EPdata.cellNames);

numFacs=length(EPdata.facNames);
numFreqs=length(EPdata.freqNames);
if numFacs==0
    numFacs=1;
end
if ~isempty(EPdata.facData)
    numCMBfacs=size(EPdata.facData,5);
else
    numCMBfacs=0;
end
numSGLfacs=numFacs-numCMBfacs;

newAreas=[];
tempData=[];
addRegFlag=0;
if ~autoPCA
    if ~isempty(EPchanGrp) && (EPmain.window.chanGrp <= length(EPchanGrp.group))
        %add regional channels if any have new names
        areaList = setdiff(unique(EPchanGrp.group(EPmain.window.chanGrp).channel),EPchanGrp.numAreas+1); %identify area numbers that were used
        newAreas=find(~ismember(EPchanGrp.group(EPmain.window.chanGrp).areaName(areaList),EPdata.chanNames));
        if ~isempty(newAreas)
            for theArea=1:length(newAreas)
                whichArea=areaList(newAreas(theArea));
                chanList=find(EPchanGrp.group(EPmain.window.chanGrp).channel == whichArea);
                if length(chanList)>1 %no point adding a regional channel if only one channel was involved.
                    tempData=ep_combineData(EPdata,'channels',{chanList,[],[],[],[],[]},[],EPchanGrp.group(EPmain.window.chanGrp).areaName{whichArea});
                    if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
                    if ~isempty(tempData)
                        EPdata=tempData;
                        addRegFlag=1;
                    end
                end
            end
        end
    end
end
if addRegFlag
    disp(['Adding windowed regional channels to ' EPdataset.dataset(EPmain.window.dataset).dataName]);
end

%add combined output cells if any have new names

newCells=find(~ismember(outCellNames,EPdata.cellNames));
if ~isempty(newCells)
    disp(['Adding windowed combined cells to ' EPdataset.dataset(EPmain.window.dataset).dataName]);
    for theOutCell=1:length(newCells)
        cellList=inputcells{newCells(theOutCell)};
        tempData=ep_combineData(EPdata,'cells',{[],[],cellList,[],[],[]},[],outCellNames{newCells(theOutCell)});
        if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
        if ~isempty(tempData)
            EPdata=tempData;
        end
    end
    if ~isempty(tempData.GAVsubs)
        for iGAV=2:size(tempData.GAVsubs,1)
            for iFac=1:size(tempData.GAVsubs,3)
                tempData.GAVsubs{iGAV,end,iFac}=tempData.GAVsubs{iGAV,1,iFac};
            end
        end
        if isempty(tempGAVsubs)
            tempGAVsubs=tempData.GAVsubs;
        else
            tempGAVsubs(1,end+1:end+size(tempData.GAVsubs,2)-1,:)=tempData.GAVsubs(:,2:end,1);
        end
    end
end

if addRegFlag || ~isempty(newCells)
    if ~exist([EPdataset.EPwork filesep 'EPwork'],'dir')
        warndlg('The work directory cannot be found.')
        return
    end

    if ~isempty(tempGAVsubs)
        %if there are virtual grand averages, convert them back to virtual form
        tempSubNames=EPdata.subNames;
        tempSubTypes=EPdata.subTypes;
        tempCellNames=EPdata.cellNames;
        tempCellTypes=EPdata.cellTypes;
        [EPdata]=ep_selectData(EPdata,{[],[],[1:numRcells],[1:numRsubs],[],[]});
        EPdata.GAVsubs=tempGAVsubs;
        EPdata.subNames=tempSubNames;
        EPdata.subTypes=tempSubTypes;
        EPdata.cellNames=tempCellNames;
        EPdata.cellTypes=tempCellTypes;
    end

    [err]=ep_checkEPfile(EPdata);
    if ~err
        theDescription=['Added regional channel corresponding to the ANOVA window.'];
        EPdata.history=ep_addHistory(EPdata.history, EPmain.preferences.records, theDescription, EPmain.handles.hMainWindow, EPmain.preferences, EPmain.preferences.EPver);
        ep_tictoc('ioStart');
        ep_saveEPdataset(EPdata,EPmain.window.dataset,'no');
        ep_tictoc('ioFinish');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ANOVAAdds(ANOVAdata, betweenCombos, whichData, GAVsubsTable)
%add new grand averages for each between group, such that the subjects that are trimmed correspond to the ones that will be dropped in the ANOVA.
%Also add trimmed averages for each cell if there is non-zero trimming.
%These will only be computed for channel areas and cells that are present in the ANOVA output.

global EPmain EPdataset EPtictoc

if ~any(ismember(ANOVAdata.cellNames,ANOVAdata.EPdata.cellNames))
    disp('No cells in the ANOVA are present in the dataset so corresponding grand averages will not be added.');
end

if ~any(ismember(ANOVAdata.areaNames,ANOVAdata.EPdata.chanNames))
    disp('No channel regions in the ANOVA are present in the dataset so corresponding grand averages will not be added.');
end

if ~isfield(ANOVAdata,'EPdata')
    return
end

numFacs=length(ANOVAdata.EPdata.facNames);
if numFacs==0
    numFacs=1;
end
if ~isempty(ANOVAdata.EPdata.facData)
    numCMBfacs=size(ANOVAdata.EPdata.facData,5);
else
    numCMBfacs=0;
end

if ~isempty(ANOVAdata.EPdata.facNames)
    facLocation=strfind(ANOVAdata.factor,':');
    facLocation2=strfind(ANOVAdata.factor,'(');
    theFactor=find(strcmp(ANOVAdata.factor(facLocation+2:facLocation2-2),ANOVAdata.EPdata.facNames));
else
    theFactor=[];
end
if isempty(theFactor)
    theFactor =1;
end
% 
% if strcmp(ANOVAdata.changrp,'autoPCA')
%     windowName='autoPCA';
% else
%     windowName=ANOVAdata.window;
% end

for iBetween=3:size(GAVsubsTable,1)
    betweenName=['[' GAVsubsTable{iBetween,1} ']'];
    if isempty(find(strcmp(betweenName,ANOVAdata.EPdata.subNames), 1))
        EPadd=[];
        EPadd.subNames{1}=betweenName;
        EPadd.subTypes{1}='GAV';
        EPadd.GAVsubs{1}=[GAVsubsTable{iBetween,2} ones(size(GAVsubsTable{iBetween,2}))];

        [ANOVAdata.EPdata]=ep_addData(ANOVAdata.EPdata,EPadd,'subjects');
        if isempty(ANOVAdata.EPdata)
            disp('Error: ANOVAadds failed.')
            return
        end
    end
end
for iWithin=3:size(GAVsubsTable,2)
    withinName=['[' GAVsubsTable{1,iWithin} ']'];
    if isempty(find(strcmp(withinName,ANOVAdata.EPdata.cellNames), 1))
        EPadd=[];
        EPadd.cellNames{1}=withinName;
        EPadd.cellTypes{1}='CMB';
        if ~isempty(ANOVAdata.inputCells)
            inputCells=[];
            for iRow=1:size(GAVsubsTable{2,iWithin},1)
                for iInCell=1:length(ANOVAdata.inputCells{GAVsubsTable{2,iWithin}(iRow,1)})
                    EEGcell=find(strcmp(ANOVAdata.inputCells{GAVsubsTable{2,iWithin}(iRow,1)}{iInCell},ANOVAdata.EPdata.cellNames));
                    if ~isempty(inputCells)
                        GAVrow=find(inputCells(:,1)==EEGcell);
                    else
                        GAVrow=[];
                    end
                    if isempty(GAVrow)
                        inputCells(end+1,1)=EEGcell;
                        inputCells(end,2)=GAVsubsTable{2,iWithin}(iRow,2);
                    else
                        inputCells(GAVrow,2)=inputCells(GAVrow,2)+GAVsubsTable{2,iWithin}(iRow,2);
                    end
                end
            end
        else
            inputCells=GAVsubsTable{2,iWithin};
        end
        EPadd.GAVsubs{1}=inputCells;
        [ANOVAdata.EPdata]=ep_addData(ANOVAdata.EPdata,EPadd,'cells');
        if isempty(ANOVAdata.EPdata)
            disp('Error: ANOVAadds failed.')
            return
        end
    end
end

numSubs=length(ANOVAdata.EPdata.subNames);
numVsubs=max(0,size(ANOVAdata.EPdata.GAVsubs,1)-1);
numRsubs=numSubs-numVsubs;
numCells=length(ANOVAdata.EPdata.cellNames);
numVcells=max(0,size(ANOVAdata.EPdata.GAVsubs,2)-1);
numRcells=numCells-numVcells;

for iBetween=3:size(GAVsubsTable,1)
    betweenName=['[' GAVsubsTable{iBetween,1} ']'];
    betweenLevel=find(strcmp(betweenName,ANOVAdata.EPdata.subNames), 1)-numRsubs+1;
    for iWithin=3:size(GAVsubsTable,2)
        withinName=['[' GAVsubsTable{1,iWithin} ']'];
        if ~strcmp(betweenName,'[all]') || ~strcmp(withinName,'[all]') %can't be [all] for both between and within
            withinLevel=find(strcmp(withinName,ANOVAdata.EPdata.cellNames), 1)-numRcells+1;
            ANOVAdata.EPdata.GAVsubs{betweenLevel,withinLevel,theFactor}=[GAVsubsTable{iBetween,iWithin}, ones(size(GAVsubsTable{iBetween,iWithin}))];
        end
    end
    for iFac=1:numFacs
        ANOVAdata.EPdata.GAVsubs{betweenLevel,1,iFac}=[[1:numRsubs]',ones(numRsubs,1)];
    end
end

if ~exist([EPdataset.EPwork filesep 'EPwork'],'dir')
    warndlg('The work directory cannot be found.')
    return
end
[err]=ep_checkEPfile(ANOVAdata.EPdata);
if ~err
    theDescription=['Added trimmed grand average corresponding to the ANOVA sample.'];
    ANOVAdata.EPdata.history=ep_addHistory(ANOVAdata.EPdata.history, EPmain.preferences.records, theDescription, EPmain.handles.hMainWindow, EPmain.preferences, EPmain.preferences.EPver);
    ep_saveEPdataset(ANOVAdata.EPdata,whichData,'no');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function quit(src,eventdata)
%quit from EP Toolkit

clearWorkingSet

close('EP Toolkit');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeWindowDataset(src,eventdata)
%change the dataset in the window function

global EPmain EPdata EPdataset

oldDataset=EPmain.window.dataset;
tempVar=get(EPmain.handles.window.dataset,'Value');
if tempVar ~=0
    EPmain.window.dataset=EPmain.window.aveData(tempVar);
end
if isempty(tempVar)
    EPmain.window.dataset=EPmain.window.aveData(tempVar);
end
EPmain.window.datasetName=EPdataset.dataset(EPmain.window.dataset).dataName;

if isempty(EPdataset.dataset(EPmain.window.dataset).timeNames)
    EPmain.window.measure=1;
end

% if strcmp(EPdataset.dataset(EPmain.window.dataset).dataType,'single_trial')
%     [u i]=unique(EPdataset.dataset(EPmain.window.dataset).cellNames,'first');
%     EPmain.window.inCells=EPdataset.dataset(EPmain.window.dataset).cellNames(sort(i));
% else
cellNames=EPdataset.dataset(EPmain.window.dataset).cellNames(find(strcmp('SGL',EPdataset.dataset(EPmain.window.dataset).cellTypes)));
[u i]=unique(cellNames,'first');
EPmain.window.inCells=cellNames(sort(i));
% end

EPmain.window.outCells=EPmain.window.inCells;

EPmain.window.specSelect=repmat(false,length(EPdataset.dataset(EPmain.window.dataset).subjectSpecNames),1);

if ~isequal(EPdataset.dataset(EPmain.window.dataset).timeNames,EPdataset.dataset(oldDataset).timeNames)
    EPmain.window.sampStart=1;
    EPmain.window.sampEnd=length(EPdataset.dataset(EPmain.window.dataset).timeNames);
end
if ~isequal(EPdataset.dataset(EPmain.window.dataset).freqNames,EPdataset.dataset(oldDataset).freqNames)
    EPmain.window.HzStart=1;
    EPmain.window.HzEnd=length(EPdataset.dataset(EPmain.window.dataset).freqNames);
end
EPmain.window.factor=1;

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeWork(src,eventdata)
%change work directory

global EPdataset EPmain

newDir = uigetdir(EPdataset.EPwork,'Select EPwork directory');
if newDir == 0
    msg{1}='No directory selected.';
    [msg]=ep_errorMsg(msg);
    return
end

seps=findstr(newDir,filesep);
if ~isempty(seps)
    theDir=newDir(seps(end)+1:end);
    newDir=newDir(1:seps(end)-1);
else
    theDir=newDir;
    newDir='.';
end

if ~strcmp(theDir,'EPwork')
    msg{1}='Directory not named EPwork.';
    [msg]=ep_errorMsg(msg);
    return
end

if ~exist([newDir filesep 'EPwork' filesep 'EPdataset.mat'],'file')
    msg{1}='EPwork directory missing EPdataset file.';
    [msg]=ep_errorMsg(msg);
    return
end

if ~exist([newDir filesep 'EPwork' filesep 'EPprefs.mat'],'file')
    msg{1}='EPwork directory missing EPprefs file.';
    [msg]=ep_errorMsg(msg);
    return
end

tempVar=load([newDir filesep 'EPwork' filesep 'EPdataset.mat']);
if isfield(tempVar,'EPdataset')
    EPdataset=tempVar.EPdataset;
end
clear tempVar;

tempVar=load([newDir filesep 'EPwork' filesep 'EPprefs.mat']);
if isfield(tempVar,'prefs')
    prefs=tempVar.prefs;
end
clear tempVar;

EPmain.preferences=prefs; 
EPdataset.EPwork=newDir;

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function createWork(src,eventdata)
%create work directory in current directory

global EPdataset EPmain

if exist([pwd filesep 'EPwork'],'dir')
    msg{1}='EPwork directory already present.';
    [msg]=ep_errorMsg(msg);
    return
end

mkdir([pwd filesep 'EPwork']);
EPdataset=[];
EPdataset.EPwork=pwd;
EPdataset.dataset=cell(0);
eval(['save ''' pwd filesep 'EPwork' filesep 'EPdataset'' EPdataset']);
prefs=EPmain.preferences;
eval(['save ''' pwd filesep 'EPwork' filesep 'EPprefs'' prefs']);

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeSegmentDataset(src,eventdata)
%update the settings in the segment function

global EPmain EPdataset

if exist('src','var')
    %if triggered by changing the dataset, then also clear out the specs.
    clearSpec(1);
end

if strcmp(EPmain.segment.subPane,'cells')
    if isfield(EPmain.segment,'dataset')
        theDataset=get(EPmain.handles.segment.dataset,'Value');
        theEvent=get(EPmain.handles.segment.eventValues,'Value');
        if EPmain.segment.contData(theDataset) ~= EPmain.segment.dataset
            EPmain.segment.cellTable=cell(1,EPmain.segment.numFixed+EPmain.segment.numSpecs*3);
            theEvent=1;
        end
    else
        theDataset=length(EPmain.segment.contData);
        theEvent=1;
    end
    
    EPmain.segment.dataset=EPmain.segment.contData(theDataset);
    EPmain.segment.event=theEvent;
end


if strcmp(EPdataset.dataset(EPmain.segment.dataset).dataType,'continuous')
    allEvents=EPdataset.dataset(EPmain.segment.dataset).events;
else
    allEvents=cell(sum(cellfun(@length,EPdataset.dataset(EPmain.segment.dataset).events)),1);
    theCount=0;
    for iTrial=1:length(EPdataset.dataset(EPmain.segment.dataset).events)
        numEvents=length(EPdataset.dataset(EPmain.segment.dataset).events{1,iTrial});
        if numEvents > 0
            allEvents{1}(theCount+1:theCount+numEvents)=EPdataset.dataset(EPmain.segment.dataset).events{1,iTrial};
            theCount=theCount+numEvents;
        end
    end
end

if ~isempty(EPdataset.dataset(EPmain.segment.dataset).events{1})
    eventTypes={allEvents{1}.type}';
else
    eventTypes=[];
end
if ~isempty(EPdataset.dataset(EPmain.segment.dataset).events{1})
    eventValues={allEvents{1}.value}';
else
    eventValues=[];
end

%[eventValues{end+1:end+length(tempValues)}]=tempValues;
EPmain.segment.eventTypes=unique(cellfun(@num2str,eventTypes(find(~cellfun(@isempty,eventTypes))),'UniformOutput',false'));
EPmain.segment.eventTypes=sort(EPmain.segment.eventTypes);
EPmain.segment.eventValues=unique(cellfun(@num2str,eventValues(find(~cellfun(@isempty,eventValues))),'UniformOutput',false'));
EPmain.segment.eventValues=sort(EPmain.segment.eventValues);

%remove from EventType if all of a Type also have the same name for the Value (e.g., 'boundary' type and 'boundary' value) and also TRSP events.
dropList=[];
for iEvent=1:length(EPmain.segment.eventTypes)
    theType=EPmain.segment.eventTypes{iEvent};
    theValues=unique(cellfun(@num2str,eventValues(find(~cellfun(@isempty,eventValues)&strcmp(theType,eventTypes))),'UniformOutput',false'));
    if strcmp(theType,'boundary') || strcmp(theType,'TRSP')
        dropList(end+1)=iEvent;
    end
end
EPmain.segment.eventTypes(dropList)=[];
EPmain.segment.eventValues(strcmp('TRSP',EPmain.segment.eventValues))=[];
EPmain.segment.eventValues(strcmp('boundary',EPmain.segment.eventValues))=[];

numEventTypes=length(EPmain.segment.eventTypes);
numEventValues=length(EPmain.segment.eventValues);
EPmain.segment.eventLocks=EPmain.segment.eventTypes;
% EPmain.segment.eventLocks(end+1:end+numEventValues)=EPmain.segment.eventValues;

EPmain.segment.trialSpecNames=cell(0); %the names of the trial spec fields (cell array)
EPmain.segment.trialSpecValues=cell(0); %the values present for each of these trial spec fields (cell of cells)

EPmain.segment.stimSpecNames=cell(numEventTypes+numEventValues,1); %the list of key field names for each event value (cell of cell arrays)
EPmain.segment.stimSpecValues=cell(numEventTypes+numEventValues,1); %the values present for each of these  fields (cell of cell of cell arrys)

EPmain.segment.critSpecNames=cell(0); %the menu of names listed for each of the six crits
EPmain.segment.critSpecNames{1}='none';
EPmain.segment.critSpecNames{2}='-precedes-';
EPmain.segment.critSpecNames{3}='-follows-';
EPmain.segment.critSpecNames{4}='-secNum-';
EPmain.segment.critSpecNames{5}='-secTrial-';
EPmain.segment.critSpecNames{6}='-secTrialFromEnd-';
EPmain.segment.critSpecNames{7}='-responseCorr-';
EPmain.segment.critSpecNames{8}='-responseErr-';

EPmain.segment.critSpecValues=cell(0); %the menu of values that would be listed for each of the names in the crit name menus.
EPmain.segment.critSpecValues{1}{1}='none';
EPmain.segment.critSpecValues{2}=EPmain.segment.eventLocks; %possible values for '-precedes-'
EPmain.segment.critSpecValues{3}=EPmain.segment.eventLocks; %possible values for '-follows-'
EPmain.segment.critSpecValues{4}={'0'};
EPmain.segment.critSpecValues{5}={'0'};
EPmain.segment.critSpecValues{6}={'0'};
EPmain.segment.critSpecValues{7}=EPmain.segment.eventLocks; %possible values for '-responseCorr-'
EPmain.segment.critSpecValues{8}=EPmain.segment.eventLocks; %possible values for '-responseErr-'

if strcmp(EPdataset.dataset(EPmain.segment.dataset).dataType,'single_trial')
    for iTRSP=1:length(EPdataset.dataset(EPmain.segment.dataset).trialSpecNames)
        EPmain.segment.trialSpecNames{end+1}=EPdataset.dataset(EPmain.segment.dataset).trialSpecNames{iTRSP};
        EPmain.segment.trialSpecValues{end+1}=cell(0);
        for iTrial=1:length(EPdataset.dataset(EPmain.segment.dataset).trialSpecs)
            if ~any(find(strcmp(num2str(EPdataset.dataset(EPmain.segment.dataset).trialSpecs{iTrial,iTRSP}),EPmain.segment.trialSpecValues{end})))
                EPmain.segment.trialSpecValues{end}{end+1}=num2str(EPdataset.dataset(EPmain.segment.dataset).trialSpecs{iTrial,iTRSP});
            end
        end
    end
else
    TRSPevents=find(strcmp('TRSP',eventValues));
    if ~isempty(TRSPevents)
        for iTRSP=1:length(TRSPevents)
            theTRSP=TRSPevents(iTRSP);
            for iKey=1:length(allEvents{1}(theTRSP).keys)
                if ~any(strcmp(['TS-' allEvents{1}(theTRSP).keys(iKey).code],EPmain.segment.trialSpecNames))
                    EPmain.segment.trialSpecNames{end+1}=['TS-' allEvents{1}(theTRSP).keys(iKey).code];
                    if strcmp(EPmain.segment.trialSpecNames{end},'TS-RT')
                        EPmain.segment.trialSpecValues{end+1}{1}='100'; %100ms always on the list of RT values.
                        if ~strcmp(num2str(allEvents{1}(theTRSP).keys(iKey).data),'100')
                            EPmain.segment.trialSpecValues{end}{2}=num2str(allEvents{1}(theTRSP).keys(iKey).data);
                        end
                    else
                        EPmain.segment.trialSpecValues{end+1}{1}=num2str(allEvents{1}(theTRSP).keys(iKey).data);
                    end
                else
                    theTrialSpec=find(strcmp(['TS-' allEvents{1}(theTRSP).keys(iKey).code],EPmain.segment.trialSpecNames));
                    if ~any(find(strcmp(num2str(allEvents{1}(theTRSP).keys(iKey).data),EPmain.segment.trialSpecValues{theTrialSpec})))
                        EPmain.segment.trialSpecValues{theTrialSpec}{end+1}=num2str(allEvents{1}(theTRSP).keys(iKey).data);
                    end
                end
            end
        end
    end
end

if EPmain.segment.spec.format==1
    specSuffix='EPM-';
else
    specSuffix='TSP-';
end
if ~isempty(EPmain.segment.spec.specTable)
    for iSpec=1:length(EPmain.segment.spec.specLabels)
        if ~isempty(EPmain.segment.spec.specLabels{iSpec}) && (EPmain.segment.spec.specField(iSpec)~=length(EPmain.segment.spec.specNames))
            EPmain.segment.trialSpecNames{end+1}=[specSuffix EPmain.segment.spec.specLabels{iSpec}];
            EPmain.segment.trialSpecValues{end+1}=unique(EPmain.segment.spec.specTable(~cellfun(@isempty,EPmain.segment.spec.specTable(:,EPmain.segment.spec.specField(iSpec))),EPmain.segment.spec.specField(iSpec)));
        end
    end
end

for iEventType=1:numEventTypes
    lockEvents=find(strcmp(EPmain.segment.eventTypes{iEventType},eventTypes));
    for iEvent=1:length(lockEvents)
        theEvent=lockEvents(iEvent);
        theValue=num2str(allEvents{1}(theEvent).value);
        if isempty(theValue)
            theValue='none';
        end
        if iEvent==1
            EPmain.segment.stimSpecNames{iEventType}{end+1}='value';
            EPmain.segment.stimSpecValues{iEventType}{end+1}{1}=theValue;
        else
            if ~any(find(strcmp(theValue,EPmain.segment.stimSpecValues{iEventType}{1})))
                EPmain.segment.stimSpecValues{iEventType}{1}{end+1}=theValue;
            end
        end
        for iKey=1:length(allEvents{1}(theEvent).keys)
            if ~isempty(allEvents{1}(theEvent).keys(iKey).code)
                if ~any(strcmp(allEvents{1}(theEvent).keys(iKey).code,EPmain.segment.stimSpecNames{iEventType}))
                    EPmain.segment.stimSpecNames{iEventType}{end+1}=allEvents{1}(theEvent).keys(iKey).code;
                    EPmain.segment.stimSpecValues{iEventType}{end+1}{1}=allEvents{1}(theEvent).keys(iKey).data;
                else
                    theTrialSpec=find(strcmp(allEvents{1}(theEvent).keys(iKey).code,EPmain.segment.stimSpecNames{iEventType}));
                    if ~any(find(strcmp(allEvents{1}(theEvent).keys(iKey).data,EPmain.segment.stimSpecValues{iEventType}{theTrialSpec})))
                        EPmain.segment.stimSpecValues{iEventType}{theTrialSpec}{end+1}=allEvents{1}(theEvent).keys(iKey).data;
                    end
                end
            end
        end
    end
end

for iEventType=1:numEventValues
    lockEvents=find(strcmp(EPmain.segment.eventValues{iEventType},eventValues));
    for iEvent=1:length(lockEvents)
        theEvent=lockEvents(iEvent);
        for iKey=1:length(allEvents{1}(theEvent).keys)
            if ~isempty(allEvents{1}(theEvent).keys(iKey).code)
                if ~any(strcmp(allEvents{1}(theEvent).keys(iKey).code,EPmain.segment.stimSpecNames{numEventTypes+iEventType}))
                    EPmain.segment.stimSpecNames{numEventTypes+iEventType}{end+1}=allEvents{1}(theEvent).keys(iKey).code;
                    EPmain.segment.stimSpecValues{numEventTypes+iEventType}{end+1}{1}=allEvents{1}(theEvent).keys(iKey).data;
                else
                    theTrialSpec=find(strcmp(allEvents{1}(theEvent).keys(iKey).code,EPmain.segment.stimSpecNames{numEventTypes+iEventType}));
                    if ~any(find(strcmp(allEvents{1}(theEvent).keys(iKey).data,EPmain.segment.stimSpecValues{numEventTypes+iEventType}{theTrialSpec})))
                        EPmain.segment.stimSpecValues{numEventTypes+iEventType}{theTrialSpec}{end+1}=allEvents{1}(theEvent).keys(iKey).data;
                    end
                end
            end
        end
    end
end

EPmain.segment.trialSpec(1:EPmain.segment.numSpecs)=1;
EPmain.segment.trialSpecRel(1:EPmain.segment.numSpecs)=1;
EPmain.segment.trialSpecVal(1:EPmain.segment.numSpecs)=1;

for iSpec=1:length(EPmain.segment.trialSpecValues)
    EPmain.segment.trialSpecValues{iSpec}=sort(EPmain.segment.trialSpecValues{iSpec});
end

EPmain.segment.critSpecNames(end+1:end+length(EPmain.segment.trialSpecNames))=EPmain.segment.trialSpecNames; %the names of the fields to be listed for a criterion that does not come after a "-follows-" (cell array)
EPmain.segment.critSpecValues(end+1:end+length(EPmain.segment.trialSpecValues))=EPmain.segment.trialSpecValues; %the contents of the fields to be listed for a criterion that does not come after a "-follows-" (array of cell arrays)
EPmain.segment.critSpecNames(end+1:end+length(EPmain.segment.stimSpecNames{EPmain.segment.event}))=EPmain.segment.stimSpecNames{EPmain.segment.event};
EPmain.segment.critSpecValues(end+1:end+length(EPmain.segment.stimSpecValues{EPmain.segment.event}))=EPmain.segment.stimSpecValues{EPmain.segment.event};

for iCrit=1:EPmain.segment.numSpecs
    EPmain.segment.critSpecItemNames{iCrit}=EPmain.segment.critSpecNames; %the actual list of field names to be listed for each of the five model criterion rows
    EPmain.segment.critSpecItemValues{iCrit}=EPmain.segment.critSpecValues{EPmain.segment.trialSpec(iCrit)}; %the possible values for each of the actual field names to be listed for each of the five model criterion rows
end
EPmain.segment.critSpec{EPmain.segment.numSpecs}=EPmain.segment.critSpecNames([1 5:end]); %drop the "-precedes-" and "-follows-" and "-rangeSamps-" options

EPmain.segment.sampStart=-floor(200/(1000/EPdataset.dataset(EPmain.segment.dataset).Fs)); %default 200 ms
EPmain.segment.sampEnd=floor(1000/(1000/EPdataset.dataset(EPmain.segment.dataset).Fs)); %default 1000 ms

if any(strcmp(EPmain.segment.trialSpecNames,'TS-ACC'))
    EPmain.segment.ACC=find(strcmp(EPmain.segment.trialSpecNames,'TS-ACC'));
end
if any(strcmp(EPmain.segment.trialSpecNames,'TS-RT'))
    EPmain.segment.RT=find(strcmp(EPmain.segment.trialSpecNames,'TS-RT'));
end

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeSegmentEvent(src,eventdata)
%change the lock event in the segment function

global EPmain

EPmain.segment.event=get(EPmain.handles.segment.eventValues,'Value');

EPmain.segment.critSpecNames=cell(0);
EPmain.segment.critSpecNames{1}='none';
EPmain.segment.critSpecNames{2}='-precedes-';
EPmain.segment.critSpecNames{3}='-follows-';
EPmain.segment.critSpecNames{4}='-secNum-';
EPmain.segment.critSpecNames{5}='-secTrial-';
EPmain.segment.critSpecNames{6}='-secTrialFromEnd-';
EPmain.segment.critSpecNames{7}='-responseCorr-';
EPmain.segment.critSpecNames{8}='-responseErr-';


EPmain.segment.critSpecValues=cell(0);
EPmain.segment.critSpecValues{1}{1}='none';
EPmain.segment.critSpecValues{2}=EPmain.segment.eventLocks;
EPmain.segment.critSpecValues{3}=EPmain.segment.eventLocks;
EPmain.segment.critSpecValues{4}={'0'};
EPmain.segment.critSpecValues{5}={'0'};
EPmain.segment.critSpecValues{6}={'0'};
EPmain.segment.critSpecValues{7}=EPmain.segment.eventLocks;
EPmain.segment.critSpecValues{8}=EPmain.segment.eventLocks;

EPmain.segment.critSpecNames(end+1:end+length(EPmain.segment.trialSpecNames))=EPmain.segment.trialSpecNames;
EPmain.segment.critSpecValues(end+1:end+length(EPmain.segment.trialSpecValues))=EPmain.segment.trialSpecValues;
EPmain.segment.critSpecNames(end+1:end+length(EPmain.segment.stimSpecNames{EPmain.segment.event}))=EPmain.segment.stimSpecNames{EPmain.segment.event};
EPmain.segment.critSpecValues(end+1:end+length(EPmain.segment.stimSpecValues{EPmain.segment.event}))=EPmain.segment.stimSpecValues{EPmain.segment.event};

EPmain.segment.trialSpec(1:EPmain.segment.numSpecs)=1; %which field
EPmain.segment.trialSpecRel(1:EPmain.segment.numSpecs)=1; %which relation
EPmain.segment.trialSpecVal(1:EPmain.segment.numSpecs)=1; %which value

for iCrit=1:EPmain.segment.numSpecs
    EPmain.segment.critSpecItemNames{iCrit}=EPmain.segment.critSpecNames;
    EPmain.segment.critSpecItemValues{iCrit}=EPmain.segment.critSpecValues{EPmain.segment.trialSpec(iCrit)};
end
EPmain.segment.critSpec{EPmain.segment.numSpecs}=EPmain.segment.critSpecNames([1 4:end]);  %drop the "-precedes-" and "-follows-" options

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeSegmentCritName(src,eventdata)
%change a criterion name in the segment function

global EPmain

theCrit=find([EPmain.handles.segment.trialSpecNames] == src);

if EPmain.segment.trialSpec(theCrit) ~= get(EPmain.handles.segment.trialSpecNames(theCrit),'Value') %don't do anything unless the setting was actually changed
    oldName=EPmain.segment.trialSpec(theCrit);
    EPmain.segment.trialSpec(theCrit)=get(EPmain.handles.segment.trialSpecNames(theCrit),'Value');
    EPmain.segment.trialSpecRel(theCrit)=1;
    EPmain.segment.trialSpecVal(theCrit)=1;
    
    %if "-precedes-" or "-follows-" options, then reconfigure the options for the following criterion item
    if any(strcmp(EPmain.segment.critSpecItemNames{theCrit}{EPmain.segment.trialSpec(theCrit)},{'-precedes-','-follows-','-responseCorr-','-responseErr-'}))
        EPmain.segment.trialSpec(theCrit+1)=1;
        EPmain.segment.trialSpecRel(theCrit+1)=1;
        EPmain.segment.trialSpecVal(theCrit+1)=1;
        EPmain.segment.critSpecItemValues{theCrit}=EPmain.segment.eventLocks;
        EPmain.segment.critSpecItemNames{theCrit+1}=EPmain.segment.stimSpecNames{EPmain.segment.trialSpecVal(theCrit)};
        EPmain.segment.critSpecItemNames{theCrit+1}{end+1}='none';
        if ~isempty(EPmain.segment.stimSpecValues{EPmain.segment.trialSpecVal(theCrit)})
            EPmain.segment.critSpecItemValues{theCrit+1}=EPmain.segment.stimSpecValues{EPmain.segment.trialSpecVal(theCrit)}{EPmain.segment.trialSpec(theCrit+1)};
        end
        EPmain.segment.critSpecItemValues{theCrit+1}{end+1}='none';
        %also change the next criterion list to include the -rangeSamps- option unless too close to the end.
        if (theCrit+2)<=EPmain.segment.numSpecs
            EPmain.segment.critSpecItemNames{theCrit+2}{end+1}='-rangeSamps-';
        end
    elseif any(strcmp(EPmain.segment.critSpecItemNames{theCrit}{oldName},{'-precedes-','-follows-','-responseCorr-','-responseErr-'}))
        %if changing away from "-follows-" option, then reconfigure the options for the following criterion item
        EPmain.segment.trialSpec(theCrit+1)=1;
        EPmain.segment.trialSpecRel(theCrit+1)=1;
        EPmain.segment.trialSpecVal(theCrit+1)=1;
        EPmain.segment.critSpecItemValues{theCrit}=EPmain.segment.critSpecValues{EPmain.segment.trialSpec(theCrit)};
        EPmain.segment.critSpecItemNames{theCrit+1}=EPmain.segment.critSpecNames;
        EPmain.segment.critSpecItemValues{theCrit+1}=EPmain.segment.critSpecValues{EPmain.segment.trialSpec(theCrit+1)};
        %also change the next criterion list to include the -rangeSamps- option unless too close to the end.
        if (theCrit+2)<=EPmain.segment.numSpecs
            if strcmp('-rangeSamps-',EPmain.segment.critSpecItemNames{theCrit+2}{EPmain.segment.trialSpec(theCrit+2)})
                EPmain.segment.critSpecItemNames{theCrit+2}(EPmain.segment.trialSpec(theCrit+2))=[];
                EPmain.segment.trialSpec(theCrit+2)=1;
                EPmain.segment.trialSpecRel(theCrit+2)=1;
                EPmain.segment.trialSpecVal(theCrit+2)=1;
                EPmain.segment.critSpecItemNames{theCrit+2}=EPmain.segment.critSpecNames;
                EPmain.segment.critSpecItemValues{theCrit+2}=EPmain.segment.critSpecValues{EPmain.segment.trialSpec(theCrit+2)};
            end
        end
    else
        if (theCrit>1) && any(strcmp(EPmain.segment.critSpecItemNames{theCrit-1}{EPmain.segment.trialSpec(theCrit-1)},{'-precedes-','-follows-','-responseCorr-','-responseErr-'}))
            %if line after a -follows- criterion
            if strcmp('none',EPmain.segment.critSpecItemNames{theCrit}{EPmain.segment.trialSpec(theCrit)})
                EPmain.segment.critSpecItemValues{theCrit}=cell(0);
                EPmain.segment.critSpecItemValues{theCrit}{1}='none';
            else
                EPmain.segment.critSpecItemValues{theCrit}=cell(0);
                EPmain.segment.critSpecItemValues{theCrit}=EPmain.segment.stimSpecValues{EPmain.segment.trialSpecVal(theCrit-1)}{EPmain.segment.trialSpec(theCrit)};
            end
        elseif strcmp('-rangeSamps-',EPmain.segment.critSpecItemNames{theCrit}{EPmain.segment.trialSpec(theCrit)})
            EPmain.segment.critSpecItemValues{theCrit}=num2cell([1:max(abs([EPmain.segment.sampStart EPmain.segment.sampEnd]))]);
            EPmain.segment.trialSpecRel(theCrit)=5;
            EPmain.segment.trialSpecVal(theCrit)=length(EPmain.segment.critSpecItemValues{theCrit});
             %don't include the starts, ends, and contains relation types
        else
            EPmain.segment.critSpecItemValues{theCrit}=cell(0);
            EPmain.segment.critSpecItemValues{theCrit}=EPmain.segment.critSpecValues{EPmain.segment.trialSpec(theCrit)};
        end
    end
end

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeSegmentCritValue(src,eventdata)
%change a criterion value in the segment function

global EPmain EPdataset

theCrit=find([EPmain.handles.segment.trialSpecVal] == src);

if EPmain.segment.trialSpecVal(theCrit) ~= get(EPmain.handles.segment.trialSpecVal(theCrit),'Value') %don't do anything unless the setting was actually changed
    EPmain.segment.trialSpecVal(theCrit)=get(EPmain.handles.segment.trialSpecVal(theCrit),'Value');
    
    %if "-precedes-" or "-follows-" options, then reconfigure the options for the following criterion item
    if any(strcmp(EPmain.segment.critSpecItemNames{theCrit}{EPmain.segment.trialSpec(theCrit)},{'-precedes-','-follows-','-responseCorr-','-responseErr-'}))
        EPmain.segment.trialSpec(theCrit+1)=1;
        EPmain.segment.trialSpecRel(theCrit+1)=1;
        EPmain.segment.trialSpecVal(theCrit+1)=1;
        theName=EPmain.segment.stimSpecNames{EPmain.segment.trialSpecVal(theCrit)};
        if (length(theName)==1) && isempty(theName{1})
            EPmain.segment.critSpecItemNames{theCrit+1}=cell(0);
            EPmain.segment.critSpecItemNames{theCrit+1}{1}='none';
        else
            EPmain.segment.critSpecItemNames{theCrit+1}=cell(0);
            EPmain.segment.critSpecItemNames{theCrit+1}=theName;
            EPmain.segment.critSpecItemNames{theCrit+1}{end+1}='none';
        end
        if ~isempty(EPmain.segment.stimSpecValues{EPmain.segment.trialSpecVal(theCrit)})
            theValue=EPmain.segment.stimSpecValues{EPmain.segment.trialSpecVal(theCrit)}{EPmain.segment.trialSpec(theCrit+1)};
        else
            theValue=cell(0);
        end
        if (length(theValue)==1) && isempty(theValue{1})
            EPmain.segment.critSpecItemValues{theCrit+1}=cell(0);
            EPmain.segment.critSpecItemValues{theCrit+1}{1}='none';
        else
            EPmain.segment.critSpecItemValues{theCrit+1}=cell(0);
            EPmain.segment.critSpecItemValues{theCrit+1}=theValue;
            EPmain.segment.critSpecItemValues{theCrit+1}{end+1}='none';
        end
    end
end

ep('start');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function segmentSampStart(src,eventdata) %change start of segment

global EPmain EPdataset

Fs=EPdataset.dataset(EPmain.segment.dataset).Fs;

if src==EPmain.handles.segment.sampStart
    sampStart=str2num(get(EPmain.handles.segment.sampStart,'String'));
    sampEnd=str2num(get(EPmain.handles.segment.sampEnd,'String'));
else
    sampStart=round((str2num(get(EPmain.handles.segment.msStart,'String'))/(1000/Fs)));
    sampEnd=round((str2num(get(EPmain.handles.segment.msEnd,'String'))/(1000/Fs)));
end

if sampStart > sampEnd
    sampEnd = sampStart;
end

set(EPmain.handles.segment.sampStart,'String',sampStart);
set(EPmain.handles.segment.sampEnd,'String',sampEnd);
set(EPmain.handles.segment.msStart,'String',(sampStart)*(1000/Fs));
set(EPmain.handles.segment.msEnd,'String',(sampEnd)*(1000/Fs));

EPmain.segment.sampStart=sampStart;
EPmain.segment.sampEnd=sampEnd;
ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segmentSampEnd(src,eventdata) %change end of segment

global EPmain EPdataset

Fs=EPdataset.dataset(EPmain.segment.dataset).Fs;

if src==EPmain.handles.segment.sampEnd
    sampStart=str2num(get(EPmain.handles.segment.sampStart,'String'));
    sampEnd=str2num(get(EPmain.handles.segment.sampEnd,'String'));
else
    sampStart=round((str2num(get(EPmain.handles.segment.msStart,'String'))/(1000/Fs)));
    sampEnd=round((str2num(get(EPmain.handles.segment.msEnd,'String'))/(1000/Fs)));
end

if sampStart > sampEnd
    sampStart = sampEnd;
end

set(EPmain.handles.segment.sampStart,'String',sampStart);
set(EPmain.handles.segment.sampEnd,'String',sampEnd);
set(EPmain.handles.segment.msStart,'String',(sampStart)*(1000/Fs));
set(EPmain.handles.segment.msEnd,'String',(sampEnd)*(1000/Fs));

EPmain.segment.sampStart=sampStart;
EPmain.segment.sampEnd=sampEnd;
ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segmentAddCell(src,eventdata) %add cell to segmenting table

global EPmain

if all(cellfun(@isempty,EPmain.segment.cellTable(end,:)))
    theCell=size(EPmain.segment.cellTable,1);
else
    theCell=size(EPmain.segment.cellTable,1)+1;
end

EPmain.segment.cellTable{theCell,1}=[];
EPmain.segment.cellTable{theCell,2}='';
EPmain.segment.cellTable{theCell,3}=EPmain.segment.eventLocks{EPmain.segment.event};
if ~EPmain.segment.flexible
    EPmain.segment.cellTable{theCell,4}=get(EPmain.handles.segment.msStart,'String');
    EPmain.segment.cellTable{theCell,5}=get(EPmain.handles.segment.msEnd,'String');
else
    EPmain.segment.cellTable{theCell,4}=EPmain.segment.eventLocks{EPmain.segment.flexEnd};
    EPmain.segment.cellTable{theCell,5}=['F' get(EPmain.handles.segment.flexLength,'String')];
end
EPmain.segment.cellTable{theCell,6}=get(EPmain.handles.segment.delay,'String');
EPmain.segment.cellTable{theCell,7}=get(EPmain.handles.segment.task,'String');
if ~isempty(EPmain.segment.trialSpecNames) || ~isempty(EPmain.segment.stimSpecNames{EPmain.segment.event})
    if isempty(EPmain.segment.trialSpecNames)
        theList=EPmain.segment.stimSpecNames{EPmain.segment.event};
    else
        theList=EPmain.segment.trialSpecNames;
    end
    EPmain.segment.cellTable{theCell,8}=theList{EPmain.segment.ACC};
    EPmain.segment.cellTable{theCell,9}=theList{EPmain.segment.RT};
end

for iSpec=1:EPmain.segment.numSpecs
    EPmain.segment.cellTable{theCell,EPmain.segment.numFixed+1+(iSpec-1)*3}=EPmain.segment.critSpecItemNames{iSpec}{EPmain.segment.trialSpec(iSpec)};
    EPmain.segment.cellTable{theCell,EPmain.segment.numFixed+2+(iSpec-1)*3}=EPmain.segment.relList{EPmain.segment.trialSpecRel(iSpec)};
    EPmain.segment.cellTable{theCell,EPmain.segment.numFixed+3+(iSpec-1)*3}=EPmain.segment.critSpecItemValues{iSpec}{EPmain.segment.trialSpecVal(iSpec)};
end

set(EPmain.handles.segment.cellTable,'Data',EPmain.segment.cellTable);
ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segmentDelCell(src,eventdata) %delete cell of segmenting table

global EPmain

if size(EPmain.segment.cellTable,1) ==1
    for i=1:size(EPmain.segment.cellTable,2)
        EPmain.segment.cellTable{i}='';
    end
else
    EPmain.segment.cellTable=EPmain.segment.cellTable(1:end-1,:);
end

 set(EPmain.handles.segment.cellTable,'Data',EPmain.segment.cellTable);
 ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segmentLoad(src,eventdata) %load segmenting table

global EPmain

[FileName,PathName,FilterIndex] = uigetfile('','Load Segmenting Table');
if FileName ~= 0
    eval(['tempVar=load( ''' PathName FileName ''');']);
    if isfield(tempVar,'EPsegmentTable')
        if size(tempVar.EPsegmentTable,2)==20 %backward compatibility for delay
            tempVar.EPsegmentTable(:,7:5*3+6)=tempVar.EPsegmentTable(:,6:5*3+5);
            [tempVar.EPsegmentTable{:,6}]=deal(0);
        end        
        if size(tempVar.EPsegmentTable,2)==21 %backward compatibility for sixth spec
            tempVar.EPsegmentTable(:,22:22+(EPmain.segment.numSpecs-5)*3-1)=deal(0);
        end
        if size(tempVar.EPsegmentTable,2)==24 %backward compatibility for task
            tempVar.EPsegmentTable(:,8:6*3+7)=tempVar.EPsegmentTable(:,7:6*3+6);
            [tempVar.EPsegmentTable{:,7}]=deal('');
        end
        if size(tempVar.EPsegmentTable,2)==25 %backward compatibility for ACC and RT
            tempVar.EPsegmentTable(:,10:6*3+9)=tempVar.EPsegmentTable(:,8:6*3+7);
            [tempVar.EPsegmentTable{:,8}]=deal('TS-ACC');
            [tempVar.EPsegmentTable{:,9}]=deal('TS-RT');
        end
        tempVar.EPsegmentTable=cellfun(@num2str,tempVar.EPsegmentTable,'UniformOutput',false); %converts number to text and table to cell string.
        EPmain.segment.cellTable=tempVar.EPsegmentTable;
        set(EPmain.handles.segment.cellTable,'Data',EPmain.segment.cellTable);
    else
        msg{1}='Error: not a segmenting table file.';
        [msg]=ep_errorMsg(msg);
    end
    ep('start');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function segmentSave(src,eventdata) %save segmenting table

global EPmain EPdataset

[fileName,PathName,FilterIndex] = uiputfile('*.mat','Save Segmenting Table',[EPdataset.dataset(EPmain.segment.dataset).dataName '_segmentTable']);

if fileName ~= 0
    EPsegmentTable=EPmain.segment.cellTable;
    eval(['save ''' PathName fileName ''' EPsegmentTable']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function convertFiles(src,eventdata) %batch convert files

global EPmain

set(EPmain.handles.save.convert,'enable','off');
set(EPmain.handles.save.done,'enable','off');

drawnow

typeNum = get(EPmain.handles.save.fileType,'value');
switch typeNum
    case 1
        dataType='continuous';
    case 2
        dataType='single_trial';
    case 3
        dataType='average';
    case 4
        dataType='grand_average';
    case 5
        dataType='factors';
end

formatNum = get(EPmain.handles.save.readFormat,'value');
[importSuffix,importFormatName,importFormat]=ep_fileFormats(dataType,EPmain.fileFormatReadList{formatNum});

if strcmp(importFormat,'ep_mat')
    dataType=''; %data type set by data file itself
end

EPmain.save.readFormat=formatNum;
EPmain.save.fileType=typeNum;

EPmain.convertMode=1;
theHandles=EPmain.handles.save;
if EPmain.save.NS4fix
    NS4fix
else
    readFiles(theHandles,importFormat,dataType);
end
set(EPmain.handles.save.convert,'enable','on');
set(EPmain.handles.save.done,'enable','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NS4fix(src,eventdata) %batch convert NS4 simple binary files

global EPmain EPdataset EPtictoc

ep_tictoc('begin');

%get list of continuous files to modify
[rawFileNames, pathname] = uigetfile('*.raw','Raw files:','MultiSelect','on');
if ~iscell(rawFileNames)
    temp=rawFileNames;
    rawFileNames=[];
    rawFileNames{1}=temp;
end
if rawFileNames{1}==0
    msg{1}='No raw filenames selected. You have to click on a name.';
    [msg]=ep_errorMsg(msg);
    return
end

theSub='';
EPdata=[];
for iFile=1:length(rawFileNames)
    ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
    disp(rawFileNames{iFile});
    [pathstr, fileroot, ~]=fileparts(rawFileNames{iFile});
    inArg{1}='file';
    inArg{2}=rawFileNames{iFile};
    inArg{3}='format';
    inArg{4}='egi_sbin';
    inArg{5}='type';
    inArg{6}='continuous';
    inArg{7}='ced';
    inArg{8}='GSN129.ced';
%     if strcmp(inputFormat,'matlab_mat') && (iFile >1)
%         inArg{end+1}='matlabDims';
%         inArg{end+1}=outInfo.matlabDims;
%     end
    
    if isempty(EPdata) %start new subject
        [EPdata, origEloc, outInfo]=ep_readData(inArg);
        if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
    else
        [newData, origEloc, outInfo]=ep_readData(inArg);
        if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
        if isempty(newData.data)
            msg{1}='Failure to read file.';
            [msg]=ep_errorMsg(msg);
            return
        end
        boundaryPoint=length(EPdata.timeNames);
        EPdata.events{1}(end+1).type='boundary';
        EPdata.events{1}(end).sample=boundaryPoint;
        EPdata.events{1}(end).value='boundary';
        EPdata.events{1}(end).duration=0;
        EPdata.events{1}(end).keys=struct('code','','data','','datatype','','description','');
        EPdata=ep_addData(EPdata,newData,'points');
    end
    if iFile<length(rawFileNames)
        newSub=rawFileNames{iFile+1}(1:2);
    else
        newSub='***';
    end
    if (iFile==length(rawFileNames))||~strcmp(theSub,newSub)
        if ~isempty(theSub)
            ep_tictoc('ioStart');
            EPdata2 = ep_readEventText(EPdata, [rawFileNames{iFile-1}(1:end-7) '.evt']);
            ep_tictoc('ioFinish');
            eventCounter=1;
            for iEvent=1:length(EPdata.events{1})
                if ~strcmp(EPdata.events{1}(iEvent).type,'boundary')
                    while strcmp(EPdata2.events{1}(eventCounter).type,'CELL')
                        eventCounter=eventCounter+1;
                    end
                    if ~strcmp(EPdata.events{1}(iEvent).value,EPdata2.events{1}(eventCounter).type)
                        msg{1}='Mismatch of event list items.';
                        [msg]=ep_errorMsg(msg);
                        return
                    end
                    EPdata.events{1}(iEvent).keys=EPdata2.events{1}(eventCounter).keys;
                    if strcmp(EPdata2.events{1}(eventCounter).type,'ssnd')
                        while strcmp(EPdata2.events{1}(eventCounter).type,'ssnd')
                            eventCounter=eventCounter+1;
                        end
                    else
                        eventCounter=eventCounter+1;
                    end
                end
            end
            if (eventCounter-1) ~= length(EPdata2.events{1})
                msg{1}='Mismatch of event list length.';
                [msg]=ep_errorMsg(msg);
                return
            end
            ep_tictoc('ioStart');
            save('-mat',[rawFileNames{iFile-1}(1:end-7) '.ept'],'EPdata');
            ep_tictoc('ioFinish');
        end
        EPdata=[];
        theSub=newSub;
    end
end

ep_tictoc('end');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function changeSampleTestDataset(src,eventdata)
%change the dataset in the sampleTest function

global EPmain EPdataset

EPmain.sampleTest.dataset=EPmain.sampleTest.datasetList(get(EPmain.handles.sampleTest.dataset,'Value'));
theDataset=EPmain.sampleTest.dataset;

EPmain.sampleTest.dataType=EPdataset.dataset(theDataset).dataType;
EPmain.sampleTest.cell1=1;
EPmain.sampleTest.cell2=2;
EPmain.sampleTest.channel=1;
EPmain.sampleTest.cellNameList=unique(EPdataset.dataset(theDataset).cellNames(find(ismember(EPdataset.dataset(theDataset).cellTypes,{'SGL','CMB'}))));

EPmain.sampleTest.cellNameList=unique(EPdataset.dataset(theDataset).cellNames(find(ismember(EPdataset.dataset(theDataset).cellTypes,{'SGL','CMB'}))));

EPmain.sampleTest.subNameList=cell(0);
EPmain.sampleTest.subList=cell(0);
EPmain.sampleTest.subNameList{1}='-all-';
EPmain.sampleTest.subList{1}=[1:length(find(strcmp(EPdataset.dataset(theDataset).subTypes,'AVG')))];
for iSpec=1:length(EPdataset.dataset(theDataset).subjectSpecNames)
    specLvLlist=cell(0);
    for iSub=1:length(EPdataset.dataset(theDataset).subjectSpecs(:,iSpec))
        theSpec=num2str(EPdataset.dataset(theDataset).subjectSpecs{iSub,iSpec});
        if ~isempty(theSpec) && ~all(cellfun(@isempty,cellstr(theSpec))) && ~any(strcmp(theSpec,specLvLlist))
            specLvLlist{end+1}=EPdataset.dataset(theDataset).subjectSpecs{iSub,iSpec};
            EPmain.sampleTest.subNameList{end+1}=[EPdataset.dataset(theDataset).subjectSpecNames{iSpec} '_' specLvLlist{end}];
            EPmain.sampleTest.subList{end+1}=find(strcmp(specLvLlist{end},EPdataset.dataset(theDataset).subjectSpecs(:,iSpec)));
        end
    end
end

if ~isempty(EPdataset.dataset(theDataset).sessNums)
    sessList=unique(EPdataset.dataset(theDataset).sessNums);
    sessList=sessList(sessList>0);
    EPmain.sampleTest.sessNameList=EPdataset.dataset(theDataset).sessNames(sessList);
else
    EPmain.sampleTest.sessNameList=cell(0);
    EPmain.sampleTest.sess2=1;
end
EPmain.sampleTest.sess1=1;
EPmain.sampleTest.sess2=min(2,length(EPmain.sampleTest.sessNameList));

EEGchans=find(strcmp('EEG',EPdataset.dataset(theDataset).chanTypes));
EPmain.sampleTest.eloc=EPdataset.dataset(EPmain.sampleTest.datasetList(EPmain.sampleTest.dataset)).eloc(EEGchans);
EPmain.sampleTest.changeFlag=1;
EPmain.sampleTest.freqFlag=~isempty(EPdataset.dataset(theDataset).freqNames);

maxRad=0.5;
[y,x] = pol2cart(([EPmain.sampleTest.eloc.theta]/360)*2*pi,[EPmain.sampleTest.eloc.radius]);  % transform electrode locations from polar to cartesian coordinates
y=-y; %flip y-coordinate so that nose is upwards.
plotrad = min(1.0,max([EPmain.sampleTest.eloc.radius])*1.02);            % default: just outside the outermost electrode location
plotrad = max(plotrad,0.5);                 % default: plot out to the 0.5 head boundary
x = x*(maxRad/plotrad);
y = y*(maxRad/plotrad);

xmin = min(-maxRad,min(x));
xmax = max(maxRad,max(x));
ymin = min(-maxRad,min(y));
ymax = max(maxRad,max(y));

EPmain.sampleTest.x=round(((x/(xmax-xmin))*EPmain.sampleTest.gridSize)+ceil(EPmain.sampleTest.gridSize/2));
EPmain.sampleTest.y=round(((y/(ymax-ymin))*EPmain.sampleTest.gridSize)+ceil(EPmain.sampleTest.gridSize/2));

EPmain.sampleTest.PCAlist=[];
for iFile=1:length(EPdataset.dataset)
    if ~isempty(EPdataset.dataset(iFile).facNames) && ~isempty(EPdataset.dataset(iFile).eloc)
        newEloc=EPdataset.dataset(iFile).eloc(find(strcmp('EEG',EPdataset.dataset(iFile).chanTypes)));
        if length(EPmain.sampleTest.eloc) == length(newEloc)
            if isequal(cellfun(@isempty,{EPmain.sampleTest.eloc.theta}),cellfun(@isempty,{newEloc.theta})) && ~any([EPmain.sampleTest.eloc.theta]-[newEloc.theta]) && ~any([EPmain.sampleTest.eloc.radius]-[newEloc.radius])
                if (~strcmp(EPmain.sampleTest.dataType,'continuous') && (length(EPdataset.dataset(iFile).timeNames) == length(EPdataset.dataset(theDataset).timeNames))) ||...
                        (strcmp(EPmain.sampleTest.dataType,'continuous') && (length(EPdataset.dataset(iFile).timeNames) <= length(EPdataset.dataset(theDataset).timeNames)))
                    if (length(EPdataset.dataset(iFile).freqNames)==length(EPdataset.dataset(theDataset).freqNames)) && all([EPdataset.dataset(iFile).freqNames] == [EPdataset.dataset(theDataset).freqNames])
                        EPmain.sampleTest.PCAlist(end+1)=iFile;
                    end
                end
            end
        end
    end
end
EPmain.sampleTest.PCAnameList={EPdataset.dataset(EPmain.sampleTest.PCAlist).dataName};

EPmain.sampleTest.AVElist=[];
for iFile=1:length(EPdataset.dataset)
    if isempty(EPdataset.dataset(iFile).facNames) && ~isempty(EPdataset.dataset(iFile).eloc) && strcmp(EPdataset.dataset(iFile).dataType,'average')
        newEloc=EPdataset.dataset(iFile).eloc(find(strcmp('EEG',EPdataset.dataset(iFile).chanTypes)));
        if length(EPmain.sampleTest.eloc) == length(newEloc)
            if isequal(cellfun(@isempty,{EPmain.sampleTest.eloc.theta}),cellfun(@isempty,{newEloc.theta})) && ~any([EPmain.sampleTest.eloc.theta]-[newEloc.theta]) && ~any([EPmain.sampleTest.eloc.radius]-[newEloc.radius])
                if (~strcmp(EPmain.sampleTest.dataType,'continuous') && (length(EPdataset.dataset(iFile).timeNames) == length(EPdataset.dataset(theDataset).timeNames))) ||...
                        (strcmp(EPmain.sampleTest.dataType,'continuous') && (length(EPdataset.dataset(iFile).timeNames) <= length(EPdataset.dataset(theDataset).timeNames)))
                    if (length(EPdataset.dataset(iFile).freqNames)==length(EPdataset.dataset(theDataset).freqNames)) && all([EPdataset.dataset(iFile).freqNames] == [EPdataset.dataset(theDataset).freqNames])
                        EPmain.sampleTest.AVElist(end+1)=iFile;
                        EPdata=ep_loadEPdataset(iFile);
                        EPmain.sampleTest.AVEdata{end+1}=EPdata.data;
                    end
                end
            end
        end
    end
end
EPmain.sampleTest.AVEnameList={EPdataset.dataset(EPmain.sampleTest.AVElist).dataName};
EPmain.sampleTest.datasetName=EPdataset.dataset(EPmain.sampleTest.dataset).dataName;
if isfield(EPdataset.dataset(EPmain.sampleTest.dataset),'lastChange') && ~isempty(EPdataset.dataset(EPmain.sampleTest.dataset).lastChange)
    EPmain.sampleTest.lastChange=EPdataset.dataset(EPmain.sampleTest.dataset).lastChange;
else
    EPmain.sampleTest.lastChange=cell(0);
end

ep('start')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function expandTransformChan(src,eventdata)
global EPmain

scrsz = EPmain.scrsz;

if isfield(EPmain.handles.transform,'expandedFigure')
    if ~isempty(EPmain.handles.transform.expandedFigure)
        if ishandle(EPmain.handles.transform.expandedFigure)
            close(EPmain.handles.transform.expandedFigure)
            EPmain.handles.transform.expandedFigure=[];
        end
    end
end

if (src == EPmain.handles.transform.timeFilter) || (src == EPmain.handles.transform.timeWaves(1)) || (src == EPmain.handles.transform.timeWaves(2))
    name='Time Domain';
    theData=EPmain.transform.timeData;
    theScale=EPmain.transform.timeScale;
    theAxis=EPmain.transform.timeAxis;
    colorOrder=get(EPmain.handles.transform.timeFilter,'ColorOrder');
else
    name='Frequency Domain';
    theData=EPmain.transform.freqData;
    theScale=EPmain.transform.freqScale;
    theAxis=EPmain.transform.freqAxis;
    colorOrder=get(EPmain.handles.transform.frequencyFilter,'ColorOrder');
end

EPmain.handles.transform.expandedFigure = figure('Name', name, 'NumberTitle', 'off', 'Position',[scrsz(3)/2 scrsz(4)/2 600 400]);

set(EPmain.handles.transform.expandedFigure,'DefaultAxesColorOrder',colorOrder);
EPmain.handles.transform.expandedAxes = axes('position',[.05 .05 .9 .9]);

plot(theScale,theData);

axis(theAxis);

line([theScale(1) theScale(end)],[0 0],'Color','black','LineWidth',1) % zero line
if src == EPmain.handles.transform.frequencyFilter
    if ~isempty(EPmain.transform.filter1)
        axis(theAxis);
        line([EPmain.transform.filter1 EPmain.transform.filter1],[theAxis(3) theAxis(4)],'Color','black','LineWidth',1) % filter line
    end
    if ~isempty(EPmain.transform.filter2)
        axis(theAxis);
        line([EPmain.transform.filter2 EPmain.transform.filter2],[theAxis(3) theAxis(4)],'Color','black','LineWidth',1) % filter line
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function readFiles(theHandles,importFormat,dataType)
%choose a set of files to read in
global EPmain EPtictoc

ep_tictoc('begin');

textPrefs.firstRow=EPmain.preferences.general.firstRow;
textPrefs.lastRow=EPmain.preferences.general.lastRow;
textPrefs.firstCol=EPmain.preferences.general.firstCol;
textPrefs.lastCol=EPmain.preferences.general.lastCol;
textPrefs.orientation=EPmain.preferences.general.orientation;
textPrefs.sampleRate=EPmain.preferences.general.sampleRate;

singleCellMode=get(theHandles.check,'value');

if ~singleCellMode
    
    [fileNames, activeDirectory]=ep_getFilesUI(importFormat);
    if isempty(activeDirectory)
        return
    end
    if isempty(fileNames) || ((length(fileNames{1})==1) && (fileNames{1}==0))
        msg{1}='No filenames selected. You have to click on a name.';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    if any(strcmp(importFormat,{'ns_avg','eeglab_set','text','brainvision_eeg'})) && (length(fileNames)>1)
        disp('Usually when importing average data of this type, you want to use Single File mode so you can combine the separate cells and subject files into a single EP Toolkit file.');
    end
    
    fileNames=sort(fileNames);
    
    ced=[];
    for iFile=1:length(fileNames)
        ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
        inArg=[];
        inArg{1}='file';
        inArg{2}=[activeDirectory fileNames{iFile}];
        inArg{3}='format';
        inArg{4}=importFormat;
        inArg{5}='screenSize';
        inArg{6}=EPmain.scrsz;
        inArg{7}='FontSize';
        inArg{8}=EPmain.fontsize;
        if ~strcmp(importFormat,'ep_mat')
            inArg{9}='type';
            inArg{10}=dataType;
        end
        inArg{end+1}='textPrefs';
        inArg{end+1}=textPrefs;
        if ~isempty(ced)
            inArg{end+1}='ced';
            inArg{end+1}=ced;
            disp('Assuming the electrode information is the same as for the last file.');
        end
        if ~strcmp(importFormat,'ep_mat')
            inArg{end+1}='montage';
            inArg{end+1}=EPmain.read.importMontage;
        end
        SMIsuffix=EPmain.preferences.general.SMIsuffix;
        if ~isempty(SMIsuffix)
            inArg{end+1}='SMIsuffix';
            inArg{end+1}=SMIsuffix;
        end
        specSuffix=EPmain.preferences.general.specSuffix;
        if ~isempty(specSuffix)
            inArg{end+1}='specSuffix';
            inArg{end+1}=specSuffix;
        end
        subjectSpecSuffix=EPmain.preferences.general.subjectSpecSuffix;
        if ~isempty(subjectSpecSuffix)
            inArg{end+1}='subjectSpecSuffix';
            inArg{end+1}=subjectSpecSuffix;
        end
        BVheader=EPmain.preferences.general.BVheader;
        if ~isempty(BVheader)
            inArg{end+1}='BVheader';
            inArg{end+1}=BVheader;
            if EPmain.save.EPheaderFix
                inArg{end}=2;
            end
        end
        noInternal=EPmain.preferences.general.noInternal;
        if ~isempty(noInternal)
            inArg{end+1}='noInternal';
            inArg{end+1}=noInternal;
        end
        if strcmp(importFormat,'matlab_mat') && (iFile >1)
            inArg{end+1}='matlabDims';
            inArg{end+1}=outInfo.matlabDims;
        end
        
        disp(['Reading: ' fileNames{iFile}])
        [EPdata, origEloc, outInfo]=ep_readData(inArg);
        if EPtictoc.stop
            EPtictoc.stop=0;
            return
        end
        if isempty(EPdata) || isempty(EPdata.data)
            beep();
            disp(['Error: Was unable to read ' inArg{2}]);
        else
            if ~isempty(EPdata.ced)
                ced=EPdata.ced;
            end
            if ~isempty(EPdata.montage)
                montage=EPdata.montage;
            end
            
            if ~strcmp(importFormat,'ep_mat')
                theDescription=['Imported the file ' fileNames{iFile} '.'];
                EPdata.history=ep_addHistory(EPdata.history, EPmain.preferences.records, theDescription, EPmain.handles.hMainWindow, EPmain.preferences, EPmain.preferences.EPver, fileNames{iFile});
            end
            
            if EPmain.convertMode
                saveFile(EPdata,inArg{2});
            else
                ep_saveEPdataset(EPdata);
            end
            if EPtictoc.stop
                EPtictoc.stop=0;
                return
            end
        end
    end
    
else %single cell mode
    
    subPosStr=get(theHandles.subject,'string');
    cellPosStr=get(theHandles.cell,'string');
    freqPosStr=get(theHandles.freq,'string');
    sessPosStr=get(theHandles.session,'string');
    subPos=str2num(subPosStr);
    cellPos=str2num(cellPosStr);
    freqPos=str2num(freqPosStr);
    sessPos=str2num(sessPosStr);
    
    if (EPmain.read.fileType==2) || strcmp(importFormat,'ep_mat')
        trialPosStr=get(theHandles.trial,'string');
        trialPos=str2num(trialPosStr);
    else
        trialPosStr='';
        trialPos=[];
    end
    
    if min(subPos) < 1
        beep
        set(theHandles.subjectLabel,'ForegroundColor','red');
        drawnow
        return;
    end
    if min(cellPos) < 1
        beep
        set(theHandles.cellLabel,'ForegroundColor','red');
        drawnow
        return;
    end
    if min(freqPos) < 1
        beep
        set(theHandles.freqLabel,'ForegroundColor','red');
        drawnow
        return;
    end
    if min(trialPos) < 1
        beep
        set(theHandles.trialLabel,'ForegroundColor','red');
        drawnow
        return;
    end
    if min(sessPos) < 1
        beep
        set(theHandles.sessLabel,'ForegroundColor','red');
        drawnow
        return;
    end
    if isempty(subPos) && isempty(cellPos) && isempty(trialPos) && isempty(sessPos) && isempty(freqPos)
        beep
        disp('None of the single file mode fields were set.')
        drawnow
        return;
    end

    if EPmain.convertMode
        EPmain.save.check=1;
        EPmain.save.subject=subPosStr;
        EPmain.save.cell=cellPosStr;
        EPmain.save.freq=freqPosStr;
        EPmain.save.trial=trialPosStr;
        EPmain.save.sess=sessPosStr;
    else
        EPmain.read.check=1;
        EPmain.read.subject=subPosStr;
        EPmain.read.cell=cellPosStr;
        EPmain.read.freq=freqPosStr;
        EPmain.read.trial=trialPosStr;
        EPmain.read.sess=sessPosStr;
    end
    
    [sessionFiles, pathname]=ep_getFilesUI(importFormat);
    if isempty(sessionFiles)
        return
    end
    if (sessionFiles{1}==0) & ~ischar(sessionFiles{1})
        msg{1}='No filenames selected. You have to click on a name.';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    sessionFiles=sort(sessionFiles);
    
    for iFile=1:length(sessionFiles)
        ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
        [pathstr, name, fileSuffix] = fileparts(sessionFiles{iFile});
        if ~isempty(subPos)
            if max(subPos) > length(name)
                msg{1}=['The file name ' sessionFiles{iFile} ' does not have ' num2str(max(subPos)) ' characters.'];
                [msg]=ep_errorMsg(msg);
                return
            end
            theSubs{iFile}=sessionFiles{iFile}(subPos);
        else
            if any(strcmp(dataType,{'grand_average','factors'}))
                theSubs{iFile}='GAV';
            else
                theSubs{iFile}='S01';
            end
        end
        if ~isempty(cellPos)
            if max(cellPos) > length(name)
                msg{1}=['The file name ' sessionFiles{iFile} ' does not have ' num2str(max(cellPos)) ' characters.'];
                [msg]=ep_errorMsg(msg);
                return
            end
            theCells{iFile}=sessionFiles{iFile}(cellPos);
        else
            theCells{iFile}='';
        end
        if ~isempty(freqPos)
            if max(freqPos) > length(name)
                msg{1}=['The file name ' sessionFiles{iFile} ' does not have ' num2str(max(freqPos)) ' characters.'];
                [msg]=ep_errorMsg(msg);
                return
            end
            theFreqs{iFile}=str2num(sessionFiles{iFile}(freqPos));
        else
            theFreqs{iFile}=[];
        end
        if ~isempty(trialPos)
            if max(trialPos) > length(name)
                msg{1}=['The file name ' sessionFiles{iFile} ' does not have ' num2str(max(trialPos)) ' characters.'];
                [msg]=ep_errorMsg(msg);
                return
            end
            trialNum=str2double(sessionFiles{iFile}(trialPos));
            if ~isnan(trialNum)
                theTrials(iFile)=trialNum;
            else
                theTrials{iFile}=0;
            end
        else
            theTrials{iFile}=0;
        end
        if ~isempty(sessPos)
            if max(sessPos) > length(name)
                msg{1}=['The file name ' sessionFiles{iFile} ' does not have ' num2str(max(sessPos)) ' characters.'];
                [msg]=ep_errorMsg(msg);
                return
            end
            theSess{iFile}=sessionFiles{iFile}(sessPos);
        else
            theSess{iFile}=[];
        end
        sessionFiles{iFile}=[pathname sessionFiles{iFile}];
        
        if iFile==1
            [path fileName ext]=fileparts(sessionFiles{iFile});
            disp(['According to your settings, the first file(' fileName ') is for:']);
            if ~isempty(subPos)
                disp(['Subject: ' fileName(subPos)]);
            end
            if ~isempty(cellPos)
                disp(['Cell: ' fileName(cellPos)]);
            end
            if ~isempty(freqPos)
                disp(['Frequency: ' fileName(freqPos)]);
            end
            if ~isempty(trialPos)
                disp(['Trial: ' fileName(trialPos)]);
            end
            if ~isempty(sessPos)
                disp(['Session: ' fileName(sessPos)]);
            end
        end
    end
    
    if ~strcmp(dataType,'continuous')
        mergeName = char(inputdlg('Name of new dataset?','Dataset name'));
        pause(1); %pause to avoid crash on Windows due to bug in Matlab 2009

        if isempty(mergeName)
            msg{1}='No filename given for the new merged file.  Read operation aborted.';
            [msg]=ep_errorMsg(msg);
            return
        end
    end
    
    eloc=[];
    ced=[];
    if strcmp(dataType,'average')
        uniqueSubs=[];
        numSubFiles=1;
    else
        uniqueSubs=unique(theSubs);
        numSubFiles=length(uniqueSubs);
    end
    
    for iSub=1:numSubFiles
        inArg=[];
        if strcmp(dataType,'average')
            theFiles=[1:length(theSubs)];
        else
            theFiles=find(strcmp(uniqueSubs(iSub),theSubs));
        end
        for iFile=1:length(theFiles)
            inArg{iFile,1}=sessionFiles{theFiles(iFile)};
            inArg{iFile,2}='format';
            inArg{iFile,3}=importFormat;
            if ~strcmp(importFormat,'ep_mat')
                inArg{iFile,4}='type';
                inArg{iFile,5}=dataType;
            end
            inArg{iFile,6}='labels';
            inArg{iFile,7}={theCells(theFiles(iFile)) theSubs(theFiles(iFile)) theSess(theFiles(iFile)) theFreqs(theFiles(iFile)) theTrials(theFiles(iFile))};
            inArg{iFile,8}='textPrefs';
            inArg{iFile,9}=textPrefs;
            inArg{iFile,12}='screenSize';
            inArg{iFile,13}=EPmain.scrsz;
            inArg{iFile,14}='FontSize';
            inArg{iFile,15}=EPmain.fontsize;
            if ~isempty(eloc)
                inArg{iFile,16}='eloc';
                inArg{iFile,17}=eloc;
            end
            if ~isempty(ced)
                inArg{iFile,18}='ced';
                inArg{iFile,19}=ced;
                disp('Assuming the electrode information is the same as for the last file.');
            end
            inArg{iFile,20}='montage';
            inArg{iFile,21}=EPmain.read.importMontage;
            SMIsuffix=EPmain.preferences.general.SMIsuffix;
            if ~isempty(SMIsuffix)
                inArg{iFile,22}='SMIsuffix';
                inArg{iFile,23}=SMIsuffix;
            end
            specSuffix=EPmain.preferences.general.specSuffix;
            if ~isempty(specSuffix)
                inArg{iFile,24}='specSuffix';
                inArg{iFile,25}=specSuffix;
            end
            subjectSpecSuffix=EPmain.preferences.general.subjectSpecSuffix;
            if ~isempty(subjectSpecSuffix)
                inArg{iFile,26}='subjectSpecSuffix';
                inArg{iFile,27}=subjectSpecSuffix;
            end
            BVheader=EPmain.preferences.general.BVheader;
            if ~isempty(BVheader)
                inArg{iFile,28}='BVheader';
                inArg{iFile,29}=BVheader;
            end
        end

        if ~strcmp(dataType,'continuous')
            if length(theFiles) ==1
                [path fullMergeName ext]=fileparts(sessionFiles{theFiles});
            elseif length(uniqueSubs) > 1 && ~strcmp(dataType,'average')
                fullMergeName=[mergeName '_' num2str(iSub)];
            else
                fullMergeName=mergeName;
            end
            [EPdata, eloc]=ep_mergeEPfiles(inArg,fullMergeName);
            ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
            if isempty(EPdata)
                beep();
                disp(['Error: Was unable to read ' fullMergeName]);
            else
                numChans=length(EPdata.chanNames);
                numPoints=length(EPdata.timeNames);
                numCells=length(unique(EPdata.cellNames));
                numTrials=length(EPdata.trialNames);
                numSubs=length(EPdata.subNames);
                numFacs=length(EPdata.facNames);
                numFreqs=length(EPdata.freqNames);
                numRels=length(EPdata.relNames);
                numSess=length(unique(EPdata.sessNames));

                disp(['The imported file ' fullMergeName])
                disp(['Channels: ' num2str(numChans)]);
                disp(['Time Points: ' num2str(numPoints)]);
                disp(['Cells: ' num2str(numCells)]);
                if strcmp(EPdata.dataType,'single_trial')
                    disp(['Trials: ' num2str(numTrials)]);
                end
                disp(['Subjects: ' num2str(numSubs)]);
                if ~isempty(EPdata.facNames)
                    disp(['Factors: ' num2str(numFacs)]);
                end
                if ~isempty(EPdata.freqNames)
                    disp(['Frequencies: ' num2str(numFreqs)]);
                end
                if ~isempty(EPdata.relNames)
                    disp(['Relations: ' num2str(numRels)]);
                end
                if ~isempty(EPdata.sessNames)
                    disp(['Sessions: ' num2str(numSess)]);
                end

                if ~isempty(EPdata.ced)
                    ced=EPdata.ced;
                end
                if ~isempty(EPdata.montage)
                    montage=EPdata.montage;
                end
                if ~strcmp(importFormat,'ep_mat')
                    theDescription=['Merged imported single mode files.'];
                    EPdata.history=ep_addHistory(EPdata.history, EPmain.preferences.records, theDescription, EPmain.handles.hMainWindow, EPmain.preferences, EPmain.preferences.EPver,sessionFiles);
                end
                if EPmain.convertMode
                    saveFile(EPdata,[pathname fullMergeName]);
                else
                    ep_saveEPdataset(EPdata);
                end
            end
        else %continuous files
            for iFile=1:length(theFiles)
                inArg2=cell(1,size(inArg,2)+1);
                inArg2(1,2:end)=inArg(iFile,:);
                inArg2{1,1}='file';
                [EPdata, origEloc, outInfo]=ep_readData(inArg2);
                ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
                if isempty(EPdata)
                    beep();
                    disp(['Error: Was unable to read ' inArg{iFile,1}]);
                else
                    numChans=length(EPdata.chanNames);
                    numPoints=length(EPdata.timeNames);
                    numCells=length(unique(EPdata.cellNames));
                    numTrials=length(EPdata.trialNames);
                    numSubs=length(EPdata.subNames);
                    numFacs=length(EPdata.facNames);
                    numFreqs=length(EPdata.freqNames);
                    numRels=length(EPdata.relNames);
                    numSess=length(unique(EPdata.sessNames));

                    disp(['The imported file ' inArg{iFile,1}])
                    disp(['Channels: ' num2str(numChans)]);
                    disp(['Time Points: ' num2str(numPoints)]);
                    disp(['Cells: ' num2str(numCells)]);
                    if ~isempty(EPdata.facNames)
                        disp(['Factors: ' num2str(numFacs)]);
                    end
                    if ~isempty(EPdata.freqNames)
                        disp(['Frequencies: ' num2str(numFreqs)]);
                    end
                    if ~isempty(EPdata.relNames)
                        disp(['Relations: ' num2str(numRels)]);
                    end
                    if ~isempty(EPdata.sessNames)
                        disp(['Sessions: ' num2str(numSess)]);
                    end

                    if ~isempty(EPdata.ced)
                        ced=EPdata.ced;
                    end
                    if ~isempty(EPdata.montage)
                        montage=EPdata.montage;
                    end
                    if ~strcmp(importFormat,'ep_mat')
                        theDescription=['Imported single mode file.'];
                        EPdata.history=ep_addHistory(EPdata.history, EPmain.preferences.records, theDescription, EPmain.handles.hMainWindow, EPmain.preferences, EPmain.preferences.EPver,sessionFiles);
                    end
                    if EPmain.convertMode
                        [path fileName ext]=fileparts(inArg{iFile,1});
                        saveFile(EPdata,[pathname fileName]);
                    else
                        ep_saveEPdataset(EPdata);
                    end
                end
            end
        end
    end
end
ep_tictoc('end');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function crossVerifyData(src, eventdata)
%choose PCA dataset for cross-verification
global EPmain

if isempty(eventdata.Indices) %if just deselecting
    return;
end

EPmain.pca.crossVerifyPCA=EPmain.pca.PCAdatasets(eventdata.Indices(1));
ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sampleTestSampStart(src,eventdata) %change start of sample test template

global EPmain

if src==EPmain.handles.sampleTest.sampStart
    EPmain.sampleTest.sampStart=str2num(get(EPmain.handles.sampleTest.sampStart,'String'))+1+EPmain.sampleTest.baseline;
else
    EPmain.sampleTest.sampStart=round(str2num(get(EPmain.handles.sampleTest.msStart,'String'))/(1000/EPmain.sampleTest.Fs))+EPmain.sampleTest.baseline+1;
end

if EPmain.sampleTest.sampStart < -EPmain.sampleTest.baseline+1
    EPmain.sampleTest.sampStart = -EPmain.sampleTest.baseline+length(EPmain.sampleTest.templateWaveform);
end

if EPmain.sampleTest.sampStart > EPmain.sampleTest.sampEnd
    EPmain.sampleTest.sampEnd = sampStart;
end

set(EPmain.handles.sampleTest.sampStart,'String',EPmain.sampleTest.sampStart-EPmain.sampleTest.baseline+1);
set(EPmain.handles.sampleTest.sampEnd,'String',EPmain.sampleTest.sampEnd-EPmain.sampleTest.baseline);
set(EPmain.handles.sampleTest.msStart,'String',(EPmain.sampleTest.sampStart-EPmain.sampleTest.baseline+1)*(1000/EPmain.sampleTest.Fs));
set(EPmain.handles.sampleTest.msEnd,'String',(EPmain.sampleTest.sampEnd-EPmain.sampleTest.baseline)*(1000/EPmain.sampleTest.Fs));

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sampleTestSampEnd(src,eventdata) 
%change end of sample test template

global EPmain

if src==EPmain.handles.sampleTest.sampEnd
    EPmain.sampleTest.sampEnd=str2num(get(EPmain.handles.sampleTest.sampEnd,'String'))+EPmain.sampleTest.baseline;
else
    EPmain.sampleTest.sampEnd=round((str2num(get(EPmain.handles.sampleTest.msEnd,'String')))/(1000/EPmain.sampleTest.Fs))+EPmain.sampleTest.baseline;
end

if EPmain.sampleTest.sampEnd > -EPmain.sampleTest.baseline+length(EPmain.sampleTest.templateWaveform)
    EPmain.sampleTest.sampEnd = -EPmain.sampleTest.baseline+length(EPmain.sampleTest.templateWaveform);
end

if EPmain.sampleTest.sampStart > EPmain.sampleTest.sampEnd
    EPmain.sampleTest.sampStart = EPmain.sampleTest.sampEnd;
end

set(EPmain.handles.sampleTest.sampStart,'String',EPmain.sampleTest.sampStart-EPmain.sampleTest.baseline-1);
set(EPmain.handles.sampleTest.sampEnd,'String',EPmain.sampleTest.sampEnd-EPmain.sampleTest.baseline);
set(EPmain.handles.sampleTest.msStart,'String',(EPmain.sampleTest.sampStart-EPmain.sampleTest.baseline-1)*(1000/EPmain.sampleTest.Fs));
set(EPmain.handles.sampleTest.msEnd,'String',(EPmain.sampleTest.sampEnd-EPmain.sampleTest.baseline)*(1000/EPmain.sampleTest.Fs));

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function clearWorkingSet(src,eventdata)
%clear working set?
global EPdataset

button = questdlg('Clear Working Set?');

switch button
    case 'Yes'
        if exist([EPdataset.EPwork filesep 'EPwork'],'dir')
            fileList=dir([EPdataset.EPwork filesep 'EPwork']);
            for iFile=1:length(fileList)
                if ~any(strcmp(fileList(iFile).name,{'.','..','EPprefs.mat'}))
                    delete([EPdataset.EPwork filesep 'EPwork' filesep fileList(iFile).name]);
                end
            end
            EPdataset.dataset=cell(0);
            ep('start');
        else
            disp('Sorry, could not find the working directory.  Is the disk drive disconnected?');
        end
    case 'No'
        
    case 'Cancel'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function specTemplate(src,eventdata)
%load example spec text file.
global EPmain

switch EPmain.segment.spec.format
    case 1
        specFormat='EPM';
    case 2
        specFormat='TSP';
    otherwise
        disp('oops')
        return
end

[experimentFieldNames, experimentData, theHeaders, trialData] = ep_loadSpecFile(specFormat,'');
if isempty(theHeaders)
    return
end

if ~isempty(trialData)
    if isempty(EPmain.segment.spec.specTable) || ~isempty(setdiff(EPmain.segment.spec.specNames(1:end-1),theHeaders)) || ~isempty(setdiff(theHeaders,EPmain.segment.spec.specNames(1:end-1)))
        if ~isempty(EPmain.segment.spec.specTable)
            disp('The example file has different columns than that used for the spec table so clearing the spec table.');
        end
        EPmain.segment.spec.specNames=theHeaders;
        EPmain.segment.spec.columnEditable=repmat(false,1,length(theHeaders));
        EPmain.segment.spec.ColumnFormat=cell(0);
        EPmain.segment.spec.ColumnWidth=cell(0);
        for iCol=1:length(theHeaders)
            EPmain.segment.spec.ColumnFormat{iCol}='char';
            EPmain.segment.spec.ColumnWidth{iCol}='auto';
        end
        trialData=cellfun(@num2str,trialData,'UniformOutput',false);
        EPmain.segment.spec.specTable=trialData;
        set(EPmain.handles.segment.specTable,'Data',EPmain.segment.spec.specTable);
        EPmain.segment.spec.specNames=theHeaders;
        EPmain.segment.spec.specNames{end+1}='none';
        EPmain.segment.spec.excludeSpec=length(EPmain.segment.spec.specNames);
        %EPmain.segment.spec.specExcludeValues=unique(EPmain.segment.spec.specTable(:,EPmain.segment.spec.excludeSpec));
        EPmain.segment.spec.specExcludeValues=cell(0);
        EPmain.segment.spec.excludeSpec2=length(EPmain.segment.spec.specNames);
        EPmain.segment.spec.specExcludeValues2=cell(0);
        EPmain.segment.spec.excludeSpec3=length(EPmain.segment.spec.specNames);
        EPmain.segment.spec.specExcludeValues3=cell(0);
        %EPmain.segment.spec.specExcludeValues2=unique(EPmain.segment.spec.specTable(:,EPmain.segment.spec.excludeSpec2));
        for iSpec=1:8
            EPmain.segment.spec.specField(iSpec)=length(EPmain.segment.spec.specNames);
            EPmain.segment.spec.specLabels{iSpec}='';
        end
    else
        EPmain.segment.spec.specTable=trialData;
    end
end
if ~isempty(EPmain.segment.contData)
    changeSegmentDataset
end
ep('start');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeSpecName(src,eventdata)
%change the name of the exclusion spec name.

global EPmain

if src==EPmain.handles.segment.excludeSpec
    newSpec=get(EPmain.handles.segment.excludeSpec,'Value');
    
    if newSpec ~= EPmain.segment.spec.excludeSpec %don't do anything unless the setting was actually changed
        EPmain.segment.spec.excludeSpec=newSpec;
        if strcmp('none',EPmain.segment.spec.specNames{EPmain.segment.spec.excludeSpec})
            EPmain.segment.spec.specExcludeValues=cell(0);
            EPmain.segment.spec.excludeValue=1;
        else
            EPmain.segment.spec.specExcludeValues=unique(EPmain.segment.spec.specTable(:,newSpec));
            EPmain.segment.spec.excludeValue=1;
        end
    end
elseif src==EPmain.handles.segment.excludeSpec2
    newSpec=get(EPmain.handles.segment.excludeSpec2,'Value');
    
    if newSpec ~= EPmain.segment.spec.excludeSpec2 %don't do anything unless the setting was actually changed
        EPmain.segment.spec.excludeSpec2=newSpec;
        if strcmp('none',EPmain.segment.spec.specNames{EPmain.segment.spec.excludeSpec2})
            EPmain.segment.spec.specExcludeValues2=cell(0);
            EPmain.segment.spec.excludeValue2=1;
        else
            EPmain.segment.spec.specExcludeValues2=unique(EPmain.segment.spec.specTable(:,newSpec));
            EPmain.segment.spec.excludeValue2=1;
        end
    end
elseif src==EPmain.handles.segment.excludeSpec3
    newSpec=get(EPmain.handles.segment.excludeSpec3,'Value');
    
    if newSpec ~= EPmain.segment.spec.excludeSpec3 %don't do anything unless the setting was actually changed
        EPmain.segment.spec.excludeSpec3=newSpec;
        if strcmp('none',EPmain.segment.spec.specNames{EPmain.segment.spec.excludeSpec3})
            EPmain.segment.spec.specExcludeValues3=cell(0);
            EPmain.segment.spec.excludeValue3=1;
        else
            EPmain.segment.spec.specExcludeValues3=unique(EPmain.segment.spec.specTable(:,newSpec));
            EPmain.segment.spec.excludeValue3=1;
        end
    end
end

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clearSpec(segmentChange)
%clear spec settings.

global EPmain

EPmain.segment.spec.format=1;
EPmain.segment.spec.specTable=cell(0);
EPmain.segment.spec.columnEditable=[false];
EPmain.segment.spec.ColumnFormat=cell(0);
EPmain.segment.spec.ColumnWidth={'auto'};
EPmain.segment.spec.specNames={'none'};
EPmain.segment.spec.excludeSpec=1;
EPmain.segment.spec.excludeValue=1;
EPmain.segment.spec.specExcludeValues=cell(0);
EPmain.segment.spec.excludeValue2=1;
EPmain.segment.spec.specExcludeValues2=cell(0);
EPmain.segment.spec.excludeValue3=1;
EPmain.segment.spec.specExcludeValues3=cell(0);

for iSpec=1:8
    EPmain.segment.spec.specField(iSpec)=length(EPmain.segment.spec.specNames);
    EPmain.segment.spec.specLabels{iSpec}='';
end

if isfield(EPmain.segment,'cellTable') %not just entering segment function
    if ~segmentChange %avoid infinite loop
        if ~isempty(EPmain.segment.contData)
            changeSegmentDataset
        end
        ep('start');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function specLoad(src,eventdata) 
%load spec table

global EPmain

[FileName,PathName,FilterIndex] = uigetfile('','Load Spec Table');
if FileName ~= 0
    eval(['tempVar=load( ''' PathName FileName ''');']);
    if isfield(tempVar,'EPspecTable')
        if ~isempty(setdiff(tempVar.EPspecTable.specNames,EPmain.segment.spec.specNames))
            EPmain.segment.spec=tempVar.EPspecTable;
            EPmain.segment.spec.specTable=[];
        end
    else
        msg{1}='Error: not a spec table file.';
        [msg]=ep_errorMsg(msg);
    end
    if ~isfield(EPmain.segment.spec,'excludeValue3') %backward compatibility
        EPmain.segment.spec.excludeValue3=1;
        EPmain.segment.spec.specExcludeValues3=cell(0);
        EPmain.segment.spec.excludeSpec3=length(EPmain.segment.spec.specNames); 
    end
    if ~isempty(EPmain.segment.contData)
        changeSegmentDataset
    end
    ep('start');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function specSave(src,eventdata)
%save spec table

global EPmain EPdataset

[fileName,PathName,FilterIndex] = uiputfile('*.mat','Save Spec Table',[EPdataset.dataset(EPmain.segment.dataset).dataName '_specTable']);

if fileName ~= 0
    EPspecTable=EPmain.segment.spec;
    eval(['save ''' PathName fileName ''' EPspecTable']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updateSpec(src,eventdata,theSpec) 
%update the values of the spec settings on the specs subpane of the Segment function.

global EPmain 

EPmain.segment.spec.specField(theSpec)=get(EPmain.handles.segment.specNames(theSpec),'Value');
EPmain.segment.spec.specLabels{theSpec}=get(EPmain.handles.segment.specLabels(theSpec),'String');
if any(src==EPmain.handles.segment.specNames) && (EPmain.segment.spec.specField(theSpec) ~= length(EPmain.segment.spec.specNames))
    theSpecField=EPmain.segment.spec.specNames{EPmain.segment.spec.specField(theSpec)};
    if (length(theSpecField)>4) && strcmp(theSpecField(end-3:end),'.ACC')
        theLabel='ACC';
    elseif (length(theSpecField)>3) && strcmp(theSpecField(end-2:end),'.RT')
        theLabel='RT';
    elseif (length(theSpecField)>5) && strcmp(theSpecField(end-4:end),'.RESP')
        theLabel='RESP';
%     elseif (length(theSpecField)>11) && strcmp(theSpecField(end-10:end),'.OnsetDelay')
%         theLabel='EPoffset';
    else
        theLabel=theSpecField;
    end
    EPmain.segment.spec.specLabels{theSpec}=theLabel;
end

changeSegmentDataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resetColors(src,eventdata)
%reset the color values of the waveform plots.

global EPmain 

for iColor=1:EPmain.maxColors
    EPmain.preferences.view.color(iColor).RGB=EPmain.defaultColor(iColor).RGB;
    EPmain.preferences.view.color(iColor).lineSize=2;
end

ep('start')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeColor(src,eventdata,theColor,theParam)
%change a color value of the waveform plots.

global EPmain 

theValue=str2double(src.String);
switch theParam
    case {1;2;3}
        if ~isnan(theValue) && (theValue>=0) && (theValue<=1)
            EPmain.preferences.view.color(theColor).RGB(theParam)=theValue;
        end
    case 4
        if ~isnan(theValue) && (theValue>=0)
            EPmain.preferences.view.color(theColor).lineSize=theValue;
        end
    case 5
        EPmain.preferences.view.color(theColor).lineStyle=src.String{src.Value};
    otherwise
        disp('Programmer Error: changeColor value not recognized.')
end

ep('start')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changeAverageDataset(src,eventdata)
%change the average function model dataset.

global EPmain EPdataset

if EPmain.average.trialSpecDataset==0
    EPmain.average.trialSpecDataset=length(EPmain.average.datasetList);
else
    EPmain.average.trialSpecDataset=get(EPmain.handles.average.trialSpecDataset,'Value');
end

if isempty(EPmain.average.datasetList)
    EPmain.average.RTspec='no RT';
    EPmain.average.ACCspec='no ACC';
elseif ~isfield(EPmain.average,'RTspec') || ~isfield(EPmain.average,'ACCspec')
    if any(strcmpi('RT',EPdataset.dataset(EPmain.average.datasetList(EPmain.average.trialSpecDataset)).trialSpecNames))
        EPmain.average.RTspec='RT';
    else
        EPmain.average.RTspec='no RT';
    end
    if any(strcmpi('ACC',EPdataset.dataset(EPmain.average.datasetList(EPmain.average.trialSpecDataset)).trialSpecNames))
        EPmain.average.ACCspec='ACC';
    else
        EPmain.average.ACCspec='no ACC';
    end
    EPmain.average.latencySpec=1;
    EPmain.average.itemSpec=1;
end

ep('start')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function loadWindowtable(src,eventdata)
%load settings in Window table.

global EPmain EPdataset

[FileName,PathName,FilterIndex] = uigetfile('','Load Window Table');
if FileName ~= 0
    try
        eval(['tempVar=load( ''' PathName FileName ''');']);
    catch
        tempVar=[];
    end
    if isfield(tempVar,'EPwindowTable')
        cellNames=EPdataset.dataset(EPmain.window.dataset).cellNames(find(strcmp('SGL',EPdataset.dataset(EPmain.window.dataset).cellTypes)));
        [uCellNames i]=unique(cellNames,'first');
        if any(~ismember(tempVar.EPwindowTable.inCells,uCellNames))
            msg{1}='Error: cell names do not correspond to this dataset.';
            [msg]=ep_errorMsg(msg);
            return
        end
        EPmain.window.undo.outCells=EPmain.window.outCells;
        EPmain.window.undo.inCells=EPmain.window.inCells;
        EPmain.window.outCells=tempVar.EPwindowTable.outCells;
        EPmain.window.inCells=tempVar.EPwindowTable.inCells;
    else
        msg{1}='Error: not an Window table file.';
        [msg]=ep_errorMsg(msg);
    end
    ep('start');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function saveWindowtable(src,eventdata)
%save settings in Window table.

global EPmain

nameStart=strfind(EPmain.window.datasetName,': ');
if isempty(nameStart)
    nameStart=1;
else
    if length(nameStart) > 1
        nameStart=nameStart(1);
    end
    nameStart=nameStart+2;
end
datasetName=EPmain.window.datasetName(nameStart:end);

[fileName,PathName,FilterIndex] = uiputfile('*.mat','Save Window Table',[datasetName '_WindowTable']);

if fileName ~= 0
    EPwindowTable.outCells=EPmain.window.outCells;
    EPwindowTable.inCells=EPmain.window.inCells;
    eval(['save ''' PathName fileName ''' EPwindowTable']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resetWindowtable(src,eventdata)
%resets settings in Window table.

global EPmain EPdataset

EPmain.window.undo.factor=EPmain.window.inCells;
EPmain.window.undo.levels=EPmain.window.outCells;

cellNames=EPdataset.dataset(EPmain.window.dataset).cellNames(find(strcmp('SGL',EPdataset.dataset(EPmain.window.dataset).cellTypes)));
[u i]=unique(cellNames,'first');
EPmain.window.inCells=cellNames(sort(i));
EPmain.window.outCells=EPmain.window.inCells;

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function undoWindowtable(src,eventdata)
%undo settings in Window table.

global EPmain

tempVar=EPmain.window;
EPmain.window.inCells=EPmain.window.undo.inCells;
EPmain.window.outCells=EPmain.window.undo.outCells;
EPmain.window.undo.inCells=tempVar.inCells;
EPmain.window.undo.outCells=tempVar.outCells;

ep('start');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resetEP(src,eventdata)
%reset button.

global EPtictoc

MATLABver=ver('MATLAB');
[MatVer, MatRel]=strtok(MATLABver.Version,'.');
MatRel=MatRel(2:end);
[STA,~]=dbstack;

if ((str2double(MatVer) >25) || ((str2double(MatVer)==25) && (str2double(MatRel) >=1))) && strcmp(STA(2).file,'UIControlController.m')
    if isempty(dbstack(5))
        EPtictoc.stop=0;
        ep('start');
    else
        EPtictoc.stop=1;
    end
else
    if isempty(STA) || (isscalar(STA) && strcmp(STA(1).name,'resetEP'))
        EPtictoc.stop=0;
        ep('start');
    else
        EPtictoc.stop=1;
    end
end
