function [EPdataOut, origEloc, outInfo]=ep_readData(varargin)
%  [EPdataOut, origEloc, outInfo]=ep_readData(varargin);
%       imports data from a variety of file formats and converts into EP file format.
%
%Inputs:
%  OPTIONAL KEYWORDS (accepts either a single cell string or a series of variables)
%  Information in keywords override contents of the file.
%  'file'       : Name of the file to be imported.  It needs to have an
%                 appropriate suffix in order for it to be recognized properly.
%                 See the tutorial file.
%                 If none is specified it will assume that it is an EGI EGIA average file.
%
%  'format'     : File format of the file.
%  'samples'    : Samples included in analysis.
%  'cells'      : Cells included in analysis (trials for single_trial data).  Refers to order of appearance in the data file.
%                 Should be specified in order desired for output.  See example.
%  'subjects'   : Subjects included in analysis.  Should be specified in
%                 order desired for output.
%  'channels'   : Channels included in analysis.
%  'montage'    : The electrode montage.  If standard 10-05 then will just be 'Generic 10-05'.
%  'name'       :  The name of the experiment.
%  'labels'      : Renames the condition, subject, session, frequency, and trial names (e.g., 'labels',{{'cell1','cell2'},{'sub1','sub2'},{},{8,10},{1,2}}).
%                  If only one cell label is specified, it will be applied to all the data.
%                  Frequency labels should be a number that indicates the Hz or an empty set if not frequency data.
%                  Trial labels should be a number that indicates which trial or an empty set if not single-trial data.
%  'type'       : The type of the data.  Should be followed by the data type.
%                 'continuous' raw data from a single subject, not segmented into epochs.
%                 'single_trial' segmented single subject not averaged.
%                 'average' file with trials averaged so there is only one waveform per condition.
%                 'grand_average' file with multiple subjects averaged together.
%  'NumberSubjects'   : The number of subjects or factors.
%  'prestim'   : The number of msec prior to the stimulus onset (prior to using the samples keyword if used).
%  'ced'      : The name of the .ced file for electrode coordinates.  Use 'none' if none and should not ask for one.  Can be invalid ced file.
%  'eloc'      : The eloc information (only if good info).
%  'silent'   : If 'on' rather than 'off', do not provide output to command line for standard messages. (default off)
%  'origReference'  : Specify the original reference channel if present explicitly in the dataset.
%                 Followed by a field for the original reference(s), no more than two, in an array.
%  'currReference'  : Specify the current reference channel if present explicitly in the dataset.
%                 Followed by a field for the current reference(s), no more than two, in an array.
%                 The field can also be the key word 'AVG' for average reference and 'CSD' for current source
%                 density.
%  'textPrefs'  : Specifies the parameters for text files.  Structured variable:
%                 .firstRow : the first row of the data.
%                 .lastRow  : the last row or the data.  Zero equals last row.
%                 .firstCol : the first column of the data.
%                 .lastCol  : the last column of the data.  Zero equals last column.
%                 .orientation : the orientation of the data.  1=channels are columns, 2=channels are rows.
%                 .sampleRate : the sample rate of the data.
%
%  'SMIsuffix'  : Specifies add any SMI files with the same file name stem and this suffix to the data.
%  'specSuffix'  : Specifies add any spec text files with the same file name stem and this suffix to the data.
%                  The first line should contain the spec names.
%  'subjectSpecSuffix'  : Specifies add any subject spec text files with the same file name stem and this suffix to the data.
%                  The first line should contain the subject spec names.
%  'matlabDims'  : Vector of how the dimensions in a Matlab matrix file format map onto EP Toolkit data dimensions.
%
%  Reads in FieldTrip I/O outputs: data, hdr, eventHdr.  Assumes 'average' data is organized
%  with the subjects grouped together by cell (e.g., C1S1 C1S2 C1S3 C2S1 C2S2 C2S3).
%  Except for simple binary average files where it is assumed to be the other way around.
%  Recommended but not necessarily required is a .ced file (EEGlab format) for electrode coordinates.
%
%  If the file is continuous and there is a file with the same name and a
%  .asf prefix in the same directory, then it will be added to the .ept file.
%
%Outputs:
%  EPdataOut         : Structured array with the data and accompanying information.
%    .data      : 7D data matrix [channels, time points, cells/trials, subjects/sessions, factors, freqs, relations].
%                 Only theSGL factors are included in .data.  REG channels for spatial PCA data is handled in facVecS, not .data, so will always equal one for spatial PCA data.
%                 Virtual subjects not included.
%    .noise     : 6D matrix [channels, time points, cells/trials, subjects, factors, freqs] mirroring .data except factors includes both SGL and CMB.  Can be empty.
%                 This is the +/- reference average (Schimmel, 1967) which provides an estimate of the noise level
%                 in averaged data by flipping every other trial.  Not applicable to spectral data.  Not from combining data other than across subjects.
%                 For grand average data, this is the subject-wise mean of the noise estimates from the first level of averaging.
%    .covAVE    : 7D matrix [channels, time points, cells/trials, subjects, factors, freqs, channels] mirroring .data except factors includes both SGL and CMB.  Can be empty.
%                 This is the trialwise channel covariance matrix from averaging across trials or for grand averages it is the average of them.
%                 The diagonals are the variance of the channel and hence its square root is the standard deviation.  These are pooled variances so deviations from the sample means.
%                 If the seventh dimension is size one then there are no covariances and the numbers are variances.
%                 Not used for relational data.  For factor data, the dimension length reflects that of the original data, not the factor reduced number.
%    .GAVsubs   : Cell array (vSubs+1, vCells+1, all factors).  This field specifies virtual grand averages that are generated as needed.  This field may be empty.
%                 The columns correspond to each of the cells.  Each cell in turn contains a two-column array with the first column being the subject row
%                 that is included in the grand average and the second column contains the weight.  For a simple grand average, the weights will all be
%                 one.  For anything else, the weights should sum to zero.
%                 The first column specifies the subjects in all the normal cells (the +1).  The rest are for the virtual cells.
%                 The first row contains the cells contributing to the virtual waveforms, hence the +1.
%                 This is set up this way to accommodate trimmed grand averages from the robust statistics ANOVAadds functionality.
%                 The virtual grand averages consist of only the subNames and subTypes and the GAVsubs entries.
%                 Virtual grand averages may not index other virtual grand averages.
%                 The virtual cells consist of only the cellNames and the cellTypes and the GAVsubs entries.
%                 The virtual cells can also be averaged cells in single-trial datasets.  In this case, the first row specifies the trials going into the average.
%    .cov        (empty with no subfields if no info)
%      .covMatrix   : Variance-covariance matrix for channels, generated during averaging as indicator of noise levels.
%             Provides original covariance matrix.  Average referenced and not updated when rereferencing.  (subject,chan,chan)
%             NaN for unknown or invalid entries.  This is provided for use by the MNE software, which uses it to whiten
%             the data during source analyses (p.123 of 2.7.3 MNE manual).  Per p.89, epochs are first baseline corrected and then
%             the matrix is actually calculated as an SSCP (sums-of-squares-cross-products) matrix.
%       .Nq     : Number of observations going into the covariance matrix computation, for weighting when combining them. (subject)
%
%    .montage   : String with the montage information.  The default is 'Generic 10-05'
%    .fileFormat: The original file format.
%    .chanNames : The channel names.
%    .timeNames : The msec of the sample onset with respect to the stimulus onset time.
%                 For time-frequency data, the msec of the middle of the .5 sec window.
%                 For flexible segments, the percentage is still left-sided (so for 5% increments, from zero to 95%).
%    .timeUnits : 'ms' for milliseconds, 'per' for flexible segments.
%    .subNames  : The subject names.  There will be one entry per session.  Should be the same name for each session if there are multiple sessions.
%    .cellNames : The cell names (once for each trial for single_trial files).
%    .trialNames: The trial number ID per cell (starting from 1). (single_trial data only)  Will be empty if not single_trial data.
%                 These should be unique numbers for each cell.
%    .facNames  : The factor names (empty for non-factor data) including both SGL and CMB intermixed in any order desired.
%    .freqNames : The frequency at the middle of the frequency bin (numeric array).
%    .relNames  : The channel names for relational data (e.g., coherence).  Mirrors the channels dimension exactly.
%    .sessNames : The names of the sessions.  Can be empty.  There can be names without corresponding entries in sessNums but not vice versa.
%    .chanTypes : The type of the channel: EEG, MGM (magnetometer MEG), MGA (axial MEG), MGP (planar MEG), ECG, ANS (autonomic), REG (EEG regional average), BSC (BOSC test output)
%                 PPL (pupil dilation), XEY (x-coordinate of eye-tracking), YEY (y-coordinate of eye-tracking)
%               : ACMx ACMy ACMz accelerometer head position.
%                 The CED file can also specify FID (fiducial) and BAD (delete) channel types,
%                 but they will not end up as a chanType.
%    .subTypes  : The type of the subject: RAW (single trial), AVG (subject average), GAV (grand average).  Including vSubs.
%    .cellTypes : The type of the cell: SGL (one cell), CMB (combination of cells), STS (sample test statistics output).
%    .facTypes  : The type of the factor: SGL (one factor), CMB (combination of factors)
%    .sessNums  : The number of the session.  The numbers correspond to the sessNames.  This field should have the same length as subNames or can be empty.  Zero means missing or not applicable label.
%    .EPver     : The EP Toolkit version information.
%    .ver       : The Matlab version information.
%    .date      : The date the file was created.
%    .Fs        : The sampling frequency in Hz at the current data resolution (not necessarily the original sampling rate).  Equals number of %age increments times 10 for flexible segments (e.g., "200" if 20 bins).
%               : For non-integer sample sizes (i.e., 1024 Hz), round to get the samples.
%               : May be empty, especially for spectral data.
%    .baseline  : The number of samples prior to the trigger event (positive number).
%                 Thus, the ms in the timeNames refers to the offset of each sample.  For 250 Hz, a baseline of 50 samples means
%                 that the epoch started 200 ms (50 samples) prior to the stimulus onset.  A baseline of 1 sample means
%                 that the epoch started 4 ms (1 sample) prior to the stimulus onset.
%                 When baseline falls within the sample, shift to left side of sample.
%    .dataName  : A descriptive name for the dataset, used to differentiate the active datasets during analysis.
%    .ename     : The name of the experiment.
%    .dataType  : The type of the data: 'continuous', 'single_trial', or 'average' (default: average)
%    .trialSpecs       : Cell array (trial,spec,subject) of specific information for trials or their averages.  Can be empty if no specs.
%    .trialSpecNames   : Cell array of the name of each trial spec type.  Can be empty if no specs.
%    .subjectSpecs     : Cell array (subject,spec) of specific information for subjects.  Can be empty if no specs.
%    .subjectSpecNames : Cell array of the name of each subject spec type.  Can be empty if no specs.
%    .taskSpecs     : Array (subject,task,measure) of task performance for subjects.  Can be empty if no specs.
%    .taskNames : Cell array of the name of each task.  Can be empty if no specs.
%    .taskMeasNames : Cell array of the name of each task measure.  Can be empty if no specs.
%    .events    : Cell array of event structured variables (subject,cell/trial).  Events cells may be blank and thus have no fields.  For continuous data, cell/trial equals one.
%      .type      = string (events of "trial" type dropped as they are redundant with cell name information.)
%                   (if there was a recording stop, start of recording is marked by an event with the value 'boundary').
%      .sample    = expressed in samples, the first sample of an epoch is 1.  Out of range would be < 1 or >= numSamps+1.
%      .value     = number or string (e.g., "stm+")
%                   (if there was a recording stop, start of recording is marked by an event with the value 'boundary').
%      .duration  = expressed in samples.  For boundary events, represents the length of the recording break.
%      .keys     = additional information (subfields are required).  All field contents are strings even if the datatype is numeric.
%           .code (name)
%           .data (data)
%           .datatype (type of variable, e.g. 'short')
%           .description (notes)
%    .avgNum    : Number of waveforms going into averages (subject,cell)
%                 0 means unknown and -1 means bad.
%    .subNum    : Number of subjects going into averages (subject,cell)
%                 0 means unknown and -1 means bad.
%    .covNum    : Adjusted number of waveforms going into averages (subject,cell)
%                 0 means unknown and -1 means bad.
%                 For use with MNE software.  When computing linear combinations of waves, the effective sample size
%                 is modified to reflect the increase in noise levels, per p. 128 of the 2.7.3 MNE manual.
%    .fileName  : Name of original file.
%    .history   : cell array of changes to the dataset
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
%       4)  preferences settings.
%       5)  input files if any, as a cell array.
%    .ced       : The name of the .ced file for electrode coordinates.
%    .eloc      : The electrode location information, one for each channel (see readlocs header).  By apparent EEGlab convention, a row vector.  Can be empty.
%                 eloc is the same length as the channels, with non-EEG channels having a blank entry.
%                 For mff files, if internal eloc are present, will be used.
%                 Must have only the following fields and in this order: 'labels' 'theta' 'radius' 'X' 'Y' 'Z' 'sph_theta' 'sph_phi' 'sph_radius' 'type' 'cX' 'cY' 'cz'
%                 The cX, cY, and cZ fields are added by EP Toolkit and contain the canonical head space (Oostenveld space) coordinates.
%    .implicit  : The electrode information for fiducial locations (see readlocs header)
%                 Must have the same fields as for .eloc.
%    .facVecT   : For temporal PCA factor files, the factor waveform.  Used to compress the data. (rows=points,cols=factors)
%    .facVecS   : For spatial PCA factor files, the factor scalp topography.  Used to compress the data.(rows=chans,cols=factors).  Includes REG channels.
%    .facVecF   : For frequency PCA factor files, the factor frequency spectrum.  Used to compress the data.(rows=frequencies,cols=factors)
%    .facData   : 7D data matrix [channels, time points, cells/trials, subjects, factors, freqs, relations] including CMB factors (only when facVecS & facVecT & facVecF are used).
%    .facVar    : The variance accounted for by factors.  Unlike the contents of .pca, changed to reflect editing.  Includes CMB factors.  (1,factors).  May be empty even for factors.
%    .facVarQ   : The unique variance accounted for by factors.  Unlike the contents of .pca, changed to reflect editing.  Includes CMB factors.  (1,factors).  May be empty even for factors.
%                 They are represented separately from .data because they can't be compacted using the facVec mechanism.
%    .reference
%        .original    : Original recording reference(s): 1-2 numbers in array
%        .current     : Current reference channel(s): 1-2 numbers in array
%        .type        : Current reference type: REG (regular), AVG (average reference, even if no longer adds to zero),
%                       CSD (current source density), PAR (PARE-corrected average reference)
%    .analysis (for continuous files, divided into one second epochs and "trials" refer to these epochs. Excess time
%                points are tacked onto final epoch)
%       .blinkTrial : Array of blink-corrected trials (subject,cell/trial)
%       .saccadeTrial : Array of saccade-corrected trials (subject,cell/trial)
%       .saccadeOnset : Array of onset in samples of detected saccades (subject,cell/trial)
%       .moveTrial : Array of movement-corrected trials (subject,cell/trial)
%       .badTrials : Array of bad trials (subject,cell/trial).  1=bad.  For averages, it is the number of bad trials that went
%                       in.  See .avgNum for which are still bad after averaging.
%       .badChans : Array of corrected bad channels (subject,cell/trial,channel). -1 in a single-trial or continuous file means still bad.
%                   Negative numbers in an average file means number of still bad channels that went into the average
%                   (or rather, were left out of the average).  NaN in an average file means still bad.
%                   If the average waveform is corrected, then the number will be set to be equal to the number of trials going into the average.
%                   For the result of averaging and contrasts, the count is the number that went into the waveform regardless of weighting.
%    .pca
%        fields from ep_doPCA and ep_doPCAst steps.  See them for documentation.
%        Not affected by editing.
%    .recTime   : Time in samples of the start of the epoch (1-based) from the start of the session (cell)
%                 For averaged data, the earliest sample of the trials going into the average.  When unknown, first will be
%                 set arbitrarily to 1 and the remaining will be in order of the data file (e.g., 1, then 251, for 250 sample epochs, etc.).
%    .stims     : cell string containing images of the screens of stimuli used in the experiment as loaded in by the imread command.  The events structure
%                 provides the timing information, via a key field where the .code field contains "stim" and the .data field contains the name of the file.
%                 When empty, it should still have the subfields.
%       .name   : name of the stimulus file including the suffix
%       .image  : the image matrix itself.
%       .AOI    : areas of interest information.
%           .Word   : word in the AOI
%           .Coords : pixel coordinates of the AOI [upper left X, upper left Y, lower right X, lower right Y]
%           .tags : string array of stimulus tags
%    .calibration (can be empty with no subfields)
%       .ET          : Calibration values for eye-tracking channels
%           .Xzero   : zero correction for XEY channel
%           .Yzero   : zero correction for YEY channel
%           .Xscale  : scale correction for XEY channel
%           .Yscale  : scale correction for YEY channel
%       .SAC         : Calibration values for saccade channels
%           .Xzero   : zero correction for Hsaccade channel
%           .Yzero   : zero correction for Vsaccade channel
%           .Xscale  : scale correction for Hsaccade channel
%           .Yscale  : scale correction for Vsaccade channel
%     .impedances
%        .channels   : impedance values of the channels (chan,subject) (can be empty)
%        .ground     : impedance values of the ground electrode (subject) (can be empty)
%     .video         : co-registered video.  Not for average files.  (trials) (can be empty)
%        .frames     : Array of images.  Number of frames does not need to line up with time points of data (can be empty)
%           .cdata   : (height x width x RGB) (must be missing if no video data due to matlab idiosynchrasy)
%           .colormap: optional colormap for the images (must be missing if no video data due to matlab idiosynchrasy)
%        .times      : row vector of the times of the images in ms including prestimulus offset (can be empty)
%  origEloc          : The raw eloc prior to any editing.
%  outInfo           : Info from the reading process.
%        .matlabDims : Vector of how the dimensions in a Matlab matrix file format map onto EP Toolkit data dimensions.
%
% With respect to the FieldTrip data structure:
%  According to Vladimir Litvak,
%  "In case of events with duration that define trials event.sample is the first sample of a trial and event.offset is
%  the offset of the trigger with respect to the trial. An offset of 0 means that the first sample of the trial
%  corresponds to the trigger. A  positive offset indicates that the first sample is later than the trigger,
%  a negative offset indicates that the trial begins before the trigger."
%
%  According to the ft_read_event header documentation,the convention is to have a "trial" event that essentially keeps track of the recording time
%  for the first sample of the epoch and thus can be used to relate the event time to the proper sample in each epoch.
%  The .offset field is used for the "trigger" event only and indicates when the trigger really happened with respect to
%  the epoch (given that the trial .sample is being used to indicate the real time start time of the epoch rather than
%  the timing of the trigger event per se.
%
%  Since this approach does not lend itself well to the EP Toolkit environment, as of 2.40 the convention will be for the .sample field to be based on epoch time.
%  A recTime field will enable translation back to recording time and file time if needed.
%
%  EEGlab files don't really have set conventions but to the extent that there is one,
%  the events should have both a 'value' field to denote the generic type of event,
%  as in 'trigger', and a 'type' field to denote the nature of this generic event,
%  as in the condition of the experiment.  EEGlab essentially ignores
%  'value' fields and seems to have them mostly for purposes of using
%  FieldTrip's fileio for importing files.
%  The unofficial norm in FieldTrip analyses seems to be to use the generic
%  'type' field to subselect the event type of interest (like 'trigger')
%  and then use the 'value' field to define the conditions, so essentially the opposite to EEGlab.
%  The EP Toolkit hews to the FieldTrip convention with regard to the
%  'value' and the 'type' fields.
%
% With respect to spectral data:
%  There are a lot of things about spectral data that are merely conventions so there doesn't seem to be a definition or central authority to refer to.
%  After canvassing the literature, I've chosen to recognize the following as being accepted practice.
%  The internal spectral data are in complex form.  This is so that the imaginary component of the FFT can be represented for its phase information.
%  Prior to read-outs or analysis outputs, the data are always converted to spectral density form.
%  If the power (pw) or dB option is chosen, then the data are converted to power.
%  If the dB option is chosen, then the data are converted to dB form.  By definition, dB is always power.
%     Smith, SW (2003) Digital Signal Processing.  Newnes: Amsterdam.  pp 263-264.
%  By convention, the spectral density conversion is applied to the power form.
%  So for example, power spectral density is computed by dividing by bin herz width and amplitude spectral density is computed by dividing by the square root of bin herz width.
%  Spectral density is computed prior to conversion to dB.
%  Internally frequency data are represented as complex numbers to preserve the phase information (both the .data and the .freqData fields).
%  However, when data are added together (as in averaging) they need to be converted to absolute amplitude form otherwise opposite phase data will cancel out.
%  So such data are represented internally as absolute amplitude.  There is no flag other than the nature of the numbers themselves (complex or real).
%  The exception is adding channels together as it is appropriate for opposite phase to cancel out to reflect reference channel effects.
%
%  Example:
%  [EPdata]=readData;
%  [EPdata]=readData('file','NT9coll.ave_af.egis','cells',[1:3]);
%  [EPdata]=readData('format','ns_avg');

%History:
%  by Joseph Dien (2/7/08)
%  jdien07@mac.com
%
%  bugfix 3/31/08 JD
%  Handles formats other than EGIS properly by adding information about number of cells, subjects, and ordering.
%
%  modified 4/1/08 JD
%  Added support for EGI's segmented simple binary format.
%
%  modified 4/30/08 JD
%  accepts EGIS session files.
%
%  modified 11/4/08 JD
%  outputs structured variable including montage and file format and name
%  information.  Also, data is now a 4D data matrix.  Also, now autodetects
%  the electrode montage in EGIS files.
%
%  modified 11/10/08 JD
%  inputs changed to keyword list approach.
%
%  modified and bugfix (1/31/09) JD
%  Added initEP.  Added version and date fields to output structure.  Fixed reversal of order for subjects and cells
%  specified with format keyword.
%
%  modified (2/3/09) JD
%  Added Channels keyword.  Samples changed to set of numbers rather than just first and last.  Format no longer
%  includes fields for numbers of cells and subjects.
%
%  modified and bugfix (3/25/09) JD
%  Added information on date, version, sampling rate, baseline, experiment name, data type, and specs.
%  Subject names for EGIS average files now correctly detected.
%  Generalized to handle session files too with addition of trial specs, trial names, and events.
%  Also handles factor output files.  Tested out with Neuroscan files.
%  Added .ced field for electrode coordinate information.
%
%  modified (4/17/09) JD
%  Groups epochs in single trial files by cell.
%  eloc and chanType and implicit fields added to the files.
%
%  modified (4/27/09) JD
%  Added subNum field.  Merged 'subject_average', 'combined_average' and 'grand_average' data types into 'average'.
%  Deblanks cell names for EGIS files.
%
%  modified and bugfix 7/22/09 JD
%  Eliminated the latency field from the event structure.  Dropping initial samples results in event offset being
%  adjusted.   Dropped events of "trial" type as they are redundant with cell name information.  Added cellTypes field.
%  Fixed some crashes when handling simple binary files.  Added factors as 5th dimension.  Added support for text files.
%  Added support for EP files.  Added analysis subfields blinkTrial, moveTrial, badTrials, and badChans.
%  Added support for factor file compression via facVecT and facVecS fields.  Added dataName and facData fields.
%  EGIS subject specs converted to strings.
%
%  bugfix 8/4/09 JD
%  Added workaround for defective EGIS average headers generated by NetStation (the LastDone field which indicates number of subjects is incorrect).
%  Also, fix for NetStation not outputing subject numbers in the EGIS average header.
%
%  bugfix 8/29/09 JD
%  For "other" format files, not distinguishing cell names that start with the same characters (e.g., "1" and "11" not
%  distinguished).  Thanks to Grega Repovs.
%  For single_trial data, was incorrectly reporting on the command line the wrong number of cells being read in.
%
%  modified 8/30/09 JD
%  Added settings for text files.
%
%  bugfix 9/5/09 JD
%  Crash when reading factor file other than EGIS format.
%  Header defective (led to subsequent crash) when reading simple binary average file.
%  Crash when setting reference channel type.
%  Unable to find ced file when already specified in the file.
%  Text files with non-number columns sometimes not being read properly.
%
%  bugfix 9/16/09 JD
%  Text file import was skipping the first timepoint.
%  Not reading correct number of text columns.
%
%  bugfix 10/12/09 JD
%  Crash when reading "other" format files in single file mode where event values are empty.
%
%  bugfix 10/18/09 JD
%  Now aborts when NaN or inf values are in the data.
%
%  modified 10/31/09 JD
%  Reads in eloc information if present in eeglab_set file.
%
%  bugfix 11/9/09 JD
%  Number of subjects for continuous simple binary files incorrectly set.
%  Crash when reading continuous simple binary files.
%
%  modified 11/15/09 JD
%  When reading simple binary files, if the cell names cannot be deduced, put all of the segments into the same single
%  cell rather than just aborting.
%
% bugfix and modified 11/20/09 JD
% Replaced "union" commands with "unique" commands because certain situations caused the "union" command to crash in
% Matlab 2007.
% When importing a text file, ignore tabs at the end of the line.
% Drops CELL and TRSP events for all simple binary files since NetStation loses the associated information when exporting simple binary files.
%
%  modified 1/15/10 JD
%  When reading in an EGIS format file, if none of the EGI montages are chosen, then it will allow for a .ced file to be chosen.
%
%  modified 1/26/10 JD
%  When reading in file formats that label the channels, they will be reordered to match the order of the ced file.  Also, if any channels
%  are missing from the data file then they will be added as bad channels.  Also, can now accommodate channels that are in the data file
%  but not in the ced file.  The type field in the .ced files is now used since the
%  EEGlab bug was apparently fixed and it is now functional.  The type field must now be present in the .ced file and assumptions
%  will no longer be made about which ones are REF or FID types.
%
%  bugfix 1/30/10
%  Fixed incorrect subject ID being extracted from EGIS session headers.
%  Fixed crash due to events with empty value fields (now uses type field instead if value field is empty).
%  Fixed crash due to empty type fields in eloc information by assuming they are EEG channels.
%
%  bugfix 2/4/10
%  Fixed crash when loading factor file with CED information.
%  Fixed loss of some factor information when loading EP format factors file.
%
%  modified & bugfix 2/27/10 JD
%  When importing data, epochs that are entirely flat across all channels are marked as bad.
%  Fixed crash when two REF channels (as in M1-M2 mean mastoid channels).
%  Analysis fields no longer optional.
%  Fixed crash when loading in an EP file with ced information.
%  Fixed crash when loading in a .set file that does not have ced information included in its header.
%  Initializes subject spec names as cell rather than empty array.
%  When importing EP file format data, correctly checks for data type even for "factors" and "grand_average".
%  Eliminated chantype field for implicit channels.
%  Now reads the nsweeps field of Neuroscan AVG files to fill in the avgNum and subNum fields.
%  Implicit reference channels will no longer be marked as bad.
%
%  bugfix 3/5/10 JD
%  One dimensional cell array fields are now standardized to be column vectors.  Some file formats were causing crashes
%  in the Edit function because some of these fields were coming out as row vectors.
%
%  modified 3/15/10 JD
%  if .type field is missing from eloc information, then add it (assume channels are "EEG").
%
%  bugfix 3/27/10 JD
%  Fixed crash when loading in .txt file format data.
%  Fixed crash when loading in data where event information is empty.
%  Fixed text files treating all channels as being implicit when used with ced files with channel names different than
%  the default channel names.
%  For file formats with fixed channel orders, use channel names from ced file if available.
%  Tried to reinstate support for .edf file format.
%
%  bugfix 4/18/10 JD
%  Fixed wrong number of channel names when reading fixed order file formats (like EGIS) and the reference channel is implicit.
%
%  bugfix 5/1/10 JD
%  Fixed crash when reading data using single file mode with .set files which contain the .eloc information.
%  Fixed crash when reading data file that was originally in .set format and then was saved in EP file format.
%
%  modified 5/12/10 JD
%  For files with ced set to "eeglab", change to either name in chaninfo.filename field if present, else "none".
%
%  bugfix 5/22/10 JD
%  Fixed crash when importing .set file or EP file with unavailable or invalid ced file named in ced field.
%  Fixed channel selection not operating on channel coordinates eloc field.
%
%  bugfix 6/17/10 JD
%  Fixed not keeping FID channels in implicit channel info for fixed channel order file formats (e.g., EGIS, text).
%  Fixed electrode information not matching the data for second file onwards when reading files in single cell file
%  mode.
%  Can now be passed eloc information so don't have to access ced file on every pass for single cell file mode.
%  Fixed cell labels and sub labels not being applied to EP file formats.
%  Fixed crash for file formats with fixed-order channels when ced file has wrong number of channels.
%  Fixed not ignoring extra tab at end of line of text files, resulting in "not-a-number" errors.
%  Fixed crash when loading EP or .set file with name of original ced file in addition to eloc information.
%  Fixed losing the electrode coordinate information (eloc and and ced) when loading in a .study file.
%
%  bugfix 8/2/10 JD
%  Fixed crash if using "trial" events to determine name of cells and the .value field is empty.
%
%  bugfix 8/25/10 JD
%  Fixed not obtaining channel names from eloc info when provided by function call, as when merging files during Single File mode.
%  Fixed errors when importing .study files where the group or the session fields were left blank.
%
%  bugfix 10/3/10 JD
%  Fixed events in the final sample of a segment being assigned to the succeeding segment and causing a crash if the
%  segment was already the last one.
%  Fixed crash when the montage keyword was followed by a blank, as when preprocessing a batch of non-EGIS files
%  containing more than one data file.
%
%  modified 10/12/10 JD
%  For continuous files, analysis fields refer to one second epochs.
%
%  modified 10/16/10 JD
%  Added support for HEOG saccade correction.
%
%  bugfix 11/3/10 JD
%  Now accepts upper case file suffixes (e.g., .EGIS).
%
%  bugfix 12/6/10 JD
%  When reading in EGIS session files with implicit ref, adds it in explicitly, thus avoiding crash in eyeblink
%  correction routine.
%
%  modified 2/15/11 JD
%  Added support for CSV text files in addition to tab-delimited text files.
%
% modified 2/26/11 JD
% Added support for EGI Matlab files
%
% modified 5/21/11 JD
% Added support for 6th dimension of frequency for time-frequency analyses
%
%  bugfix 6/26/11 JD
%  Eliminated crash when loading EP format file with both electrode coordinates and a regional channel.
%
% bugfix 12/19/11 JD
% Fixed not assigning cell names to .set files generated by ERPlab.
%
% modified 1/26/12 JD
% Eliminated REF channel type and added reference field to keep track of original and current reference scheme.
%
% modified 1/29/12 JD
% When loading older EP files, if data had implicit reference and was EGI data, then converted to explicit reference.
%
% modified 2/22/12 JD
% Noise and std fields no longer optional.  Now set to empty if not used.
%
% modified 3/25/12 JD
% Added support for freq label when reading in text files.
%
% bugfix 6/6/12 JD
% Fixed crash when loading in EP format file with no implicit channels.
%
% bugfix 9/4/12 JD
% Fixed crash when loading older EP format files with factor data.
%
% modified 9/8/12 JD
% Improved ability to figure out the cell names of EEGlab files.
%
% modified 10/16/12 JD
% Added option to set last row to be imported in text files.
%
% bugfix 10/18/12 JD
% Fixed subNames field not necessarily being a column vector.
%
% modified 1/10/13 JD
% Added power field to the EP format.
%
% bugfix 1/17/13 JD
% Events assigned to wrong subject (off by one).
%
% modified 1/17/13 JD
% Added support for .set files generated by Widmann's pop_grandaverage function.
% Added support for ERPlab .erp file format.
%
% bugfix 1/18/13 JD
% Fixed erroneous "labels" error message when trying to load .study file.
%
% bugfix 1/23/13 JD
% Fixed problem that loading EP files with frequency PCA data results in damaged data file and error messages.
%
% bugfix 2/7/13 JD
% Fixed crash when loading in .study average file where the .set files have no prestimulus period.
%
% bugfix 2/20/13 JD
% Fixed crash when data file has no event information.
%
% bugfix 4/24/13 JD
% Better handles ced files where there are channel types other than EEG, FID, and REF.
%
% bugfix 5/8/13 JD
% Fixed waveforms in simple binary average files getting scrambled!  It no longer makes assumptions about the ordering
% of the cells.
% Fixed bug where if baseline was zero ms then instead it was being recorded as being 4 ms.
% Fixed not handling Neuroscan files with two physically linked explicit reference sites.
%
% modified 5/9/13 JD
% Added support for EGI's epoch-marked simple binary format for session files.
%
% bugfix 5/20/13 JD
% Fixed crashing when merging fixed channel files (like Neuroscan) where there are channels in the data that are not in the CED file.
%
% modified 5/20/13 JD
% Added option to eliminate unwanted channels, as in a GFP channel, in the data by marking the channel as BAD in the CED file.
% Added original eloc output.
%
% modified 5/26/13 JD
% Text files with multiple delimiters between values (as in space-space) now treated as a single delimiter.
%
% bugfix 6/24/13 JD
% Fixed crash when segmented simple binary session file is incorrectly specified to be an average file by the user.
%
% bugfix 9/18/13 JD
% Fixed incorrect inference of reference scheme when a single ref channel is designated.
%
% bugfix 9/20/13 JD
% Fixed channel type not changed from REF to EEG for file formats with flexible order channels.
% Fixed REF channels assumed to be last for fixed order channel file formats.
% General overhaul of CED code to better support MEG and ANS chan types, BAD CED code, and added ECG chan type.
%
% modified 9/22/13 JD
% Added support for mff files with recording stops and for ep_recordingStart
% events.
%
% bugfix 9/24/13 JD
% Workaround for EEGlab issue where if a CED file has just the label and the type filled out, the type info migrates over to the theta column for some reason.
%
% modified 9/26/13 JD
% Improved detection of Simple Binary files with scrambled cells.
%
% bugfix 10/2/13 JD
% Fixed crash if cell names are a mix of numbers and strings.
% Fixed crash if type field from CED file contains numbers for some reason.
%
% modified 10/10/13 JD
% Added recTime field.
% No longer rearranging single trial data to group by cell.
% Eliminated offset field from events structure.
% Fixed hitory field getting overwritten.
%
% modified 10/16/13 JD
% Epochs with boundary events are marked as bad.
%
% modified 10/21/13 JD
% Added support for reading rejected fields from EEGlab files.
% Added support for reading nTrials and chanlocs fields from ERPlab files.
%
% bugfix 10/21/13 JD
% Ensures that power comes after analysis field so order of fields is always the same.
%
% bugfix 10/29/13 JD
% Fixed crash when reading file with TRSP events that was not in mff format.
%
% modified 10/30/13 JD
% "boundary" in .type as well as .value fields for boundary events.
% Added .keys field to events.
% Added subject specs support for mff files.
%
% bugfix 11/14/13 JD
% Shifted data selection code to ep_selectData to update implementation of data selection protocols.
%
% bugfix 11/22/13 JD
% No longer redo eloc when loading in EP format files and ced is specified to be the same as what it already is.
% No longer extract trial spec names from mff continuous files.
%
% bugfix 12/24/13 JD
% Fixed crash when keyword list includes empty cells.
%
% bugfix 1/9/14 JD
% Fixed crash when reading mff file with more than one subject field.
%
% modified 1/19/14 JD
% When reading in boundary events in mff files, duration field contains the length of the recording pause.
%
% modified 2/26/14 JD
% pca no longer optional field.  No fields are optional.
%
% bugfix 2/27/14 JD
% Fixed crash for EGIS average files with custom cell header lengths.
%
% bugfix 3/13/14 JD
% Fixed crash when reading mff file with subject field where the field was left blank.
%
% modified 3/13/14 JD
% Added support for reading FIFF files.
% ced label for electrode coordinates provided by file (e.g., eeglab, MFF, FIFF formats) is "internal".
%
% modified 3/14/14 JD
% Eliminated file type check for EGIS files since NetStation generates average EGIS headers that are incorrectly marked as being session files.
%
% modified 3/16/14 JD
% Uses internal electrode coordinates provided by MFF and FIFF files.  Added elecPrefs.
%
% bufix 3/19/14 JD
% Fixed recTime field not including space for the 'all' cell in PCA output, resulting in crashes when edited.
%
% modified 3/24/14 JD
% Added .cov field.  Added bad channel info for FIFF files.
%
% bufix 4/8/14 JD
% Fixed keys field of events not being added when missing, resulting in EP files created by older versions of Toolkit
% not being usable.
% Fixed not putting factor variance information in correct location when loading PCA .ept files,
% resulting in "There are 0 factors that meet the minimum variance criterion" error messages when trying to autoPCA them.
%
% modified 4/9/14 JD
% Added relNames dimension to data to handle coherence and phase-locking measure.
%
% modified 4/14/14 JD
% Added .covNum field.
%
% bufix 4/25/14 JD
% Added conversion of REF channel type to EEG for older EP files.
%
% bufix 4/28/14 JD
% Fixed not able to load in EP files with frequency PCA data.
%
% bufix 5/21/14 JD
% Fixed crash when loading an ept file with no theta values for the electrode coordinates.
%
% bufix 6/12/14 JD
% Fixed blank keys field of events being produced without .key (e.g., .keys.keyCode instead of .keys.key.keyCode)
%
% bufix 6/29/14 JD
% Added fix for crash when trying to read mff files that erroneously label their COM channel as being a reference channel.
%
% modified 7/16/14 JD
% Simplified keys code field structure.
%
% bufix 7/31/14 JD
% Fixed .cov.Nq not being updated when stripping off subject adds.
%
% modified 8/24/14 JD
% Workaround for NetStation bug where EGIS average files have number of trials in cell in wrong spot in the cell header.
%
% modified 10/2/14 JD
% Added support for reading in event text files when present.
% Added support for bdf recording stop events (by marking them as 'boundary' events).
%
% bufix 3/20/15 JD
% Fixed crash when FontSize not provided as in reading .study files.
%
% modified 5/29/15 JD
% Added support for reading edf files with channels of varying sampling rates.
%
% bufix 7/5/15 JD
% Fixed crash when reading data with an event prior to first sample of data.
%
% bufix 9/25/15 JD
% Fixed crash when preferences set to rotate mff/eeglab electrode coordinates 180
% or 270 degrees.
%
% bufix 10/9/15 JD
% Fixed crash when mff average file has only one subject.
%
% modified 10/13/15 JD
% Standardized electrode labeling (e.g., E10) when reading mff average files.
%
% bufix 10/19/15 JD
% Fixed crash when continuous set file has a single eventHdr event and it
% has an NaN or numeric .value.
% Handles the NaN values that EEGlab seems to set channels to when it
% automatically edits them to being bad data.
%
% bufix 12/10/15 JD
% Fixed crash when cell name was undefined and was therefore set to "cell001"
%
% modified 12/18/15 JD
% SMI word files now expected to end in .txt rather than _smi.txt
%
% bufix 2/3/16 JD
% Check to see if SMI info has already been added.
%
% bufix 2/17/16 JD
% Fixed couldn't read text files if lastrow equalled zero and firstrow is larger than 1.
%
% modified 3/8/16 JD
% Eliminated support for .power field.
%
% modified 9/24/16 JD
% NaN values no longer set to zero and bad data flags set for non-EEG channels.
%
% modified 10/15/16 JD
% Added stims field.
%
% modified 11/5/16 JD
% Added support for reading subject spec text files.
%
% modified 11/13/16 JD
% Added support for .calibration field
%
% modified 3/16/17 JD
% Added support for New Segment events in brainvision files.
%
% bufix 5/19/17 JD
% Adds .stims and .calibration fields to older .ept PCA data files to avoid crashes in two-step PCA.
% When reading NS5 continuous mff files, now recognizes reference channel type correctly.
% Added option to flip the electrode locations for .mff and .fiff files.
% Fixed crash when rotating electrode coordinates 180 or 270 degrees for .mff and .fiff files.
%
% modified 5/30/17 JD
% When importing an average file with multiple cells with the same name, modifies the names to be unique rather than assuming single-file type.
%
% modified & bugfix 6/19/17 JD
% Added support for flexible segments with addition of the .timeUnits field.
% Presence of NaN in EEG channels no longer zeroed as bad data if original file was an .ept file (as fix for zeroing out plv transform).
%
% modified & bugfix  11/28/17 JD
% Eliminated restrictions on location of CED files.
% Added support for impedances field.
% Fixed not recognizing ECG channels in mff files.
%
% bufix 12/7/17 JD
% Fixed boundary events not being handled correctly in BrainVision eeg files.
% Fixed odd behavior due to Matlab bug in which str2num executes word strings which are function names.
%
% bugfix 3/1/18 JD
% Fixed unable to use eloc information from a batch of .ept files when they originally used different ced files.
%
% bugfix 3/18/18 JD
% Now handles situation where one of the directories has a space after its name.
% Fixed crash when .ept file has event with empty key.
%
% modified 5/13/18 JD
% No longer adds every read event to the .history field.
%
% modified 5/25/18 JD
% Added preference option to turn off BV header encoding.
%
% bugfix 6/12/18 JD
% Fixed not ensuring subject spec names and trial spec names are column vectors.
%
% modified & bugfix 6/22/18 JD
% Now detects when fif files have average reference projection and applies average reference to the channels.
% Fixed covNum and avgNum fields being set to cell arrays for .fif files, resulting in crashes later on.
%
% bugfix & modified 12/14/18 JD
% Added support for trial names to single file mode.
% Added support for specifying sample rate for text data.
% Fixed not reading subject specs when present as secondary file.
%
% modified 12/21/18 JD
% Shifted mff import procedure to FieldTrip's v3.
%
% bugfix 12/24/18 JD
% Handles multi-subject average mff files under new v3 mff code by making some assumptions.
%
% modified 1/8/19 JD
% FacVar and FacVarQ now include CMB factors.
%
% bugfix 1/13/19 JD
% Fixed problem that could keep average .fiff file from being read.
%
% bugfix 1/18/19 JD
% Fixed .fiff data not being scaled correctly whenever the eloc data were not initially provided.
%
% bugfix 2/17/19 JD
% Fixed mff not identifying cell names correctly, although using a heuristic that may not work for all files.
%
% modified 3/30/19 JD
% Added support for task level performance measures.
%
% bugfix 4/23/19 JD
% Implemented workaround for limitations in reading averaged mff files.
%
% modified 5/01/19 JD
% Added option to fix truncated EP headers.
%
% modified 5/13/19 JD
% Can also fix missing EP headers.
%
% bugfix 11/4/19 JD
% Will no longer keep asking for ced file in single file mode if user chooses cancel to first request.
% Updated reading of mff trial specs to work with the new mffmatlabio routine.
%
% modified 11/10/20 JD
% Added sessNums sessNames fields.
%
% modified & bugfix 12/30/19 JD
% Enabled reading of text files with a greater range of variations in character encoding, end-of-line markers, and field separation markers.
% Now allowing NaN in EEG channels.  I imagine it'll break some things but better than treating them as zeros I think.
% Upgraded support of std information by adding .covAVE and .GAVsubs fields and eliminating .std and .stdCM fields.
% Reworked initialization of new EP files.
% Fixed crash when applying Transform to a file with channel adds, as in preprocessed data.
%
% bugfix 1/20/20 JD
% Fixed crash when reading in text files.
% Fixed failure to read mff combined subject average files.
% Fixed crash when reading in frequency data with single file mode.
%
% modified 3/13/20 JD
% When there is no eloc information, automatically add canonical coordinates based on channel names, including EOG channels.
% Fixed under some conditions mff files will be read in as single precision, resulting in crashes.
% Improved support for edf files with varying sampling rates.
% Added mffmatlabio to external directory since FieldTrip does not have the latest version.
%
% bugfix & modified 3/19/20 JD
% Added support for reading fieldtrip and PTB files.
% Restored support for fiducial coordinates for mff files.
% Fixed crashes when reading EGIS average files and simple binary files.
%
% bugfix 4/9/20 JD
% Fixed fiducial information not being formatted correctly in the implicit field, resulting in crashes down the line.
% Fixed crash when loading EEGlab .study dataset.
%
% bugfix & modified 5/3/20 JD
% Made further improvements to PTB file format support.
% Fixed assumption that channel names derived from hdr.labels is a column-vector, resulting in failure to read when they are not.
%
% bugfix 5/25/20 JD
% Fixed crash when reading a .study file and the session field is empty.
% Fixed and updated support for EEGlab .study session, group, and run fields.
%
% modified 7/23/20 JD
% Removed preference option to rotate electrode coordinates when importing mff or eeglab files that have internal coordinates as no longer needed.
% Fixed crash when importing data from a data file with no event information.
% No longer drops TRSP events when reading in simple binary files as they can still be useful even though they are empty.
% Fixed error when reading simple binary files with multiple subjects.
%
% bugfix 10/14/20 JD
% Fixed crash when reading a .set average file.
% Fixed crash when a file has a mix of numeric and string values.
%
% bugfix 11/12/20 JD
% Fixed crash when using BV headers and the recording was cut short so it ends on an incomplete TRSP.
%
% bugfix & modified 12/29/20 JD
% Fixed crash when reading a continuous file with a boundary event that falls in the extra points that are appended to the last epoch.
% noInternal preference option now applies to .ept files as well.
%
% bugfix 2/4/21 JD
% Fixed crash when importing data file with no internal electrode coordinates and no CED file has been selected.
%
% modified 3/9/21 JD
% Eliminated legacy 'factors' file type from Read function.
% PCA data can now be represented without factor vector compression.
% Updated file suffix checking.
% Disabled automatic correction of file type to "continuous" from "single-trial" if there is only one trial.
%
% bugfix & modified 8/8/21 JD
% Fixed error message when reading in a data file where the file type has been designated as grand_average.
% Deprecated old SLAY code and EGI montage code.
% Added support for canonical 10-05 coordinates in eloc structure.
%
% bugfix 2/10/22 JD
% Fixed crash when loading continuous mff files.
%
% bugfix & modified 6/14/22 JD
% Fixed crash when merging files while applying a CED file with BAD channels.
% ChanTypes in .fif files are now ignored as they are apparently specified by users rather than being standardized.  Instead CED files will need to be used.
% Will now accommodate case where .fif file erroneously codes the EOG channels as being fiducial points.
%
% modified 7/24/22 JD
% Added support for reading Matlab .mat files.
%
% bugfix & modified 10/13/22 JD
% Fixed crash when loading continuous files in single-file mode and the Cell field has been left empty.
% Added support for Neuroelectrics .easy files.
%
% modified 1/23/23 JD
% Improved support for fieldtrip .mat files.
%
% modified 11/5/23 JD
% Added ability to fix BV Headers when recorded using both Response and Stimulus events rather than just Stimulus events.
%
% modified 11/17/23 JD
% For EGI montages, if not otherwise specified (as in mff files using internal electrode information rather than external ced file), assumes last channel is the original reference.
%
% bugfix 11/29/23 JD
% Minor fix to handling of Response/Stimulus coding for BV Headers to avoid error message.
%
% bugfix 3/8/24 JD
% When reading in an mff file, will only consider EEG channels when assuming the last is the reference.
%
% modified 4/18/25 JD
% Added support for reading in supplemental .asf video file.
% Added support for reading in supplemental .edf EyeLink eye-tracker file.
% Added support for reading in supplemental .sfp file.
% Added support for virtual grand averages.
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

global EPtictoc EPmain

EPdataOut=[];

EPdata=ep_newFile;
origEloc=[];
outInfo.matlabDims=[];

samples=[];
cells=[];
subjects=[];
channels=[];
eventHdr=[];
origRefChan=[];
currRefChan=[];

numSubs=[];
numRsubs=[];
numVsubs=0;
numFacs=1;
numFreqs=1;
numCells=[];
numWaves=[];
numRels=1;
cellLabels=[];
subLabels=[];
sessLabels=[];
freqLabels=[];
trialLabels=[];
silent='off';
refChan=[];
textPrefs.firstRow=1;
textPrefs.lastRow=0;
textPrefs.firstCol=1;
textPrefs.lastCol=0;
textPrefs.orientation=1;
textPrefs.sampleRate=250;
subjectsGrouped=0; %assume average data normally grouped by cells unless specified otherwise
SMIsuffix='';
specSuffix='';
subjectSpecSuffix='';
scrsz=[];
FontSize=[];
BVheader=0;
sevenDdata=[];
noInternal=0;
matlabDims=[];

if ~isempty(varargin)
    if isa(varargin{1},'cell') && nargin==1 %if keywords were input as a single cell string
        inputSet=varargin{1};
    else
        inputSet=varargin;
    end
else
    inputSet=[];
end

argNum=length(inputSet);
if mod(argNum,2) ~= 0
    msg{1}='The keywords need to all be in pairs, with a keyword followed by the keyword information.';
    [msg]=ep_errorMsg(msg);
    return
end

argCount=1;
while argCount <= argNum
    if isempty(inputSet{argCount})
        argCount=argCount+1;
    else
        switch inputSet{argCount}
            case 'file'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''file'' keyword must be followed by a file name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''file'' keyword must be followed by a file name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                EPdata.fileName=inputSet{argCount};
                argCount=argCount+1;
            case 'format'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''format'' keyword must be followed by a format name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''format'' keyword must be followed by a format name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                EPdata.fileFormat=inputSet{argCount};
                argCount=argCount+1;
            case 'samples'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''samples'' keyword must be followed by a list of samples.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~isnumeric(inputSet{argCount})
                    msg{1}='The ''samples'' keyword must be followed by a set of numbers (e.g., [1:250]).';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                samples=inputSet{argCount};
                argCount=argCount+1;
            case 'cells'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''cells'' keyword must be followed by a list of cells.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~isnumeric(inputSet{argCount})
                    msg{1}='The ''cells'' keyword must be followed by a set of numbers (e.g., [1:3]).';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                cells=inputSet{argCount};
                argCount=argCount+1;
            case 'subjects'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''subjects'' keyword must be followed by a list of subjects.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~isnumeric(inputSet{argCount})
                    msg{1}='The ''subjects'' keyword must be followed by a set of numbers (e.g., [1:3]).';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                subjects=inputSet{argCount};
                argCount=argCount+1;
            case 'channels'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''channels'' keyword must be followed by a list of channels.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~isnumeric(inputSet{argCount})
                    msg{1}='The ''channels'' keyword must be followed by a set of numbers (e.g., [1:100 105:129]).';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                channels=inputSet{argCount};
                argCount=argCount+1;
            case 'montage'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''montage'' keyword must be followed by a montage name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~ischar(inputSet{argCount}) && ~isempty(inputSet{argCount})
                    msg{1}='The ''montage'' keyword must be followed by a montage name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                EPdata.montage=inputSet{argCount};
                argCount=argCount+1;
            case 'name'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''name'' keyword must be followed by an experiment name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''name'' keyword must be followed by an experiment name.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                EPdata.ename=inputSet{argCount};
                argCount=argCount+1;
            case 'labels'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''labels'' keyword must be followed by a set of condition, subject, session, freq, and trial labels, all within a set.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~iscell(inputSet{argCount})
                    msg{1}='The ''labels'' keyword must be followed by a set of condition, subject, session, freq, and trial labels, all within a set.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if length(inputSet{argCount}) ~=5
                    msg{1}='The ''labels'' keyword must be followed by a set of condition, subject, session, freq, and trial labels, all within a set.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                cellLabels=inputSet{argCount}{1};
                subLabels=inputSet{argCount}{2};
                sessLabels=inputSet{argCount}{3};
                freqLabels=inputSet{argCount}{4};
                trialLabels=inputSet{argCount}{5};
                if ~iscell(cellLabels)
                    cellLabels=cellstr(cellLabels);
                end
                if ~iscell(subLabels)
                    subLabels=cellstr(subLabels);
                end
                if ~iscell(sessLabels)
                    sessLabels=cellstr(sessLabels);
                end
                if ~iscell(freqLabels)
                    freqLabels=cellstr(freqLabels);
                end
                if ~iscell(trialLabels)
                    trialLabels=cellstr(trialLabels);
                end
                argCount=argCount+1;
            case 'type'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''type'' keyword must be followed by a data type.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''type'' keyword must be followed by a data type.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~any(strcmp(inputSet{argCount},{'continuous','single_trial','average','grand_average', 'factors'}))
                    msg{1}='The data type must be one of the following: continuous, single_trial, average, grand_average, factors.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                EPdata.dataType=inputSet{argCount};
                argCount=argCount+1;
            case 'NumberSubjects'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''NumberSubjects'' keyword must be followed by number of subjects or factors.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~isnumeric(inputSet{argCount})
                    msg{1}='The ''NumberSubjects'' keyword must be followed by number of subjects or factors.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                numSubs=inputSet{argCount};
                argCount=argCount+1;
            case 'prestim'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''prestim'' keyword must be followed by the number of msec prior to the event onset.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~isnumeric(inputSet{argCount})
                    msg{1}='The ''prestim'' keyword must be followed by the number of msec prior to the event onset.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                EPdata.baseline=inputSet{argCount};
                argCount=argCount+1;
            case 'ced'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''ced'' keyword must be followed by the name of the .ced file.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''ced'' keyword must be followed by the name of the .ced file.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                EPdata.ced=inputSet{argCount};
                argCount=argCount+1;
            case 'eloc'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''eloc'' keyword must be followed by the electrode information.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~isempty(inputSet{argCount}) && ~isstruct(inputSet{argCount})
                    msg{1}='The ''eloc'' keyword must be followed by the electrode information.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                EPdata.eloc=inputSet{argCount};
                if ~isempty(EPdata.eloc)
                    EPdata.implicit=EPdata.eloc(1); %setup up implicit to have the same structure as eloc.
                    EPdata.implicit(1)=[];
                end
                argCount=argCount+1;
            case 'silent'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''silent'' keyword must be followed by on or off.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''silent'' keyword must be followed by on or off.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~any(strcmp(inputSet{argCount},{'on','off'}))
                    msg{1}='The ''silent'' keyword must be followed by on or off.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                silent=inputSet{argCount};
                argCount=argCount+1;
            case 'origReference'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''origReference'' keyword must be followed by the number of the original recording reference channel(s).';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                origRefChan=inputSet{argCount};
                argCount=argCount+1;
            case 'currReference'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''currReference'' keyword must be followed by the number of the current recording reference channel(s).';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                
                currRefChan=inputSet{argCount};
                argCount=argCount+1;
            case 'textPrefs'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''textPrefs'' keyword must be followed by a structured variable with the preferences.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~isstruct(inputSet{argCount})
                    msg{1}='The ''textPrefs'' keyword must be followed by a structured variable with the preferences.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~isfield(inputSet{argCount},'firstRow') || ~isfield(inputSet{argCount},'lastRow') || ~isfield(inputSet{argCount},'firstCol') || ~isfield(inputSet{argCount},'lastCol') || ~isfield(inputSet{argCount},'orientation') || ~isfield(inputSet{argCount},'sampleRate')
                    msg{1}='The structured variable must have all six fields included.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                textPrefs=inputSet{argCount};
                argCount=argCount+1;
            case 'SMIsuffix'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''SMIsuffix'' keyword must be followed by a string with the file name suffix.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''SMIsuffix'' keyword must be followed by a string with the file name suffix.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                SMIsuffix=inputSet{argCount};
                argCount=argCount+1;
            case 'specSuffix'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''specSuffix'' keyword must be followed by a string with the file name suffix.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''specSuffix'' keyword must be followed by a string with the file name suffix.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                specSuffix=inputSet{argCount};
                argCount=argCount+1;
            case 'subjectSpecSuffix'
                argCount=argCount+1;
                if argCount > argNum
                    msg{1}='The ''subjectSpecSuffix'' keyword must be followed by a string with the file name suffix.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                if ~ischar(inputSet{argCount})
                    msg{1}='The ''subjectSpecSuffix'' keyword must be followed by a string with the file name suffix.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                subjectSpecSuffix=inputSet{argCount};
                argCount=argCount+1;
            case 'screenSize'
                argCount=argCount+1;
                if (argCount > argNum) || ~isnumeric(inputSet{argCount}) || length(inputSet{argCount}) ~= 4
                    msg{1}='The ''screenSize'' keyword must be followed by a vector with the four screen size values in pixels.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                scrsz=inputSet{argCount};
                argCount=argCount+1;
            case 'FontSize'
                argCount=argCount+1;
                if (argCount > argNum) || ~isnumeric(inputSet{argCount}) || length(inputSet{argCount}) ~= 1
                    msg{1}='The ''FontSize'' keyword must be followed by a keyword must be followed by a number.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                FontSize=inputSet{argCount};
                argCount=argCount+1;
            case 'BVheader'
                argCount=argCount+1;
                if (argCount > argNum) || ~isnumeric(inputSet{argCount}) || length(inputSet{argCount}) ~= 1 || ~ismember(inputSet{argCount},[0 1 2])
                    msg{1}='The ''BVheader'' keyword must be followed by either zero or one or two.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                BVheader=inputSet{argCount};
                argCount=argCount+1;
            case 'noInternal'
                argCount=argCount+1;
                if (argCount > argNum) || ~isnumeric(inputSet{argCount}) || length(inputSet{argCount}) ~= 1 || ~ismember(inputSet{argCount},[0 1])
                    msg{1}='The ''noInternal'' keyword must be followed by either zero or one.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                noInternal=inputSet{argCount};
                argCount=argCount+1;
            case 'matlabDims'
                argCount=argCount+1;
                if (argCount > argNum) || ~isnumeric(inputSet{argCount}) || length(inputSet{argCount}) ~= 7 || ~all(ismember(inputSet{argCount},[1:7]))
                    msg{1}='The ''matlabDims'' keyword must be followed by a vector of the numbers one through seven.';
                    [msg]=ep_errorMsg(msg);
                    return
                end
                matlabDims=inputSet{argCount};
                argCount=argCount+1;
            otherwise
                keyword=inputSet{argCount};
                if isnumeric(inputSet{argCount})
                    keyword=num2str(inputSet{argCount});
                end
                msg{1}=[keyword ' is not a keyword.'];
                [msg]=ep_errorMsg(msg);
                return
        end
    end
end

fileFormat=EPdata.fileFormat;

if isempty(scrsz)
    scrsz = get(0,'ScreenSize');
    if length(scrsz) < 2
        scrsz=[1 1 800 600];
    end
end

if isempty(FontSize)
    FontSize = 10;
end

origEloc=EPdata.eloc;

if isempty(EPdata.fileName)
    [EPdata.fileName, pathname] = uigetfile('*.*','Open:');
    if ~ischar(EPdata.fileName)
        EPdata.fileName=[];
    end
    
    if (isempty(EPdata.fileName))
        msg{1}='No filename selected. You have to click on a name';
        [msg]=ep_errorMsg(msg);
        return
    end
    EPdata.fileName = [pathname EPdata.fileName];
else
    if ~exist(EPdata.fileName,'file')
        [pathstr, name, fileSuffix] = fileparts(EPdata.fileName);
        if exist([EPdata.fileName ' '],'file')
            msg{1}=['The file ' EPdata.fileName ' has a space after its name which needs to be deleted.'];
        elseif ~exist([pathstr],'dir')
            msg{1}=['One of the directories containing ' EPdata.fileName ' may have a space after its name which needs to be deleted.'];
        else
            msg{1}=['The file ' EPdata.fileName ' does not exist.'];
        end
        [msg]=ep_errorMsg(msg);
        return
    end
end

if exist(EPdata.fileName,'file') && ~exist(deblank(EPdata.fileName),'file')
    msg{1}=['The file ' EPdata.fileName ' has a space after its name which needs to be deleted.'];
    [msg]=ep_errorMsg(msg);
    return
end

%determine file type
[pathstr, name, fileSuffix] = fileparts(EPdata.fileName);
[pathstr2, name2, fileSuffix2] = fileparts(name);

[formatSuffix,formatName]=ep_fileExtensions(EPdata.fileFormat);

if ~strcmpi(EPdata.fileFormat,'stacked')
    if ~any(strcmpi(fileSuffix,formatSuffix))
        if length(formatSuffix)==1
            msg{1}=['The file format ' formatName ' must end in: ' formatSuffix{1}];
        else
            msg{1}=['The file format ' formatName ' must end in one of the following:'];
            for iSuffix=1:length(formatSuffix)
                msg{iSuffix+1}=formatSuffix{iSuffix};
            end
        end
        if any(strcmpi(fileSuffix2,formatSuffix)) && ~isempty(fileSuffix2)
            msg{end+1}=['It appears that you may have a hidden suffix at the end of your file name: ' fileSuffix '.'];
        end
        [msg]=ep_errorMsg(msg);
        return
    end
end

% if any(strcmp(EPdata.fileFormat,{'egi_egia','egi_egis'})) && ~any(strcmpi(fileSuffix,{'.egis','.ave','.gav','.raw','.ses'}))
%     msg{1}='EGIS file names must end in .egis, .ave, .gav, .raw, or .ses to be recognized.';
%     if any(strcmpi(fileSuffix2,{'.egis','.ave','.gav','.raw','.ses'}))
%         msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
%     end
%     [msg]=ep_errorMsg(msg);
%     return
% end
% 
% if strcmp(EPdata.fileFormat,'eeglab_erp') && ~strcmpi(fileSuffix,'.erp')
%     msg{1}='ERPlab file names must end in .erp to be recognized.';
%     if strcmpi(fileSuffix,'.erp')
%         msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
%     end
%     [msg]=ep_errorMsg(msg);
%     return
% end
% 
% if strcmp(EPdata.fileFormat,'ns_cnt') && ~strcmpi(fileSuffix,'.cnt')
%     msg{1}='Neuroscan continuous file names must end in .cnt to be recognized.';
%     if strcmpi(fileSuffix,'.cnt')
%         msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
%     end
%     [msg]=ep_errorMsg(msg);
%     return
% end
% 
% if strcmp(EPdata.fileFormat,'ns_eeg') && ~strcmpi(fileSuffix,'.eeg')
%     msg{1}='Neuroscan single-trial file names must end in .cnt to be recognized.';
%     if strcmpi(fileSuffix,'.eeg')
%         msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
%     end
%     [msg]=ep_errorMsg(msg);
%     return
% end
% 
% if strcmp(EPdata.fileFormat,'ns_avg') && ~strcmpi(fileSuffix,'.avg')
%     msg{1}='Neuroscan average file names must end in .avg to be recognized.';
%     if strcmpi(fileSuffix,'.avg')
%         msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
%     end
%     [msg]=ep_errorMsg(msg);
%     return
% end
% 
% if strcmp(EPdata.fileFormat,'egi_sbin') && ~any(strcmpi(fileSuffix,{'.sbin','.raw'}))
%     msg{1}='SBIN file names must end in .raw or .sbin to be recognized.';
%     if any(strcmp(fileSuffix2,{'.sbin','.raw'}))
%         msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
%     end
%     [msg]=ep_errorMsg(msg);
%     return
% end
% 
% if strcmp(EPdata.fileFormat,'ep_mat') && ~strcmpi(fileSuffix,'.ept')
%     msg{1}='EP file names must end in .ept to be recognized.';
%     if strcmp(fileSuffix2,'.ept')
%         msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
%     end
%     [msg]=ep_errorMsg(msg);
%     return
% end
% 
% if strcmp(EPdata.fileFormat,'egi_mff') && ~strcmpi(fileSuffix,'.mff')
%     msg{1}='mff file names must end in .mff to be recognized.';
%     if strcmp(fileSuffix2,'.mff')
%         msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
%     end
%     [msg]=ep_errorMsg(msg);
%     return
% end
% 
% if any(strcmp(EPdata.fileFormat,'edf')) && ~strcmpi(fileSuffix,'.edf')
%     msg{1}='edf file names must end in .edf to be recognized.';
%     if strcmp(fileSuffix2,'.edf')
%         msg{2}='It appears that you may have a hidden suffix at the end of your file name.';
%     end
%     [msg]=ep_errorMsg(msg);
%     return
% end
% 
% if isempty(EPdata.fileFormat)
%     switch fileSuffix
%         case '.sbin'
%             EPdata.fileFormat='egi_sbin';
%         case '.egis'
%             EPdata.fileFormat='egi_egia';
%         case '.mat'
%             EPdata.fileFormat='ep_mat';
%         otherwise
%             if strcmp(silent,'off')
%                 disp('Assuming file format is EGIS format.');
%             end
%             EPdata.fileFormat='egi_egia';
%     end
% end

if strcmp(EPdata.fileFormat,'egi_egis') && ~strcmp(EPdata.dataType,'single_trial')
    EPdata.dataType='single_trial';
    if strcmp(silent,'off')
        disp('EGIS session files are always single trial format.');
    end
end

if strcmp(EPdata.fileFormat,'egi_egia') && any(strcmp(EPdata.dataType,{'single_trial','continuous'}))
    EPdata.dataType='average';
    if strcmp(silent,'off')
        disp('EGIS average files are never single trial or continuous data types.  Defaulting to assume the file is an average file.');
    end
end

if strcmp(EPdata.fileFormat,'eeglab_erp') && any(strcmp(EPdata.dataType,{'single_trial','continuous'}))
    EPdata.dataType='average';
    if strcmp(silent,'off')
        disp('ERPlab files are never single trial or continuous data types.  Defaulting to assume the file is an average file.');
    end
end

if strcmp(EPdata.fileFormat,'biosemi_bdf') && any(strcmp(EPdata.dataType,{'average','grand_average'}))
    EPdata.dataType='continuous';
    if strcmp(silent,'off')
        disp('Only continuous and single trial bdf files are supported.  Defaulting to assume the file is a continuous file.');
    end
end

if strcmp(EPdata.fileFormat,'edf') && ~strcmp(EPdata.dataType,'continuous')
    EPdata.dataType='continuous';
    if strcmp(silent,'off')
        disp('Only continuous edf files are supported.  Defaulting to assume the file is a continuous file.');
    end
end

EPdata2=EPdata;
if strcmp(EPdata.fileFormat,'ep_mat')
    try
        ep_tictoc('ioStart');
        tempVar=load('-mat', EPdata.fileName);
        ep_tictoc('ioFinish');
        
        if isfield(tempVar,'EPdata')
            EPdata=tempVar.EPdata;
        else
            msg{1}=['The file ' EPdata.fileName 'did not contain EP Toolkit data in it.'];
            [msg]=ep_errorMsg(msg);
            return
        end
        clear tempVar;
    catch ME
        msg{1}=['The attempt to load in the file ' EPdata.fileName ' resulted in the error:' ME.identifier];
        msg{2}=ME.message;
        [msg]=ep_errorMsg(msg);
        return
    end
    try
        EPver=ver('EP_Toolkit');
        if str2num(EPdata.EPver.Version) > str2num(EPver.Version)
            msg{1}=['The file ' EPdata.fileName 'was created by a more recent version of the EP Toolkit (' EPdata.EPver.Version ') than is running on this computer (' EPver.Version ') and therefore cannot be read.'];
            [msg]=ep_errorMsg(msg);
            return
        end
    catch
        EPver='unavailable'; %workaround for bug in earlier version of Matlab
    end
    
    %add in information provided via keywords
    if ~isempty(EPdata2.fileName)
        EPdata.fileName=EPdata2.fileName;
    end
    if ~isempty(EPdata2.fileFormat)
        EPdata.fileFormat=EPdata2.fileFormat;
    end
    if ~isempty(EPdata2.montage)
        EPdata.montage=EPdata2.montage;
    end
    if ~isempty(EPdata2.ename)
        EPdata.ename=EPdata2.ename;
    end
    if ~isempty(EPdata2.dataType)
        EPdata.dataType=EPdata2.dataType;
    end
    if ~isempty(EPdata2.baseline)
        EPdata.baseline=EPdata2.baseline;
    end
    if ~isempty(EPdata2.ced)
        EPdata.ced=EPdata2.ced;
    end
    if ~isempty(EPdata2.eloc)
        EPdata.eloc=EPdata2.eloc;
        EPdata.implicit=EPdata2.implicit;
    end
    
    %backward compatibility conversions
    EPdata=ep_updateEPfile(EPdata);
    if isempty(EPdata)
        return
    end
    
    %extract contents of the EP file.
    
    sevenDdata=EPdata.data;
    numSubs=length(EPdata.subNames);
    numVsubs=max(0,size(EPdata.GAVsubs,1)-1);
    numRsubs=numSubs-numVsubs;
    numWaves=length(EPdata.cellNames);
    numCells=length(unique(EPdata.cellNames));
    numVcells=max(0,size(EPdata.GAVsubs,2)-1);
    numRcells=numCells-numVcells;
    numFreqs=length(EPdata.freqNames);
    numRels=length(EPdata.relNames);
    
    ep_tictoc;if EPtictoc.stop;return;end
    if isfield(EPdata,'implicit')
        if ~isempty(EPdata.implicit)
            if isfield(EPdata.implicit,'chantype')
                EPdata.implicit=rmfield(EPdata.implicit,'chantype'); %backward compatibility for older files with a chantype field.
            end
            implicitTypes={EPdata.implicit.type};
            implicitChan=find(strcmp('REF',implicitTypes));
            if ~isempty(EPdata.montage)
                if (~isempty(strfind(EPdata.montage,'GSN')) || ~isempty(strfind(EPdata.montage,'Hydrocel'))) && (isscalar(implicitChan)) %if implicit reference and was EGI data, make explicit
                    
                    sevenDdata(end+1,:,:,:,:)=zeros(size(sevenDdata,2),size(sevenDdata,3),size(sevenDdata,4),size(sevenDdata,5),size(sevenDdata,6),size(sevenDdata,7));
                    if ~isempty(EPdata.facVecS)
                        EPdata.facVecS(end+1,:)=zeros(size(EPdata.facVecS,2));
                    end
                    
                    if ~isempty(EPdata.facData)
                        EPdata.facData(end+1,:,:,:,:)=zeros(size(EPdata.facData,2),size(EPdata.facData,3),size(EPdata.facData,4),size(EPdata.facData,5),size(EPdata.facData,6),size(EPdata.facData,7));
                    end
                    
                    EPdata.analysis.badChans(:,:,end+1)=zeros(size(EPdata.analysis.badChans,1),size(EPdata.analysis.badChans,2));
                    
                    if ~isempty(EPdata.noise)
                        EPdata.noise(end+1,:,:,:,:)=zeros(size(EPdata.noise,2),size(EPdata.noise,3),size(EPdata.noise,4),size(EPdata.noise,5));
                    end
                    if ~isempty(EPdata.covAVE)
                        if size(EPdata.covAVE,7)==1
                            EPdata.covAVE(end+1,:,:,:,:,:,:)=zeros(1,size(EPdata.covAVE,2),size(EPdata.covAVE,3),size(EPdata.covAVE,4),size(EPdata.covAVE,5),size(EPdata.covAVE,6),1);
                        else
                            EPdata.covAVE(end+1,:,:,:,:,:,:)=zeros(1,size(EPdata.covAVE,2),size(EPdata.covAVE,3),size(EPdata.covAVE,4),size(EPdata.covAVE,5),size(EPdata.covAVE,6),size(EPdata.covAVE,7)+1);
                        end
                    end
                    
                    EPdata.chanTypes{end+1}='EEG';
                    EPdata.chanNames{end+1}=EPdata.implicit(implicitChan).labels;
                    EPdata.eloc(end+1)=EPdata.implicit(implicitChan);
                    EPdata.implicit(implicitChan)=[];
                    EPdata.reference.current(end+1)=length(EPdata.eloc);
                    EPdata.reference.original(end+1)=length(EPdata.eloc);
                end
            end
        end
    end
    
elseif (strcmp(EPdata.fileFormat,'eeglab_set') && strcmp(fileSuffix,'.study'))
    try
        ep_tictoc('ioStart');if EPtictoc.stop;return;end
        eval(['load(''-mat'',  ''' EPdata.fileName ''');']);
        ep_tictoc('ioFinish');if EPtictoc.stop;return;end
    catch ME
        msg{1}=['The attempt to load in the file ' EPdata.fileName 'resulted in the error:' ME.identifier];
        msg{2}=ME.message;
        [msg]=ep_errorMsg(msg);
        return
    end
    
    fileList=cell(length(STUDY.datasetinfo),1);
    
    if (length(fileList) > 1) && strcmp(EPdata.dataType,'continuous')
        msg{1}='Multiple files cannot be combined into a single continuous data file.  You will need to load them in separately.';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    if (length(unique({STUDY.datasetinfo.subject})) > 1) && strcmp(EPdata.dataType,'single_trial')
        msg{1}='Multiple subjects cannot be combined into a single single-trial data file.  You will need to load them in separately.';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    if (length(unique([STUDY.datasetinfo.session])) > 1) && strcmp(EPdata.dataType,'single_trial')
        msg{1}='Multiple sessions cannot be combined into a single single-trial data file.  You will need to load them in separately.';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    for setFile=1:length(STUDY.datasetinfo)
        theFile=[STUDY.datasetinfo(setFile).filepath filesep STUDY.datasetinfo(setFile).filename];
        if ~exist(theFile,'file')
            theFile=STUDY.datasetinfo(setFile).filename;
            if ~exist(theFile,'file')
                msg{1}=['The file ' theFile ' is missing.'];
                [msg]=ep_errorMsg(msg);
                return
            end
        end
        fileList{setFile}=theFile;
    end
    mergeName=STUDY.name;
    if isempty(mergeName)
        mergeName='mergedData';
    end
    
    ep_tictoc;if EPtictoc.stop;return;end
    mergeArg=[];
    for file=1:length(fileList)
        mergeArg{file,1}=fileList{file};
        mergeArg{file,2}='format';
        mergeArg{file,3}='eeglab_set';
        mergeArg{file,4}='type';
        mergeArg{file,5}=EPdata.dataType;
        mergeArg{file,6}='labels';
        
        theSession=num2str(STUDY.datasetinfo(file).session);
        if isempty(theSession)
            theSession='';
        end
        if isfield(STUDY.datasetinfo(file),'run')
            theRun=num2str(STUDY.datasetinfo(file).run);
            if ~isempty(theRun)
                if ~isempty(theSession)
                    theSession=[theSession '-' theRun];
                else
                    theSession=theRun;
                end
            end
        end
        theCondition=STUDY.datasetinfo(file).condition;
        if isempty(theCondition)
            theCondition='Cond';
        end
        theSubject=STUDY.datasetinfo(file).subject;
        if isempty(theSubject)
            theSubject='1';
        end
        if isfield(STUDY.datasetinfo(file),'group')
            theGroup=STUDY.datasetinfo(file).group;
            if ~isempty(theGroup)
                theSubject=[theSubject '-' theGroup];
            end
        end
        
        mergeArg{file,7}={theCondition theSubject theSession cell(0) cell(0)};
        mergeArg{file,8}='FontSize';
        mergeArg{file,9}=FontSize;
    end
    
    EPdata2=EPdata;
    [EPdata]=ep_mergeEPfiles(mergeArg,mergeName);
    if isempty(EPdata)
        return;
    end
    
    %add in information provided via keywords
    if ~isempty(EPdata2.fileName)
        EPdata.fileName=EPdata2.fileName;
    end
    if ~isempty(EPdata2.fileFormat)
        EPdata.fileFormat=EPdata2.fileFormat;
    end
    if ~isempty(EPdata2.montage)
        EPdata.montage=EPdata2.montage;
    end
    if ~isempty(EPdata2.ename)
        EPdata.ename=EPdata2.ename;
    end
    if ~isempty(EPdata2.dataType)
        EPdata.dataType=EPdata2.dataType;
    end
    if ~isempty(EPdata2.baseline)
        EPdata.baseline=EPdata2.baseline;
    end
    if ~isempty(EPdata2.ced)
        EPdata.ced=EPdata2.ced;
    end
    if ~isempty(EPdata2.eloc)
        EPdata.eloc=EPdata2.eloc;
        EPdata.implicit=EPdata2.implicit;
    end
    
    %update variables
    hdr=[];
    sevenDdata=EPdata.data;
    numSubs=length(EPdata.subNames);
    numWaves=length(EPdata.cellNames);
    numCells=length(unique(EPdata.cellNames));
    numFreqs=length(EPdata.freqNames);
    numRels=length(EPdata.relNames);
    
elseif strcmp(EPdata.fileFormat,'PTB_mat')
    try
        ep_tictoc('ioStart');
        tempVar=load(EPdata.fileName);
        ep_tictoc('ioFinish');
    catch ME
        msg{1}=['The attempt to load in the file ' EPdata.fileName ' resulted in the error:' ME.identifier];
        msg{2}=ME.message;
        [msg]=ep_errorMsg(msg);
        return
    end
    
    if isfield(tempVar,'erp')
        if isfield(tempVar.erp,'subs')
            if isfield(tempVar.erp.subs,'name')
                EPdata.subNames=tempVar.erp.subs.name(:);
            end
        end
        if isfield(tempVar.erp,'stim')
            if isfield(tempVar.erp.stim,'study')
                sessList=unique(tempVar.erp.stim.study);
                for iSess=1:length(sessList)
                    EPdata.sessNames{iSess,1}=num2str(sessList(iSess));
                end
            end
        end
        EPdata.chanNames=tempVar.erp.elecnames(:);
        EPdata.Fs=tempVar.erp.samplerate;
        EPdata.baseline=tempVar.erp.tbin-1;
        subList=unique(tempVar.erp.subnum);
        subList=setdiff(subList,0);
        numSubs=length(subList);
        if (numSubs > 1) && ~strcmp(EPdata.dataType,'average')
            msg{1}=['EP Toolkit does not currently support data files with multiple subjects other than average files.'];
            [msg]=ep_errorMsg(msg);
            return
        end
        numWaves=max(tempVar.erp.sweep);
        numChan=length(EPdata.chanNames);
        numSamples=size(tempVar.erp.data,2);
        EPdata.trialSpecNames{1}='ttype';
        EPdata.trialSpecNames{2}='ACC';
        EPdata.trialSpecNames{3}='RT';
        EPdata.trialSpecNames{4}='resp';
        EPdata.trialSpecs=cell(numWaves,4,numSubs);
        EPdata.events=cell(numSubs,numWaves);
        EPdata.events{1,1}=struct('value','','sample',[],'type','','duration',[],'keys',struct('code','','data','','datatype','','description',''));
        EPdata.timeNames(:,1)=(((0:numSamples-1)-EPdata.baseline)*(1000/EPdata.Fs));
        EPdata.recTime=[1:EPdata.Fs:EPdata.Fs*(numWaves-1)+1]';
        EPdata.subNum=zeros(numSubs,numWaves);
        EPdata.covNum=zeros(numSubs,numWaves);
        for iCell=1:numWaves
            EPdata.cellNames{iCell,1}='missing';
        end
        
        %group subjects together by cell.
        sevenDdata=nan(numChan,numSamples,numWaves,numSubs);
        dataCount=zeros(numChan,numWaves,numSubs);
        if strcmp(EPdata.dataType,'average')
            EPdata.analysis.badChans=nan(numSubs,numWaves,numChan);
            EPdata.analysis.badTrials=zeros(numSubs,numWaves);
            EPdata.avgNum=ones(numSubs,numWaves)*(-1);
        else
            EPdata.analysis.badChans=ones(numSubs,numWaves,numChan)*(-1);
            EPdata.analysis.badTrials=ones(numSubs,numWaves);
        end
        
        for iTrial=1:size(tempVar.erp.data,1)
            theChan=tempVar.erp.elec(iTrial);
            theWave=tempVar.erp.sweep(iTrial);
            if theWave>0
                theSub=find(ismember(subList,tempVar.erp.subnum(iTrial)));
                switch tempVar.erp.domain
                    case 'time'
                        sevenDdata(theChan,:,theWave,theSub)=tempVar.erp.data(iTrial,:);
                    case 'freq-amplitude'
                        sevenDdata(theChan,1,theWave,theSub,1,:)=tempVar.erp.data(iTrial,:);
                    otherwise
                        msg{1}=['The file ' EPdata.fileName ' has an unsupported domain: ' tempVar.erp.domain '.'];
                        [msg]=ep_errorMsg(msg);
                        return
                end
                if dataCount(theChan,theWave,theSub)==1
                    fprintf('Warning: there were duplicate data entries for channel %d, subject %d, sweep %d.\n',theChan,theWave,theSub);
                end
                dataCount(theChan,theWave,theSub)=dataCount(theChan,theWave,theSub)+1;
                switch EPdata.dataType
                    case {'continuous','single_trial'}
                        EPdata.analysis.badChans(theSub,theWave,theChan)=0;
                        EPdata.analysis.badTrials(theSub,theWave)=0;
                    case 'average'
                        EPdata.analysis.badChans(theSub,theWave,theChan)=0;
                        EPdata.avgNum(theSub,theWave)=0;
                end
                %if there are conflicting meta information, will just take the last one given.
                if tempVar.erp.accept(iTrial)==0
                    switch EPdata.dataType
                        case {'continuous','single_trial'}
                            EPdata.analysis.badChans(theSub,theWave,theChan)=1;
                        case 'average'
                            EPdata.analysis.badChans(theSub,theWave,theChan)=NaN;
                    end
                end
                if isfield(tempVar.erp,'stim')
                    if isfield(tempVar.erp.stim,'catcodes')
                        theCellName=tempVar.erp.stim.catcodes(iTrial);
                        if isempty(theCellName)
                            theCellName='missing';
                        end
                        EPdata.cellNames{theWave,1}=theCellName;
                    end
                end
                if isfield(tempVar.erp,'stim')
                    if isfield(tempVar.erp.stim,'study')
                        EPdata.sessNums(theSub,1)=find(strcmp(num2str(tempVar.erp.stim.study(iTrial)),EPdata.sessNames));
                    end
                end
                EPdata.trialSpecs{theWave,1,theSub}=tempVar.erp.ttype(iTrial);
                EPdata.trialSpecs{theWave,2,theSub}=tempVar.erp.correct(iTrial);
                EPdata.trialSpecs{theWave,3,theSub}=tempVar.erp.rt(iTrial);
                EPdata.trialSpecs{theWave,4,theSub}=tempVar.erp.response(iTrial);
            end
        end
    else
        msg{1}=['The file ' EPdata.fileName ' did not contain Psychophysiology Toolbox data in it.'];
        [msg]=ep_errorMsg(msg);
        return
    end
    clear tempVar;

else
    %data imported via FieldTrip data structures or otherwise starting with partial structure so convert first into FieldTrip format to standardize including 3D theData matrix
    
    ep_tictoc;if EPtictoc.stop;return;end
    %determine data type
    if strcmp(EPdata.fileFormat,'egi_egis') && isempty(EPdata.dataType)
        EPdata.dataType='single_trial';
    end
    
    if strcmp(EPdata.fileFormat,'egi_egia') && isempty(EPdata.dataType)
        EPdata.dataType='average';
    end
    
    if strcmp(EPdata.fileFormat,'ns_avg') && isempty(EPdata.dataType)
        EPdata.dataType='average';
    end
    
    if strcmp(EPdata.fileFormat,'ns_eeg') && isempty(EPdata.dataType)
        EPdata.dataType='single_trial';
    end
    
    if strcmp(EPdata.fileFormat,'ns_cnt') && isempty(EPdata.dataType)
        EPdata.dataType='continuous';
    end
    
    if strcmp(EPdata.fileFormat,'biosemi_bdf')
        EPdata.dataType='continuous';
    end
    
    if isempty(EPdata.dataType)
        EPdata.dataType='single_trial';
        if strcmp(silent,'off')
            disp('Defaulting to assuming the file is single trial segmented data from one subject.');
        end
    end
    
    %read in data, generating theData, hdr, and eventHdr
    switch EPdata.fileFormat
        case 'text'
            ep_tictoc('ioStart');if EPtictoc.stop;return;end
            [theHeader, theRawData, theDelim] = ep_textScan(EPdata.fileName,textPrefs.firstRow,textPrefs.lastRow,textPrefs.firstCol,textPrefs.lastCol);
            ep_tictoc('ioFinish');if EPtictoc.stop;return;end
            if isempty(theRawData)
                msg{1}='No data were read.  Did you check the EP Toolkit preference settings for reading text files to make sure they were set correctly?';
                [msg]=ep_errorMsg(msg);
                return
            end
            theData=cellfun(@str2double,theRawData);
            
            if textPrefs.orientation ==1
                theData=theData';
            end
            
            hdr.nChans=size(theData,1);
            hdr.nSamples=size(theData,2);
            hdr.nSamplesPre=0;
            hdr.nTrials=0;
            hdr.Fs=textPrefs.sampleRate;
            if ~isfield(hdr,'label')
                for i = 1:hdr.nChans
                    hdr.label{i}  = ['e' num2str(i)];
                end
            end
            if ~isempty(freqLabels) && ~isempty(freqLabels{1})
                if iscell(freqLabels) && ~isempty(freqLabels{1})
                    EPdata.freqNames(1)=freqLabels{1};
                else
                    EPdata.freqNames=freqLabels;
                end
            end
        case 'stacked'
            ep_tictoc('ioStart');if EPtictoc.stop;return;end
            [theHeader, theRawData, theDelim] = ep_textScan(EPdata.fileName,2,0,1,0);
            ep_tictoc('ioFinish');if EPtictoc.stop;return;end
            if isempty(theRawData)
                msg{1}='No data were read.  Did you check the EP Toolkit preference settings for reading text files to make sure they were set correctly?';
                [msg]=ep_errorMsg(msg);
                return
            end
            if size(theRawData,2) < 9
                msg{1}=['There are too few columns of data (' num2str(size(theRawData,2)) ') for this to be the stacked format (at least nine).'];
                [msg]=ep_errorMsg(msg);
                return
            end
            theRawMatrix=cellfun(@str2double,theRawData(:,9:end));
            numObs=size(theRawMatrix,1);
            
            chanList=unique(theRawData(:,1),'stable');
            cellList=unique(theRawData(:,2),'stable');
            trialList=unique(str2double(theRawData(:,3)),'stable');
            subList=unique(theRawData(:,4),'stable');
            facList=unique(theRawData(:,5),'stable');
            freqList=unique(theRawData(:,6),'stable');
            relList=unique(theRawData(:,7),'stable');
            sessList=unique(theRawData(:,8),'stable');
            waveIndex=zeros(numObs,2);
            sessIndex=zeros(numObs,2);
            for iObs=1:numObs
                waveIndex(iObs,1)=find(strcmp(theRawData{iObs,2},cellList));
                waveIndex(iObs,2)=find(trialList==str2double(theRawData{iObs,3}));
                sessIndex(iObs,1)=find(strcmp(theRawData{iObs,4},subList));
                sessIndex(iObs,2)=find(strcmp(theRawData{iObs,8},sessList));
            end
            waveList=unique(waveIndex,'rows','stable');
            subSessList=unique(sessIndex,'rows','stable');
                        
            numChans=length(chanList);
            numPoints=size(theRawMatrix,2);
            numCells=length(cellList);
            numTrials=size(waveList,1);
            numSubs=size(subSessList,1);
            numFacs=length(facList);
            numFreqs=length(freqList);
            numRels=length(relList);
            numSess=length(sessList);
            
            theData=nan(numChans,numPoints,numTrials,numSubs,numFacs,numFreqs,numRels);
            multCount=0;
            for iRow=1:size(theRawMatrix,1)
                theChan=find(strcmp(theRawData(iRow,1),chanList));
                theWave=find(ismember(waveList,waveIndex(iRow,:),'rows'));
                theSub=find(ismember(subSessList,sessIndex(iRow,:),'rows'));
                theFac=find(strcmp(theRawData(iRow,5),facList));
                theFreq=find(strcmp(theRawData(iRow,6),freqList));
                theRel=find(strcmp(theRawData(iRow,7),relList));
                if ~all(isnan(theData(theChan,:,theWave,theSub,theFac,theFreq,theRel)))
                    multCount=multCount+1;
                end
                theData(theChan,:,theWave,theSub,theFac,theFreq,theRel)=theRawMatrix(iRow,:);
            end
            ep_tictoc('ioFinish');if EPtictoc.stop;return;end
            theData=reshape(theData,numChans,numPoints,[]);
            if multCount>0
                disp(['Warning: there are ' num2str(multCount) ' rows of data with the same labels.  Only the last of such sets of duplicates will be used.']);
            end
            
            EPdata.chanNames=chanList;
            EPdata.timeNames=cellfun(@str2double,theHeader{1}(9:end));
            if isscalar(EPdata.timeNames)
                EPdata.timeNames=[];
            end
            for iWave=1:numTrials
                EPdata.cellNames{iWave,1}=cellList{waveList(iWave,1)};
                EPdata.trialNames(iWave,1)=trialList(waveList(iWave,2));
            end
            if all(EPdata.trialNames==EPdata.trialNames(1))
                EPdata.trialNames=[];
            end
            if isscalar(facList)
                EPdata.facNames=cell(0);
            else
                EPdata.facNames=facList;
            end
            if isscalar(freqList)
                EPdata.freqNames=[];
            else
                EPdata.freqNames=str2double(freqList);
            end
            if isscalar(relList)
                EPdata.relNames=cell(0);
            else
                EPdata.relNames=relList;
            end
            
            if (isscalar(sessList)) && strcmp(sessList{1},'-')
                EPdata.sessNums=[];
                EPdata.sessNames=cell(0);
                EPdata.subNames=subList;
            else
                for iSub=1:numSubs
                    theSub=find(ismember(subSessList,sessIndex(iRow,:),'rows'));
                    EPdata.subNames{iSub}=subList{subSessList(theSub,1)};
                    EPdata.sessNums(iSub)=subSessList(theSub,2);
                    EPdata.sessNames=sessList;
                end
            end
            
            hdr.nChans=numChans;
            hdr.nSamples=numPoints;
            if ~isempty(EPdata.timeNames)
                sampLength=mean(diff(EPdata.timeNames));
                hdr.nSamplesPre=-EPdata.timeNames(1)/sampLength;
                hdr.Fs=1000/sampLength;
            else
                hdr.nSamplesPre=0;
                hdr.Fs=[];
            end
            hdr.nTrials=numObs/numChans;
            hdr.labels=chanList;
        case 'ns_mat'
            ep_tictoc('ioStart');if EPtictoc.stop;return;end
            [theNSFData,theEvents]=ep_readNSF(EPdata.fileName);
            ep_tictoc('ioFinish');if EPtictoc.stop;return;end
            if isempty(theNSFData)
                msg{1}=['The attempt to load in the file ' EPdata.fileName 'resulted in no data'];
                [msg]=ep_errorMsg(msg);
                return
            end
            
            theData=zeros(size(theNSFData{1,1},1),size(theNSFData{1,1},2),size(theNSFData,2));
            switch EPdata.dataType
                case {'continuous','single_trial'}
                    EPdata.cellNames=theNSFData(2,:)';
                    for i=1:length(EPdata.cellNames)
                        theIndex=strfind(EPdata.cellNames{i},'_Segment');
                        EPdata.cellNames{i,1}=EPdata.cellNames{i,1}(1:theIndex(end)-1);
                        theData(:,:,i)=theNSFData{1,i};
                    end
                case {'average','grand_average','factors'}
                    EPdata.subNames=theNSFData(2,:)';
                    EPdata.cellNames=theNSFData(2,:)';
                    for i=1:length(EPdata.cellNames)
                        theIndex=strfind(EPdata.cellNames{i},'_Average');
                        EPdata.cellNames{i,1}=EPdata.cellNames{i,1}(1:theIndex(end)-1);
                        EPdata.subNames{i,1}=EPdata.subNames{i,1}(theIndex(end)+8:end);
                        theData(:,:,i)=theNSFData{1,i};
                    end
            end
            hdr.nChans=size(theData,1);
            hdr.nSamples=size(theData,2);
            hdr.nSamplesPre=0;
            hdr.nTrials=0;
            if ~isempty(EPdata.Fs)
                hdr.Fs=EPdata.Fs;
            else
                hdr.Fs=250;
                disp('Assuming default sample rate of 250Hz.');
            end
            
            for i = 1:hdr.nChans
                hdr.label{i}  = ['e' num2str(i)];
            end
            
            eventHdr=[];
            if (theEvents{2,end}/(1000/hdr.Fs)) < (size(theData,2)*size(theData,3))
                for i = 1:size(theEvents,2)
                    eventHdr(i).value=theEvents{1,i};
                    eventHdr(i).type=[];
                    eventHdr(i).sample=theEvents{2,i}/(1000/hdr.Fs); %convert latency to samples from milliseconds.
                end
            else
                disp('Some events fall outside the range of data points.  One possible explanation is that when NetStation is set to drop bad segments, it only drops');
                disp('them from the voltage data, not the event data.  If this was done, then it is no longer possible');
                disp('to determine which event data go with which voltage data.  In any case, the event data will therefore be ignored.');
            end
            
        case 'egi_mff'
            hdr = mff_fileio_read_header(EPdata.fileName);
            eventHdr = mff_fileio_read_event(EPdata.fileName);
            ep_tictoc('ioStart');if EPtictoc.stop;return;end
            theData = mff_fileio_read_data(EPdata.fileName, 'header', hdr);
            ep_tictoc('ioFinish');if EPtictoc.stop;return;end
%             if ~any(strcmp(EPdata.dataType,{'average','grand_average'}))
%                 [hdr2]=ft_read_header(EPdata.fileName,'dataFormat','egi_mff_v1', 'headerformat','egi_mff_v1'); %v2 and v3 code doesn't provide subject specs and electrode coordinates yet
%                 hdr.origv2=hdr2.orig;
%             end

            %         try
            %            EPdata.events = read_mff_event(EPdata.fileName, hdr);
            %              if strcmp(EPdata.fileFormat,'egi_mff_v2')
            %                  EPdata.fileFormat='egi_mff_v1';
            %                 [tempHdr]=ft_read_header(EPdata.fileName,'dataFormat',EPdata.fileFormat, 'headerformat','egi_mff_v1'); %EGI's v2 routines don't provide all the header information yet.
            %                 hdr.orig=tempHdr.orig;
            %              end
            
        case 'edf'
            %currently only set up for continuous files.  I'd be happy to support other data types but would need an example to work with.
            ep_tictoc('ioStart');if EPtictoc.stop;return;end
            [hdr]=ft_read_header(EPdata.fileName,'dataFormat',EPdata.fileFormat, 'headerformat',EPdata.fileFormat,'chanindx',[1]);
            ep_tictoc('ioFinish');if EPtictoc.stop;return;end
            chanList=setdiff([1:hdr.orig.NS],hdr.orig.annotation);
            chanList=setdiff(chanList,find(strcmp(cellstr(hdr.orig.Label), 'ESUTimestamp')));
            chanList=setdiff(chanList,find(strcmp(cellstr(hdr.orig.Label), 'SystemTimestamp')));
            chanList=setdiff(chanList,find(strcmp(cellstr(hdr.orig.Label), 'Battery')));
            if isempty(chanList)
                hdr=[];
                eventHdr=[];
            else
                hdr2=hdr;
                modeFs=mode(hdr.orig.SampleRate);
                if ~all(hdr.orig.SampleRate==hdr.orig.SampleRate(1))
                    disp('EDF file has channels with different sampling rates.  Interpolating all to that of the mode rate.');
                    disp('May result in edge artifacts at end of epoch, which will be zeroed out.');
                end
                modeSampLength=1000/modeFs;
                newTimes=[modeSampLength:modeSampLength:hdr.orig.NRec*hdr.orig.Dur*1000];
                theData=zeros(length(chanList),length(newTimes));
                hdr.Fs=modeFs;
                hdr.nChans=length(chanList);
                for iChan=1:length(chanList)
                    theChan=chanList(iChan);
                    sampLength=1000/hdr.orig.SampleRate(theChan);
                    chanHdr = ft_read_header(EPdata.fileName, 'chanindx', theChan);
                    chanData=ft_read_data(EPdata.fileName,'dataformat',EPdata.fileFormat, 'headerformat',EPdata.fileFormat,'chanindx',1,'header',chanHdr);
                    chanTimes=[sampLength:sampLength:hdr.orig.NRec*hdr.orig.Dur*1000];
                    [theData(iChan,:)]=interp1(chanTimes(1:length(chanData)),chanData,newTimes);
                    hdr.label{iChan,1}=chanHdr.label{1};
                    hdr.chantype{iChan,1}=chanHdr.chantype{1};
                    hdr.chanunit{iChan,1}=chanHdr.chanunit{1};
                end
                ep_tictoc('ioStart');if EPtictoc.stop;return;end
                [eventHdr]=ft_read_event(EPdata.fileName,'dataformat',EPdata.fileFormat, 'headerformat',EPdata.fileFormat,'eventformat',EPdata.fileFormat,'chanindx',1,'header',hdr2);
                ep_tictoc('ioFinish');if EPtictoc.stop;return;end
            end
        case 'fieldtrip'
            try
                ep_tictoc('ioStart');
                tempVar=load(EPdata.fileName);
                ep_tictoc('ioFinish');
            catch ME
                msg{1}=['The attempt to load in the file ' EPdata.fileName ' resulted in the error:' ME.identifier];
                msg{2}=ME.message;
                [msg]=ep_errorMsg(msg);
                return
            end

            if isscalar(fieldnames(tempVar))
                fieldList=fieldnames(tempVar);
                evalc(['tempVar=tempVar.' fieldList{1} ';']);
            end

            if isfield(tempVar,'trial') && ~isfield(tempVar,'dat')
                tempVar.dat.trial=tempVar.trial;
                tempVar=rmfield(tempVar,'trial');
            end
            
            if (~isfield(tempVar,'hdr') || ~isfield(tempVar,'dat')) && (~isfield(tempVar,'freq'))
                msg{1}=['The file ' EPdata.fileName ' did not contain valid FieldTrip data in it.  EP Toolkit can only import a dat structure or a freq structure, as well as optionally an eventHdr structure.'];
                [msg]=ep_errorMsg(msg);
                return
            end
            
            if isfield(tempVar,'dat')
                dat=tempVar.dat;
                if isfield(tempVar.dat,'hdr')
                    hdr=tempVar.dat.hdr;
                else
                    hdr=[];
                end
            else
                dat=[];
            end
            if isfield(tempVar,'eventHdr')
                eventHdr=tempVar.eventHdr;
            else
                eventHdr=[];
            end
            if isfield(tempVar,'freq')
                freq=tempVar.freq;
            else
                freq=[];
            end
            if ~isempty(freq)
                if isfield(freq,'powspctrm')
                    msg{1}=['The file ' EPdata.fileName ' has frequency data with power units.  Only complex fourier values can be imported by EP Toolkit.'];
                    [msg]=ep_errorMsg(msg);
                    return
                elseif isfield(freq,'fourierspctrm')
                    switch freq.dimord
                        case 'rpt_chan_freq_time'
                            sevenDdata=ipermute(freq.fourierspctrm,[3 1 6 2 4 5]);
                        case 'rpt_chan_freq'
                            sevenDdata=ipermute(freq.fourierspctrm,[3 1 6 2 4 5]);
                        case 'chan_freq_time'
                            sevenDdata=ipermute(freq.fourierspctrm,[1 6 2 3 4 5]);
                        case 'chan_freq'
                            sevenDdata=ipermute(freq.fourierspctrm,[1 6 2 3 4 5]);
                        case 'rpttap_chan_freq_time'
                            sevenDdata=ipermute(freq.fourierspctrm,[3 1 6 2 4 5]);
                        case 'rpttap_chan_freq'
                            sevenDdata=ipermute(freq.fourierspctrm,[3 1 6 4 5 2]);
                        otherwise
                            msg{1}=['The file ' EPdata.fileName ' has a data structure with an unsupported dimord: ' freq.dimord '.'];
                            [msg]=ep_errorMsg(msg);
                            return
                    end
                else
                    msg{1}=['The file ' EPdata.fileName ' has a data structure with an unsupported data field.'];
                    [msg]=ep_errorMsg(msg);
                    return
                end
                EPdata.chanNames=freq.label(:);
                EPdata.freqNames=freq.freq(:);
                EPdata.timeNames=freq.time(:);
            else
                %                 if ~isfield(dat,'dimord')
                %                     msg{1}=['The file ' EPdata.fileName ' has a data structure with no dimord field so it cannot be interpreted.'];
                %                     [msg]=ep_errorMsg(msg);
                %                     return
                %                 end

                if isfield(dat,'avg')
                    if ~isfield(dat,'dimord')
                        if size(dat.avg) ==2
                            dat.dimord='chan_time';
                            disp('No dimord field so defaulting to chan_time.')
                        elseif size(dat.avg) ==3
                            dat.dimord='rpt_chan_time';
                            disp('No dimord field so defaulting to rpt_chan_time.')
                        end
                    end
                    switch dat.dimord
                        case 'chan_time'
                            sevenDdata=dat.avg;
                        case 'rpt_chan_time'
                            sevenDdata=shiftdim(dat.avg,1);
                        otherwise
                            msg{1}=['The file ' EPdata.fileName ' has a data structure with an unsupported dimord: ' dat.dimord '.'];
                            [msg]=ep_errorMsg(msg);
                            return
                    end
                    if ~strcmp(EPdata.dataType,'average')
                        disp('The fieldtrip data structure indicates it is actually an average file.')
                        EPdata.dataType='average';
                    end
                elseif isfield(dat,'trial')
                    if ~isfield(dat,'dimord')
                        if iscell(dat.trial)
                            dat.dimord='{rpt}_chan_time';
                            disp('No dimord field so defaulting to {rpt}_chan_time.')
                        else
                            dat.dimord='rpt_chan_time';
                            disp('No dimord field so defaulting to rpt_chan_time.')
                        end
                    end
                    switch dat.dimord
                        case '{rpt}_chan_time'
                            if ~iscell(dat.trial)
                                msg{1}=['The file ' EPdata.fileName ' has a dimord that is inconsistent with the data: ' dat.dimord '.'];
                                [msg]=ep_errorMsg(msg);
                                return
                            end
                            sevenDdata=zeros(size(dat.trial{1},1),size(dat.trial{1},2),length(dat.trial));
                            for iTrial=1:length(dat.trial)
                                sevenDdata(:,:,iTrial)=dat.trial{iTrial};
                            end
                        case 'rpt_chan_time'
                            sevenDdata=shiftdim(dat.trial,1);
                        otherwise
                            msg{1}=['The file ' EPdata.fileName ' has a data structure with an unsupported dimord: ' dat.dimord '.'];
                            [msg]=ep_errorMsg(msg);
                            return
                    end
                    if ~strcmp(EPdata.dataType,'single_trial')
                        disp('The fieldtrip data structure indicates it is actually a single-trial file.')
                        EPdata.dataType='single_trial';
                    end
                else
                    msg{1}=['The file ' EPdata.fileName ' has a data structure with an unsupported data field.'];
                    [msg]=ep_errorMsg(msg);
                    return
                end
            end
            theData=sevenDdata;
            hdr=tempVar.hdr;
            if length(hdr.label) > size(theData,1)
                eegList=find(strcmpi('eeg',hdr.chantype));
                if length(eegList) == size(theData,1)
                    hdr.label=hdr.label(eegList);
                    hdr.nChans=length(eegList);
                    hdr.chantype=hdr.chantype(eegList);
                    hdr.chanunit=hdr.chanunit(eegList);
                else
                    msg{1}=['The file ' EPdata.fileName ' has a hdr structure with a different number of channels than the data.'];
                    [msg]=ep_errorMsg(msg);
                    return
                end
            end

        case 'matlab_mat'

            tempVar=load('-mat', EPdata.fileName);
            tempFields=fieldnames(tempVar);
            if length(tempFields) > 1
                msg{1}=['The .mat file contains more than one variable.'];
                [msg]=ep_errorMsg(msg);
                return
            end

            eval(['matlabMatrix=tempVar.' tempFields{1} ';']);
            if ~isnumeric(matlabMatrix)
                msg{1}=['The .mat file does not contain a numeric variable.'];
                [msg]=ep_errorMsg(msg);
                return
            end

            numMatlabDims=length(size(matlabMatrix));
            if numMatlabDims > 7
                msg{1}=['The .mat file contains a matrix with too many dimensions.'];
                [msg]=ep_errorMsg(msg);
                return
            end

            if isempty(matlabDims)
                theList={'Channels';'TimePoints';'Cells/Trials';'Subjects';'Factors';'Frequencies';'Relations'};
                numDims=length(theList);
                EPmain.read.matlabDims=[1:numDims];
                matlabPaneCoords=get(EPmain.handles.hMainWindow,'Position');
                matlabPaneCoords(1)=matlabPaneCoords(1)+EPmain.panesize(1);
                handles.read.matlab = figure('Name', 'Matlab Matrix Dimensions', 'NumberTitle', 'off', 'Position',matlabPaneCoords, 'MenuBar', 'none');
                colormap jet;

                uicontrol('Style','text','FontSize',EPmain.fontsize,...
                    'String','All seven dimensions must be set with no repetitions.','Position',[5 20*numDims+20 150 40]);

                for iDim=1:numDims
                    EPmain.handles.read.matlabDims(iDim) = uicontrol('Style','popupmenu','FontSize',EPmain.fontsize,...
                        'String',theList,...
                        'Value',EPmain.read.matlabDims(iDim),'Position',[0 20*numDims-(iDim*20) 200 20], 'Callback', ['EPmain.read.matlabDims(' num2str(iDim) ')=get(EPmain.handles.read.matlabDims(' num2str(iDim) '),''Value'');']);
                end

                handles.read.SMItableDone = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',FontSize,...
                    'Position', [ 0 0 60 20], 'Callback', 'uiresume(gcbf)');

                uiwait(handles.read.matlab);
                close(handles.read.matlab);

                if length(EPmain.read.matlabDims) ~= length(unique(EPmain.read.matlabDims))
                    msg{1}=['Duplicate dimensions specified.'];
                    [msg]=ep_errorMsg(msg);
                    return
                end
                matlabDims=EPmain.read.matlabDims;
                EPmain.read.matlabDims=[];
            end

            [newDims, matlabOrder]=sort(matlabDims);
            theData=ipermute(matlabMatrix,matlabOrder);

            numChans=size(theData,1);
            numPoints=size(theData,2);
            numCells=size(theData,3);
            numWaves=numCells;
            numSubs=size(theData,4);
            numFacs=size(theData,5);
            numFreqs=size(theData,6);
            numRels=size(theData,7);

            theData=reshape(theData,numChans,numPoints,[]);

            hdr.nChans=numChans;
            hdr.nSamples=numPoints;
            hdr.nSamplesPre=0;
            hdr.nTrials=numCells;
            hdr.Fs=textPrefs.sampleRate;
            if ~isfield(hdr,'label')
                for i = 1:hdr.nChans
                    hdr.label{i}  = ['e' num2str(i)];
                end
            end
            if ~isempty(freqLabels) && ~isempty(freqLabels{1})
                if iscell(freqLabels) && ~isempty(freqLabels{1})
                    EPdata.freqNames(1)=freqLabels{1};
                else
                    EPdata.freqNames=freqLabels;
                end
            end

        case 'Neuroelectrics_easy'
            ep_tictoc('ioStart');if EPtictoc.stop;return;end
            [theHeader, theRawData, theDelim] = ep_textScan(EPdata.fileName,1,0,1,0);
            ep_tictoc('ioFinish');if EPtictoc.stop;return;end
            if isempty(theRawData)
                msg{1}='No data were read.';
                [msg]=ep_errorMsg(msg);
                return
            end
            theData=cellfun(@str2double,theRawData)';
            
            hdr.nChans=size(theData,1)-2;
            hdr.nSamples=size(theData,2);
            hdr.nSamplesPre=0;
            hdr.nTrials=0;

            infoName=[EPdata.fileName(1:end-5) '.info'];
            if ~exist(infoName,'file')
                msg{1}='No corresponding .info file.';
                [msg]=ep_errorMsg(msg);
                return
            end
            fid=fopen(infoName,'r');
            if fid==-1
                msg{1}='Unable to open .info file.';
                [msg]=ep_errorMsg(msg);
                return
            end
            infoData=textscan(fid,'%s','Delimiter',char(9));
            fclose(fid);
            if ~isfield(hdr,'label')
                chanCounter=1;
            else
                chanCounter=[]; %channels already named.
            end

            for iRow=1:size(infoData{1},1)
                if strfind(infoData{1}{iRow},'EEG sampling rate: ')==1
                    hdr.Fs=str2double(infoData{1}{iRow}(20:end-15));
                end
                if ~isempty(chanCounter) && ~isempty(strfind(infoData{1}{iRow},'Channel ')) && (strfind(infoData{1}{iRow},'Channel ')==1)
                    [token,remain] = strtok(infoData{1}{iRow},':');
                    hdr.label{chanCounter,1}  = remain(3:end);
                    chanCounter=chanCounter+1;
                end
            end
            hdr.label=[hdr.label;'ACMx'; 'ACMy'; 'ACMz'];
            EPdata.chanTypes=[cellstr(repmat('EEG',size(theData,1)-5,1));'ACMx'; 'ACMy'; 'ACMz'];

            eventHdr=struct('type',{},'sample',{},'value',{},'duration',{},'keys',struct('code','','data','','datatype','','description',''));
            eventHdr(length(find(theData(end-1,:)))).type='';
            eventCounter=1;
            for iSamp=1:size(theData,2)
                if theData(end-1,iSamp)
                    eventHdr(eventCounter).type='trigger';
                    eventHdr(eventCounter).sample=iSamp;
                    eventHdr(eventCounter).value=theData(end-1,iSamp);
                    eventHdr(eventCounter).duration=0;
                    eventCounter=eventCounter+1;
                end
            end
            theData=theData(1:end-2,:); %drop events and clock info
            theData(1:end-3)=theData(1:end-3)/1000; %convert from nanovolts to microvolts.

        otherwise
            ep_tictoc('ioStart');if EPtictoc.stop;return;end
            try
                [eventHdr]=ft_read_event(EPdata.fileName,'dataformat',EPdata.fileFormat, 'headerformat',EPdata.fileFormat,'eventformat',EPdata.fileFormat);
                [hdr]=ft_read_header(EPdata.fileName,'dataFormat',EPdata.fileFormat, 'headerformat',EPdata.fileFormat);
                %                 end
                %                 if strcmp(EPdata.fileFormat,'egi_mff_v1')
                %                     COMchan=find(strcmp('COM',hdr.label));
                %                     if ~isempty(COMchan) && (hdr.nChans ~= length(hdr.label))
                %                         %fix bug present in some mff files, resulting in the COM channel being identified in the header as an EEG channel.
                %                         hdr.label(COMchan)=[];
                %                         hdr.chantype(COMchan)=[];
                %                         hdr.chanunit(COMchan)=[];
                %                     end
                %                 end
                
                %                 if strcmp(EPdata.fileFormat,'egi_mff_v2') || strcmp(EPdata.fileFormat,'egi_mff_v1')
                %                     for iChan=1:length(hdr.label) %standardize mff channel labels
                %                         if ~isempty(str2num(hdr.label{iChan}))
                %                             newName=sprintf('E%s', hdr.label{iChan});
                %                             if isempty(find(strcmp(newName,hdr.label)))
                %                                 hdr.label{iChan}=newName;
                %                             end
                %                         elseif any(strcmp(hdr.label{iChan},{'Cz','REF','VREF'}))
                %                             newName=sprintf('E%s', num2str(iChan));
                %                             if isempty(find(strcmp(newName,hdr.label)))
                %                                 hdr.label{iChan}=newName;
                %                             end
                %                         end
                %                     end
                %                 end
                
                %                 if strcmp(EPdata.fileFormat,'egi_mff_v1') && strcmp(EPdata.dataType,'continuous') %workaround for continuous mff files with recording stops not being read
                %                     if isfield(hdr.orig,'epochdef')
                %                         theData=zeros(hdr.nChans,hdr.nSamples);
                %                         if isempty(eventHdr)
                %                             eventHdr=struct('type',{},'sample',{},'value',{},'duration',{},'keys',struct('code','','data','','datatype','','description',''));
                %                         end
                %                         if (size(hdr.orig.epochdef,1) > 2) && ~any(diff(hdr.orig.epochdef(:,2)-hdr.orig.epochdef(:,1)))
                %                             disp('Looks like this might actually be single trial data.');
                %                             [theData]=ft_read_data(EPdata.fileName,'dataformat',EPdata.fileFormat, 'headerformat',EPdata.fileFormat);
                %                             EPdata.dataType='single_trial';
                %                         else
                %                             for i=1:size(hdr.orig.epochdef,1)
                %                                 [tempData]=ft_read_data(EPdata.fileName,'dataformat',EPdata.fileFormat, 'headerformat',EPdata.fileFormat,'begsample',hdr.orig.epochdef(i,1),'endsample',hdr.orig.epochdef(i,2));
                %                                 theData(:,hdr.orig.epochdef(i,1):hdr.orig.epochdef(i,2))=tempData;
                %                                 if i > 1
                %                                     startEvents=find([eventHdr.sample] >= hdr.orig.epochdef(i,1));
                %                                     if ~isempty(startEvents)
                %                                         newEvent=startEvents(1);
                %                                         eventHdr=[eventHdr(1:newEvent-1) eventHdr(1) eventHdr(newEvent:end)];
                %                                     else
                %                                         newEvent=length(eventHdr)+1;
                %                                         eventHdr(end+1)=eventHdr(1);
                %                                     end
                %                                     eventHdr(newEvent).type='boundary';
                %                                     eventHdr(newEvent).sample=hdr.orig.epochdef(i,1);
                %                                     eventHdr(newEvent).duration=hdr.orig.epochdef(i,3); %duration field holds length of recording pause in samples
                %                                     eventHdr(newEvent).value='boundary';
                %                                     eventHdr(newEvent).keys=struct('key',struct('code','','data','','datatype','','description',''));
                %                                 end
                %                             end
                %                         end
                %                     else
                %                         [theData]=ft_read_data(EPdata.fileName,'dataformat',EPdata.fileFormat, 'headerformat',EPdata.fileFormat);
                %                     end
                %                 else
                [theData]=double(ft_read_data(EPdata.fileName,'dataformat',EPdata.fileFormat, 'headerformat',EPdata.fileFormat));
            catch ME
                msg{1}=['No data were read.  The error message was:' ME.identifier];
                msg{2}=ME.message;
                if strcmp(EPdata.fileFormat,'egi_egis')
                    msg{3}='Is it possible this is actually an average file? Perhaps you chose the wrong option in NetStation?';
                end
                if strcmp(EPdata.fileFormat,'egi_egia')
                    msg{3}='Is it possible this is actually a single-trial file? Perhaps you chose the wrong option in NetStation?';
                end
                if strcmp(EPdata.fileFormat,'eeglab_set')
                    msg{3}='EEGlab allows you to manually change the name of the .set file but not the .fdt file that goes with it.  Is that what happened here?';
                end
                [msg]=ep_errorMsg(msg);
                ep_tictoc('ioFinish');
                return
            end
            ep_tictoc('ioFinish');if EPtictoc.stop;return;end
            
            if strcmp(EPdata.fileFormat,'brainvision_eeg')
                segEvents=find(strcmp('New Segment',{eventHdr.type}));
                if ~isempty(segEvents)
                    boundaryCounter=0;
                    for iEvent=1:length(segEvents)
                        theEvent=segEvents(iEvent);
                        if eventHdr(theEvent).sample > 1
                            boundaryCounter=boundaryCounter+1;
                            eventHdr(theEvent).type = 'boundary';
                            eventHdr(theEvent).duration=0;
                            eventHdr(theEvent).value='boundary';
                        end
                    end
                    if boundaryCounter > 0
                        disp(['There were ' num2str(boundaryCounter) ' recording discontinuities.']);
                    end
                end
            end
    end

    if ~any(theData)
        msg{1}='Error: The data contained nothing but zeroes.';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    if any(any(any(isinf(theData))))
        msg{1}='Error: The data contained infinite numbers.';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    numChans=size(theData,1);
    numPoints=size(theData,2);
    
    if ~isempty(eventHdr)
        if isfield(eventHdr,'sample')
            if any([eventHdr.sample] < 0)
                if ~strcmp(EPdata.fileFormat,'egi_mff_v2') || any([eventHdr.sample] < -2)
                    disp('Warning: The events data contained negative samples.  These events will be dropped.');
                end
                if strcmp(EPdata.fileFormat,'egi_mff_v2') && any([eventHdr.sample] == -1)
                    disp(['There were ' num2str(length(find([eventHdr.sample] == -1))) ' events falling between epochs, which will be dropped.']);
                end
                if strcmp(EPdata.fileFormat,'egi_mff_v2') && any([eventHdr.sample] == -2)
                    disp(['There were ' num2str(length(find([eventHdr.sample] == -2))) ' samples falling after the last epoch, which will be dropped.']);
                end
                eventHdr=eventHdr(find([eventHdr.sample] > 0));
            end
        else
            eventHdr.sample=[];
        end

        if ~isfield(eventHdr,'value')
            eventHdr.value=[];
        end

        if ~isfield(eventHdr,'type')
            eventHdr.type=[];
        end

        if ~isfield(eventHdr,'duration')
            eventHdr.duration=[];
        end
    else
        eventHdr= struct('value','','sample',[],'type','','duration',[]);
        eventHdr(1)=[];
    end

    eventValues=cell(0);
    for event=1:length(eventHdr)
        if ~isempty(eventHdr(event).value)
            eventValues{event}=num2str(eventHdr(event).value);
        else
            eventValues{event}=num2str(eventHdr(event).type);
        end
    end
    
    if isempty(theData) && strcmp(EPdata.fileFormat,'egi_egia') %if it is actually  egi_egis, may return empty data
        EPdata.fileFormat='egi_egis';
        EPdata.dataType='single_trial';
        [theData]=ft_read_data(EPdata.fileName,'dataformat',EPdata.fileFormat, 'headerformat',EPdata.fileFormat);
    end
    
    if isempty(theData)
        msg{1}='No data were read.  Perhaps file format was misspecified?';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    % commented out.  Unfortunately, NetStation does not do this part of the EGIS header correctly.
    %     if (strcmp(EPdata.fileFormat,'egi_egia') && hdr.orig.fhdr(2) == 3)
    %         disp('The file indicates that this is actually an EGIS session file.');
    %         EPdata.fileFormat='egi_egis';
    %         EPdata.dataType = 'single_trial';
    %         numSubs=1;
    %         [theData]=ft_read_data(EPdata.fileName,'dataformat',EPdata.fileFormat, 'headerformat',EPdata.fileFormat);
    %     end
    
    if (strcmp(EPdata.fileFormat,'egi_egis') && hdr.orig.fhdr(2) == -1)
        disp('The file indicates that this is actually an EGIS average file.');
        EPdata.fileFormat='egi_egia';
        if strcmp(EPdata.dataType,'single_trial')
            EPdata.dataType = 'average';
        end
    end
    
    %     if strcmp(EPdata.fileFormat,'egi_sbin') && any(strcmp(EPdata.dataType,{'average','grand_average'})) && (hdr.orig.header_array(14) == 0)
    %         disp('The file indicates that this is actually a continuous file.');
    %         EPdata.dataType = 'continuous';
    %     end
    %
    if strcmp(EPdata.fileFormat,'biosemi_bdf') && ~isempty(strcmp('Status',hdr.label))
        statusChan=find(strcmp('Status',hdr.label)); %remove status chan from bdf files
        hdr.label(statusChan)=[];
        theData(statusChan,:)=[];
        hdr.chantype(statusChan)=[];
        hdr.chanunit(statusChan)=[];
        hdr.nChans=hdr.nChans-1;
        numChans=numChans-1;
    end
    
    if strcmp(EPdata.fileFormat,'neuromag_fif')
        if ~any(strcmp(EPdata.dataType,{'average','grand_average'})) && hdr.orig.isaverage
            disp('The file indicates that this is actually an average file.');
            EPdata.dataType = 'average';
        end
        if ~strcmp(EPdata.dataType,'continuous') && hdr.orig.iscontinuous && ~hdr.orig.isepoched
            disp('The file indicates that this is actually a continuous file.');
            EPdata.dataType = 'continuous';
        end
        if ~strcmp(EPdata.dataType,'single_trial') && hdr.orig.isepoched && ~hdr.orig.isaverage
            disp('The file indicates that this is actually a single-trial file.');
            EPdata.dataType = 'single_trial';
        end
        
        %From what I understand, .fif channel types are specified by users rather than being standardized, so there is no use to using them.
%         if isfield(hdr,'chantype')
%             EPdata.chanTypes=cellstr(repmat('EEG',size(theData,1),1));
%             MGMchans=find(strcmp('megmag',hdr.chantype));
%             if ~isempty(MGMchans)
%                 [EPdata.chanTypes{MGMchans}]=deal('MGM');
%             end
%             MGAchans=find(strcmp('megaxial',hdr.chantype));
%             if ~isempty(MGAchans)
%                 [EPdata.chanTypes{MGAchans}]=deal('MGA');
%             end
%             MGPchans=find(strcmp('megplanar',hdr.chantype));
%             if ~isempty(MGPchans)
%                 [EPdata.chanTypes{MGPchans}]=deal('MGP');
%             end
%             EEGchans=find(strcmp('eeg',hdr.chantype));
%             if ~isempty(EEGchans)
%                 [EPdata.chanTypes{EEGchans}]=deal('EEG');
%             end
%             EOGchans=find(strcmp('eog',hdr.chantype));
%             if ~isempty(EOGchans)
%                 [EPdata.chanTypes{EOGchans}]=deal('REG'); %assumed to be bipolar array
%             end
%             ATchans=find(strcmp('analog trigger',hdr.chantype));
%             if ~isempty(ATchans)
%                 [EPdata.chanTypes{ATchans}]=deal('BAD');
%             end
%             DTchans=find(strcmp('digital trigger',hdr.chantype));
%             if ~isempty(DTchans)
%                 [EPdata.chanTypes{DTchans}]=deal('BAD');
%             end
%             OTchans=find(strcmp('other trigger',hdr.chantype));
%             if ~isempty(OTchans)
%                 [EPdata.chanTypes{OTchans}]=deal('BAD');
%             end
%         end
        if isfield(hdr.orig,'projs')
            if isfield(hdr.orig.projs,'desc') && ~isempty(hdr.orig.projs)
                disp('I am not exactly sure how projections should be implemented so am ignoring them other than average reference.');
                if any(strcmp('Average EEG reference',{hdr.orig.projs(:).desc}))
                    disp('Average reference projection being applied to the data.')
                    EPdata.reference.type='AVG';
                    for iTrial=1:size(theData,3)
                        for iPoint=1:size(theData,2)
                            theData(:,iPoint,iTrial)=theData(:,iPoint,iTrial)-mean(theData(:,iPoint,iTrial));
                        end
                    end
                end
            end
        end
    end
    
    if strcmp(EPdata.fileFormat,'brainvision_eeg') && strcmp(EPdata.dataType,'continuous') && BVheader
        %look for subject header using EP convention
        ep_tictoc;if EPtictoc.stop;return;end
        iEvent=1;
        [badFlag,badBitFlag,headerFlag,hdrStart,hdrEnd,subjectSpecNames2,subjectSpecs2,trialSpecNames2]=readBVheader(eventHdr,iEvent,EPdata.subjectSpecNames,EPdata.subjectSpecs);

        if badBitFlag
            disp('Attempting to fix bad bit setting problem.')
            eventHdr2=eventHdr(1);
            eventHdr2(1)=[];
            iFix=1;
            while iFix <= (length(eventHdr)-1)
                if ~isempty(eventHdr(iFix).value)
                    if strcmp(eventHdr(iFix).value(1),'R') && strcmp(eventHdr(iFix+1).value(1),'S') && (eventHdr(iFix).sample==eventHdr(iFix+1).sample)
                        eventHdr2(end+1)=eventHdr(iFix);
                        eventHdr2(end).value=['S' pad(num2str(str2double(eventHdr(iFix+1).value(2:4)) + bitshift(str2double(eventHdr(iFix).value(2:4)),4)),3,'left')];
                        eventHdr2(end).type=eventHdr(iFix+1).type;
                        iFix=iFix+1;
                    elseif strcmp(eventHdr(iFix).value(1),'R')
                        eventHdr2(end+1).value=['S' pad(num2str(bitshift(str2double(eventHdr(iFix).value(2:4)),4)),3,'left')];
                        eventHdr2(end).type='Stimulus';
                    else
                        eventHdr2(end+1)=eventHdr(iFix);
                    end
                end
                iFix=iFix+1;
            end
            if iFix < length(eventHdr)
                eventHdr2(end+1)=eventHdr(end);
            end
            eventHdr=eventHdr2;
            [badFlag,badBitFlag,headerFlag,hdrStart,hdrEnd,subjectSpecNames2,subjectSpecs2,trialSpecNames2]=readBVheader(eventHdr,iEvent,EPdata.subjectSpecNames,EPdata.subjectSpecs);
        end
        partialHeader=0;
        if badFlag || ~headerFlag
            if BVheader==2 %ask user for file with replacement EP header.
                iEvent=1;
                while iEvent < length(eventHdr)-7 %minimum size of a valid subject header is eight
                    if strcmp(eventHdr(iEvent).value,'S114') && strcmp(eventHdr(iEvent+1).value,'S100') && strcmp(eventHdr(iEvent+2).value,'S104')
                        disp('Partial header found.  Please select file with intact EP header that can be used to replace it.');
                        partialHeader=1;
                        break
                    end
                    iEvent=iEvent+1;
                end
                if ~partialHeader
                    disp('No header found.  Please select file with intact EP header that can be used to replace it.');
                end
                [fixFileName, fixFilePathname] = uigetfile('*.*','Open Replacement File:','MultiSelect','off');
                if fixFileName ~= 0
                    [eventHdr2]=ft_read_event([fixFilePathname filesep fixFileName],'dataformat',EPdata.fileFormat, 'headerformat',EPdata.fileFormat,'eventformat',EPdata.fileFormat);
                    if ~isempty(eventHdr2)
                        [badFlag2,badBitFlag2,headerFlag2,hdrStart2,hdrEnd2,subjectSpecNames2,subjectSpecs2,trialSpecNames2]=readBVheader(eventHdr2,1,EPdata.subjectSpecNames,EPdata.subjectSpecs);
                        if ~headerFlag2
                            badFlag2=1;
                        end
                        if badFlag2
                            disp('***Replacement file also had problems with its EP header.  Will proceed without header information.***')
                        else
                            headerFlag=1;
                            badFlag=0;
                            if partialHeader
                                hdrEnd=iEvent+2;
                            else
                                hdrEnd=0;
                            end
                        end
                    end
                end
            end
        end
        
        if ~headerFlag
            disp(['Note: No EP Toolkit header detected in ' EPdata.fileName '.']);
            badFlag=1;
        end
        if ~badFlag
            EPdata.subjectSpecNames=subjectSpecNames2;
            EPdata.subjectSpecs=subjectSpecs2;
            EPdata.trialSpecNames=trialSpecNames2;
            if hdrEnd
                if hdrStart
                    eventHdr(hdrStart:hdrEnd)=[];
                else
                    eventHdr(1:hdrEnd)=[];
                end
            end
            iEvent=1;
            newEventHdr=[];
            TRSPcount=0;
            numKeys=length(EPdata.trialSpecNames);
            while iEvent <= length(eventHdr)
                ep_tictoc;if EPtictoc.stop;return;end
                if iEvent <= (length(eventHdr)+2)
                    if strcmp(eventHdr(iEvent).value,'S104') && strcmp(eventHdr(iEvent+1).value,'S100') && strcmp(eventHdr(iEvent+2).value,'S114')
                        disp(['***Warning: there is more than one header present.  You may want to trim this file before segmenting.***']);
                    end
                end
                if strcmp(eventHdr(iEvent).value,'S255')
                    TRSPcount=TRSPcount+1;
                    hdrStart=iEvent;
                    iEvent=iEvent+1;
                    if (iEvent > length(eventHdr)) || (length(EPdata.trialSpecNames) ~= str2double(eventHdr(iEvent).value(2:end)))
                        disp(['TRSP error.  Aborting effort to read TRSP #' num2str(TRSPcount) '.']);
                        badFlag=1;
                        continue
                    end
                    iEvent=iEvent+1; %skip over field indicating number of keys
                    theTRSP=[];
                    theTRSP.type='TRSP';
                    theTRSP.sample=eventHdr(iEvent).sample;
                    theTRSP.value='TRSP';
                    theTRSP.duration=0;
                    for iKey=1:numKeys
                        theTRSP.keys(iKey)=struct('code','','data','','datatype','','description','');
                        theTRSP.keys(iKey).code=EPdata.trialSpecNames{iKey};
                        switch eventHdr(iEvent).value
                            case {'S  1','S  4'}
                                theTRSP.keys(iKey).datatype='short';
                            case {'S  2','S  5','S  6','S  7'}
                                theTRSP.keys(iKey).datatype='long';
                            case {'S  3','S  8'}
                                theTRSP.keys(iKey).datatype='text';
                            otherwise
                                disp(['TRSP anomaly in TRSP #' num2str(TRSPcount) '.']);
                                badFlag=1;
                        end
                        [iEvent,fieldOut]=readBVheaderField(eventHdr,iEvent);
                        if isnan(fieldOut)
                            if strcmp(theTRSP.keys(1).data,'257') && (iKey==numKeys) %error due to several 1's in a row in the trigger line not being distinguished by PyCorder
                                disp(['Attempting fix in TRSP #' num2str(TRSPcount) '.']);
                                badFlag=0;
                                fieldOut=theTRSP.keys(numKeys-1).data;
                                for keyFix=numKeys:-1:3
                                    theTRSP.keys(keyFix).data=theTRSP.keys(keyFix-1).data;
                                    theTRSP.keys(keyFix).datatype=theTRSP.keys(keyFix-1).datatype;
                                end
                                theTRSP.keys(2).data=1;
                                theTRSP.keys(2).datatype='short';
                            else
                                disp(['TRSP error.  Aborting effort to read TRSP #' num2str(TRSPcount) '.']);
                                badFlag=1;
                                continue
                            end
                        end
                        theTRSP.keys(iKey).data=fieldOut;
                    end
                    if ~strcmp(eventHdr(iEvent).value,'S255')
                        disp(['TRSP error.  Aborting effort to read TRSP #' num2str(TRSPcount) '.']);
                        badFlag=1;
                    end
                    if badFlag
                        %                         newEventHdr(end+1:end+(iEvent-hdrStart+1))=eventHdr(hdrStart:iEvent); %add back into new eventHdr
                        badFlag=0;
                    else
                        newEventHdr(end+1).type=theTRSP.type;
                        newEventHdr(end).sample=theTRSP.sample;
                        newEventHdr(end).value=theTRSP.value;
                        newEventHdr(end).duration=theTRSP.duration;
                        newEventHdr(end).keys=theTRSP.keys;
                    end
                    iEvent=iEvent+1;
                else
                    %not part of TRSP so include in new eventHdr
                    newEventHdr(end+1).type=eventHdr(iEvent).type;
                    newEventHdr(end).sample=eventHdr(iEvent).sample;
                    newEventHdr(end).value=eventHdr(iEvent).value;
                    newEventHdr(end).duration=eventHdr(iEvent).duration;
                    if isfield(eventHdr,'keys')
                        newEventHdr(end).keys=eventHdr(iEvent).keys;
                    end
                    iEvent=iEvent+1;
                end
            end
            eventHdr=newEventHdr;
        end
        
        %add impedance info
        if isfield(hdr,'orig') && isfield(hdr.orig,'impedances')
            if isfield(hdr.orig.impedances,'channels')
                EPdata.impedances.channels=hdr.orig.impedances.channels;
            end
            if isfield(hdr.orig.impedances,'ground')
                EPdata.impedances.ground=hdr.orig.impedances.ground;
            end
        end
    end
    
    if strcmp(EPdata.fileFormat,'egi_mff') && strcmp(EPdata.dataType,'continuous') && isfield(eventHdr,'mffkeysbackup')
        for iEvent=1:length(eventHdr)
            theKeys=[];
            if ~isempty(eventHdr(iEvent).mffkeysbackup)
                keyData=eventHdr(iEvent).mffkeysbackup;
                [a keyData]=strtok(keyData,char(39));
                keyData=keyData(2:end);
                while length(keyData)>3
                    if strcmp(keyData(1),char(39)) %since strtok ignores the delimiter if in the first position
                        theName='';
                        keyData=keyData(2:end);
                    else
                        [theName keyData]=strtok(keyData,[char(39) ',']);
                    end
                    keyData=keyData(5:end);
                    keyCount=0;
                    while ~strcmp(keyData(1),'}')
                        if strcmp(keyData(1),char(39)) %since strtok ignores the delimiter if in the first position
                            theValue='';
                        else
                            [theValue keyData]=strtok(keyData,char(39));
                        end
                        keyCount=keyCount+1;
                        eval(['theKeys(' num2str(keyCount) ').' theName '=''' theValue ''';']);
                        keyData=keyData(2:end);
                        if strcmp(keyData(1:2),[',' char(39)])
                            keyData=keyData(3:end);
                        end
                    end
                    keyData=keyData(4:end);
                end
                
            end
            if strcmp(eventHdr(iEvent).code,'SESS') %if there are multiple, then will end up using the last one.
                if isfield(eventHdr,'label')
                    EPdata.ename=eventHdr(iEvent).label;
                end
                for iKey=1:length(theKeys)
                    if isfield(theKeys,'code') && isfield(theKeys,'data')
                        EPdata.subjectSpecNames{iKey,1}=theKeys(iKey).code;
                        EPdata.subjectSpecs{1,iKey}=theKeys(iKey).data;
                        if strcmp(theKeys(iKey).code,'subj')
                            EPdata.subNames{1}=theKeys(iKey).data;
                        end
                    end
                end
            else
                eventHdr(iEvent).keys=theKeys;
            end
        end
        fieldList=fieldnames(eventHdr);
        for iField=1:length(fieldList)
            if (length(fieldList{iField})>5) && strcmp(fieldList{iField}(1:6),'mffkey')
                eventHdr=rmfield(eventHdr,fieldList{iField});
            end
        end
    end
    
    %Determine number of subjects and factors
    ep_tictoc;if EPtictoc.stop;return;end
    if strcmp(EPdata.fileFormat,'egi_egia') && isempty(numSubs)
        numSubs = hdr.orig.chdr(1,2); %since NetStation does not properly set the fhdr(11) field, use the number of subjects from the chdr instead
    end
    
    if strcmp(EPdata.fileFormat,'egi_sbin') && isempty(numSubs)
        if ~any(strcmp(EPdata.dataType,{'single_trial','continuous'}))
            numSubs = size(theData,3)/hdr.orig.header_array(14);
            if floor(numSubs) ~= numSubs
                disp(['The number of segments (' num2str(size(theData,3)) ') is not an even multiple of the number of categories (' num2str(hdr.orig.header_array(14)) '), which suggests this is really a single trial file.']);
                EPdata.dataType='single_trial';
                numSubs=1;
            end
        else
            numSubs=1;
        end
    end
    
    if strcmp(EPdata.dataType,'factors') && ~isempty(numSubs)
        numFacs=numSubs; %non-EP files will normally have only subjects OR factors.
        numSubs=1;
        if strcmp(EPdata.fileFormat,'egi_egia')
            if numFacs == hdr.orig.chdr(1,2)-1
                numFacs=numFacs+1; %it is assumed that factor files always have an additional summed factors waveform
            end
        end
    end
    
    if strcmp(EPdata.fileFormat,'egi_egia')
        if (numSubs ~= hdr.orig.chdr(1,2)) && (hdr.orig.chdr(1,2)*hdr.orig.fhdr(18)==size(theData,3))
            numSubs = hdr.orig.chdr(1,2);
            disp(['According to the file, the number of subjects is actually: ' num2str(hdr.orig.chdr(1,2))]);
        end
    end
    
    if isempty(numSubs)
        numSubs=1;
    end

    if strcmp(silent,'off')
        if strcmp(EPdata.fileFormat,'egi_egis')
            theMeanStd=mean(mean(std(reshape(theData,numChans,[]),'omitnan')));
            if theMeanStd > 200
                disp(['Warning: The mean standard deviation of the data is excessively high ( ' num2str(theMeanStd) ' ) and so it is possible that NetStation screwed up the EGIS data scaling again, although perhaps it is due to not filtering.']);
            end
        end
    end
    
    if strcmp(EPdata.dataType,'continuous') && (size(theData,3) > 1)
        disp('A file that has multiple ''trials'' is actually either a ''single-trial'' or ''average'' file.  It will be assumed this is a ''single-trial'' file.');
        EPdata.dataType='single_trial';
    end
    
    if strcmp(EPdata.dataType,'continuous')
        if numSubs > 1
            disp(['For continuous data, the number of subjects must be one, rather than ' num2str(numSubs) '.']);
        end
        numSubs=1;
    end
    
    if strcmp(EPdata.dataType,'single_trial')
        if numSubs > 1
            disp(['For single trial data, the number of subjects must be one, rather than ' num2str(numSubs) '.']);
        end
        numSubs=1;
    end
    
    %Determine other general information.
    if isfield(hdr,'Fs') && ~isempty(hdr.Fs)
        EPdata.Fs=hdr.Fs;
    else
        EPdata.Fs=250;
        disp('Assuming default sample rate of 250Hz.');
    end
    
    if isempty(EPdata.chanNames) && isfield(hdr,'label')
        EPdata.chanNames=hdr.label;
        EPdata.chanNames=EPdata.chanNames(:);
    end
    
    if length(EPdata.chanNames) ~= length(unique(EPdata.chanNames))
        msg{1}='Error: Channel names must be unique.';
        [msg]=ep_errorMsg(msg);
        return
    end
    
    if isempty(EPdata.baseline)
        if isfield(hdr,'nSamplesPre')
            EPdata.baseline=hdr.nSamplesPre;
        end
        if isempty(EPdata.baseline)
            EPdata.baseline=0;
        end
    else
        if mod(EPdata.baseline,(1000/EPdata.Fs)) >0
            EPdata.baseline=EPdata.baseline-mod(EPdata.baseline,(1000/EPdata.Fs));
            disp(['The specified baseline was not evenly divided by the sample lengths of ' num2str(1000/EPdata.Fs) ' msec.']);
            disp(['The baseline has been changed to ' num2str(EPdata.baseline) ' msec.']);
        end
        EPdata.baseline=EPdata.baseline/(1000/EPdata.Fs); %change to samples
    end
    
    numSamples=size(theData,2);
    if isempty(EPdata.timeNames) && isfield(EPdata,'baseline')
        EPdata.timeNames(:,1)=(((0:numSamples-1)-EPdata.baseline)*(1000/EPdata.Fs));
    end
    
    %determine cellnames and so forth from what is available in the file format.
    
    %EGIS session file
    ep_tictoc;if EPtictoc.stop;return;end
    if strcmp(EPdata.fileFormat,'egi_egis')
        EPdata.trialSpecs=cell(hdr.orig.fhdr(18),(max(hdr.orig.chdr(:,5))/2));
        trialCounter=0;
        for theCell = 1:hdr.orig.fhdr(18)
            numSpecs=(hdr.orig.chdr(theCell,5)/2);
            for trial = 1:hdr.orig.chdr(theCell,2)
                trialCounter=trialCounter+1;
                EPdata.cellNames{trialCounter,1} = deblank(hdr.orig.cnames{theCell});
                EPdata.trialNames(trialCounter,1)=trial;
                for spec=1:numSpecs
                    EPdata.trialSpecs{trialCounter,spec} = hdr.orig.chdr(theCell,5+(trial-1)*numSpecs+spec);
                end
            end
        end
        numWaves=trialCounter;
        numCells=hdr.orig.fhdr(18);
        
        EPdata.subNames{1,1}=['sub' sprintf('%03.0f',hdr.orig.fhdr(11))];
        
        for i=1:size(EPdata.trialSpecs,2)
            EPdata.trialSpecNames{i,1}=['spec' sprintf('%03.0f',i)];
        end
        EPdata.trialSpecNames{1,1}='edit';
        
        EPdata.subjectSpecs=cell(1,13);
        EPdata.subjectSpecNames{1,1}='RunDateMo';
        EPdata.subjectSpecNames{2,1}='RunDateDay';
        EPdata.subjectSpecNames{3,1}='RunDateYr';
        EPdata.subjectSpecNames{4,1}='RunTimeHr';
        EPdata.subjectSpecNames{5,1}='RunTimeMin';
        EPdata.subjectSpecNames{6,1}='RunTimeSec';
        EPdata.subjectSpecNames{7,1}='SubjID';
        EPdata.subjectSpecNames{8,1}='Handed';
        EPdata.subjectSpecNames{9,1}='Sex';
        EPdata.subjectSpecNames{10,1}='Age';
        EPdata.subjectSpecNames{11,1}='ExperID';
        EPdata.subjectSpecNames{12,1}='Comment';
        EPdata.subjectSpecNames{13,1}='Text';
        
        for i=1:11
            EPdata.subjectSpecs{i}=num2str(hdr.orig.fhdr(i+4));
        end
        EPdata.subjectSpecs{12}=char(hdr.orig.fcom);
        EPdata.subjectSpecs{13}=char(hdr.orig.ftext);
        
        %EGIS average file
    elseif strcmp(EPdata.fileFormat,'egi_egia')
        numCells = hdr.orig.fhdr(18);
        numWaves=numCells; %EGIS average files are never session files
        
        cellHeaderInfoSize=size(hdr.orig.chdr,2)-numSubs*hdr.orig.chdr(1,5)/2; %normally 11 fields but in some Dien datasets there were less to save on limited header space
        
        EPdata.cellNames=deblank(hdr.orig.cnames);
        EPdata.cellNames=EPdata.cellNames(:);
        if strcmp(EPdata.dataType,'factors')
            for i=1:numFacs
                EPdata.facNames{i,1}=['fac' sprintf('%03.0f',i)];
                EPdata.facTypes{i,1}='SGL';
            end
            EPdata.facNames{numFacs,1}='summed factors'; %assume the last one is a summed factor since that is what the EP Toolkit generates by default
            EPdata.facTypes{numFacs,1}='CMB';
            EPdata.subNames{1,1}='summed subjects';
        else
            if ~any(hdr.orig.chdr(1,cellHeaderInfoSize+(1:numSubs)*hdr.orig.chdr(1,5)/2))
                for i=1:numSubs
                    EPdata.subNames{i,1}=['sub' sprintf('%03.0f',i)]; %correct for NetStation bug of not outputing subject numbers
                end
            else
                for i=1:numSubs
                    EPdata.subNames{i,1}=['sub' sprintf('%03.0f',hdr.orig.chdr(1,cellHeaderInfoSize+(i-1)*hdr.orig.chdr(1,5)/2))];
                end
            end
        end
        
        EPdata.subjectSpecs=cell(numSubs,cellHeaderInfoSize+2);
        numSpecs=(hdr.orig.chdr(1,5)/2); %assume subject information is the same for each cell
        for theSub = 1:numSubs
            for spec=1:cellHeaderInfoSize
                EPdata.subjectSpecs{theSub,spec} = num2str(hdr.orig.chdr(1,(theSub-1)*numSpecs+spec));
            end
            EPdata.subjectSpecs{theSub,cellHeaderInfoSize+1}=char(hdr.orig.fcom);
            EPdata.subjectSpecs{theSub,cellHeaderInfoSize+2}=char(hdr.orig.ftext);
        end
        
        if cellHeaderInfoSize ==11
            EPdata.subjectSpecNames{1,1}='RunDateMo';
            EPdata.subjectSpecNames{2,1}='RunDateDay';
            EPdata.subjectSpecNames{3,1}='RunDateYr';
            EPdata.subjectSpecNames{4,1}='RunTimeHr';
            EPdata.subjectSpecNames{5,1}='RunTimeMin';
            EPdata.subjectSpecNames{6,1}='RunTimeSec';
            EPdata.subjectSpecNames{7,1}='SubjID';
            EPdata.subjectSpecNames{8,1}='Handed';
            EPdata.subjectSpecNames{9,1}='Sex';
            EPdata.subjectSpecNames{10,1}='Age';
            EPdata.subjectSpecNames{11,1}='ExperID';
            EPdata.subjectSpecNames{12,1}='Comment';
            EPdata.subjectSpecNames{13,1}='Text';
        else
            for i=1:cellHeaderInfoSize
                EPdata.subjectSpecNames{i,1}=sprintf('spec%02.0f',i);
            end
            EPdata.subjectSpecNames{cellHeaderInfoSize+1,1}='Comment';
            EPdata.subjectSpecNames{cellHeaderInfoSize+2,1}='Text';
        end
        
        for theCell=1:numCells
            for theSub=1:numSubs
                numSpecs=(hdr.orig.chdr(theCell,5)/2);
                if strcmp(EPdata.dataType,'factors') || numSubs == 1
                    EPdata.subNum(theSub,theCell)=hdr.orig.chdr(theCell,5+(theSub-1)*numSpecs+1);
                    EPdata.avgNum(theSub,theCell)=0;
                    EPdata.covNum(theSub,theCell)=0;
                else
                    theNum=hdr.orig.chdr(theCell,5+(theSub-1)*numSpecs+1);
                    if theNum==-1
                        theNum=hdr.orig.chdr(theCell,5+(theSub-1)*numSpecs+2); %workaround for NetStation bug
                    end
                    EPdata.avgNum(theSub,theCell)=theNum;
                    EPdata.covNum(theSub,theCell)=theNum;
                    EPdata.subNum(theSub,theCell)=1;
                end
            end
        end
        
        %EGI Simple Binary
    elseif strcmp(EPdata.fileFormat,'egi_sbin')
        disp('Reminder: I do not recommend using simple binary format if it can be avoided.');
        numCells = hdr.orig.header_array(14);
        for i=1:numSubs
            EPdata.subNames{i,1}=['sub' sprintf('%03.0f',i)];
        end
        numWaves=size(theData,3)/numSubs;
        
        if any(strcmp(EPdata.dataType,{'average','grand_average'}))
            EPdata.subjectSpecs=cell(numSubs,7);
            for spec=1:7
                EPdata.subjectSpecs{numSubs,spec} = hdr.orig.header_array(spec+1);
            end
            EPdata.subjectSpecNames{1,1}='RunDateYr';
            EPdata.subjectSpecNames{2,1}='RunDateMo';
            EPdata.subjectSpecNames{3,1}='RunDateDay';
            EPdata.subjectSpecNames{4,1}='RunTimeHr';
            EPdata.subjectSpecNames{5,1}='RunTimeMin';
            EPdata.subjectSpecNames{6,1}='RunTimeSec';
            EPdata.subjectSpecNames{7,1}='RunTimeMsec';
        end
        
        droppedEvent=0;
        eventIndex=ones(length(eventValues),1);
        for i=1:length(eventValues)
            if isempty(eventValues{i})
                eventIndex(i)=0;
            elseif strcmp(eventValues{i},'CELL') %the information associated with these two events are lost when NS exports the file.
                eventIndex(i)=0;
                droppedEvent=1;
            end
        end
        
        eventValues=eventValues(find(eventIndex));
        eventHdr=eventHdr(find(eventIndex));
        
        if droppedEvent
            disp('Dropping CELL events since NetStation loses the associated information when it exports a simple binary file.');
        end
        
        subjectsGrouped=1;
        if strcmp(EPdata.dataType,'continuous')
            disp('Note that NetStation, at least as of version 4.1.2, does not import the event information associated with');
            disp('continuous simple binary files, so if you have problems with this, this is not the Toolkit''s fault.');
            numCells=1; %continuous data only has one cell
            EPdata.cellNames{1,1}='cell01';
        else
            if ~isempty(cellLabels)
                numCells=length(unique(cellLabels));
                if isscalar(cellLabels)
                    EPdata.cellNames=repmat(cellLabels,numWaves,1);
                else
                    EPdata.cellNames=cellLabels;
                end
            elseif length(eventHdr)==numWaves
                numCells=length(unique(eventValues));
                EPdata.cellNames=eventValues;
            elseif sum(strcmp('trial',{eventHdr.type})) ==numWaves*numSubs
                if (hdr.orig.header_array(14))==0 && (hdr.orig.header_array(15) > 1) %epoch-marked simple binary file format
                    numCells=1; %epoch-marked doesn't allow for separate cells
                    for i=1:numWaves
                        EPdata.cellNames{i,1}='cell01';
                    end
                    disp('Epoch-marked simple binary is not able to separate the cells so putting all the segments into the same cell.');
                    disp('You will need to use the Edit function to provide them with their correct cell names.');
                    disp('It will not be possible to use this format for combined subject average files.');
                else
                    EPdata.cellNames={eventHdr(strcmp('trial',{eventHdr.type})).value}; %the category names in the SegHdrs have been stored here.
                    numCells=length(unique(EPdata.cellNames));
                end
            else
                theTrial=ceil([eventHdr.sample]/numSamples);
                theOffset=mod([eventHdr.sample],numSamples);
                if (length(find(theOffset==1))==numWaves) && length(theTrial(unique(find(theOffset==1))))==numWaves
                    numCells=length(unique(eventValues(find(theOffset==1)))); %an event is at the offset every time
                    EPdata.cellNames=eventValues(find(theOffset==1));
                else
                    numCells=1; %give up trying to separate the cells and put them all into the same one.
                    for i=1:numWaves
                        EPdata.cellNames{i,1}='cell01';
                    end
                    disp('Can''t figure out the cell names so putting all the segments into the same cell.');
                end
            end
            
            uniqueCellNames=unique(EPdata.cellNames);
            if any(strcmp(EPdata.dataType,{'average','grand_average'}))
                subjectsGrouped=-1; %whether a subject is grouped such that all its cells are consecutive in the file
                for i=1:length(uniqueCellNames)
                    whichCells=find(strcmp(uniqueCellNames{i},EPdata.cellNames));
                    if length(whichCells) ~= numSubs
                        disp('The subjects do not all have one of each cell so it is not possible to decrypt which waveform goes with which subject.');
                        return
                    end
                    if length(whichCells) > 1
                        if all(diff(whichCells) == 1) %the waveforms are grouped by cell
                            if subjectsGrouped==-1
                                subjectsGrouped=0;
                            elseif subjectsGrouped==1
                                disp('The waveforms are not consistently organized in the file so it is not possible to decrypt which waveform goes with which subject.');
                                return
                            end
                        elseif ~any(diff(diff(whichCells)))
                            %the waveforms are grouped by subject
                            if subjectsGrouped==-1
                                subjectsGrouped=1;
                            elseif subjectsGrouped==0
                                disp('The waveforms are not consistently organized in the file so it is not possible to decrypt which waveform goes with which subject.');
                                return
                            end
                        else
                            disp('The waveforms are not consistently organized in the file so it is not possible to decrypt which waveform goes with which subject.');
                        end
                    end
                end
                if subjectsGrouped
                    EPdata.cellNames=EPdata.cellNames(1:length(uniqueCellNames));
                else
                    EPdata.cellNames=EPdata.cellNames(1:numSubs:length(EPdata.cellNames));
                end
            else
                if ~isempty(trialLabels)
                    EPdata.trialNames=trialLabels;
                else
                    trialCounter=zeros(length(uniqueCellNames),1);
                    for theCell=1:length(EPdata.cellNames)
                        newCell=strcmp(EPdata.cellNames(theCell),uniqueCellNames);
                        trialCounter(newCell)=trialCounter(newCell)+1;
                        EPdata.trialNames(theCell,1)=trialCounter(newCell);
                    end
                end
                EPdata.trialSpecs=cell(numWaves,0);
            end
        end
        EPdata.cellNames=EPdata.cellNames(:);
        EPdata.trialNames=EPdata.trialNames(:);
        
        %any other sort of file other than EGIS or Simple Binary
    else
        
        if isscalar(eventHdr)
            if strcmp(eventHdr.type,'average') && ~any(strcmp(EPdata.dataType,{'average','grand_average'}))
                EPdata.dataType='average';
                disp('The file indicates that it is actually an average file.');
            end
        end
        
        if isempty(EPdata.subNames)
            for iSub=1:numSubs
                EPdata.subNames{iSub,1}=['sub' sprintf('%03.0f',iSub)];
            end
        end
        numWaves=size(theData,3)/(numSubs*numFacs*numFreqs*numRels);
        
        if strcmp(EPdata.fileFormat,'ns_avg')
            if strcmp(EPdata.dataType,'grand_average')
                EPdata.subNum=hdr.orig.nsweeps;
            else
                EPdata.avgNum=hdr.orig.nsweeps;
                EPdata.covNum=hdr.orig.nsweeps;
            end
        end
        
        if strcmp(EPdata.fileFormat,'egi_mff_v2') && any(strcmp(EPdata.dataType,{'average','grand_average'}))
            if isempty(hdr.orig.epochSubjects)
                EPdata.subNames{1,1}='sub001';
                numSubs=1;
            else
                if length(unique(hdr.orig.multiSubj)) >1
                    numSubs=length(hdr.orig.epochSubjects);
                elseif length(unique(hdr.orig.epochFilenames)) >1
                    numSubs=length(unique(hdr.orig.epochFilenames));
                else
                    disp('Unable to determine number of subjects.');
                end
                for i=1:numSubs
                    if ~isempty(hdr.orig.epochSubjects{i})
                        EPdata.subNames{i,1}=hdr.orig.epochSubjects{i};
                    elseif ~isempty(hdr.orig.epochFilenames{i})
                        EPdata.subNames{i,1}=hdr.orig.epochFilenames{i};
                    else
                        EPdata.subNames{i,1}=['sub' sprintf('%03.0f',i)];
                    end
                end
            end
            numWaves=size(theData,3)/numSubs;
            EPdata.cellNames=hdr.orig.epochLabels;
            EPdata.cellNames=EPdata.cellNames(1:numSubs:length(EPdata.cellNames)); %averages are grouped by cell
            numCells=length(EPdata.cellNames);
        end
        
        if strcmp(EPdata.fileFormat,'egi_mff') && any(strcmp(EPdata.dataType,{'average','grand_average'}))
            %in egi_mff_v3,the cell name labels appear in events with zero duration.  Unfortunately, Arno is giving me only partial information about
            %how this information is meant to be extracted from his code's output so I'm having to make a best guess.
            %if there are any other events with zero duration, then this approach will not work.
            %There is also no subject information so I can only make assumptions about how to tell the subjects apart.
            
            disp('Warning: EGI''s mffmatlabio code provides no clear information on the identity of averaged segments so the EP Toolkit will have to make a best guess as to cell and subjects.');
            disp('Make sure to double check the waveforms to verify that they were correctly apportioned.  You will also need to restore the subject names manually.');
            disp('If you are not happy about the current state of affairs, please feel free to inform EGI.');
            
            EPdata.cellNames=eventValues(find([eventHdr.duration]==0))';
            if length(EPdata.cellNames) ~= numWaves
                EPdata.cellNames=cell(0);
                disp('EP Toolkit was unable to determine the cell names using mffmatlabio.');
                disp('If you contact me, I am willing to see if I can come up with a solution for you.');
                disp('Otherwise, you will have to seek further help from EGI support or Arno Delorme.')
            else
                uniqueCellNames=unique(EPdata.cellNames,'stable');
                if length(EPdata.cellNames) ~= length(uniqueCellNames)
                    [~, ~, cellCount]=unique(EPdata.cellNames);
                    occ = histc(cellCount, 1:numel(uniqueCellNames));
                    if any(diff(occ))
                        disp('EP Toolkit was unable to determine the number of subjects using mffmatlabio.');
                        disp('If you contact me, I am willing to see if I can come up with a solution for you.');
                        disp('Otherwise, you will have to seek further help from EGI support or Arno Delorme.')
                    else
                        numSubs=numWaves/length(uniqueCellNames);
                        numWaves=size(theData,3)/numSubs;
                        EPdata.cellNames=uniqueCellNames;
                        for i=1:numSubs
                            EPdata.subNames{i,1}=['sub' sprintf('%03.0f',i)];
                        end
                        subjectsGrouped=1; %averages are assumed to be grouped by subject
                    end
                    numCells=length(EPdata.cellNames);
                end
            end
        end
        
        %         if strcmp(EPdata.fileFormat,'egi_mff') && any(strcmp(EPdata.dataType,{'average','grand_average'}))
        %             labelList=unique({hdr.orig.epoch.eventtype});
        %             countList=zeros(length(labelList),1);
        %             for iLabel=1:length(labelList)
        %                 countList(iLabel)=length(find(strcmp(labelList{iLabel},{hdr.orig.epoch.eventtype})));
        %             end
        %             if ~any(diff(countList))
        %                 numSubs=countList(1);
        %                 numWaves=size(theData,3)/numSubs;
        %                 EPdata.cellNames={hdr.orig.epoch.eventtype};
        %                 EPdata.cellNames=EPdata.cellNames(1:numSubs:length(EPdata.cellNames)); %averages are grouped by cell
        %                 numCells=length(EPdata.cellNames);
        %                 for i=1:numSubs
        %                     EPdata.subNames{i}=['sub' sprintf('%03.0f',i)];
        %                 end
        %
        %             else
        %                 disp('Unable to determine number of subjects.');
        %             end
        %         end
        
        if any(strcmp(EPdata.fileFormat,{'egi_mff_v1';'egi_mff_v2'})) && any(strcmp(EPdata.dataType,{'continuous','single_trial'}))
            %not sure yet how this structure works so may need to be revised
            for iSpec=1:length(hdr.orig.xml.subject.fields)
                EPdata.subjectSpecNames{end+1,1}=hdr.orig.xml.subject.fields(iSpec).field.name;
                if isfield(hdr.orig.xml.subject.fields(iSpec).field.data,'data')
                    EPdata.subjectSpecs{1,end+1}=hdr.orig.xml.subject.fields(iSpec).field.data.data;
                else
                    EPdata.subjectSpecs{1,end+1}=[];
                end
            end
            %             if ~isempty(find(strcmp('SESS',eventValues))) && strcmp(EPdata.dataType,'continuous')
            %                 theSess=find(strcmp('SESS',eventValues));
            %                 theSess=theSess(end); %if there are multiple SESS events, assume that the earlier ones were false starts.
            %                 for iKey=1:length(eventHdr(theSess).orig.keys)
            %                     EPdata.subjectSpecNames{end+1}=eventHdr(theSess).orig.keys(iKey).key.keyCode;
            %                     EPdata.subjectSpecs{1,end+1}=eventHdr(theSess).orig.keys(iKey).key.data.data;
            %                 end
            %             end
        end
        
        ep_tictoc;if EPtictoc.stop;return;end
        if strcmp(EPdata.fileFormat,'eeglab_set')
            if isfield(hdr.orig.event,'setname')
                if ~strcmp(EPdata.dataType,'average')
                    EPdata.dataType='average';
                    disp('The header indicates this is actually an average file generated by Widmann''s pop_grandaverage function.')
                end
                numSubs=numWaves; %third dimension of data is subjects.  one cell per file.
                numWaves=1;
                numCells=1;
                if isempty(EPdata.cellNames)
                    EPdata.cellNames{1,1}=hdr.orig.event(1).setname;
                end
                EPdata.avgNum=shiftdim(EPdata.avgNum);
                EPdata.covNum=shiftdim(EPdata.covNum);
                EPdata.subNum=ones(numSubs,numCells);
                for i=1:numSubs
                    EPdata.subNames{i,1}=['sub' sprintf('%03.0f',i)];
                    EPdata.avgNum(i,1)=hdr.orig.epoch(i).eventtrials;
                    EPdata.covNum(i,1)=hdr.orig.epoch(i).eventtrials;
                end
            elseif strcmp(EPdata.dataType,'single_trial') && length(hdr.orig.event) == numWaves
                EPdata.recTime=cell2mat({hdr.orig.event(:).latency});
            end
            if strcmp(EPdata.dataType,'single_trial')
                if isfield(hdr.orig,'reject')
                    if isfield(hdr.orig.reject,'manualE')
                        EPdata.analysis.badChans(1,:,:)=hdr.orig.reject.manualE';
                    end
                    if isfield(hdr.orig.reject,'manual')
                        EPdata.analysis.badTrials(1,:)=hdr.orig.reject.manual';
                    end
                end
            end
        end
        
        if strcmp(EPdata.fileFormat,'eeglab_erp')
            if isempty(EPdata.subNames)
                EPdata.subNames{1}=hdr.orig.erpname; %erplab files have one subject/multiple cells whereas Widmann eeglab variant has multiple subjects/one cell
            end
            for i=1:numWaves
                EPdata.avgNum(1,i)=hdr.orig.ntrials.accepted(i);
                EPdata.covNum(1,i)=hdr.orig.ntrials.accepted(i);
            end
            if length(hdr.orig.EVENTLIST)>1 %grand average file
                EPdata.subNum=ones(1,numWaves)*length(hdr.orig.EVENTLIST);
                if ~strcmp(EPdata.dataType,'grand_average')
                    disp('Header indicates this is actually a grand average file.');
                    EPdata.dataType='grand_average';
                end
            end
            EPdata.avgNum(1,:)=hdr.orig.ntrials.accepted;
            EPdata.covNum(1,:)=hdr.orig.ntrials.accepted;
            EPdata.analysis.badTrials(1,:)=hdr.orig.ntrials.rejected;
        end
        
        if strcmp(EPdata.fileFormat,'neuromag_fif')
            if any(strcmp(EPdata.dataType,{'average','grand_average','single_trial'}))
                if isfield(hdr.orig,'evoked')
                    theNames=fieldnames(hdr.orig.evoked);
                    tempVar=struct2cell(hdr.orig.evoked);
                    theAspects=cell2mat(squeeze(tempVar(strcmp('aspect_kind',theNames),1,:)));
                    theEpochs=find(ismember(theAspects,[100 102 103 104]));
                    if isempty(theEpochs)
                        disp('Error: No averages in this dataset.');
                        return
                    else
                        theData=theData(:,:,theEpochs);
                        numWaves=length(theEpochs);
                        numCells=numWaves;
                        %                         EPdata.cellTypes=EPdata.cellTypes(theEpochs);
                        %                         EPdata.events=EPdata.events(theEpochs);
                    end
                    if any(cell2mat(squeeze(tempVar(strcmp('is_smsh',theNames),1,:))))
                        disp('Warning: Data collected using MaxShield active shielding needs to be processed with MaxFilter(tm) to produce reliable results.');
                    end
                    EPdata.avgNum=cell2mat(squeeze(tempVar(strcmp('nave',theNames),1,:)));
                    EPdata.avgNum=double(EPdata.avgNum(:)');
                    EPdata.covNum=EPdata.avgNum;
                    if isempty(EPdata.cellNames)
                        EPdata.cellNames=squeeze(tempVar(strcmp('comment',theNames),1,:));
                    end
                end
            end
        end
        
        ep_tictoc;if EPtictoc.stop;return;end
        if ~strcmp(EPdata.dataType,'single_trial')
            if isempty(numCells)
                numCells = size(theData,3); %number of cells is not ambiguous if not a single_trial file.
            end
            if isempty(EPdata.cellNames)
                if strcmp(EPdata.dataType,'continuous')
                    EPdata.cellNames{1}='cell001';
                else
                    if ~isempty(eventHdr)
                        if sum(strcmp('trial',{eventHdr.type})) ==numCells
                            EPdata.cellNames={eventHdr(strcmp('trial',{eventHdr.type})).value};
                        elseif sum(strcmp('trigger',{eventHdr.type})) ==numCells
                            EPdata.cellNames={eventHdr(strcmp('trigger',{eventHdr.type})).value};
                        else
                            theTrial=ceil([eventHdr.sample]/numSamples);
                            theOffset=mod([eventHdr.sample],numSamples);
                            if (length(find(theOffset==1))==numWaves) && length(theTrial(unique(find(theOffset==1))))==numWaves && ~any(cellfun(@isempty,eventValues(find(theOffset==1))))
                                numCells=length(unique(eventValues(theOffset==1))); %an event is at the offset every time
                                EPdata.cellNames=eventValues(theOffset==1);
                            else
                                for iCell=1:numCells
                                    EPdata.cellNames{iCell,1}=['cell' sprintf('%03.0f',i)];
                                end
                            end
                        end
                        EPdata.cellNames=EPdata.cellNames(:);
                        for iCell=1:length(EPdata.cellNames)
                            if isempty(EPdata.cellNames{iCell})
                                EPdata.cellNames{iCell,1}=['cell' sprintf('%03.0f',iCell)];
                            elseif isnumeric(EPdata.cellNames{iCell})
                                if isnan(EPdata.cellNames{iCell})
                                    EPdata.cellNames{iCell,1}=['cell' sprintf('%03.0f',iCell)];
                                else
                                    EPdata.cellNames{iCell,1}=num2str(EPdata.cellNames{iCell});
                                end
                            end
                        end
                    else
                        for iCell=1:numCells
                            EPdata.cellNames{iCell,1}=['cell' sprintf('%03.0f',iCell)];
                        end
                    end
                end
            end
            
            if length(unique(EPdata.cellNames)) ~= length(EPdata.cellNames)
                cellList=unique(EPdata.cellNames);
                for iCell=1:length(cellList)
                    cellIndex=find(strcmp(cellList{iCell},EPdata.cellNames));
                    if length(cellIndex) > 1
                        for iList=1:length(cellIndex)
                            EPdata.cellNames{cellIndex(iList)}=[EPdata.cellNames{cellIndex(iList)} '-' sprintf('%03d',iList)];
                        end
                    end
                end
                disp('EP Toolkit requires each average condition to have a unique name.  The names are being modified to make them unique.');
            elseif numSubs==1
                if all(strncmp(EPdata.cellNames,'Sub',3)) %if all the cellNames start with 'Sub' then may be a combined subject average file
                    tempVar=char(EPdata.cellNames);
                    subNums=str2double(tempVar(:,4:6));
                    tempCellNames=cellstr(deblank(tempVar(:,7:end)));
                    if all(subNums) %all of them must also have three digit numbers after the Sub
                        if all(hist(subNums,unique(subNums))==length(unique(tempCellNames))) %if the putative cells and subject names are fully crossed
                            EPdata.cellNames=unique(tempCellNames);
                            tempSubNames=cellstr(num2str(subNums));
                            EPdata.subNames=unique(tempSubNames);
                            numSubs=length(EPdata.subNames);
                            numCells=length(EPdata.cellNames);
                            numWaves=numCells;
                            theData2=zeros(length(EPdata.chanNames),numSamples,numCells,numSubs);
                            for iCell=1:numCells
                                theData2(:,:,strcmp(tempCellNames{iCell},EPdata.cellNames),strcmp(tempSubNames{iCell},EPdata.subNames))=theData(:,:,iCell);
                            end
                            for iSub=1:numSubs
                                EPdata.subNames{iSub,1}=['sub' strtrim(EPdata.subNames{iSub})];
                            end
                        end
                    end
                end
            end
        end
        
        ep_tictoc;if EPtictoc.stop;return;end
        if strcmp(EPdata.dataType,'single_trial')
            
            if isempty(EPdata.cellNames)
                if ~isempty(cellLabels)
                    if isscalar(cellLabels)
                        EPdata.cellNames=repmat(cellLabels,numWaves,1);
                    else
                        EPdata.cellNames=cellLabels;
                    end
                elseif strcmp(EPdata.fileFormat,'egi_mff')
                    %in egi_mff_v3,the cell name labels appear in events with zero duration.  Unfortunately, Arno is declining to tell me
                    %how this information is meant to be extracted from his code's output so I'm having to make a best guess.
                    %if there are any other events with zero duration, then this approach will not work.
                    EPdata.cellNames=eventValues(find([eventHdr.duration]==0))';
                    if length(EPdata.cellNames) ~= numWaves
                        EPdata.cellNames=cell(0);
                        disp('EP Toolkit was unable to determine the cell names using mffmatlabio.');
                        disp('If you contact me, I am willing to see if I can come up with a solution for you.');
                        disp('Otherwise, you will have to seek further help from EGI support or Arno Delorme.')
                    end
                elseif strcmp(EPdata.fileFormat,'ns_eeg')
                    for i=1:numWaves
                        EPdata.cellNames{i,1}='cell001';
                    end
                elseif isscalar(eventHdr) && strcmp(eventHdr.type,'average')
                    for i=1:numWaves
                        EPdata.cellNames{i,1}='cell001';
                    end
                elseif length(find(~strcmp(eventValues,'trial')))==numWaves
                    EPdata.cellNames=eventValues(~strcmp(eventValues,'trial'))';
                elseif strcmp(EPdata.fileFormat,'egi_mff') && (length(hdr.orig.epoch)==numWaves)
                    for iWave=1:numWaves
                        EPdata.cellNames{iWave,1}=hdr.orig.epoch(iWave).eventlabel{1};
                    end
                elseif length(find(strcmp({eventHdr.type},'trial')))==numWaves
                    EPdata.cellNames=eventValues(strcmp({eventHdr.type},'trial'))';
                    if strcmp(EPdata.fileFormat,'egi_mff_v1')
                        EPdata.recTime=[hdr.orig.epochdef(:,3)];
                    end
                else
                    if ~isempty(eventHdr)
                        theTrial=ceil([eventHdr.sample]/numSamples);
                        theOffset=mod([eventHdr.sample],numSamples);
                        if (length(find(theOffset==1))==numWaves) && length(theTrial(unique(find(theOffset==1))))==numWaves && all(cellfun(@isempty,eventValues))
                            EPdata.cellNames=eventValues(theOffset==1);
                        end
                        n = hist(theOffset,length(unique(theOffset)));
                        if isscalar(find(n==numWaves)) %if there is one and only one sample for which the number of events equals the number of trials
                            uniqueEventSamples=unique(theOffset);
                            for i=1:length(uniqueEventSamples)
                                whichEvents=find(theOffset==uniqueEventSamples(i));
                                if length(whichEvents) == numWaves
                                    if strcmp(eventHdr(whichEvents(1)).type,'trigger')
                                        if ~isempty(eventHdr(whichEvents(1)).value)
                                            EPdata.cellNames=eventValues(whichEvents);
                                        end
                                        
                                    end
                                    continue %no need to go through entire loop since found the one
                                end
                            end
                        end
                    end
                end
            end
            
            if isempty(EPdata.cellNames) %if was unable to determine cell names assume just one cell
                disp('Was unable to identify distinct cell names from the event header so assuming all just one condition.');
                for i=1:numWaves
                    EPdata.cellNames{i,1}='cell01';
                end
            end
            
            for i=1:length(EPdata.cellNames)
                if isnumeric(EPdata.cellNames{i})
                    EPdata.cellNames{i,1}=num2str(EPdata.cellNames{i});
                end
            end
            
            uniqueCellNames=unique(EPdata.cellNames);
            numCells=length(uniqueCellNames);
            if ~isempty(trialLabels) && ~any(cell2mat(trialLabels)==0)
                EPdata.trialNames=trialLabels;
            else
                trialCounter=zeros(length(uniqueCellNames),1);
                for theCell=1:length(EPdata.cellNames)
                    newCell=strmatch(EPdata.cellNames(theCell),uniqueCellNames,'exact');
                    trialCounter(newCell)=trialCounter(newCell)+1;
                    EPdata.trialNames(theCell,1)=trialCounter(newCell);
                end
            end
            EPdata.trialSpecs=cell(numWaves,0);
        end
    end

    
    numRsubs=numSubs-numVsubs;
    
    ep_tictoc;if EPtictoc.stop;return;end
    if ~isempty(EPdata.freqNames) && isscalar(EPdata.timeNames)
        EPdata.timeNames=[]; %may assume that data is spectral data so no time points
    end
    
    if isempty(EPdata.subNum)
        EPdata.subNum=ones(numRsubs,numCells);
    end
    
    if isempty(EPdata.avgNum)
        if any(strcmp(EPdata.dataType,{'continuous','single_trial'}))
            EPdata.avgNum=ones(numRsubs,numCells);
        else
            EPdata.avgNum=zeros(numRsubs,numCells); %number in average is unknown
        end
    end
    
    if isempty(EPdata.covNum)
        if any(strcmp(EPdata.dataType,{'continuous','single_trial'}))
            EPdata.covNum=ones(numRsubs,numCells);
        else
            EPdata.covNum=zeros(numRsubs,numCells); %number in average is unknown
        end
    end
    
    
    if strcmp(EPdata.dataType,'continuous')
        EPdata.avgNum=ones(1,1);
        EPdata.covNum=ones(1,1);
        EPdata.subNum=ones(1,1);
    end
    
    if strcmp(EPdata.dataType,'single_trial')
        EPdata.avgNum=ones(1,length(EPdata.trialNames));
        EPdata.covNum=ones(1,length(EPdata.trialNames));
        EPdata.subNum=ones(1,length(EPdata.trialNames));
    end
    
    if ~isempty(sevenDdata) && (mod(size(theData,3),numRsubs) ~= 0) && strcmp(EPdata.dataType,'average')
        msg{1}=['Number of subjects (' num2str(numRsubs) ') does not divide evenly into the number of epochs (' num2str(size(theData,3)) ').'];
        [msg]=ep_errorMsg(msg);
        return
    end
    
    if isempty(EPdata.ename)
        if any(strcmp(EPdata.fileFormat,{'egi_egia','egi_egis'}))
            EPdata.ename = hdr.orig.ename;
        else
            EPdata.ename = [];
        end
    end
    
    if size(EPdata.ename,1) > size(EPdata.ename,2)
        EPdata.ename=EPdata.ename';
        if size(EPdata.ename,1) >1
            EPdata.ename=[];
        end
    end
    
    if ~isempty(EPdata.ename)
        EPdata.ename=strtrim(EPdata.ename);
    end
    
    if sum(subjects)==0
        subjects=[];
    end
    
    if isempty(EPdata.recTime)
        EPdata.recTime=[1:EPdata.Fs:EPdata.Fs*(numWaves-1)+1]';
    end
    
    if ~isempty(cells) && (max(cells) > numCells)
        msg{1}=['Largest cell to be kept (' num2str(max(cells)) ') is larger than number of cells (' num2str(numCells) ') in dataset.'];
        [msg]=ep_errorMsg(msg);
        return
    end
    
    if (numWaves*numRsubs*numFacs*numFreqs*numRels ~= size(theData,3))
        msg{1}=['Number of epochs (' num2str(numWaves) ') times subjects (' num2str(numSubs) ') times factors (' num2str(numFacs) ') times frequencies (' num2str(numFreqs) ') times relations (' num2str(numRels) ') does not equal number of observations (' num2str(size(theData,3)) ').'];
        [msg]=ep_errorMsg(msg);
        return
    end
    
    if numChans==1
        chanStr='channel';
    else
        chanStr='channels';
    end
    if numPoints==1
        pointStr='point';
    else
        pointStr='points';
    end
    if isscalar(EPdata.trialNames)
        trialStr='trial';
    else
        trialStr='trials';
    end
    if numCells==1
        cellStr='cell';
    else
        cellStr='cells';
    end
    if numSubs==1
        subStr='subject';
    else
        subStr='subjects';
    end
    if numFacs==1
        facStr='factor';
    else
        facStr='factors';
    end
    if numFreqs==1
        freqStr='frequency';
    else
        freqStr='frequencies';
    end
    if numRels==1
        relStr='relation';
    else
        relStr='relations';
    end

    ep_tictoc;if EPtictoc.stop;return;end
    
    if strcmp(silent,'off')
        disp(['Read in ' num2str(numChans) ' ' chanStr '.'])
        if ~isempty(EPdata.timeNames)
            disp(['Read in ' num2str(numPoints) ' ' pointStr '.'])
        end
        if strcmp(EPdata.dataType,'single_trial')
            disp(['Read in ' num2str(length(unique(EPdata.cellNames))) ' ' cellStr ' with a total of ' num2str(length(EPdata.trialNames)) ' ' trialStr '.']);
        elseif strcmp(EPdata.dataType,'average')
            disp(['Read in ' num2str(numCells) ' ' cellStr '.']);
        end
        disp(['Read in ' num2str(numSubs) ' ' subStr '.'])
        if ~isempty(EPdata.facNames)
            disp(['Read in ' num2str(numFacs) ' ' facStr '.'])
        end
        if ~isempty(EPdata.freqNames)
            disp(['Read in ' num2str(numFreqs) ' ' freqStr '.'])
        end
        if ~isempty(EPdata.relNames)
            disp(['Read in ' num2str(numRels) ' ' relStr '.'])
        end
    end
    
    if strcmp(EPdata.dataType,{'continuous'}) && ~isempty(cells)
        disp('Cannot select subset of cells in a continuous raw data file.');
        cells=[];
    end
    
    if ~isempty(samples)
        if max(samples) > size(theData,2)
            msg{1}='Largest sample is larger than number of timepoints in data.';
            [msg]=ep_errorMsg(msg);
            return
        end
        if min(samples) < 0
            error('No negative samples.');
        end
        
    end
    
    if ~isempty(channels)
        if max(channels) > size(theData,1)
            msg{1}='Largest channel is larger than number of channels in data.';
            [msg]=ep_errorMsg(msg);
            return
        end
        
        if min(channels) < 0
            error('No negative channels.');
        end
    end
    
    if ~isempty(subjects)
        if min(subjects) < 0
            error('No negative subjects.');
        end
        if ~isempty(subjects) && (max(subjects) > numSubs)
            msg{1}=['Largest subject to be kept (' num2str(max(subjects)) ') is larger than number of subjects (' num2str(numSubs) ') indicated for dataset.'];
            [msg]=ep_errorMsg(msg);
            return
        end
    end
    
    %Reorganize waveforms and events and extract events from EventHdr
    %mff is file time whereas everything else from FieldTrip is real time
    
    ep_tictoc;if EPtictoc.stop;return;end
    if ~isempty(eventHdr)
        EPdata.events=cell(numRsubs,numWaves);
        %EPdata.events{1,1}=struct('value','','sample',[],'type','','duration',[],'keys',struct('code','','data','','datatype','','description',''));
        samples_trials = [eventHdr(find(strcmp('trial', {eventHdr.type}))).sample]; %find "trial" events
        if any(strcmp('TRSP',{eventHdr.value})) && isfield(eventHdr,'orig') && ~strcmp(EPdata.dataType,'continuous')
            for i=1:length(eventHdr(min(find(strcmp('TRSP',{eventHdr.value})))).orig.keys) %assume TRSP keys same for all trials
                EPdata.trialSpecNames{end+1}=eventHdr(min(find(strcmp('TRSP',{eventHdr.value})))).orig.keys(i).code;
            end
        end
        for theEvent=1:length(eventHdr)
            ep_tictoc;if EPtictoc.stop;return;end
            if ~isempty(eventHdr(theEvent).sample) && ~isnan(eventHdr(theEvent).sample)
                if isempty(samples_trials)
                    epoch=floor((eventHdr(theEvent).sample-1)/numSamples)+1;
                    sample=mod((eventHdr(theEvent).sample-1),numSamples)+1;
                else
                    epoch=max(find(samples_trials<=eventHdr(theEvent).sample));
                    sample=eventHdr(theEvent).sample-samples_trials(epoch)+1;
                end
                if any(strcmp(EPdata.fileFormat,{'egi_sbin','egi_mff'})) && subjectsGrouped
                    %cells are grouped together by subject.
                    subject=floor((epoch-1)/numWaves)+1;
                    cellTrial=mod((epoch-1),numWaves)+1;
                else
                    %assume that subjects are grouped together by cell.
                    cellTrial=floor((epoch-1)/numRsubs)+1;
                    subject=mod((epoch-1),numRsubs)+1;
                end
                if (epoch > numWaves*numRsubs) || epoch < 1
                    disp('Warning: event falls outside recorded data.  Discarded.');
                else
                    if strcmp(eventHdr(theEvent).value,'TRSP') && strcmp(EPdata.dataType,'single_trial') && isfield(eventHdr,'orig')
                        for trsp=1:length(eventHdr(theEvent).orig.keys)
                            EPdata.trialSpecs{cellTrial,find(strcmp(eventHdr(theEvent).orig.keys(trsp).code,EPdata.trialSpecNames))}=eventHdr(theEvent).orig.keys(trsp).data;
                        end
                    elseif ~isempty(EPdata.trialSpecNames)
                        oneEvent=[];
                        oneEvent.type=eventHdr(theEvent).type;
                        oneEvent.sample=sample;
                        oneEvent.value=eventHdr(theEvent).value;
                        oneEvent.duration=eventHdr(theEvent).duration;
                        if isfield(eventHdr,'keys')
                            oneEvent.keys=eventHdr(theEvent).keys;
                        else
                            oneEvent.keys=struct('code','','data','','datatype','','description','');
                        end
                        EPdata.events{subject,cellTrial}=[EPdata.events{subject,cellTrial} oneEvent];
                    elseif any(strcmp(EPdata.fileFormat,{'eeglab_set','eeglab_erp'}))
                        oneEvent=[];
                        oneEvent.type=eventHdr(theEvent).type;
                        oneEvent.sample=sample;
                        oneEvent.value=eventHdr(theEvent).value;
                        oneEvent.duration=eventHdr(theEvent).duration;
                        oneEvent.keys=struct('code','','data','','datatype','','description','');
                        nameList=fieldnames(eventHdr(theEvent));
                        nameList=setdiff(nameList,{'type','value','sample','offset','duration'});
                        for iKey=1:length(nameList)
                            eval(['fieldEmpty=isempty(eventHdr(theEvent).' nameList{iKey} ');']);
                            if ~fieldEmpty
                                fieldLabel=nameList{iKey};
                                if ~isempty(strfind(fieldLabel,'hash_'))
                                    fieldLabel=strrep(fieldLabel,'hash_','#');
                                end
                                if ~isempty(strfind(fieldLabel,'plus_'))
                                    fieldLabel=strrep(fieldLabel,'plus_','+');
                                end
                                eval(['oneEvent.keys(end+1).code=''' fieldLabel ''';']);
                                eval(['oneEvent.keys(end).data=eventHdr(theEvent).' nameList{iKey} ';']);
                            end
                        end
                        EPdata.events{subject,cellTrial}=[EPdata.events{subject,cellTrial} oneEvent];
                    else
                        oneEvent=[];
                        oneEvent.type=eventHdr(theEvent).type;
                        oneEvent.sample=sample;
                        oneEvent.value=eventHdr(theEvent).value;
                        oneEvent.duration=eventHdr(theEvent).duration;
                        if isempty(oneEvent.value)
                            oneEvent.value=oneEvent.type;
                            oneEvent.type='trigger';
                        end
                        if strcmp(EPdata.fileFormat,'egi_mff_v2')
                            if strcmp(oneEvent.type,'break cnt') || strcmp(oneEvent.value,'break cnt')
                                oneEvent.type='boundary';
                                oneEvent.value='boundary';
                            end
                        end
                        if strcmp(EPdata.fileFormat,'biosemi_bdf')
                            if strcmp(oneEvent.type,'trigger') && strcmp(oneEvent.value,'Epoch')
                                oneEvent.type='boundary';
                                oneEvent.value='boundary';
                            end
                        end
                        if isfield(eventHdr,'orig')
                            if isfield(eventHdr(theEvent).orig,'keys')
                                if strcmp(EPdata.fileFormat,'egi_mff_v1')
                                    if isfield(eventHdr(theEvent).orig.keys,'key')
                                        for iKey=1:length(eventHdr(theEvent).orig.keys)
                                            oneEvent.keys(iKey).code=eventHdr(theEvent).orig.keys(iKey).key.keyCode;
                                            oneEvent.keys(iKey).data=eventHdr(theEvent).orig.keys(iKey).key.data.data;
                                            oneEvent.keys(iKey).datatype=eventHdr(theEvent).orig.keys(iKey).key.data.dataType;
                                            oneEvent.keys(iKey).description='';
                                        end
                                    else
                                        oneEvent.keys=struct('code','','data','','datatype','','description','');
                                    end
                                else
                                    oneEvent.keys=cell2mat(eventHdr(theEvent).orig.keys);
                                end
                            else
                                oneEvent.keys=struct('code','','data','','datatype','','description','');
                            end
                        elseif isfield(eventHdr,'keys') && ~isempty(eventHdr(theEvent).keys)
                            oneEvent.keys=eventHdr(theEvent).keys;
                        else
                            oneEvent.keys=struct('code','','data','','datatype','','description','');
                        end
                        EPdata.events{subject,cellTrial}=[EPdata.events{subject,cellTrial} oneEvent];
                    end
                end
            end
        end
    else
        EPdata.events=cell(numRsubs,length(EPdata.cellNames));
        %EPdata.events{1,1}=struct('value','','sample',[],'type','','duration',[],'keys',struct('code','','data','','datatype','','description',''));
    end
    
    %reorganize data into 6D array
    ep_tictoc;if EPtictoc.stop;return;end
    if isempty(sevenDdata)
        if strcmp(EPdata.fileFormat,'egi_sbin') && subjectsGrouped
            sevenDdata=nan(size(theData,1),size(theData,2),numWaves,numRsubs,numFacs,numFreqs,numRels);
            %cells are grouped together by subject.
            for sub=1:numRsubs
                sevenDdata(:,:,:,sub,1,1,1)=theData(:,:,1+(sub-1)*numWaves:sub*numWaves);
            end
        else
            sevenDdata=reshape(theData,size(theData,1),size(theData,2),numWaves,numRsubs,numFacs,numFreqs,numRels);
        end
    end
end

numVcells=max(0,size(EPdata.GAVsubs,2)-1);
numRcells=numCells-numVcells;

EPdata.data=sevenDdata;

numChan=length(EPdata.chanNames);

% apply cell and sub names if separately specified
ep_tictoc;if EPtictoc.stop;return;end
if ~isempty(cellLabels) && ~all(cellfun(@isempty,cellLabels))
    if isscalar(cellLabels)
        if length(unique(EPdata.trialNames)) ~= length(EPdata.trialNames)
            disp('This file appears to have multiple cell names in it and therefore the cell field of the Single File Mode will be ignored.');
        else
            if ~strcmp(EPdata.dataType,'average') || (isscalar(EPdata.cellNames))
                EPdata.cellNames=repmat(cellLabels,numWaves,1);
            end
        end
    elseif numWaves == length(cellLabels)
        EPdata.cellNames=cellLabels;
    else
        disp(['The number of epochs (' num2str(numWaves) ') did not match the number of condition labels (' num2str(length(cellLabels)) ').']);
    end
end

if ~isempty(subLabels)
    if isscalar(subLabels)
        for iSub=1:numSubs
            EPdata.subNames{iSub}=subLabels{1};
        end
    elseif numSubs == length(subLabels)
        EPdata.subNames=subLabels;
    else
        disp(['The number of subjects (' num2str(numSubs) ') did not match the number of subject labels (' num2str(length(subLabels)) ').']);
    end
end

if ~isempty(sessLabels) && ~((isscalar(sessLabels)) && isempty(sessLabels{1}))
    if isscalar(sessLabels)
        EPdata.sessNames=sessLabels;
        for iSub=1:numRsubs
            EPdata.sessNums(iSub)=1;
        end
    elseif numSess == length(sessLabels)
        EPdata.sessNames=unique(sessLabels);
        for iSub=1:numRsubs
            theSess=find(strcmp(sessLabels{iSub},EPdata.sessNames));
            if isempty(theSess)
                EPdata.sessNums(iSub)=0;
            else
                EPdata.sessNums(iSub)=theSess;
            end
        end
    else
        disp(['The number of sessions (' num2str(numSess) ') did not match the number of session labels (' num2str(length(sessLabels)) ').']);
    end
end

numEpochs=numWaves;
if strcmp(EPdata.dataType,'continuous')
    numEpochs=floor(size(EPdata.data,2)/ceil(EPdata.Fs)); %excess time points are tacked onto final epoch
    if numEpochs == 0
        numEpochs =1;
    end
end

%add Type information
if isempty(EPdata.chanTypes)
    EPdata.chanTypes=cellstr(repmat('EEG',numChan,1));
end

if isempty(EPdata.cellTypes)
    EPdata.cellTypes=cellstr(repmat('SGL',numWaves,1));
end

if isempty(EPdata.facTypes) && ~isempty(EPdata.facNames)
    EPdata.facTypes=cellstr(repmat('SGL',numFacs,1));
end

if isempty(EPdata.subTypes)
    switch EPdata.dataType
        case {'single_trial', 'continuous'}
            EPdata.subTypes=cellstr(repmat('RAW',numSubs,1));
        case 'average'
            EPdata.subTypes=cellstr(repmat('AVG',numSubs,1));
        case 'grand_average'
            EPdata.subTypes=cellstr(repmat('GAV',numSubs,1));
            EPdata.dataType='average';
    end
end

%add time unit information
if isempty(EPdata.timeUnits)
    EPdata.timeUnits='ms';
end

%add preprocessing information
if isempty(EPdata.analysis.blinkTrial)
    EPdata.analysis.blinkTrial=zeros(numRsubs,numEpochs);
end
if isempty(EPdata.analysis.saccadeTrial)
    EPdata.analysis.saccadeTrial=zeros(numRsubs,numEpochs);
end
if isempty(EPdata.analysis.saccadeOnset)
    EPdata.analysis.saccadeOnset=zeros(numRsubs,numEpochs);
end
if isempty(EPdata.analysis.moveTrial)
    EPdata.analysis.moveTrial=zeros(numRsubs,numEpochs);
end
if isempty(EPdata.analysis.badTrials)
    EPdata.analysis.badTrials=zeros(numRsubs,numEpochs);
end
if isempty(EPdata.analysis.badChans)
    EPdata.analysis.badChans=zeros(numRsubs,numEpochs,length(EPdata.chanNames));
end

ep_tictoc;if EPtictoc.stop;return;end
if ~strcmp(EPdata.fileFormat,'ep_mat')
    %if an EEG channel has any NaN values, as from EEGlab automatic editing, it
    %will be set to bad data.
    if any(any(any(any(any(any(any(isnan(EPdata.data))))))))
        for iSub=1:size(EPdata.data,4)
            for iCell=1:size(EPdata.data,3)
                for iChan=1:size(EPdata.data,1)
                    if any(strcmp(EPdata.chanTypes{iChan},{'EEG','REG'}))
                        if any(any(any(any(isnan(squeeze(EPdata.data(iChan,:,iCell,iSub)))))))
                            %                             EPdata.data(iChan,:,iCell,iSub,:,:,:)=0;
                            if any(strcmp(EPdata.dataType,{'continuous','single_trial'}))
                                EPdata.analysis.badChans(iSub,iCell,iChan)=-1;
                            elseif strcmp(EPdata.dataType,'average')
                                EPdata.analysis.badChans(iSub,iCell,iChan)=NaN;
                            end
                        end
                    end
                end
            end
        end
    end
end

%channel coordinate information
if any(strcmp(EPdata.fileFormat,{'eeglab_set','eeglab_erp'})) && ~noInternal
    if isfield(hdr,'orig')
        if isfield(hdr.orig,'chanlocs')
            if isfield(hdr.orig,'chaninfo')
                if isfield(hdr.orig.chaninfo,'filename')
                    [pathstr, name, ext] = fileparts(hdr.orig.chaninfo.filename);
                    EPdata.ced=[name ext];
                else
                    EPdata.ced='internal';
                end
            else
                EPdata.ced='internal';
            end
            EPdata.eloc=hdr.orig.chanlocs;
            if isempty(EPdata.eloc)
                EPdata.ced=[];
            else
                EPdata=ep_updateEPfile(EPdata);
            end
        end
    end
end
if any(strcmp(EPdata.fileFormat,{'eeglab_set','eeglab_erp','ep_mat'})) && noInternal
    if ~isempty(EPdata.eloc)
        disp('Ignoring internal electrode coordinates information per noInternal preferences setting.')
    end
    EPdata.eloc=[];
    if any(strcmp(EPdata.ced,{'internal';'none'}))
        EPdata.ced=[];
    end
end

ep_tictoc;if EPtictoc.stop;return;end
if strcmp(EPdata.fileFormat,'neuromag_fif')
    %the .chs coordinates will always be Neuromag, meaning in meters, with +Y through the nose and +X through the right pre-auricular fiducial.
    %but since the .chs coordinates don't include the fiducials, which are important, the .dig coordinates will be used
    %instead unless they are not available to keep all the coordinates consistent with each other.
    %since the .dig coordinates cannot be assumed to have any specific coordinate system or orientation, the head rotation preferences
    %will be applied.
    
    if isempty(EPdata.eloc) && isempty(EPdata.ced)
        implicitCoords=[];
        if isfield(hdr.orig,'dig')
            if isfield(hdr.orig.dig,'kind')
                FIDchans=find([hdr.orig.dig.kind]==1);
                if length(FIDchans)==3
                    if isempty(EPdata.implicit)
                        implicitCoords(1).labels='LPA';
                        implicitCoords(1).X=hdr.orig.dig(FIDchans(1)).r(1);
                        implicitCoords(1).Y=hdr.orig.dig(FIDchans(1)).r(2);
                        implicitCoords(1).Z=hdr.orig.dig(FIDchans(1)).r(3);
                        implicitCoords(2).labels='nasion';
                        implicitCoords(2).X=hdr.orig.dig(FIDchans(2)).r(1);
                        implicitCoords(2).Y=hdr.orig.dig(FIDchans(2)).r(2);
                        implicitCoords(2).Z=hdr.orig.dig(FIDchans(2)).r(3);
                        implicitCoords(3).labels='RPA';
                        implicitCoords(3).X=hdr.orig.dig(FIDchans(3)).r(1);
                        implicitCoords(3).Y=hdr.orig.dig(FIDchans(3)).r(2);
                        implicitCoords(3).Z=hdr.orig.dig(FIDchans(3)).r(3);
                    end
                end
                %I don't understand MEG well enough yet to implement reading its coordinate information so will be set to [0 0 0].
                EEGchans=find(strcmp(EPdata.chanTypes,'EEG'));
                badDigFlag=0;
                kindCode=3;
                if (length(EEGchans)+1) == length(find([hdr.orig.dig.kind]==3))
                    %there is an implicit reference so add it to the data as an explicit channel.
                    theRef=intersect(find([hdr.orig.dig.ident]==0),find([hdr.orig.dig.kind]==3));
                    addRefFlag=0;
                    if length(theRef)==1
                        if ~isempty(EPdata.facVecS)
                            EPdata.facVecS(end+1,:)=zeros(size(EPdata.facVecS,2),1);
                        end
                        EPdata.data(end+1,:,:,:,:,:,:)=0;
                        if ~isempty(EPdata.noise)
                            EPdata.noise(end+1,:,:,:,:)=zeros(1,size(EPdata.noise,2),size(EPdata.noise,3),size(EPdata.noise,4),size(EPdata.noise,5));
                        end
                        if ~isempty(EPdata.covAVE)
                            if size(EPdata.covAVE,7)==1
                                EPdata.covAVE(end+1,:,:,:,:,:,:)=zeros(1,size(EPdata.covAVE,2),size(EPdata.covAVE,3),size(EPdata.covAVE,4),size(EPdata.covAVE,5),size(EPdata.covAVE,6),1);
                            else
                                EPdata.covAVE(end+1,:,:,:,:,:,:)=zeros(1,size(EPdata.covAVE,2),size(EPdata.covAVE,3),size(EPdata.covAVE,4),size(EPdata.covAVE,5),size(EPdata.covAVE,6),size(EPdata.covAVE,7)+1);
                            end
                        end
                        EPdata.analysis.badChans(:,:,end+1)=zeros(size(EPdata.analysis.badChans,1),size(EPdata.analysis.badChans,2),1); %implicit reference channels are not bad
                        EPdata.chanNames{end+1}='REF';
                        EPdata.chanTypes{end+1}='EEG';
                        EPdata.reference.original=length(EPdata.chanNames);
                        EPdata.reference.current=EPdata.reference.original;
                        addRefFlag=1;
                    else
                        disp('Was unable to determine which channel was the implicit reference.');
                        badDigFlag=1;
                    end
                elseif length(EEGchans) ~= length(find([hdr.orig.dig.kind]==3))
                    if length(EEGchans) == (length(find([hdr.orig.dig.kind]==1))+length(find([hdr.orig.dig.kind]==3)))
                        disp('The number of electrode coordinates provided by the fif files do not match up with the number of EEG channels,');
                        disp('but did add up when fiducial points were also included, so on the assumption of user error, all will be included as EEG channels.');
                        kindCode=[1; 3];
                    else
                        badDigFlag=1;
                    end
                end
                if ~badDigFlag
                    digEEGchans=find(ismember(kindCode,[hdr.orig.dig.kind]));
                    tempEloc=struct('labels',EPdata.chanNames);
                    for iChan=1:length(EEGchans)
                        theChan=EEGchans(iChan);
                        tempEloc(theChan).X=hdr.orig.dig(digEEGchans(iChan)).r(1);
                        tempEloc(theChan).Y=hdr.orig.dig(digEEGchans(iChan)).r(2);
                        tempEloc(theChan).Z=hdr.orig.dig(digEEGchans(iChan)).r(3);
                    end
                    if addRefFlag
                        tempEloc(end).X=hdr.orig.dig(theRef).r(1);
                        tempEloc(end).Y=hdr.orig.dig(theRef).r(2);
                        tempEloc(end).Z=hdr.orig.dig(theRef).r(3);
                    end
                    %rotate head so it is pointing upwards
                    EPdata.eloc=tempEloc;
                end
            end
        end
        if isempty(EPdata.eloc) && isfield(hdr.orig,'chs') %if it didn't work out with .dig, try .chs for electrode coordinates
            if isfield(hdr.orig.chs,'eeg_loc')
                tempEloc=struct('labels',hdr.label);
                refCoord=[];
                for iChan=1:length(hdr.label)
                    tempVar=hdr.orig.chs(iChan).eeg_loc;
                    if ~isempty(tempVar)
                        tempEloc(iChan).X=tempVar(1,1);
                        tempEloc(iChan).Y=tempVar(2,1);
                        tempEloc(iChan).Z=tempVar(3,1);
                        if isempty(refCoord) && (size(tempVar,2) == 2)
                            refCoord.X=tempVar(1,2);
                            refCoord.Y=tempVar(2,2);
                            refCoord.Z=tempVar(3,2);
                        end
                    end
                end
                
                %add the reference channel, which is implicit if coordinates are given.  Assume a common reference site.
                if ~isempty(refCoord) && any([refCoord.X refCoord.Y refCoord.Z]) %if no reference, then the coords will be [0,0,0] in FIFF files
                    M1=length(hdr.label)+1;
                    EPdata.reference.original=M1;
                    EPdata.reference.current=EPdata.reference.original; %assume that if there are implicit references then they were the original reference and still are
                    tempEloc(M1).X=refCoord.X;
                    tempEloc(M1).Y=refCoord.Y;
                    tempEloc(M1).Z=refCoord.Z;
                    tempEloc(M1).labels='REF';
                    EPdata.chanTypes{M1}='EEG';
                    EPdata.data(end+1,:,:,:,:,:,:)=0;
                end
            end
        end
        if ~isempty(EPdata.eloc)
            EPdata.ced='internal';
            try
                EPdata.eloc = convertlocs(EPdata.eloc, 'cart2all');
                [EPdata.eloc.type] = EPdata.chanTypes{:};
            catch
                disp('Warning: Reading of electrode coordinates from file has failed.');
                EPdata.eloc=[];
            end
            if isempty(EPdata.implicit) && ~isempty(implicitCoords)
                try
                    EPdata.implicit = convertlocs(implicitCoords, 'cart2all');
                    [EPdata.implicit.type]=deal('FID');
                catch
                    disp('Warning: Reading of fiducial coordinates from file has failed.');
                end
            end
        end
    end
    %convert to microvolts
    for iChan=1:length(hdr.orig.chs)
        if hdr.orig.chs(iChan).kind==2 %EEG channels
            EPdata.data(iChan,:,:,:,:,:,:)=EPdata.data(iChan,:,:,:,:,:,:)*10^double(6+hdr.orig.chs(iChan).unit_mul);
        end
    end
    %not sure if cov matrices are normally included in the fif average file but just in case.
    try
        [tempVar] = mne_read_noise_cov(EPdata.fileName);
        tempVar=tempVar*10^12; %convert to microvolts
        EPdata.cov.covMatrix(1,:,:)=tempVar;
        EPdata.cov.Nq=NaN;
    catch
        %no cov file
    end
    for iBad=1:length(hdr.orig.bads)
        if any(strcmp(EPdata.dataType,{'average','grand_average'}))
            EPdata.analysis.badChans(:,:,find(strcmp(hdr.orig.bads{iBad},EPdata.chanNames)))=NaN;
        else
            EPdata.analysis.badChans(:,:,find(strcmp(hdr.orig.bads{iBad},EPdata.chanNames)))=-1;
        end
    end
end

ep_tictoc;if EPtictoc.stop;return;end
if strcmp(EPdata.fileFormat,'egi_mff')
    listcolformat = ep_elocFormat;
    if isfield(hdr,'orig') && isfield(hdr.orig,'chaninfo') && isfield(hdr.orig.chaninfo,'nodatchans')  %fiducials
        for iChan=1:length(hdr.orig.chaninfo.nodatchans)
            for iField=1:length(listcolformat)
                if isfield(hdr.orig.chaninfo.nodatchans,listcolformat{iField})
                    eval(['theField=hdr.orig.chaninfo.nodatchans(' num2str(iChan) ').' listcolformat{iField} ';'])
                else
                    theField='';
                end
                if isnumeric(theField)
                    if isempty(theField)
                        eval(['EPdata.implicit(' num2str(iChan) ').' listcolformat{iField} '=[];'])
                    else
                        eval(['EPdata.implicit(' num2str(iChan) ').' listcolformat{iField} '=' num2str(theField) ';'])
                    end
                else
                    eval(['EPdata.implicit(' num2str(iChan) ').' listcolformat{iField} '=''' theField ''';'])
                end
            end
            if isfield(hdr.orig.chaninfo.nodatchans,'description')
                switch hdr.orig.chaninfo.nodatchans(iChan).description
                    case 'Nasion'
                        EPdata.implicit(iChan).labels='Nz';
                    case 'Left periauricular point'
                        EPdata.implicit(iChan).labels='LPA';
                    case 'Right periauricular point'
                        EPdata.implicit(iChan).labels='RPA';
                end
            end
        end
    end
    
    if isfield(hdr.orig,'chanlocs') && ~isempty(hdr.orig.chanlocs)
        for iChan=1:length(hdr.orig.chanlocs)
            if isfield(hdr.orig.chanlocs,'type')
                switch hdr.orig.chanlocs(iChan).type
                    case {'eeg','EEG'}
                        EPdata.chanTypes{iChan}='EEG';
                    case {'ecg','ECG'}
                        EPdata.chanTypes{iChan}='ECG';
                    case {'pns','PNS'}
                        if strcmp(hdr.orig.chanlocs(iChan).labels,'ECG')
                            EPdata.chanTypes{iChan}='ECG';
                        else
                            EPdata.chanTypes{iChan}='ANS';
                        end
                    otherwise
                        EPdata.chanTypes{iChan}='ANS'; %default to ANS if unknown or unrecognized
                end
            end
            for iField=1:length(listcolformat)
                if isfield(hdr.orig.chanlocs,listcolformat{iField})
                    eval(['theField=hdr.orig.chanlocs(' num2str(iChan) ').' listcolformat{iField} ';'])
                else
                    theField='';
                end
                if isnumeric(theField)
                    if isempty(theField)
                        eval(['EPdata.eloc(' num2str(iChan) ').' listcolformat{iField} '=[];'])
                    else
                        eval(['EPdata.eloc(' num2str(iChan) ').' listcolformat{iField} '=' num2str(theField) ';'])
                    end
                else
                    eval(['EPdata.eloc(' num2str(iChan) ').' listcolformat{iField} '=''' theField ''';'])
                end
            end
        end
        if isfield(hdr.orig.chaninfo,'nodatchans')
        end
        %rotate head so it is pointing upwards
        tempEloc=EPdata.eloc;
        tempFID=EPdata.implicit;
        EPdata.ced='internal';
    end
end

% if any(strcmp(EPdata.fileFormat,{'egi_mff_v1','egi_mff_v2'}))
%     for iChan=1:length(EPdata.chanTypes)
%         switch hdr.chantype{iChan}
%             case 'eeg'
%                 EPdata.chanTypes{iChan}='EEG';
%             case {'ecg','ECG'}
%                 EPdata.chanTypes{iChan}='ECG';
%             otherwise
%                 EPdata.chanTypes{iChan}='ANS'; %default to ANS if unknown or unrecognized
%         end
%         if any(strcmp(hdr.label{iChan},{'REF','VREF','Cz','vertex reference'}))
%             EPdata.chanTypes{iChan}='EEG';
%         end
%     end
% end

if strcmp(EPdata.fileFormat,'edf')
    for iChan=1:numChan
        if strcmp(EPdata.chanNames{iChan},'ECG')
            EPdata.chanTypes{iChan}='ECG';
        end
        %         if strcmp(EPdata.chanNames{iChan},'Tilt X')
        %             EPdata.chanTypes{iChan}='ACMx';
        %         end
        %         if strcmp(EPdata.chanNames{iChan},'Tilt Y')
        %             EPdata.chanTypes{iChan}='ACMy';
        %         end
        %         if strcmp(EPdata.chanNames{iChan},'Tilt Z')
        %             EPdata.chanTypes{iChan}='ACMz';
        %         end
    end
end

%test out the ced file and set to null if it's no good if EP file gives name of ced file but does not actually have eloc info.  This is legacy code.
ep_tictoc;if EPtictoc.stop;return;end
if ~isempty(EPdata.ced) && ~any(strcmp(EPdata.ced,{'none','internal'})) && isempty(EPdata.eloc)
    if exist(EPdata.ced,'file')
        whichCED=EPdata.ced;
        [pathstr, name, fileSuffix] = fileparts(whichCED);
        EPdata.ced=[name fileSuffix];
    else
        whichCED=which(EPdata.ced);
    end
    if isempty(whichCED)
        disp(['Could not find the file ' EPdata.ced '.  Please put it either in the electrodes folder of the EP Toolkit or in the active directory.']);
        EPdata.ced = [];
        EPdata.eloc = [];
        EPdata.implicit=[];
    else
        try
            eval('EPdata.eloc = ep_readlocsWrapper([whichCED],''filetype'',''chanedit'');');
            if EPtictoc.stop;return;end
            EPdata.implicit=EPdata.eloc(1); %setup up implicit to have the same structure as eloc.
            EPdata.implicit(1)=[];
            EPdata.eloc=[];
        catch
            disp(['The ced file ' EPdata.ced ' did not work for some reason.  The error message was:']);
            disp(lasterr)
            EPdata.ced = [];
            EPdata.eloc = [];
            EPdata.implicit=[];
        end
    end
end

if strcmp(EPdata.ced,'none')
    cedFlag=1;
else
    cedFlag=0;
end

if any(strcmp(EPdata.ced,{'none','internal'})) && isempty(EPdata.eloc)
    EPdata.ced=[];
end

if isempty(EPdata.ced) && ~cedFlag
    if (~isempty(strfind(EPdata.montage,'GSN')) || ~isempty(strfind(EPdata.montage,'Hydrocel')))
        [EPdata.ced]=ep_whichEGIced(EPdata.montage);
    end
    whichCED=EPdata.ced;
    if isempty(EPdata.ced)
        temp=dir(pwd);
        temp2={temp.name};
        if ~all(cellfun(@isempty,regexp(temp2,'.*ced$')))
            pathstr=pwd;
        else
            thisFile = mfilename('fullpath');
            [pathstr, name, ext] = fileparts(thisFile);
            theSeps=findstr(pathstr,filesep);
            pathstr=[pathstr(1:theSeps(end)-1) filesep 'electrodes'];
        end
        [EPdata.ced, pathstr] = ep_getuifile('*.ced',['Electrode Coordinate file (' num2str(numChan) ' channels):'],pathstr);
        whichCED=[pathstr EPdata.ced];
    end
    if isnumeric(EPdata.ced)
        EPdata.ced='none';
        whichCED=EPdata.ced;
    end
else
    whichCED=EPdata.ced;
end

if (~isempty(EPdata.ced) && ~any(strcmp(EPdata.ced,{'none','internal'})) && ~strcmp(EPdata.fileFormat,'ep_mat')) || (~isempty(EPdata.ced) && (isempty(EPdata.eloc) || all(isempty([EPdata.eloc.radius])) || all(isempty([EPdata.eloc.theta]))) && ~strcmp('none',EPdata.ced))
    if ~exist('hdr','var')
        hdr=[];
    end
    if ~isempty(EPdata.eloc)
        if  all(isempty([EPdata.eloc.theta])) && all(~isempty([EPdata.eloc.radius]))
            disp('Warning: There may be something wrong with your electrode coordinate information.');
            disp('All your theta information is missing.');
            disp('Check your ced file, if that is what you used, to verify that it is properly formed.');
        end
        EPdata.eloc(cellfun(@isempty,{EPdata.eloc.labels}))=[]; %if there are blank extra elocs added by something like SMI remove them
    end
    [EPdata]=ep_addEloc(whichCED,EPdata.eloc,EPdata.fileFormat,EPdata,silent);
    if EPtictoc.stop;return;end
    if isempty(EPdata)
       return 
    end
end

if isempty(EPdata.eloc)
    %if no electrode coordinates are available, use the canonical 10-05 ones.
    %for CED Y is +left -right, Z is anterior+ posterior-, X is dorsal+ ventral -.  The sfp file needs to map X to Y, Y to X, and Z to Z.
    EPdata.ced='none';
    tempCED='Standard-10-5-Cap385-VEOG.ced';
    eloc = ep_readlocsWrapper(tempCED,'filetype','chanedit');
    if EPtictoc.stop;return;end
    if isempty(eloc)
        return
    end
    EPdata.eloc=ep_elocFormat('initialize');
    for iChan=1:length(EPdata.chanNames)
        theEloc=[];
        if strcmp(EPdata.chanTypes{iChan},'EEG')
            switch EPdata.chanNames{iChan}
                case 'LHEOG'
                    theEloc=find(strcmpi('AFp9',{eloc.labels}));
                case 'RHEOG'
                    theEloc=find(strcmpi('AFp10',{eloc.labels}));
                otherwise
                    theEloc=find(strcmpi(EPdata.chanNames{iChan},{eloc.labels}));
            end
        end
        EPdata.eloc(iChan).labels=EPdata.chanNames{iChan};
        if ~isempty(theEloc)
            EPdata.eloc(iChan).theta=eloc(theEloc).theta;
            EPdata.eloc(iChan).radius=eloc(theEloc).radius;
            EPdata.eloc(iChan).X=eloc(theEloc).X;
            EPdata.eloc(iChan).Y=eloc(theEloc).Y;
            EPdata.eloc(iChan).Z=eloc(theEloc).Z;
            EPdata.eloc(iChan).sph_theta=eloc(theEloc).sph_theta;
            EPdata.eloc(iChan).sph_phi=eloc(theEloc).sph_phi;
            EPdata.eloc(iChan).sph_radius=eloc(theEloc).sph_radius;
        else
            EPdata.eloc(iChan).theta=[];
            EPdata.eloc(iChan).radius=[];
            EPdata.eloc(iChan).X=[];
            EPdata.eloc(iChan).Y=[];
            EPdata.eloc(iChan).Z=[];
            EPdata.eloc(iChan).sph_theta=[];
            EPdata.eloc(iChan).sph_phi=[];
            EPdata.eloc(iChan).sph_radius=[];
        end
        EPdata.eloc(iChan).type=EPdata.chanTypes{iChan};
    end
end

if ~isempty(EPdata.eloc) && ~any([EPdata.eloc.cX])
    %add canonical 10-05 electrode coordinates
    [outEloc outFid] = ep_transformEloc(EPdata.eloc, EPdata.implicit, '', '', EPdata.chanNames, '', EPdata.montage, '');
    if EPtictoc.stop;return;end
    if isempty(outEloc)
        return
    end
    if ~isempty(outEloc)
        for iChan=1:length(EPdata.eloc)
            EPdata.eloc(iChan).cX=outEloc(iChan).X;
            EPdata.eloc(iChan).cY=outEloc(iChan).Y;
            EPdata.eloc(iChan).cZ=outEloc(iChan).Z;
        end
        for iChan=1:length(EPdata.implicit)
            EPdata.implicit(iChan).cX=outFid(iChan).X;
            EPdata.implicit(iChan).cY=outFid(iChan).Y;
            EPdata.implicit(iChan).cZ=outFid(iChan).Z;
        end
    end
end

% if strcmp(EPdata.ced,'internal') && any(strcmp(EPdata.fileFormat,{'egi_mff_v1';'egi_mff_v2'}))
%     %since mff files have internal electrode coordinates, will not normally be using ced files to designate the channels types.
%     %mff file format is not sufficiently well documented yet for me to implement a more generalized fix.
%     if any(strcmp('ECG',EPdata.chanNames))
%         EPdata.chanTypes{find(strcmp('ECG',EPdata.chanNames))}='ECG';
%     end
%     if any(strcmp('s2_unknown259',EPdata.chanNames))
%         EPdata.chanTypes{find(strcmp('s2_unknown259',EPdata.chanNames))}='BAD';
%     end
% end

%input settings for function override the data file
if ~isempty(origRefChan)
    EPdata.reference.original=origRefChan;
end
if ~isempty(currRefChan)
    if ischar(currRefChan)
        EPdata.reference.type=currRefChan;
    else
        EPdata.reference.current=currRefChan;
        EPdata.reference.type='REG';
    end
end
if isempty(EPdata.reference.original) && strcmp(EPdata.ced,'internal')
    if (~isempty(strfind(EPdata.montage,'GSN')) || ~isempty(strfind(EPdata.montage,'Hydrocel')))
    %now that mff defaults to using internal ced info, is no longer providing default reference information.
        EEGchans=find(strcmp(EPdata.chanTypes,'EEG'));
        EPdata.reference.original=EEGchans(end);
        disp('Assuming the last EEG channel is reference.')
    end
end

%Detect bad epochs and mark them
if ~strcmp(EPdata.fileFormat,'ep_mat')
    for theSub=1:numRsubs
        for theWave=1:numWaves
            if ~any(any(any(EPdata.data(:,:,theWave,theSub,:,:,:),1),2),5) %if an epoch is entirely zeroes and NaNs, then mark epoch as bad.
                if any(strcmp(EPdata.dataType,{'single_trial','continuous'}))
                    EPdata.analysis.badTrials(theSub,theWave)=1;
                else
                    EPdata.avgNum(theSub,theWave)=-1;
                    EPdata.covNum(theSub,theWave)=-1;
                    EPdata.subNum(theSub,theWave)=-1;
                end
            end
        end
    end
end

if ~any(EPdata.data,'all')
    disp('Warning: All data values are zero.');
end

if all(isnan(EPdata.data),'all')
    disp('Warning: All data values are NaN.');
end

%Display summary of results

ep_tictoc;if EPtictoc.stop;return;end
if strcmp(silent,'off')
    disp(['The name of the experiment is: ' EPdata.ename]);
    if ~isempty(EPdata.baseline) && ~isempty(EPdata.Fs)
        disp(['The pre-stimulus period of the data is: ' num2str(EPdata.baseline * (1000/EPdata.Fs)) ' msec.']);
        if min(samples)>1
            disp('after dropping samples as instructed.');
        end
        if any(strcmp(EPdata.fileFormat,{'egi_egis','text'}))
            disp(['Note that the data file format (' EPdata.fileFormat ') is unable to specify the baseline period so it may be in error, in which case you will need to manually fix it using the Edit function.']);
        end
        if strcmp(EPdata.fileFormat,'egi_sbin')
            if ~(hdr.orig.header_array(14))==0 && (hdr.orig.header_array(15) > 1)
                disp(['Note that the data file format (' EPdata.fileFormat ') is unable to specify the baseline period so it may be in error, in which case you will need to manually fix it using the Edit function.']);
            end
        end
    end
end

%Final construction of output file
if ~isempty(EPdata.fileName)
    [pathstr, name, ext] = fileparts(EPdata.fileName);
    EPdata.dataName=name;
end

%Add secondary files to the data file if present

if ~isempty(specSuffix) %add in trial events information if desired
    [pathstr, name, fileSuffix] = fileparts(EPdata.fileName);
    specFileName=[pathstr filesep name specSuffix];
    if exist(specFileName,'file')
        outData = ep_readEventText(EPdata, specFileName);
        disp(['Adding trial events information from ' specFileName])
        if ~isempty(outData)
            EPdata=outData;
        end
    end
end

if ~isempty(subjectSpecSuffix) %add in subject specs information if desired
    [pathstr, name, fileSuffix] = fileparts(EPdata.fileName);
    specFileName=[pathstr filesep name subjectSpecSuffix];
    if exist(specFileName,'file')
       outData = ep_readSubjectSpecText(EPdata, specFileName);
       disp(['Adding subject specs information from ' specFileName])
        if ~isempty(outData)
            EPdata=outData;
        end
    end
end

%add in EyeLink eye-tracker information if present and not already added
if strcmp(EPdata.dataType,'continuous') && ~all(ismember({'pupil','x-eye','y-eye'},EPdata.chanNames))
    %  EyeLink edf files aren't actually edf file format.  To read them, one needs either their file reader or one developed by a Zurich lab.
    %  The Zurich one needs Wine to be set up for non-Windows computers.  The SR one is simpler to set up, but does need the SR Developers kit to be installed.
    %
    [pathstr, name, fileSuffix] = fileparts(EPdata.fileName);
    EyeLinkFileName=[pathstr filesep name '.edf'];
    if exist(EyeLinkFileName,'file')
        disp(['Adding eye-tracking from ' pathstr filesep name '.edf'])
        edfStruct=[];
        try
            edfStruct = edfmex(EyeLinkFileName); %SR file reader
        catch ME
            disp('Unable to read the file.  Did you install the SR Research EyeLink Developers Kit?')
            disp(ME.message);
        end
        EyeLinkMatchFileName=[pathstr filesep 'EyeLinkmatch.txt'];
        if exist(EyeLinkMatchFileName,'file')
            fid=fopen(EyeLinkMatchFileName);
            if fid ~=-1
                theData=textscan(fid, '%s%s','Delimiter','\t');
                matchTable=cell(length(theData{1}),length(theData));
                for iCol=1:length(theData)
                    for iRow=1:length(theData{iCol})
                        matchTable{iRow,iCol}=theData{iCol}{iRow};
                    end
                end
                fclose(fid);
            end
        else
            disp('No match table.')
            edfStruct=[];
        end
        if ~isempty(edfStruct)
            [EPdata msgLog]=ep_readEyelink(EPdata,edfStruct,matchTable);
        end
    end
end

%add in SMI eye-tracker information if desired and not already added
ep_tictoc;if EPtictoc.stop;return;end
if ~isempty(SMIsuffix) && ~all(ismember({'pupil','x-eye','y-eye'},EPdata.chanNames))
    if strcmp(EPdata.dataType,'continuous')
        SMItableData=[];
        [pathstr, name, fileSuffix] = fileparts(EPdata.fileName);
        SMIfileName=[pathstr filesep name SMIsuffix];
        SMImatchFileName=[pathstr filesep name '_SMImatch.txt'];
        if exist(SMIfileName,'file')
            %if already previously compiled word files then don't need to do so again.
        elseif exist([pathstr filesep name '_word01.txt'],'file') || exist([pathstr filesep name '_word01' SMIsuffix],'file')
            fileCounter=1;
            wordFid=0;
            theWordData=[];
            while wordFid ~= -1
                wordFile=[pathstr filesep name '_word' sprintf('%02d',fileCounter) '.txt'];
                wordFile2=[pathstr filesep name '_word' sprintf('%02d',fileCounter) SMIsuffix];
                if exist(wordFile,'file')
                    wordFid=fopen(wordFile);
                elseif exist(wordFile2,'file')
                    wordFid=fopen(wordFile2);
                else
                    wordFid=-1;
                end
                if wordFid ~= -1
                    if fileCounter == 1
                        theData=textscan(wordFid,'%s','Delimiter','\b'); %backspace as delimiter as don't want real delimiter
                        fclose(wordFid);
                    else
                        commentLine=1;
                        numComments=0;
                        while commentLine
                            tempLine=fgetl(wordFid);
                            if strcmp(tempLine(1:2),'##')
                                numComments=numComments+1;
                            else
                                commentLine=0;
                            end
                        end
                        
                        frewind(wordFid);
                        for i=1:numComments+1 %comments plus the header line
                            tempLine=fgetl(wordFid);
                        end
                        
                        theData=textscan(wordFid, '%s','Delimiter','\b');
                        fclose(wordFid);
                    end
                    theWordData=[theWordData;theData{1}];
                end
                fileCounter=fileCounter+1;
            end
            writetable(cell2table(theWordData),SMIfileName,'Delimiter','\t','WriteVariableNames',false,'QuoteStrings',false);
        else
            if exist(SMImatchFileName,'file')
                disp(['Error: No SMI files available.']);
            end
            SMIfileName=[];
        end
        if ~isempty(SMIfileName)
            disp(['Adding SMI eye-tracking data to: ' EPdata.fileName]);
            SMIfid=fopen(SMIfileName);
            if SMIfid==-1
                disp(['Error: Unable to open the SMI file: ' SMIfileName]);
            else
                if exist(SMImatchFileName,'file')
                    fid=fopen(SMImatchFileName);
                    if fid ~=-1
                        theData=textscan(fid, '%s%s','Delimiter','\t');
                        SMItableData=cell(length(theData{1}),length(theData));
                        for iCol=1:length(theData)
                            for iRow=1:length(theData{iCol})
                                SMItableData{iRow,iCol}=theData{iCol}{iRow};
                            end
                        end
                        fclose(fid);
                    end
                end
                if isempty(SMItableData) %if there is no match file, then present GUI
                    commentLine=1;
                    numComments=0;
                    while commentLine
                        tempLine=fgetl(SMIfid);
                        if strcmp(tempLine(1:2),'##')
                            numComments=numComments+1;
                        else
                            commentLine=0;
                        end
                    end
                    
                    delim='\t';
                    numcols=length(regexp(tempLine,delim))+1; %determine number of columns based on number of delimiters
                    
                    frewind(SMIfid);
                    for i=1:numComments
                        tempLine=fgetl(SMIfid);
                    end
                    
                    theData=textscan(SMIfid, repmat('%s',1,numcols),'Delimiter',delim, 'MultipleDelimsAsOne', 1);
                    fclose(SMIfid);
                    MSGlines=find(strcmp('MSG',theData{2}));
                    MSGtypes=unique(theData{4}(MSGlines));
                    EEGevents={EPdata.events{1}.value};
                    EEGtypes=unique(cellfun(@num2str,EEGevents(find(~cellfun(@isempty,EEGevents))),'UniformOutput',false'));
                    
                    for i=1:length(EEGtypes)
                        SMItableData{i,1}=EEGtypes{i};
                        SMItableData{i,2}='none';
                    end
                    
                    tableNames{1}='EEG';
                    tableNames{2}='SMI';
                    columnEditable =  [false true];
                    ColumnFormat{1}=[];
                    ColumnFormat{2}=[MSGtypes', 'none'];
                    
                    handles.read.SMI = figure('Name', 'Match SMI events', 'NumberTitle', 'off', 'Position',[scrsz(1) scrsz(4)-1100+scrsz(2) 200 500], 'MenuBar', 'none');
                    colormap jet;
                    
                    handles.read.SMItable = uitable('Data',SMItableData,'ColumnName',tableNames,'FontSize',FontSize,...
                        'ColumnEditable', columnEditable, 'ColumnFormat', ColumnFormat,'Position',[0 50 200 450]);
                    
                    handles.read.SMItableDone = uicontrol('Style', 'pushbutton', 'String', 'Done','FontSize',FontSize,...
                        'Position', [ 0 0 60 20], 'Callback', 'uiresume(gcbf)');
                    
                    uiwait(handles.read.SMI);
                    SMItableData=get(handles.read.SMItable,'Data');
                    close(handles.read.SMI);
                    drawnow;
                end
                
                if isempty(SMItableData)
                    disp(['Error: Addition of SMI eye-tracker information failed.']);
                else
                    [EPdataSMI logMsg]=ep_readSMI(EPdata,SMIfileName,SMItableData);
                    ep_tictoc;if EPtictoc.stop;return;end
                    if isempty(EPdataSMI)
                        disp(['Error: Addition of SMI eye-tracker information failed.']);
                    else
                        EPdata=EPdataSMI;
                    end
                end
            end
        end
    end
end

%add .asf video file
if strcmp(EPdata.dataType,'continuous') && isempty(EPdata.video)
    [pathstr, name, fileSuffix] = fileparts(EPdata.fileName);
    videoFileName=[pathstr filesep name '.asf'];
    if exist(videoFileName,'file')
        disp(['Reading video from ' pathstr filesep name '.asf'])
        outData=[];
        try
            outData = mmread(videoFileName);
        catch ME
            disp(['The attempt to load in the file ' videoFileName ' resulted in the error:' ME.identifier]);
            disp(ME.message);
        end
        if ~isempty(outData)
            for iVideo=1:length(outData)
                if outData(iVideo).nrFramesTotal > 0
                    break
                end
            end
            if outData(iVideo).nrFramesTotal > 0
                disp(['Adding video from ' pathstr filesep name '.asf'])
                EPdata.video(1).frames=outData(iVideo).frames;
                %rate and totalDuration fields are not reliable
                %times are assumed to be right side of video samples since they don't start with zero.
                EPdata.video(1).times=((outData(iVideo).times-median(diff(outData(iVideo).times)))*1000)-EPdata.baseline*(1000/EPdata.Fs);
            end
        end
    end
end

%add .sfp electrode coordinates file
if strcmp(EPdata.dataType,'continuous')
    [pathstr, name, fileSuffix] = fileparts(EPdata.fileName);
    fileList=dir(pathstr);
    matchList=cell(0);
    for iFile=1:length(fileList)
        theFile=fileList(iFile).name;
        if (length(theFile)>4) && (strcmp(theFile(end-3:end),'.sfp')) && (length(theFile(1:end-3)) <= length(name))
            if strcmp(theFile(1:end-4),name(1:length(theFile(1:end-4))))
                matchList{end+1}=theFile;
            end
        end
    end
    if length(matchList)>1
        disp('Multiple .sfp files match the file, so not using them.')
    elseif ~isempty(matchList)
        sfpFileName=[pathstr filesep matchList{1}];
        if exist(sfpFileName,'file')
            disp(['Reading electrode coordinates from ' pathstr filesep name '.sfp'])
            newElocs=ep_readlocsWrapper(sfpFileName,'filetype','sfp');
            if ~isempty(newElocs)
                disp(['Adding electrode coordinates from ' pathstr filesep name '.sfp'])
                FIDchans=find(strcmp({newElocs.type},'FID'));
                if isempty(EPdata.implicit)
                    EPdata.implicit=newElocs(FIDchans);
                else
                    for iFID=1:length(FIDchans)
                        if ~any(strcmp(EPdata.implicits.labels,newElocs(FIDchans(iFID)).labels))
                            EPdata.implicits(end+1)=newElocs(FIDchans(iFID));
                        else
                            EPdata.implicits(strcmp(EPdata.implicits.labels,newElocs(FIDchans(iFID)).labels))=newElocs(FIDchans(iFID));
                        end
                    end
                end
                newElocs(FIDchans)=[];
                for iChan=1:length(EPdata.eloc)
                    matchChan=find(strcmp({newElocs.labels},EPdata.eloc(iChan).labels));
                    if length(matchChan)>1
                        matchChan=matchChan(1);
                        disp('There are multiple channels with the same name.')
                    end
                    theCx=EPdata.eloc(iChan).cX;
                    theCy=EPdata.eloc(iChan).cY;
                    theCz=EPdata.eloc(iChan).cZ;
                    if ~isempty(matchChan)
                        EPdata.eloc(iChan)=newElocs(matchChan);
                    else
                        theName=EPdata.eloc(iChan).labels;
                        theType=EPdata.eloc(iChan).type;
                        tempEloc=ep_elocFormat('initialize');
                        tempEloc(1).labels=theName;
                        tempEloc(1).type=theType;
                        EPdata.eloc(iChan)=tempEloc;
                    end
                    EPdata.eloc(iChan).cX=theCx;
                    EPdata.eloc(iChan).cY=theCy;
                    EPdata.eloc(iChan).cZ=theCz;
                end
            else
                disp('Programming error: not reading sfp file.')
            end
        end
    end
end

%drop events of "trial" type as they are redundant with cell name information
for iSub=1:numRsubs
    for iWave=1:numRcells
        eventList=[];
        for event=1:length(EPdata.events{iSub,iWave})
            if ~strcmp(EPdata.events{iSub,iWave}(event).type,'trial')
                eventList=[eventList event];
            end
        end
        EPdata.events{iSub,iWave}=EPdata.events{iSub,iWave}(eventList);
    end
end

%mark epochs with a boundary event as being bad
if ~isempty(EPdata.events)
    if strcmp(EPdata.dataType,'continuous')
        if ~isempty(EPdata.events{1})
            epochSamples=[EPdata.events{1}(:).sample];
            boundaryList=ceil(epochSamples(find(strcmp('boundary',{EPdata.events{1}(:).value})))/EPdata.Fs);
            boundaryList(boundaryList>numEpochs)=numEpochs; %if a boundary event occurs in the extra points appended to the last epoch, assign it to the last epoch.
            EPdata.analysis.badTrials(1,unique(boundaryList))=1;
        end
    else
        for iSub=1:numRsubs
            for iWave=1:numRcells
                if ~isempty(EPdata.events{iSub,iWave})
                    if any(strcmp('boundary',{EPdata.events{iSub,iWave}(:).value}))
                        EPdata.analysis.badTrials(iSub,iWave)=1;
                    end
                end
            end
        end
    end
end

%Select subsets of data

ep_tictoc;if EPtictoc.stop;return;end
badCEDchans=find(strcmp('BAD',EPdata.chanTypes));
if ~isempty(badCEDchans)
    if isempty(channels)
        channels=[1:numChan];
    end
    channels=setdiff(channels,badCEDchans); %drop "channels" that the CED file indicates should be deleted entirely.
    for i=1:length(badCEDchans)
        msgString=['Dropping channels: '];
        for i=1:length(badCEDchans)
            msgString=[msgString EPdata.chanNames{badCEDchans(i)} ';'];
        end
    end
    disp(msgString);
end

EPdata2=ep_selectData(EPdata,{channels,samples,cells,subjects,[],[]});
if isempty(EPdata2)
    msg{1}='Defective file will not be loaded.';
    if strcmp(fileFormat,'ep_mat')
        msg{2}='File can be fixed by loading manually into Matlab and editing, as in load(''filename.ept'',''-mat''); and then save(''filename.ept'',''-mat'');';
    end
    [msg]=ep_errorMsg(msg);
    EPdata=[];
    return
else
    EPdata=EPdata2;
    EPdata2=[];
end

if ~isempty(cells) || ~isempty(subjects)
    if isscalar(EPdata.trialNames)
        trialStr='trial';
    else
        trialStr='trials';
    end
    if isscalar(EPdata.cellNames)
        cellStr='cell';
    else
        cellStr='cells';
    end
    if numSubs==1
        subStr='subject';
    else
        subStr='subjects';
    end
    if strcmp(silent,'off')
        if strcmp(EPdata.dataType,'single_trial')
            disp(['After selecting, there are ' num2str(numCells) ' ' cellStr ' with a total of ' num2str(length(EPdata.trialNames)) ' ' trialStr '.']);
        else
            disp(['After selecting, there are ' num2str(numCells) ' ' cellStr ' and ' num2str(numSubs) ' ' subStr '.']);
        end
    end
end

EPdataOut=EPdata;
outInfo.matlabDims=matlabDims;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [indexOut fieldOut]=readBVheaderField(eventHdr,indexIn)
%translate a field of BV header info
fieldType=str2double(eventHdr(indexIn).value(2:end));
indexOut=indexIn;
if isempty(fieldType) || isnan(fieldType) || (fieldType < 1) || (round(fieldType) ~= fieldType)
    fieldOut=NaN;
    return
end
switch fieldType
    case 1 %uint
        if (indexIn+1) > length(eventHdr)
            fieldOut=NaN;
            return
        end
        fieldOut=num2str(str2double(eventHdr(indexIn+1).value(2:end)));
        indexOut=indexOut+2;
    case 2 %slong
        if (indexIn+2) > length(eventHdr)
            fieldOut=NaN;
            return
        end
        fieldOut=num2str(str2double(eventHdr(indexIn+1).value(2:end))+bitshift(str2double(eventHdr(indexIn+2).value(2:end)),8));
        indexOut=indexOut+3;
    case 3 %text
        if (indexIn+1) > length(eventHdr)
            fieldOut=NaN;
            return
        end
        textLength=str2double(eventHdr(indexIn+1).value(2:end));
        if isnan(textLength) || isempty(textLength) || (textLength < 1) || (round(textLength) ~= textLength)
            fieldOut=NaN;
            return
        end
        theChar='';
        indexIn=indexIn+2;
        if indexIn > length(eventHdr)
            disp('Ran out of events.  Perhaps the recording was cut short?')
            fieldOut=NaN;
            return
        end
        for iChar=1:textLength
            [indexIn,charOut]=readBVheaderField(eventHdr,indexIn);
            if ~isnan(charOut)
                theChar(iChar)=char(str2double(charOut));
            else
                theChar(iChar)=' ';
            end
        end
        indexOut=indexIn;
        fieldOut=theChar;
    case 4 %zero unint
        fieldOut='0';
        indexOut=indexOut+1;
    case 5 %low byte slong zero
        if (indexIn+1) > length(eventHdr)
            fieldOut=NaN;
            return
        end
        fieldOut=num2str(bitshift(str2double(eventHdr(indexIn+1).value(2:end)),8));
        indexOut=indexOut+2;
    case 6 %hi byte slong zero
        if (indexIn+1) > length(eventHdr)
            fieldOut=NaN;
            return
        end
        fieldOut=num2str(str2double(eventHdr(indexIn+1).value(2:end)));
        indexOut=indexOut+2;
    case 7 %slong zero
        fieldOut='0';
        indexOut=indexOut+1;
    case 8 %empty text
        fieldOut='';
        indexOut=indexOut+1;
    otherwise
        fieldOut=NaN;
        return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [badFlag,badBitFlag,headerFlag,hdrStart,hdrEnd,subjectSpecNames,subjectSpecs,trialSpecNames]=readBVheader(eventHdr,iEvent,subjectSpecNames,subjectSpecs)
%reads BVheader
hdrStart=0;
hdrEnd=0;
badFlag=0;
trialSpecNames=cell(0);
headerFlag=0;
badBitFlag=0;
while iEvent < length(eventHdr)-7 %minimum size of a valid subject header is eight
    if strcmp(eventHdr(iEvent).value,'S104') && strcmp(eventHdr(iEvent+1).value,'S100') && strcmp(eventHdr(iEvent+2).value,'S114') && (hdrStart==0)
        headerFlag=1;
        hdrStart=iEvent;
        iEvent=iEvent+3;
        numFields=str2double(eventHdr(iEvent).value(2:end));
        iEvent=iEvent+1;
        if isnan(numFields) || isempty(numFields) || (numFields < 1) || (round(numFields) ~= numFields)
            disp('Header error.  Aborting effort to read it.');
            badFlag=1;
            continue
        end
        for iField=1:numFields
            if ~strcmp(eventHdr(iEvent).value(2:end),'  3')
                disp('Header error.  Aborting effort to read it.');
                continue
            end
            [iEvent,charOut]=readBVheaderField(eventHdr,iEvent);
            if ~ischar(charOut) || isempty(charOut)
                disp('Header error.  Aborting effort to read it.');
                badFlag=1;
                continue
            else
                subjectSpecNames{end+1,1}=charOut;
            end
            [iEvent,fieldOut]=readBVheaderField(eventHdr,iEvent);
            if isnan(fieldOut)
                disp('Header error.  Aborting effort to read it.');
                badFlag=1;
                continue
            else
                subjectSpecs{1,end+1}=fieldOut;
            end
        end
        if badFlag
            continue
        end
        numSpecs=str2double(eventHdr(iEvent).value(2:end));
        iEvent=iEvent+1;
        if isnan(numSpecs) || isempty(numSpecs) || (numSpecs < 1) || (round(numSpecs) ~= numSpecs)
            disp('Header error.  Aborting effort to read it.');
            badFlag=1;
            continue
        end
        for iSpec=1:numSpecs
            if ~strcmp(eventHdr(iEvent).value(2:end),'  3') %should be a 3 to indicate text coming up
                disp('Header error.  Aborting effort to read it.');
                badFlag=1;
                continue
            end
            [iEvent,charOut]=readBVheaderField(eventHdr,iEvent);
            if ~ischar(charOut) || isempty(charOut)
                disp('Header error.  Aborting effort to read it.');
                badFlag=1;
                continue
            else
                trialSpecNames{end+1,1}=charOut;
            end
        end
        if badFlag
            continue
        end
        if ~strcmp(eventHdr(iEvent).value,'S114') || ~strcmp(eventHdr(iEvent+1).value,'S100') || ~strcmp(eventHdr(iEvent+2).value,'S104')
            disp('Header error.  Aborting effort to read it.');
            badFlag=1;
            continue
        else
            hdrEnd=iEvent+2;
        end
        continue
    elseif strcmp(eventHdr(iEvent).value,'R  6') && strcmp(eventHdr(iEvent+1).value,'S  8') && strcmp(eventHdr(iEvent+2).value,'R  6') && (hdrStart==0)
        badBitFlag=1;
        hdrStart=iEvent;
        iEvent=iEvent+3;
        disp('Warning: It appears that the trigger bits setting on Pycorder/Recorder was set incorrectly (see tutorial).')
        continue
    else
        iEvent=iEvent+1;
    end
end
