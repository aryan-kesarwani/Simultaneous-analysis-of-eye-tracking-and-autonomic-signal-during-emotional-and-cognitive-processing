function [subtractedBlinks, subtractedSaccades, rawAverage, corrAverage, corRawData, logFile, blink, saccade] = EMCP(rawData, binData, parameters, logFileName)
%EMCP - [subtractedBlinks, subtractedSaccades, rawAverage, corrAverage, corRawData, logFile] = EMCP(rawData, binData, parameters)
% Eye movement correction program for correcting blink and saccade related artifacts.
%
%Input:
%  rawData: array of segmented microvolt data (trials, points, chans)
%  binData: array of accompanying bin numbers (trials).  Assumed to be members of a set of consecutive numbers (e.g., 1 to 5).
%  parameters:
%      .verbose                 0=no output to screen (for background jobs), 1=print progress to screen
%      .beckman                 0=not from Neuroscan lab, 1=from Neuroscan lab
%      .digitizingRate          Digitizing period (ms)
%      .VEOG_CHANNEL            position of VEOG in rawData */
%      .HEOG_CHANNEL            position of HEOG in rawData */
%      .FZ_CHANNEL              position of Fz channel in rawData.  Removed assumption that the first EEG channel is Fz.  JD
%      .eegChannels             Number of EEG channels.  It is assumed non-EEG channels are all positioned after the EEG channels.
%      .criteria                Microvolt number that is deemed out-of-range.
%  logFileName                  If a string, then the name of the file to print the log out to.  If empty, then will instead pass it on to the theLog output.
%
%Output:
%  subtractedBlinks: array of data subtracted from rawData as blink artifact (trials, points, chans)
%  subtractedSaccades: array of data subtracted from rawData as saccade artifact (trials, points, chans)
%  rawAverage: array of raw averaged data (bins, points, chans)
%  corrAverage: array of corrected averaged data (bins, points, chans)
%  corRawData: array of corrected segmented microvolt data (trials, points, chans)
%  logFile    : cell array with the log messages.
%  blink      : structure with blink propagation factors and so forth.
%  saccade    : structure with saccade propagation factors and so forth.
%
%History
%  by Joseph Dien (5/14/20)
%  jdien07@mac.com
%
% bugfix 8/27/20 JD
% Fixed translation bug causing some channels not to be corrected.
%
% bugfix 11/7/21 JD
% Fixed readTrial function not rejecting trials that contained bad channels unless it was the last channel.
%
% This is a Matlab port of the CPL C++ code (1998 version).  I have tried to cleave as closely as possible to the original code
% to avoid inadvertant changes to its functionality.  It is always possible that differences between C++ and Matlab may cause
% some changes in behavior.  I have dropped the file I/O code as no longer being needed.
% While modern amplifiers are unlikely to saturate, it still makes sense to specify an out-of-scale value so I have added
% a new .criteria parameter to replace that computed from the original .ADBits and .EOGSensitivity parameters.
% This also does not implement the scard system, making the assumption that the segments were classified during segmentation.
% There is therefore no trash bin as it is assumed all trials are already classified.
%
% The primary citation for this procedure is:
% Gratton, G., Coles, M. G. H., & Donchin, E. (1983).
% A new method for off-line removal of ocular artifact.
% Electroencephalography and Clinical Neurophysiology, 55, 468–484.
%
% %***************************************************************************/
% %                                                                          */
% % Originally written in Fortran by Gabriele Gratton : June 7, 1983.        */
% %      Revised : May 22, 1986                                              */
% %      Adapted to Microsoft Fortran : April 5, 1990. University of Illinois*/
% %      Dept. of Psychology, Cognitive Psychophysiology Lab.                */
% %                                                                          */
% % Current version written by James M. Turner : August 4, 1993.             */
% %      jturner@p300.cpl.uiuc.edu                                           */
% %      Neuroscience, Washington & Lee University                           */
% %      for the Dept. of Psychology, Cognitive Psychophysiology Lab.        */
% %      Champaign, Illinois.                                                */
% %                                                                          */
% % Thanks to Marten Scheffers for helping with comprehending and commenting */
% % the code.  Thanks to Brian Foote for lots of help with programming.      */
% % Thanks to Leun Otten for helping with the data formats and blink         */
% % criteria.  Thanks to Greg Miller and Fran Graham for explaining          */
% % critical parts of emcp.                                                  */
% %***************************************************************************/

if ~isempty(logFileName)
    logFile=fopen(logFileName,'w');
    if fid==-1
        terminate('ERROR',1,logFile,'Unable to open log file.');
        return
    end
else
    logFile=cell(0);
end

% msTen=[];                  % # samples (points) in 10 ms
% thirdOfBlink=[];           % # pts in 1/3 of blink window
% middleOfBlink=[];          % pt in middle of blink window
% initBlinkScan=[];          % 1st pt to scan for blinks
% endBlinkScan=[];           % last pt to scan for blinks
% blinkCriteria=[];          % criteria for blink
windVariance=2.0;            % window variance
% lengthOfWind=[];           % # points in blink template

rejectionCount=[0,0,0,0,0];  % why trials are rejected
% note: there are five possible tags for each trial which
% correspond to rejectionCount(0...4).
%       case 0:  trial accepted
%       case 1:  AD point out-of-scale in an EEG channel
%       case 2:  10 consecutive AD points out-of-scale
%       case 3:  more than 10 VEOG points out-of-scale
%       case 4:  unable to recover VEOG

% defined in <files.c> -- checks input string */
% reads in the +v xx.crd xx.dat xx.cal */
parameters.numTrials=size(rawData,1);
parameters.trialLength=size(rawData,2);
parameters.totalChannels=size(rawData,3);
parameters.numStorageBins=length(unique(binData));
parameters.numEEGChannels=parameters.eegChannels-2; %not including VEOG and HEOG channels

for iTrial=1:parameters.numTrials
    trialInfo(iTrial).status=0;
    trialInfo(iTrial).classification=0;
end
trialsInBin=zeros(max(binData),1);

markBlink=zeros(parameters.numTrials,parameters.trialLength);
% declare matrix of ints for blink analysis -- a 1 where there is a */
% blink present so initialize to all 0's later in checkBlink[]*/

trialBaseline=zeros(parameters.totalChannels,1);
% the baseline for a trial is stored here  */

binBaseline=zeros(parameters.totalChannels,parameters.numStorageBins);

blink.pointsInBin=zeros(parameters.trialLength,parameters.numStorageBins);
% allocate space used to determine correlation */
blink.sum=zeros(parameters.trialLength,parameters.totalChannels,parameters.numStorageBins);
blink.raw.variance=zeros(parameters.totalChannels,1);
blink.raw.veogCovar=zeros(parameters.totalChannels,1);
blink.raw.heogCovar=zeros(parameters.totalChannels,1);
% all array above used to calculate variance for blinks */

saccade.pointsInBin=zeros(parameters.trialLength,parameters.numStorageBins);
% allocate space used to determine correlation */
saccade.sum=zeros(parameters.trialLength,parameters.totalChannels,parameters.numStorageBins);
saccade.raw.variance=zeros(parameters.totalChannels,1);
saccade.raw.veogCovar=zeros(parameters.totalChannels,1);
saccade.raw.heogCovar=zeros(parameters.totalChannels,1);
% all array above used to calculate variance for saccades */

%*** Set up some preliminary variables for accept/reject ****/
%                                                           */
%************************************************************/
blinkCriteria = 14.0;
% determine criteria for blinks */
% value of 14 determined by gmiller@s.psych.uiuc.edu */
% in file emf11p.txt based on Gabrielle Gratton's original */
% code */
msTen = floor((10.0/parameters.digitizingRate)+0.5);
% # points = 10 ms */
thirdOfBlink = msTen*7;
if ~mod(thirdOfBlink,2)
    thirdOfBlink=thirdOfBlink+1;
end
% # points in 1/3 of blink template, make odd */
middleOfBlink = floor((thirdOfBlink+1)/2);
% point in middle of blink window */
lengthOfWind = thirdOfBlink*3.0;
% # points in blink template */
initBlinkScan = thirdOfBlink+middleOfBlink;
% first point to scan for blinks */
endBlinkScan = parameters.trialLength-initBlinkScan-thirdOfBlink+1;
% last point to scan for blink */

%*** Begin loop through trials -- accept/reject ****/
%                                                  */
% NOTES:  In this section of the program each      */
% trial is read into the arrays pointed to by      */
% rawData and rawIDs.  The VEOG channel (1) is     */
% checked to see if there are:                     */
%  1.  Ten consecutive points out-of-scale         */
%      rawData(t)(p)(VEOG_CHANNEL)>criteria        */
%  2.  If < 10 then program recovers epoch that    */
%      is out-of-scale.                            */
%                                                  */
%  The remaining channels are checked to see if:   */
%  1.  If 1 AD points is out-of-scale.             */
%  2.  If 10 consecutive equal data points are     */
%      in a channel.                               */
%  This rejection/acception loop is optimized for  */
%  speed -- not memory conservation (as memory is  */
%  in abundance these days).  Therefore these four */
%  criteria are checked for simultaneously,        */
%  which inherently makes the code more difficult  */
%  to read.                                        */
%  The criteria are checked in <files.c> using the */
%  routine checkData[] which in turn checks all    */
%  four criteria.                                  */
%  The values returned are used to assess the      */
%  the status of the trial.                        */
%***************************************************/

for curTrial=1:parameters.numTrials
    [trialInfo(curTrial).status, logFile]= readTrial(rawData, logFile, parameters);
    % readTrial from <files.c> second parameter '0' means to store trials */
    % in data arrays this time through. */
    % returns a # from 0...5 where 0 means trial accepted, 1..4 means */
    % trial rejected and 5 means EOF has been reached. */
    rejectionCount(trialInfo(curTrial).status+1)=rejectionCount(trialInfo(curTrial).status+1)+1;
    % update rejection count (0...4) as mentioned in prev comment */
    if (trialInfo(curTrial).status==0)
        % if trial accepted */
        [blinkReturn, trialInfo]=classified(curTrial, trialInfo, binData);
        if (blinkReturn==0)
            % if classification says to check blinks (==0) */
            % and single trial opt is not averaging only */
            % if returns 1 for classified then means not to check for */
            % blinks */
            markBlink=checkBlink(curTrial,windVariance,blinkCriteria,markBlink,lengthOfWind, initBlinkScan, endBlinkScan, thirdOfBlink, middleOfBlink, rawData, parameters);
            % check for blinks */
            [rawData, binBaseline, trialBaseline]= computeTrialBaseline('REMOVE', curTrial, rawData, binBaseline, trialBaseline, trialInfo, parameters);
            % compute && subtract average of each channel */
            % from that channel */
            [blink, saccade]= sumTrial(curTrial, blink, saccade, markBlink, rawData, trialInfo, parameters);
            % sums for correlation values */
            [rawData, binBaseline, trialBaseline]= computeTrialBaseline('RESTORE', curTrial, rawData, binBaseline, trialBaseline, trialInfo, parameters);
            % restores the mean of the channel to the */
            % channel */
        end
    else % if trial not accepted */

    end % end {if (trialInfo(curTrial).status==0)} */
    trialsInBin(trialInfo(curTrial).classification)=trialsInBin(trialInfo(curTrial).classification)+1;
    % increment counter for # trials in each bin */
    outputBins(curTrial, trialsInBin, parameters);
    % provide status of each bin */
end % end {while (readTrial(++curTrial,0,criteria))} */
%***    End of accept/rejection pass   *****/
%       and preliminary summming           */
%*******************************************/

%***   Status of accept/rejection pass  ****/
%                                          */
%*******************************************/
temp=sprintf('# of trials accepted: %d',rejectionCount(1));
logFile=display(logFile,temp,parameters);
temp=sprintf('# of trials with EEG out-of-scale: %d',rejectionCount(2));
logFile=display(logFile,temp,parameters);
temp=sprintf('# of trials EEG flat: %d',rejectionCount(3));
logFile=display(logFile,temp,parameters);
temp=sprintf('# of trials 10 EOG points out-of-scale: %d',rejectionCount(4));
logFile=display(logFile,temp,parameters);
temp=sprintf('# of trials EOG recovery failed: %d',rejectionCount(5));
logFile=display(logFile,temp,parameters);

%*******************************************/
%                                          */
%        Compute the average blink         */
%****    and saccade activity by ch    *****/
[blink.channelMean, blink.ptsInBase]=computeBaseline(blink.pointsInBin, blink.sum, parameters);
% compute average blink activity for each channel blink.channelMean(# ch) */

[saccade.channelMean, saccade.ptsInBase]=computeBaseline(saccade.pointsInBin, saccade.sum, parameters);
% compute average saccade activity for each channel saccade.channelMean (# ch) */

totalPoints=saccade.ptsInBase+blink.ptsInBase;
% saccade.ptsInBase + blink.ptsInBase = total # of points */

logFile=outputProportions(totalPoints, parameters, blink, saccade, logFile);
% prints out info on % blink data % saccade data */

%*******************************************/
%                                          */
%****    Compute the raw average blink *****/

[rawAverage, binBaseline]=computeBinAvg(blink.sum, saccade.sum, trialsInBin, binBaseline, parameters);
% computes the average wave and baseline in each bin by combining */
% blink and saccade components to form total sum of raw waveforms */
% and then divide by number of trials in each bin                 */

rawAverage=outputBinAvg(rawAverage, binBaseline, parameters);
% prepares for output by removing binBaseline

%*******************************************/
%                                          */
%***     Compute correction factors     ****/

blink.raw=computeVarForTrialWaves(blink.channelMean, blink.ptsInBase, blink.raw, parameters);
% compute the total variance and covariance for blink data */
% across trials                                            */

saccade.raw=computeVarForTrialWaves(saccade.channelMean, saccade.ptsInBase, saccade.raw, parameters);
% compute the total variance and covariance for saccade data */
% across trials                                              */

[blink.adjust]=computeVarForRawWaves(blink.sum, rawAverage, blink.ptsInBase, blink.channelMean, blink.pointsInBin, parameters);
% adjustments to variance due to the averages */
% #ifdef NEVER
% {
% 	int iii;
% for (iii=1 ; iii<=parameters.eegChannels ; iii++)
% {
%     fprintf(stderr,'%d veogc %12.9f heogc %12.9f\n',
% 	iii,
% 	blink.adjust.veogCovar(iii),
% 	blink.adjust.heogCovar(iii)) ; % Ehehehehehehehe */
% }
% }
% #endif

[saccade.adjust]=computeVarForRawWaves(saccade.sum, rawAverage, saccade.ptsInBase, saccade.channelMean, saccade.pointsInBin, parameters);
% adjustments to variance due to the averages */

[blink.residual, blink.stanDev]=computeResidualVariance(blink.raw, blink.adjust, parameters);
% #ifdef NEVER
% {
% 	int jjj ;
% for (jjj=1 ; jjj<=parameters.eegChannels ; jjj++)
% {
% 	fprintf(stderr,'standev %d %12.9f \n',
% 		jjj,
% 		blink.stanDev(jjj)) ;
% }}
% #endif
[saccade.residual, saccade.stanDev]=computeResidualVariance(saccade.raw, saccade.adjust, parameters);

% computes the difference between the results */
% of the previous two functions */

outputVarianceTables(blink, saccade, logFile, parameters);
% outputs variance/covariance tables to file */

[blink.veogCorr, blink.heogCorr]=computeCorrWithEOG(blink.stanDev, blink.residual,  parameters);
% #ifdef NEVER
% {
% 	int iii;
% fprintf(stderr,'after corr with eog\n') ;
% for (iii=1 ; iii<=parameters.eegChannels ; iii++)
% {
%     fprintf(stderr,'%d veogc %12.9f heogc %12.9f\n',
% 	iii,
% 	blink.veogCorr(iii),
% 	blink.heogCorr(iii)) ; % Ehehehehehehehe */
% }
% }
% #endif
% computes corr between residual variance and EOG */

[saccade.veogCorr, saccade.heogCorr]=computeCorrWithEOG(saccade.stanDev, saccade.residual,  parameters);
% computes corr between residual variance and EOG */

[blink.heogFactor, blink.veogFactor]=computePropogation(blink.veogCorr, blink.heogCorr, blink.stanDev, parameters);
% #ifdef NEVER
% {
% 	int iii;
% fprintf(stderr,'after prop\n') ;
% for (iii=1 ; iii<=parameters.eegChannels ; iii++)
% {
%     fprintf(stderr,'%d veogc %12.9f heogc %12.9f\n',
% 	iii,
% 	blink.veogCorr(iii),
% 	blink.heogCorr(iii)) ; % Ehehehehehehehe */
% 	blink.veogCorr(iii)=iii/10.0 ;
% }
% }
% #endif
% computes the propogation factor for correction */

[saccade.heogFactor, saccade.veogFactor]=computePropogation(saccade.veogCorr, saccade.heogCorr, saccade.stanDev, parameters);
% computes the propogation factor for correction */

% #ifdef NEVER
% { int iii;
% fprintf(stderr,'after prop\n') ;
% for (iii=1 ; iii<=parameters.eegChannels ; iii++)
% {
%     fprintf(stderr,'%d veogc %12.9f heogc %12.9f\n',
% 	iii,
% 	blink.veogCorr(iii),
% 	blink.heogCorr(iii)) ; % Ehehehehehehehe */
% 	blink.veogCorr(iii)=iii/10.0 ;
% }
% }
% #endif

outputPropogation(blink, saccade, logFile, parameters);
% #ifdef NEVER
% { int  iii;
% fprintf(stderr,'after prop\n') ;
% for (iii=1 ; iii<=parameters.eegChannels ; iii++)
% {
%     fprintf(stderr,'%d veogc %12.9f heogc %12.9f\n',
% 	iii,
% 	blink.veogCorr(iii),
% 	blink.heogCorr(iii)) ; % Ehehehehehehehe */
% 	blink.veogCorr(iii)=iii/10.0 ;
% }
% }
% #endif
% outputs the progogation factors to file */

[corrAverage, ~, ~]=computeCorrectAvgs(blink, saccade, trialsInBin, rawAverage, binBaseline, parameters);
% corrects the averages for eye movement */

[corRawData, subtractedBlinks, subtractedSaccades]=computeCorrectedSingleTrials(rawData, blink, saccade, binBaseline, trialBaseline, windVariance, lengthOfWind, initBlinkScan, endBlinkScan, thirdOfBlink, middleOfBlink, trialInfo, blinkCriteria, markBlink, parameters);
% computes the corrected single trials */

%*******************************************/
%                                          */
%***      Single trial corrections      ****/


%***  Gracefully close out the program  ****/
%                                          */
%*******************************************/

terminate('NO_ERROR',parameters, logFile, '**** EMCP completed. ****');
return
% terminate program <emcp.c> */

%***********************************************************************/
%  checkBlink  : subroutine checks to see if there is a blink        */
%***********************************************************************/
function markBlink=checkBlink(trial,windVariance,blinkCriteria,markBlink,lengthOfWind, initBlinkScan, endBlinkScan, thirdOfBlink, middleOfBlink, rawData, parameters)

% extern  int     msTen,                  % # samples (points) in 10 ms */
% thirdOfBlink,           % # pts in 1/3 of blink window */
% middleOfBlink,          % pt in middle of blink window */
% initBlinkScan,          % 1st pt to scan for blinks */
% endBlinkScan;           % last pt to scan for blinks */
% extern double   blinkCriteria,          % criteria for blink */
% windVariance,           % window variance */
% lengthOfWind;          % # points in blink template */
%
% int     refPoint,                       % center of blink window */
% checkPoint,                     % point in blink window */
% beginWind,                      % beginning of blink window */
% endWind,                        % end of blink window */
% posInTemplate,                  % position in window */
% startMark,                      % start marking blink here */
% endMark,                        % stop marking as blink here */
% point;                          % just a counter */
% double  covar,                          % covariance in window */
% slope;                          % slope in window */

for point=1:parameters.trialLength
    markBlink(trial,point)=2;
end
% clear blink record */

for refPoint=initBlinkScan:endBlinkScan
    % define scan template around refPoint */
    beginWind=refPoint-(initBlinkScan-1);
    % beginning of scan template */
    endWind=refPoint+(initBlinkScan-1);
    % end of scan template */
    covar=0;
    % init covariance */
    
    for checkPoint=beginWind:endWind
        posInTemplate=checkPoint-beginWind+1;
        % current location in window */
        if(posInTemplate<=(thirdOfBlink+1)||posInTemplate>(2*thirdOfBlink+1))
            covar=covar-rawData(trial,checkPoint,parameters.VEOG_CHANNEL);
        else
            covar=covar+2.0*rawData(trial,checkPoint,parameters.VEOG_CHANNEL);
        end
    end % end for (checkPoint) */
    covar=covar/lengthOfWind;
    % determine covariance */
    slope=covar/(windVariance*windVariance);
    % determine slope */
    
    if (abs(slope)>blinkCriteria)
        % if the slope exceeds the maximum rate of range mark as blink */
        startMark = refPoint-middleOfBlink+1;
        % this range is a blink */
        endMark = refPoint+middleOfBlink-1;
        % stop marking here */
        for checkPoint=startMark:endMark
            markBlink(trial,checkPoint)=1;
            % mark blink at checkPoint */
        end
    end % end if slope > blinkCriteria */
end % end for (refPoint) */


%***********************************************************************/
%  terminate[]  : write log about error                                */
%***********************************************************************/
function terminate(error, parameters, logFile, message)
switch error
    case 'NO_ERROR'
        logFile=display(logFile, message, parameters);
    case 'ERROR'
        logFile=display(logFile, message, parameters);
        logFile=display(logFile, '****Terminating EMCP****', parameters);
end
if ~iscell(logFile)
    fclose(logFile);
end
return
% close log file <files.c> */

%************************************************************************/
%  display  : function that prints text to screen and file as designated*/
%************************************************************************/
function logFile=display(logFile,text,parameters)
if parameters.verbose
    fprintf([text '\n']);
end
if ~iscell(logFile)
    err=fprintf(logFile,text);
    if err==-1
        disp('Error printing to log file.');
        terminate('ERROR', parameters.verbose, logFile, 'Error printing to log file.')
    end
else
    logFile{end+1}=sprintf(text);
end

%**********************************************************************/
%  classified[]  : classifies the trials based on their ids           */
%**********************************************************************/
% return 0: check blinks */
% return 1: don't check for blinks */
function [blinkReturn, trialInfo] = classified(curTrial, trialInfo, binData)

trialInfo(curTrial).classification=binData(curTrial);

blinkReturn=0;

%***************************************************************************/
%                                                                          */
% FILES.c: Handles input and output of data.  Checks for artifacts in data */
%          during the input phase.                                         */
%                                                                          */
% Current version written by James M. Turner : August 4, 1993.             */
%      jturner@p300.cpl.uiuc.edu                                           */
%      Neuroscience, Washington & Lee University                           */
%      for the Dept. of Psychology, Cognitive Psychophysiology Lab.        */
%      Champaign, Illinois.                                                */
%                                                                          */
% NOTE:  Portions of this code are direct copies of code from previous     */
%        version.                                                          */
%***************************************************************************/

% #include <stdio.h>
% #include <math.h>
% #include <stdlib.h>
% #include <string.h>
% #include 'files.h'
% #include 'allocate.h'
% #include 'selectc.h'
% #include 'emcp.h'
% 
% %** Variables ***/
% extern FILE    *parameterFile,                % file containing params */
%                *dataFile,                     % raw data file */
%                *calibrationFile,              % calibration values for AD conv */
%                *correctedSngTrlFile,          % corrected single trial file */
%                *logFile;                      % log of activity */
% 
% 
% extern struct  Parameters parameters;         % parameters for input file */
% 
% extern float   ***rawData,                    % raw data array */
%                ***rawIDs;                     % raw IDs array     */
% 
% extern int     *trialsInBin,                  % number of trials in each bin  */
%                 **markBlink;                  % blink present at [point] */
%                                               % in curTrial   */
% extern double  *curIDs,                       % current IDs */
%                *trialBaseline,                % buffer for baseline of trial */
%                **binBaseline;                 % buffer for baseline of bins */
% 
% extern struct  TrialInfo *trialInfo;          % trial info array */
%             
% extern int     verbose;                       % display info on screen, default to yes */
% 
% %** Structs for correlation determinations ***/
% extern struct  Corr    blink,                 % used to organize arrays  */
%                        saccade;               % associated with calculation */
%                                               % of correlations among */
%                                               % residual.heogCovar & HEOG channels */
%                                               % and Fz */
% 
% %** Variables for blink detection  ***/
% extern int     msTen,                         % # samples (points) in 10 ms */
%                thirdOfBlink,                  % # pts in 1/3 of blink window */
%                initBlinkScan,                 % 1st pt to scan for blinks */
%                endBlinkScan,                  % last pt to scan for blinks */
%                middleOfBlink;                 % pt in middle of blink window */
% extern double  blinkCriteria,                 % criteria for blink */
%                windVariance=2.0,              % window variance */
%                lengthOfWind;                  % # points in blink template */
% 
% %** Variables for EEG rejection ***/
% extern double  criteria;                      % criteria for out-of-scale */
% 
% %** Variables for artifact removal ***/
% extern double  totalPoints;                   % total # data points */
% extern float   ***rawAverage;                 % average wave in each bin */
% extern double  *calibration,                  % calibration vals for AD converters */
%                **meanIDs;                     % average of ID values */
% 
% int     numVEOGoffScale=0,  % num consecutive VEOG points of scale */
%         firstVEOGoffScale=0;% used for recovery of consecutive off scale VEOG */
%                             % defined globally because can't be reinitialized */
%                             % at each call of checkData(); */
% 
% char *temp;

%*********************************************************************/
%  readTrial() :function that reads data from data file by one trial */
%*********************************************************************/    
% if return 0 then good trial */
% if return 1 then bad trial */
function [accept, logFile]= readTrial(rawData, logFile, parameters)

firstVEOGoffScale=0;
numVEOGoffScale=0;
accept=0;

for j=1:parameters.totalChannels
    for k=1:parameters.trialLength %changed to trialLength from numIDs JD
        if accept~=0
            [accept, firstVEOGoffScale, numVEOGoffScale, logFile]=checkData(rawData, parameters.numTrials, k, j, firstVEOGoffScale, numVEOGoffScale, logFile, parameters);
        end
    end
end

%*********************************************************************/
%  checkData()  : checks to see if incoming data in channels is okay */
%*********************************************************************/
function [accept, firstVEOGoffScale, numVEOGoffScale, logFile]=checkData(rawData, curTrial, point, channel, firstVEOGoffScale, numVEOGoffScale, logFile, parameters)

%      *    return 0 means data point okay                          *
%      *    return 1 means EEG data point out-of-range              *
%      *    return 2 means ten consecutive equal EEG data points    *
%      *    return 3 means 10 VEOG points out-of-range              *
%      *    return 4 means VEOG points unrecoverable                */

i=9;                        % counters */
if (0 == parameters.criteria)
    accept=0;
    return;                       % if told not to reject trials
    % then exit nicely */
else
    if (parameters.VEOG_CHANNEL == channel)
        % if VEOG channel */
        if ((point<=2)||(point>=(parameters.trialLength-1)))
            % if first two or last */
            if (abs(rawData(curTrial,point,parameters.VEOG_CHANNEL)) >= parameters.criteria)
                % two points in VEOG >=parameters.criteria then exit */
                accept=4;
                return;
            end
            % VEOG unrecoverable */
            if (numVEOGoffScale>0)
                [accept, numVEOGoffScale, firstVEOGoffScale, logFile]=recoverVEOG(numVEOGoffScale, firstVEOGoffScale, curTrial, rawData, logFile, parameters);
                return ;
                % if at 2 points then time to correct   */
                % if there are off-scale points */
            end % this was forgotten in first version */
            accept=0;
            return;
            % then point is okay */
        else %to{if((point<=2)||(point>=(parameters.trialLength-1)))} */
            %********************************************************/
            %  NOTE:  This section attempts to recover VEOG points  */
            %         if a series of continuous equal points        */
            %         exists that is > 1 but < 10 points long.      */
            %********************************************************/
            if ((abs(rawData(curTrial,point,parameters.VEOG_CHANNEL)) >= parameters.criteria) && (numVEOGoffScale<10))
                % if VEOG points>=parameters.criteria save them */
                numVEOGoffScale=numVEOGoffScale+1;
                if (numVEOGoffScale==1)
                    firstVEOGoffScale=point;
                    % if first point off scale then mark  */
                    % it for reference */
                    accept=0;
                    return;
                end
                % continue reading points */
            elseif (numVEOGoffScale==10)
                % stop if > 10 VEOG points off scale */
                accept=3;
                return;
            elseif ((abs(rawData(curTrial,point,parameters.VEOG_CHANNEL))<parameters.criteria)&&(numVEOGoffScale>0))
                % if VEOG points>=parameters.criteria && if > 0 points */
                % off scale try to  */
                % recover them */
                [accept, numVEOGoffScale, firstVEOGoffScale, logFile]=recoverVEOG(numVEOGoffScale, firstVEOGoffScale, curTrial, rawData, logFile, parameters);
                return
                % return whether recovery successful */
            end
        end % end if {if ((point<=2)||(point>=(parameters.trialLength-1)))} */
        % KS (2/9/97): Modified this if statement so that checkData wouldn't apply the rejection criteria to non-EEG channels, such as an event code channel. */
        %else if (parameters.VEOG_CHANNEL ~= channel) */
    elseif ((channel ~= parameters.VEOG_CHANNEL) && (channel <= 2+parameters.numEEGChannels))
        % else to {if (parameters.VEOG_CHANNEL == channel)} */
        if (abs(rawData(curTrial,point,channel))>=parameters.criteria)
            accept=1;
            return;
        end
        % if data point out-of-range && is not in VEOG */
        if (point>9 && channel<=parameters.eegChannels)
            % if read in 10 points and EEG data channel */
            if (parameters.HEOG_CHANNEL~=channel)
                % can't be VEOG and don't want HEOG because they are frequently */
                % flat when there is not eye movement */
                while (rawData(curTrial,point-i,channel)==rawData(curTrial,point,channel))
                    % check backwards for 10 consecutive equal data points */
                    i=i-1;
                    if (i==0)
                        break;
                        % if i==0 stop checking */
                    end
                end
                if (i==0)
                    accept=2;
                    return;
                end
            end
            % if counter reaches 0 then must be 9 equal pts -- return 2 */
        end
    end % end {if (channel ~= parameters.VEOG_CHANNEL etc)} */
end % end {if (0 == parameters.criteria)} */
accept=0;
return;
% rawData[curTrial][point][channel] is okay -- return 0 */


%*********************************************************************/
%  recoverVEOG()  : attempts to recover out-of-scale VEOG points     */
%*********************************************************************/
function [accept, number, startPoint, logFile]=recoverVEOG(number,startPoint, curTrial, rawData, logFile, parameters)

% point=[];              % counter last point */
% status=[];             % used to assess results of cov f(x) */
% scalingFactor=[];      % scaling factor based on length of out-of-scale epoch */
% deltaUp=0;          % rate of change before epoch */
% deltaDown=0;        % rate of change after epoch */
% incrementUp=0;      % increment up */
% incrementDown=0;    % increment down */
% constantUp=[];         % * constant based on dist from start */
% constantDown=[];       % * constant based on dist from end */
% rSq=[];                % squared correlation value */
% slope=[];              % slope value based on covariance/variance of VEOG to Fz */
% changeUp=[];          % storage for est. up component of new points */
% changeDown=[];        % storage for est. down component of new points */
% temp='';          % buffer for output string */

endPoint=number+startPoint;
changeUp=zeros(endPoint,1); %points 1:startPoint-1 simply not used.
changeDown=zeros(endPoint,1); %points 1:startPoint-1 simply not used.
scalingFactor=sqrt((number+1)/2);
% determine scaling factor */
deltaUp=rawData(curTrial,startPoint-1,parameters.VEOG_CHANNEL)-rawData(curTrial,endPoint-2,parameters.VEOG_CHANNEL);
% difference between preceding 2 points */
deltaDown=rawData(curTrial,endPoint+1,parameters.VEOG_CHANNEL)-rawData(curTrial,endPoint+2,parameters.VEOG_CHANNEL);
% difference between following 2 points */

for point=startPoint:endPoint
    constantUp=1+sin((point+1-startPoint)*pi/number);
    incrementUp=deltaUp*constantUp*scalingFactor;
    changeUp(point)=rawData(curTrial,startPoint-1,parameters.VEOG_CHANNEL)+incrementUp;
    
    constantDown=1+sin((endPoint+1-point)*pi/number);
    incrementDown=deltaDown*constantDown*scalingFactor;
    changeDown(point)=rawData(curTrial,endPoint+1,parameters.VEOG_CHANNEL)+incrementDown;
    
    rawData(curTrial,point,parameters.VEOG_CHANNEL)=(changeUp(point)+changeDown(point))/2;
    if rawData(curTrial,point,parameters.VEOG_CHANNEL) < parameters.criteria
        % take average for estimatedVEOG and substitute in to VEOG channel */
        number=0;
        startPoint=0;
        % these variables */
        accept=4;
        return;
        % return rejection */
    end
end % end {for(point=startPoint;point<=end;point++)} */
[status, rSq, slope] = covariance(curTrial,parameters.VEOG_CHANNEL,parameters.FZ_CHANNEL,startPoint,endPoint);
% get rSq and slope based on VEOG and Fz correlations */
if (1==status)
    accept=4;
    return;
end
% if xVariance*yVariance<0 unrecoverable */
slope=slope/parameters.EOGSensitivity;
% convert slope */
if ((rSq<.81)||(abs(slope)>.40)||(abs(slope)<.10))
    accept=4;
    return;
end
% recovered successfully? */
temp=sprintf('Recovered VEOG-- trial %d, point %d to %d.',curTrial,startPoint,endPoint);
logFile=display(logFile,temp,parameters);
temp=sprintf('R-square= %6.4f, slope= %6.4f',rSq,slope);
logFile=display(logFile,temp,parameters);
for point=startPoint:endPoint
    % output info to log */
    fprintf(logFile,' %7.1f\n',rawData(curTrial,point,parameters.VEOG_CHANNEL));
    if (parameters.verbose)
        fprintf(' %7.1f\n',rawData(curTrial,point,parameters.VEOG_CHANNEL));
    end
end
number = 0;                         % recovery attempt over */
startPoint = 0;                          % therefore need to rest #s */
accept=0;
return;                           % return that it worked */


%********************************************************************/
%  covariance()  : computers covariance between two channels        */
%********************************************************************/
function [status, rSq, slope] = covariance(curTrial, x, y, startPoint, endPoint)

status=0;       % status */

sumX=0;         % sum of x */
sumY=0;         % sum of y */
sumXX=0;        % sum of x*x */
sumYY=0;        % sum of y*y */
sumXY=0;        % sum of x*y */

for point=startPoint:endPoint
    % sum points up */
    sumX=sumX+rawData(curTrial,point,x);
    sumY=sumY+rawData(curTrial,point,y);
    sumXX=sumXX+rawData(curTrial,point,x)*rawData(curTrial,point,x);
    sumYY=sumYY+rawData(curTrial,point,y)*rawData(curTrial,point,y);
    sumXY=sumXY+rawData(curTrial,point,x)*rawData(curTrial,point,y);
end

numPoints=startPoint-endPoint+1;    % number of points for average */
xMean=sumX/numPoints;               % calculate x average */
yMean=sumY/numPoints;               % calculate y average */

xVariance=(sumXX/numPoints)-(xMean*xMean);% calculate xVariance */
yVariance=(sumYY/numPoints)-(yMean*yMean);% calculate yVariance */
xyCovariance=(sumXY/numPoints)-(xMean*yMean);% calc xyCovariance */

if ((xVariance*yVariance)<=0)
    status=1;
end
% conciliation to previous version */
rSq=(xyCovariance*xyCovariance)/(xVariance*yVariance);
% calculate r squared */
slope=xyCovariance/xVariance;
% calculate slope */

return;
  
%***********************************************************************/
%  outputPropogation()  : outputs the progogation factors              */
%***********************************************************************/
function outputPropogation(blink, saccade, logFile, parameters)
    
logFile=display(logFile, sprintf('\n\n\n>>>>Correction factor for blinks \n'), parameters);

% Write table to correlation file */
logFile=display(logFile, sprintf('\n                  '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf( '%6d ',channel), parameters);
end

logFile=display(logFile, sprintf('\n>>Standard Dev:  '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf( '%6.0f ',blink.stanDev(channel)), parameters);
end

logFile=display(logFile, sprintf( '\n>>Beta weights:\n      Vertical_bb    '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf( '%6.2f ',blink.veogCorr(channel)), parameters);
end
logFile=display(logFile, sprintf( '\n      Horizontal_bb  '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf( '%6.2f ',blink.heogCorr(channel)), parameters);
end

logFile=display(logFile, sprintf('\n>>Propagation factors:\n      Vertical_bp    '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf( '%6.2f ',blink.veogFactor(channel)), parameters);
end
logFile=display(logFile, sprintf( '\n      Horizontal_bp  '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf( '%6.2f ',blink.heogFactor(channel)), parameters);
end

logFile=display(logFile, sprintf('\n\n>>>>Correction factor for saccades \n'), parameters);

% Write table to correlation file */
logFile=display(logFile, sprintf('\n                  '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf( '%6d ',channel), parameters);
end

logFile=display(logFile, sprintf('\n>>Standard Dev:  '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf( '%6.0f ',saccade.stanDev(channel)), parameters);
end

logFile=display(logFile, sprintf( '\n>>Beta weights:\n      Vertical_sb    '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf( '%6.2f ',saccade.veogCorr(channel)), parameters);
end
logFile=display(logFile, sprintf( '\n      Horizontal_sb  '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf( '%6.2f ',saccade.heogCorr(channel)), parameters);
end

logFile=display(logFile, sprintf('\n>>Propagation factors:\n      Vertical_sp    '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf( '%6.2f ',saccade.veogFactor(channel)), parameters);
end
logFile=display(logFile, sprintf( '\n      Horizontal_sp  '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf( '%6.2f ',saccade.heogFactor(channel)), parameters);
end

%***********************************************************************/
%  outputBinAvg()  : converts to microv, removes binBaseline, to file  */
%***********************************************************************/
function rawAverage=outputBinAvg(rawAverage, binBaseline, parameters)

%** Prepare for output ***/
for point=1:parameters.trialLength
    for bin=1:parameters.numStorageBins
        for channel=1:parameters.totalChannels
            rawAverage(point,channel,bin)=rawAverage(point,channel,bin)+binBaseline(channel,bin);
            % take out the baseline for bin */
        end % end channel */
    end % end bin */
end % end point */

%***********************************************************************/
%  outputProportions()  :  output info on % blinks % saccades in data  */
%                                                                      */
%  Dependent on:                                                       */
%     totalPoints, .ptsInBase                                          */
%***********************************************************************/
function logFile=outputProportions(totalPoints, parameters, blink, saccade, logFile)
%     double  blinkRatio,                     % % blink data points */
%             saccadeRatio;                   % % saccade data points */

if (totalPoints)
    blinkRatio=blink.ptsInBase/totalPoints;
    % # blink data points over total points */
    saccadeRatio=saccade.ptsInBase/totalPoints;
    % # saccade data points over total points */
end
% write proportions to correlation file */
logFile=display(logFile, sprintf('\n    Sample size for regression computation: \n'), parameters);
logFile=display(logFile, sprintf('\nREGRESSION     # OF POINTS  PROPORTION\n'), parameters);
logFile=display(logFile, sprintf('\nSACCADES       %10.0f  %10.5f', saccade.ptsInBase, saccadeRatio), parameters);
logFile=display(logFile, sprintf('\nBLINKS         %10.0f  %10.5f', blink.ptsInBase, blinkRatio), parameters);
if (parameters.verbose)
    fprintf('\n>>    Sample size for regression computation: \n');
    fprintf('\n>> REGRESSION     # OF POINTS  PROPORTION\n');
    fprintf('\n>> SACCADES       %10.0f  %10.5f', saccade.ptsInBase, saccadeRatio);
    fprintf('\n>> BLINKS         %10.0f  %10.5f', blink.ptsInBase, blinkRatio);
end

%***********************************************************************/
%  outputBins()  : output status                                       */
%***********************************************************************/
function outputBins(trial, trialsInBin, parameters)
    
    % prints out some standard info that looks likes this:
    
    %+ 1 #/bin 0 0 0 0 0 0 0 0 0 ... num of bins
    
    if (parameters.verbose)
        % if we are display to screen then do it */
        fprintf('+ trl# %4d  #/bin:',trial);
        for bin=1:parameters.numStorageBins
            fprintf('%3d ',trialsInBin(bin));
        end
        fprintf('\n');
    end

%***********************************************************************/
%  outputVariance()  : writes out variance-covariance tables           */
%***********************************************************************/
function outputVarianceTables(blink, saccade, logFile, parameters)

%**  Total variance ***/
logFile=display(logFile, sprintf( '\n\n>>Total variance-covariance table '), parameters);

xmax = 0.0;                 % determine scaling factor */
for channel=1:parameters.eegChannels
    if blink.raw.variance(channel) > xmax
        xmax = blink.raw.variance(channel);
    end
    if blink.raw.veogCovar(channel) > xmax
        xmax = blink.raw.veogCovar(channel);
    end
    if blink.raw.heogCovar(channel) > xmax
        xmax = blink.raw.heogCovar(channel);
    end
    if saccade.raw.variance(channel) >xmax
        xmax = saccade.raw.variance(channel);
    end
    if saccade.raw.veogCovar(channel) >xmax
        xmax = saccade.raw.veogCovar(channel);
    end
    if saccade.raw.heogCovar(channel) >xmax
        xmax = saccade.raw.heogCovar(channel);
    end
end

scal = 1.0;

if (xmax >= 1000000.0)
    scal = 10000.0;
elseif (xmax >= 100000.0)
    scal = 1000.0;
elseif (xmax >= 10000.0)
    scal = 100.0;
elseif (xmax >= 1000.0)
    scal = 10.0;
end

logFile=display(logFile, sprintf('\n>>Scaling factor = %6.0f\n',scal), parameters);

% blinks */
logFile=display(logFile, sprintf('\n Blinks\n                 '), parameters);
for channel=1:parameters.eegChannels    % write channel ind */
    logFile=display(logFile, sprintf('%6d ',channel), parameters);
end

logFile=display(logFile, sprintf('\n Variance        '), parameters);
for channel=1:parameters.eegChannels    % write variance for blinks */
    logFile=display(logFile, sprintf('%6.0f ', blink.raw.variance(channel)/scal), parameters);
end

logFile=display(logFile, sprintf('\n covar chn-Veog  '), parameters);
for channel=1:parameters.eegChannels    % write covariance for blinks */
    logFile=display(logFile, sprintf('%6.0f ', blink.raw.veogCovar(channel)/scal), parameters);
end
logFile=display(logFile, sprintf('\n covar chn-Heog  '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf('%6.0f ', blink.raw.heogCovar(channel)/scal), parameters);
end

% Saccades */
logFile=display(logFile, sprintf('\n\n Saccades\n                 '), parameters);
for channel=1:parameters.eegChannels    % write channel index */
    logFile=display(logFile, sprintf('%6d ',channel), parameters);
end

logFile=display(logFile, sprintf('\n Variance        '), parameters);
for channel=1:parameters.eegChannels    % write variance for saccades */
    logFile=display(logFile, sprintf('%6.0f ', saccade.raw.variance(channel)/scal), parameters);
end

logFile=display(logFile, sprintf('\n covar chn-Veog  '), parameters);
for channel=1:parameters.eegChannels    % write covariance for saccades */
    logFile=display(logFile, sprintf('%6.0f ', saccade.raw.veogCovar(channel)/scal), parameters);
end
logFile=display(logFile, sprintf('\n covar chn-Heog  '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf('%6.0f ', saccade.raw.heogCovar(channel)/scal), parameters);
end

%****************************/
%  Adjustments to variance  */
logFile=display(logFile, sprintf( '\n\n\n>>Adjustments to the variance-covariance table due to the variance subtraction\n'), parameters);

xmax = 0.0;                 % determine scaling factor */
for channel=1:parameters.eegChannels
    if (blink.adjust.variance(channel) > xmax)
        xmax = blink.adjust.variance(channel);
    end
    if (blink.adjust.veogCovar(channel) > xmax)
        xmax = blink.adjust.veogCovar(channel);
    end
    if (blink.adjust.heogCovar(channel) > xmax)
        xmax = blink.adjust.heogCovar(channel);
    end
    if (saccade.adjust.variance(channel) > xmax)
        xmax = saccade.adjust.variance(channel);
    end
    if (saccade.adjust.veogCovar(channel) > xmax)
        xmax = saccade.adjust.veogCovar(channel);
    end
    if (saccade.adjust.heogCovar(channel) > xmax)
        xmax = saccade.adjust.heogCovar(channel);
    end
end

scal = 1.0;

if (xmax >= 1000000.0)
    scal = 10000.0;
elseif (xmax >= 100000.0)
    scal = 1000.0;
elseif (xmax >= 10000.0)
    scal = 100.0;
elseif (xmax >= 1000.0)
    scal = 10.0;
end

logFile=display(logFile, sprintf('\n Scaling factor = %6.0f\n',scal), parameters);

% blinks */
logFile=display(logFile, sprintf('\n Blinks\n                 '), parameters);
for channel=1:parameters.eegChannels    % write channel index */
    logFile=display(logFile, sprintf('%6d ',channel), parameters);
end

logFile=display(logFile, sprintf('\n Variance        '), parameters);
for channel=1:parameters.eegChannels    % write variance for blinks */
    logFile=display(logFile, sprintf('%6.0f ', blink.adjust.variance(channel)/scal), parameters);
end

logFile=display(logFile, sprintf('\n covar chn-Veog  '), parameters);
for channel=1:parameters.eegChannels    % write covariance for blinks */
    logFile=display(logFile, sprintf('%6.0f ', blink.adjust.veogCovar(channel)/scal), parameters);
end
logFile=display(logFile, sprintf('\n covar chn-Heog  '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf('%6.0f ', blink.adjust.heogCovar(channel)/scal), parameters);
end

% Saccades */
logFile=display(logFile, sprintf('\n\n Saccades\n                 '), parameters);
for channel=1:parameters.eegChannels    % write channel index */
    logFile=display(logFile, sprintf('%6d ',channel), parameters);
end

logFile=display(logFile, sprintf('\n Variance        '), parameters);
for channel=1:parameters.eegChannels    % write variance for saccades */
    logFile=display(logFile, sprintf('%6.0f ', saccade.adjust.variance(channel)/scal), parameters);
end

logFile=display(logFile, sprintf('\n covar chn-Veog  '), parameters);
for channel=1:parameters.eegChannels    % write covariance for saccades */
    logFile=display(logFile, sprintf('%6.0f ', saccade.adjust.veogCovar(channel)/scal), parameters);
end
logFile=display(logFile, sprintf('\n covar chn-Heog  '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf('%6.0f ', saccade.adjust.heogCovar(channel)/scal), parameters);
end

%***********************/
%  Residual variance   */
logFile=display(logFile, sprintf( '\n\n Residual variance-covariance table '), parameters);

xmax = 0.0;                 % determine scaling factor */
for channel=1:parameters.eegChannels
    if (blink.residual.variance(channel) > xmax)
        xmax = blink.residual.variance(channel);
    end
    if (blink.residual.veogCovar(channel) > xmax)
        xmax = blink.residual.veogCovar(channel);
    end
    if (blink.residual.heogCovar(channel) > xmax)
        xmax = blink.residual.heogCovar(channel);
    end
    if (saccade.residual.variance(channel) > xmax)
        xmax = saccade.residual.variance(channel);
    end
    if (saccade.residual.veogCovar(channel) > xmax)
        xmax = saccade.residual.veogCovar(channel);
    end
    if (saccade.residual.heogCovar(channel) > xmax)
        xmax = saccade.residual.heogCovar(channel);
    end
end

scal = 1.0;

if (xmax >= 1000000.0)
    scal = 10000.0;
elseif (xmax >= 100000.0)
    scal = 1000.0;
elseif (xmax >= 10000.0)
    scal = 100.0;
elseif (xmax >= 1000.0)
    scal = 10.0;
end

logFile=display(logFile, sprintf('\n Scaling factor = %6.0f\n',scal), parameters);

% blinks */
logFile=display(logFile, sprintf('\n Blinks\n                 '), parameters);
for channel=1:parameters.eegChannels    % write channel index */
    logFile=display(logFile, sprintf('%6d ',channel), parameters);
end

logFile=display(logFile, sprintf('\n Variance        '), parameters);
for channel=1:parameters.eegChannels    % write variance for blinks */
    logFile=display(logFile, sprintf('%6.0f ', blink.residual.variance(channel)/scal), parameters);
end

logFile=display(logFile, sprintf('\n covar chn-Veog  '), parameters);
for channel=1:parameters.eegChannels    % write covariance for blinks */
    logFile=display(logFile, sprintf('%6.0f ', blink.residual.veogCovar(channel)/scal), parameters);
end
logFile=display(logFile, sprintf('\n covar chn-Heog  '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf('%6.0f ', blink.residual.heogCovar(channel)/scal), parameters);
end

% Saccades */
logFile=display(logFile, sprintf('\n\n Saccades\n                 '), parameters);
for channel=1:parameters.eegChannels    % write channel index */
    logFile=display(logFile, sprintf('%6d ',channel), parameters);
end

logFile=display(logFile, sprintf('\n Variance        '), parameters);
for channel=1:parameters.eegChannels    % write variance for saccades */
    logFile=display(logFile, sprintf('%6.0f ', saccade.residual.variance(channel)/scal), parameters);
end

logFile=display(logFile, sprintf('\n covar chn-Veog  '), parameters);
for channel=1:parameters.eegChannels    % write covariance for saccades */
    logFile=display(logFile, sprintf('%6.0f ', saccade.residual.veogCovar(channel)/scal), parameters);
end
logFile=display(logFile, sprintf('\n covar chn-Heog  '), parameters);
for channel=1:parameters.eegChannels
    logFile=display(logFile, sprintf('%6.0f ', saccade.residual.heogCovar(channel)/scal), parameters);
end

%***************************************************************************/
% compute.c : handles all of the computational subroutines necessary for   */
%             the regression performed (and described in the routines      */
%                                                                          */
% Current version written by James M. Turner : August 4, 1993.             */
%      jturner@p300.cpl.uiuc.edu                                           */
%      Neuroscience, Washington & Lee University                           */
%      for the Dept. of Psychology, Cognitive Psychophysiology Lab.        */
%      Champaign, Illinois.                                                */
%***************************************************************************/
% #include <stdio.h>
% #include <stdlib.h>
% #include <math.h>
% #include 'emcp.h'
% #include 'compute.h'
% 
% %** Variables ***/
% extern FILE    *parameterFile,                % file containing params */
%                *dataFile,                     % raw data file */
%                *calibrationFile,              % calibration values for AD conv */
%                *correctedSngTrlFile,          % corrected single trial file */
%                *logFile;                      % log of activity */
% 
% 
% extern struct  Parameters parameters;         % parameters for input file */
% 
% extern float   ***rawData,                    % raw data array */
%                ***rawIDs;                     % raw IDs array     */
% 
% extern int     *trialsInBin,                  % number of trials in each bin  */
%                 **markBlink;                  % blink present at [point] */
%                                               % in curTrial   */
% extern double  *curIDs,                       % current IDs */
%                *trialBaseline,                % buffer for baseline of trial */
%                **binBaseline;                 % buffer for baseline of bins */
% 
% extern struct  TrialInfo *trialInfo;          % trial info array */
%             
% extern int     verbose;                       % display info on screen, default to yes */
% 
% %** Structs for correlation determinations ***/
% extern struct  Corr    blink,                 % used to organize arrays  */
%                        saccade;               % associated with calculation */
%                                               % of correlations among */
%                                               % residual.heogCovar & HEOG channels */
%                                               % and Fz */
% 
% %** Variables for blink detection  ***/
% extern int     msTen,                         % # samples (points) in 10 ms */
%                thirdOfBlink,                  % # pts in 1/3 of blink window */
%                initBlinkScan,                 % 1st pt to scan for blinks */
%                endBlinkScan,                  % last pt to scan for blinks */
%                middleOfBlink;                 % pt in middle of blink window */
% extern double  blinkCriteria,                 % criteria for blink */
%                windVariance=2.0,              % window variance */
%                lengthOfWind;                  % # points in blink template */
% 
% %** Variables for EEG rejection ***/
% extern double  criteria;                      % criteria for out-of-scale */
% 
% %** Variables for artifact removal ***/
% extern double  totalPoints;                   % total # data points */
% extern float   ***rawAverage;                 % average wave in each bin */
% extern double  *calibration,                  % calibration vals for AD converters */
%                **meanIDs;                     % average of ID values */


%***********************************************************************/
%  computeCorrectedSingleTrials()  : corrects the trials               */
%                                                                      */
%  Follows same procedure as computeCorrectedAverageTrials expect      */
%  loops through trials instead of through bins.                       */
%***********************************************************************/
function [corRawData, subtractedBlinks, subtractedSaccades]=computeCorrectedSingleTrials(rawData, blink, saccade, binBaseline, trialBaseline, windVariance, lengthOfWind, initBlinkScan, endBlinkScan, thirdOfBlink, middleOfBlink, trialInfo, blinkCriteria, markBlink, parameters)

subtractedBlinks=zeros(size(rawData)); %JD
subtractedSaccades=zeros(size(rawData)); %JD

for trial=1:parameters.numTrials
    % loop through all trials */
    if trialInfo(trial).status==0
        % if this is an accepted trial continue */
        [rawData, binBaseline, trialBaseline]= computeTrialBaseline('REMOVE', trial, rawData, binBaseline, trialBaseline, trialInfo, parameters);
        % remove average activity from trial */
        markBlink=checkBlink(trial,windVariance,blinkCriteria,markBlink,lengthOfWind, initBlinkScan, endBlinkScan, thirdOfBlink, middleOfBlink, rawData, parameters);
        % check for blinks sans trial baseline */
        for point=1:parameters.trialLength
            % loop through length of trial */
            if markBlink(trial,point)==1.0
                % if blink at this point*/
                deltaBlinkVEOG=rawData(trial,point,parameters.VEOG_CHANNEL)-blink.channelMean(parameters.VEOG_CHANNEL);
                % difference between average blink activity at bin and total average activity */
                deltaBlinkHEOG=rawData(trial,point,parameters.HEOG_CHANNEL)-blink.channelMean(parameters.HEOG_CHANNEL);
                % difference between average saccade activity in bin and total average activity */
                for channel=1:parameters.eegChannels %eliminated assumption that Fz is the first channel JD
                    % loop through just the EEG channels */
                    avgAdjustment=blink.veogFactor(channel)*deltaBlinkVEOG+blink.heogFactor(channel)*deltaBlinkHEOG;
                    rawData(trial,point,channel)=rawData(trial,point,channel)-(avgAdjustment+blink.channelMean(channel));
                    subtractedBlinks(trial,point,channel)=(avgAdjustment+blink.channelMean(channel)); %JD
                end % end channel */
            else
                % must be saccade */
                deltaSaccadeVEOG=rawData(trial,point,parameters.VEOG_CHANNEL)-saccade.channelMean(parameters.VEOG_CHANNEL);
                % difference between average saccade activity at bin and total average activity */
                deltaSaccadeHEOG=rawData(trial,point,parameters.HEOG_CHANNEL)-saccade.channelMean(parameters.HEOG_CHANNEL);
                % difference between average saccade activity in bin and total average activity */
                for channel=1:parameters.eegChannels %eliminated assumption that Fz is the first channel JD
                    % loop through the EEG channels */
                    avgAdjustment=saccade.veogFactor(channel)*deltaSaccadeVEOG+saccade.heogFactor(channel)*deltaSaccadeHEOG;
                    rawData(trial,point,channel)=rawData(trial,point,channel)-(avgAdjustment+saccade.channelMean(channel));
                    subtractedSaccades(trial,point,channel)=(avgAdjustment+saccade.channelMean(channel)); %JD
                end % end channel */
            end % end if markBlink */
        end % end point */
        [rawData, binBaseline, trialBaseline]= computeTrialBaseline('RESTORE', trial, rawData, binBaseline, trialBaseline, trialInfo, parameters);
        % restore average activity in trial */
    end %end if accepted */
end % end trial */
corRawData=rawData;

%***********************************************************************/
%  computeCorrectAvgs()  : corrects the averages for eye movements     */
%                                                                      */
%  In this function we compute the raw average for blink points and    */
%  saccade points separately (based on blink.sum and saccade.sum)      */
%  and then after subtracting the channel mean (blink.channelMean and  */
%  saccade.channelMean) we regress every chan on VEOG and HEOG.        */
%  Note, the regression model does not have the constant term.         */
%  (intercept).                                                        */
%***********************************************************************/
function [corrAverage, subtractedBlinks, subtractedSaccades]=computeCorrectAvgs(blink, saccade, trialsInBin, rawAverage, binBaseline, parameters)

%   int     bin,                       % counters */
%           point,
%           channel;
% 
%   double  deltaBlinkVEOG,
%           deltaBlinkHEOG,
%           deltaSaccadeVEOG,
%           deltaSaccadeHEOG,
%           avgAdjustment;

subtractedBlinks=zeros(size(rawAverage)); %JD
subtractedSaccades=zeros(size(rawAverage)); %JD

for bin=1:parameters.numStorageBins
    % loop through each bin */
    if trialsInBin(bin)
        % if there are trials in storage bin continue */
        for point=1:parameters.trialLength
            % loop through points */
            deltaBlinkVEOG=0;
            deltaBlinkHEOG=0;
            deltaSaccadeVEOG=0;
            deltaSaccadeHEOG=0;
            blink.weight=0;
            saccade.weight=0;
            if blink.pointsInBin(point,bin)>=1.0
                for channel=1:parameters.totalChannels
                    % loop through all channels */
                    blink.sum(point,channel,bin)=blink.sum(point,channel,bin)/blink.pointsInBin(point,bin);
                    % divide to give average activity */
                end
                
                deltaBlinkVEOG=blink.sum(point,parameters.VEOG_CHANNEL,bin)-blink.channelMean(parameters.VEOG_CHANNEL);
                % difference between average blink activity at bin and total average activity */
                deltaBlinkHEOG=blink.sum(point,parameters.HEOG_CHANNEL,bin)-blink.channelMean(parameters.HEOG_CHANNEL);
                % difference between average saccade activity in bin and total average activity */
                
                blink.weight=blink.pointsInBin(point,bin)/trialsInBin(bin);
                % determine weight of blink at time point */
            end % end if (blink.pointsInBin) */
            
            if saccade.pointsInBin(point,bin)>=1.0
                for channel=1:parameters.totalChannels
                    % loop through all channels EEG & nonEEG */
                    saccade.sum(point,channel,bin)=saccade.sum(point,channel,bin)/saccade.pointsInBin(point,bin);
                    % divide to give average saccade activity */
                end
                
                deltaSaccadeVEOG=saccade.sum(point,parameters.VEOG_CHANNEL,bin)-saccade.channelMean(parameters.VEOG_CHANNEL);
                % difference between average saccade activity at bin and total average activity */
                deltaSaccadeHEOG=saccade.sum(point,parameters.HEOG_CHANNEL,bin)-saccade.channelMean(parameters.HEOG_CHANNEL);
                % difference between average saccade activity in bin and total average activity */
                
                saccade.weight=saccade.pointsInBin(point,bin)/trialsInBin(bin);
                % determine weight of saccade at time point */
            end
            
            for channel=1:parameters.eegChannels %removed assumption that the first EEG channel is Fz. JD
                % loop through all EEG channels */
                blink.adjustment=blink.weight*(blink.veogFactor(channel)*deltaBlinkVEOG+blink.heogFactor(channel)*deltaBlinkHEOG);
                % calculate the adjustment for blinks */
                saccade.adjustment=saccade.weight*(saccade.veogFactor(channel)*deltaSaccadeVEOG+saccade.heogFactor(channel)*deltaSaccadeHEOG);
                % calculate the adjustment for saccades */
                avgAdjustment=blink.adjustment+saccade.adjustment;
                % combine the two */
                
                rawAverage(point,channel,bin)=rawAverage(point,channel,bin)-(avgAdjustment+blink.weight*blink.channelMean(channel)+saccade.weight*saccade.channelMean(channel));
                subtractedBlinks(point,channel,bin)=(blink.weight*blink.channelMean(channel))+blink.adjustment; %JD
                subtractedSaccades(point,channel,bin)=(saccade.weight*saccade.channelMean(channel))+saccade.adjustment; %JD
                % make the correction */
            end% end channel */
            for channel=1:parameters.totalChannels
                % loop through all channels */
                rawAverage(point,channel,bin)=rawAverage(point,channel,bin)+binBaseline(channel,bin);
                % add in the average activity for the bin (ERP wave) */
            end % end channel */
        end % end points */
    end % end if */
end % end bin */
corrAverage=rawAverage;

%***********************************************************************/
%  computePropogation()  : computes contribution of EOG to EEG chans   */
%                                                                      */
%  This routine returns veogFactors and heogFactors which are beta     */
%  weights for a multivariate linear regression of every EEG chan      */
%  on VEOG and HEOG.  The standardized beta weights are computed using */
%  computational formulas based on the correlations between EEG and    */
%  VEOG and HEOG and between VEOG and HEOG.  The formula is:           */
%                                        Ry1-Ry2*R12                   */
%  standardized B weights = xeogFactor= --------------                 */
%                                        1 - R12*R12                   */
%                                                                      */
%  After this the standardized B weights are converted into            */
%  unstandardized B weights which are propogation factors used for     */
%  correcting the data.                                                */
%***********************************************************************/
function [heogFactor, veogFactor]=computePropogation(veogCorr, heogCorr, stanDev, parameters)

% 						% veog corr w/ all channels from computeCorrEOG */
%                         double *heogCorr,
%                         % heog corr w/ all channels from computeCorrEOG */
%                         double *stanDev,
%                         % standard deviation from computeResVar */
%                         double **veogFactor,
%                         % pass pnt to array of correction factors */
%                         double **heogFactor
%                         % pass pnt to array of correction factors */
%                         )
%     double  *rawVEOGCorr,                       % temp storage for raw */
%             *rawHEOGCorr;                       % correlations */

veogFactor=zeros(parameters.eegChannels,1);
heogFactor=zeros(parameters.eegChannels,1);
rawVEOGCorr=zeros(parameters.eegChannels,1);
rawHEOGCorr=zeros(parameters.eegChannels,1);
% allocate memory for temporary storage */

%** for blinks ***/
for channel=1:parameters.eegChannels
    rawVEOGCorr(channel)=veogCorr(channel);
    % copy over blink corr to temp variable */
    rawHEOGCorr(channel)=heogCorr(channel);
    % copy over blink corr to temp variable */
end

%* Brian Foote  3/9/95 -- Added the else parts to the computations below...

for channel=1:parameters.eegChannels
    if rawVEOGCorr(parameters.HEOG_CHANNEL)~=1
        % compute the standardized Beta Weights */
        veogCorr(channel)=(rawVEOGCorr(channel)-rawHEOGCorr(channel)*rawVEOGCorr(parameters.HEOG_CHANNEL))/(1.0-(rawVEOGCorr(parameters.HEOG_CHANNEL)*rawVEOGCorr(parameters.HEOG_CHANNEL)));
    else
        veogCorr(channel)=0.0 ; % Is this reasonable? */
    end
    
    if(rawHEOGCorr(parameters.VEOG_CHANNEL)~=1)
        heogCorr(channel)=(rawHEOGCorr(channel)-rawVEOGCorr(channel)*rawHEOGCorr(parameters.VEOG_CHANNEL))/(1.0-(rawHEOGCorr(parameters.VEOG_CHANNEL)*rawHEOGCorr(parameters.VEOG_CHANNEL)));
    else
        heogCorr(channel)=0.0 ; % Is this reasonable? */
    end
    
    % compute the propogation factors (unstandardized B weights) */
    if(stanDev(parameters.VEOG_CHANNEL))
        veogFactor(channel)=veogCorr(channel)*stanDev(channel)/stanDev(parameters.VEOG_CHANNEL);
    else
        veogFactor(channel)=0.0 ;
    end
    
    if(stanDev(parameters.HEOG_CHANNEL))
        heogFactor(channel)=heogCorr(channel)*stanDev(channel)/stanDev(parameters.HEOG_CHANNEL);
    else
        heogFactor(channel)=0.0 ;
    end
    
end % end channel */
                    
%***********************************************************************/
%  computeResidualVariance() : computes residual variance and stan dev */
%                                                                      */
%  Returns the difference between the variance over trials and the     */
%  variance over raw averages (computed in computeVarForTrialWaves and */
%  computeVarForRawWaves).                                             */
%  Note also that variance refers to covariance for heog and veog.     */
%***********************************************************************/
function [res, stanDev]=computeResidualVariance(raw, adjust, parameters)
% 										% trial based co/variance */
%                                         struct Variance adjust,
%                                         % average based co/variance */
%                                         double **stanDev)
%                                         % pass pnt to array of stan dev */

stanDev=zeros(parameters.eegChannels,1);

res.variance=zeros(parameters.eegChannels,1);
res.veogCovar=zeros(parameters.eegChannels,1);
res.heogCovar=zeros(parameters.eegChannels,1);

for channel=1:parameters.eegChannels
    res.variance(channel)=raw.variance(channel)-adjust.variance(channel);
    res.veogCovar(channel)=raw.veogCovar(channel)-adjust.veogCovar(channel);
    res.heogCovar(channel)=raw.heogCovar(channel)-adjust.heogCovar(channel);
    % take the differences to get the residual variance */
    
    if res.variance(channel)<0
        res.variance(channel)=0;
        % if negative then make positive */
    end
    stanDev(channel)=sqrt(res.variance(channel));
    % take sqroot for standard deviation */
end % end channel */
return;


%***********************************************************************/
%  computeCorrWithEOG() : computes corr between every chan & EOG       */
%                                                                      */
%  Compute all pair-wise correlations between the EEG channels and     */
%  VEOG and HEOG.  These correlations are based on residual covars     */
%  and standard deviations which were both computed in                 */
%  computeResidualVar()                                                */
%***********************************************************************/
function [veogCorr, heogCorr]=computeCorrWithEOG(stanDev, res,  parameters)
% 						% stan dev from computeResVar */
%                         struct Variance res,
%                         % residual variance from computeResVar */
%                         double **veogCorr,
%                         % pointer to array of correlations for return */
%                         double **heogCorr)
%                         % pointer to array of correlations for return */
veogCorr=zeros(parameters.eegChannels,1);
heogCorr=zeros(parameters.eegChannels,1);
% allocate memory */

for channel=1:parameters.eegChannels
    % loop through EOG and EEG channels */
    % computes correlation by taking covariance/(stan dev*stan dev) */
    % if no negatives in denom then determine corr w/ HEOG */
    % Added the else parts <*BF*> 27 March 1995 */
    if (stanDev(channel)) && (stanDev(parameters.VEOG_CHANNEL))
        veogCorr(channel)=res.veogCovar(channel)/(stanDev(channel)*stanDev(parameters.VEOG_CHANNEL));
    else
        veogCorr(channel) = 0.0 ;
    end
    % if no negatives in denom then determine corr w/ VEOG */
    if (stanDev(channel)) && (stanDev(parameters.HEOG_CHANNEL))
        heogCorr(channel)=res.heogCovar(channel)/(stanDev(channel)*stanDev(parameters.HEOG_CHANNEL));
    else
        heogCorr(channel) = 0.0 ; %changed from veogCorr to heogCorr JD
    end
    % if no negatives in denom then determine corr w/ HEOG */
end % end channel */

%***********************************************************************/
%  computeVarForRawWaves()  : computes adjustments due to averages     */
%                                                                      */
%  You compute the variance around an estimate of the systematic       */
%  activity in the data.  Where systematic activity is the raw         */
%  average wave form computed across raw single trials, for all bins   */
%  and channels.  The variance is computed across all bins, points     */
%  for channels.                                                       */
%                                                                      */
%***********************************************************************/
function [adjust]=computeVarForRawWaves(sum, rawAverage, ptsInBase, baseline, pointsInBin, parameters)
% 									  % blink/sac component of rawAvg */
% 									  % from sumTrials			       */
%                                       float  ***rawAverage,
%                                       % raw average for each bin      */
%                                       double ptsInBase,
%                                       % number of blink/saccade pts   */
%                                       double *baseline,
%                                       % mean total activity by channel */
%                                       int    **pointsInBin)
%                                       % # of blink/sac pts by pts by bin */

weightedAvg=zeros(parameters.eegChannels,1);
adjust.variance=zeros(parameters.eegChannels,1);
adjust.veogCovar=zeros(parameters.eegChannels,1);
adjust.heogCovar=zeros(parameters.eegChannels,1);

for bin=1:parameters.numStorageBins
    % loop through each bin */
    for point=1:parameters.trialLength
        % loop through all points */
        if pointsInBin(point,bin)>=1.0
            % if there are points for correction */
            for channel=1:parameters.eegChannels
                % loop through VEOG,HEOG and EEG channels */
                adjust.variance(channel)=adjust.variance(channel)+rawAverage(point,channel,bin)*((2.0*sum(point,channel,bin))-(pointsInBin(point,bin)*rawAverage(point,channel,bin)));
                adjust.veogCovar(channel)=adjust.veogCovar(channel)+(rawAverage(point,channel,bin)*sum(point,parameters.VEOG_CHANNEL,bin))+(rawAverage(point,parameters.VEOG_CHANNEL,bin)*sum(point,channel,bin))-(pointsInBin(point,bin)*rawAverage(point,channel,bin)*rawAverage(point,parameters.VEOG_CHANNEL,bin));
                adjust.heogCovar(channel)=adjust.heogCovar(channel)+(rawAverage(point,channel,bin)*sum(point,parameters.HEOG_CHANNEL,bin))+(rawAverage(point,parameters.HEOG_CHANNEL,bin)*sum(point,channel,bin))-(pointsInBin(point,bin)*rawAverage(point,channel,bin)*rawAverage(point,parameters.HEOG_CHANNEL,bin));
                weightedAvg(channel)=weightedAvg(channel)+(pointsInBin(point,bin)*rawAverage(point,channel,bin));
                % weighted average based on number of saccade points used later */
            end % end channel */
        end % end if points */
    end % end points */
end % end bin */

%** Compute covariance with total mean ***/
for channel=1:parameters.eegChannels
    % loop through VEOG,HEOG and EEG channels */
    % Added the else part...  <*BF*> 27 March 1995 */
    if ptsInBase
        weightedAvg(channel)=weightedAvg(channel)/ptsInBase;
    else
        weightedAvg(channel) = 0;
    end
    % make average based on total blink points */
end % end channel */

%** Compute actual adjustment values of the variance and ***/
%** covariance due to the averages                       ***/
for channel=1:parameters.eegChannels
    % loop through VEOG,HEOG and EEG channels */
    if ptsInBase
        adjust.variance(channel)=(adjust.variance(channel)/ptsInBase)-weightedAvg(channel)*((2.0*baseline(channel))-weightedAvg(channel));
        adjust.veogCovar(channel)=(adjust.veogCovar(channel)/ptsInBase)+(weightedAvg(channel)*weightedAvg(parameters.VEOG_CHANNEL))-(weightedAvg(channel)*baseline(parameters.VEOG_CHANNEL))-(weightedAvg(parameters.VEOG_CHANNEL)*baseline(channel));
        adjust.heogCovar(channel)=(adjust.heogCovar(channel)/ptsInBase)+(weightedAvg(channel)*weightedAvg(parameters.HEOG_CHANNEL))-(weightedAvg(channel)*baseline(parameters.HEOG_CHANNEL))-(weightedAvg(parameters.HEOG_CHANNEL)*baseline(channel));
    end
end % end channel */

return;
    
%***********************************************************************/
%  computeVarForTrialWaves()  : computes the total variance in raw data        */
%  That is,                                                            */
%  variance across trials, points and bins by channel.  Note the       */
%  distinction between computeVarForTrialWaves() and computeAdjustVar() is     */
%  that this is a measure of total variance in the uncorrected trial   */
%  whereas computeAdjustVar() is the measure of total variance based   */
%  on the rawAverage array.                                            */
%***********************************************************************/
function raw=computeVarForTrialWaves(baseline, ptsInBase, raw, parameters)
%                              % [# ch] mean activity in each channel */
%                              double ptsInBase,
%                              % total # of blink or saccade points in data set */
%                              struct Variance *raw)
%                              % ptr to trial based co/variance for return */

total.variance=zeros(parameters.eegChannels,1);
total.veogCovar=zeros(parameters.eegChannels,1);
total.heogCovar=zeros(parameters.eegChannels,1);
% allocate memory for total data */

for channel=1:parameters.eegChannels
    total.variance(channel)=baseline(channel)*baseline(channel);
    total.veogCovar(channel)=baseline(channel)*baseline(parameters.VEOG_CHANNEL);
    total.heogCovar(channel)=baseline(channel)*baseline(parameters.HEOG_CHANNEL);
    % baseline from baseline(), average blink activity by channel */
end

% computer total variance and covariance using the formula */
% variance=sqrt[(SUM(i=1..N)(Xi*Xi/N)-(Xbar*Xbar))] */
for channel=1:parameters.eegChannels
    if (ptsInBase)
        % make sure we don't divide by zero */
        raw.variance(channel)=raw.variance(channel)/ptsInBase-total.variance(channel);
        raw.veogCovar(channel)=raw.veogCovar(channel)/ptsInBase-total.veogCovar(channel);
        raw.heogCovar(channel)=raw.heogCovar(channel)/ptsInBase-total.heogCovar(channel);
    else % <*BF*> 27 March 1995 */
        raw.variance(channel)= 0.0;
        raw.veogCovar(channel)=0.0;
        raw.heogCovar(channel)=0.0;
    end
end % end channel */

%***********************************************************************/
%  computeBinAvg()  : produces the average wave form in each bin       */
%                                                                      */
%  Computes the raw average for each bin.  The mean activity in each   */
%  trial has been removed since based on .sum.                         */
%                                                                      */
%  Dependent on:                                                       */
%***********************************************************************/
function [rawAverage, binBaseline]=computeBinAvg(bComponent, sComponent, trialsInBin, binBaseline, parameters)

%     int     point,                          % counters */
%             channel,
%             bin;
%     
%     float   ***rawAverage;

rawAverage=zeros(parameters.trialLength,parameters.totalChannels,parameters.numStorageBins);
% array to hold raw average wave for each bin */

for bin=1:parameters.numStorageBins
    for channel=1:parameters.totalChannels
        for point=1:parameters.trialLength
            if (trialsInBin(bin))
                rawAverage(point,channel,bin)=bComponent(point,channel,bin)+sComponent(point,channel,bin);
                rawAverage(point,channel,bin)=rawAverage(point,channel,bin)/trialsInBin(bin);
            else % <*BF*> 27 March 1995 */
                rawAverage(point,channel,bin)=0.0 ;
            end
        end % end point */
        binBaseline(channel,bin)=binBaseline(channel,bin)/trialsInBin(bin);
        % compute avg baseline activity in bin also */
    end % end channel */
end % end bin */
return;

%***********************************************************************/
%  computeBaseline()  :  average blink/saccade activity in each channel*/
%                                                                      */
%  Get mean blink and saccade activity collapsing over bins and        */
%  points.  The average activity of each trial has been removed,       */
%  since this is based on .sum which is determined in sumTrial.        */
%                                                                      */
%  Dependent on:                                                       */
%***********************************************************************/
function [baseline, ptsInBase]=computeBaseline(pointsInBin, sum, parameters)
%                  % # pts used in calc *baseline */
%                  int     **pointsInBin, 
%                  % [# pts][# bins] # pts in each bin over all pts */
%                  % in trial */
%                  float   ***sum 
%                  % [# pts][# ch][# bins] contribution to raw average */
%                  % of blink or saccade data */
%                 )                    

%   double   *baseline;   % [channels] stores mean activity */

baseline=zeros(parameters.totalChannels,1);

ptsInBase=0;
% zero denominators */

for bin=1:parameters.numStorageBins
    % go through each storage bin */
    for point=1:parameters.trialLength
        % go through each point in a trial */
        if (pointsInBin(point,bin))
            % determined in sumTrial */
            % if there are points in this bin @ this point do...*/
            ptsInBase=ptsInBase+pointsInBin(point,bin);
            % add # pts to total for denominator */
            for channel=1:parameters.totalChannels
                % go through each channel */
                baseline(channel)=baseline(channel)+sum(point,channel,bin);
                % add activity to baseline */
                % sum determined in sumTrial */
            end
        end
        %this seems to be an error check for the case where there are no channels but it doesn't parse right and is not going to come up normally anyway.  JD
%         else % <*BF*> 27 March 1995 */
%             baseline(channel) = 0 ;
%         end
    end % end if points */
end % end for storage bins */
for channel=1:parameters.totalChannels
    % go through each channel */
    if ptsInBase>0
        baseline(channel)=baseline(channel)/ptsInBase;
        % gives average activity in each channel */
    end
end % end for channel */
return;

%***********************************************************************/
%  computeTrialBaseline()  : removes or adds baseline to trial         */
%                                                                      */
%  This routine takes average activity in each channel for the trial   */
%  and subtracts it from each point.  It also restores this data when  */
%  remove = 0.                                                         */
%***********************************************************************/
function [rawData, binBaseline, trialBaseline]= computeTrialBaseline(remove, trial, rawData, binBaseline, trialBaseline, trialInfo, parameters)

if strcmp(remove,'REMOVE')
    % if we are removing the baseline then do: */
    for channel=1:parameters.totalChannels % jt */
        trialBaseline(channel)=0;
    end
    
    % clear baseline buffer */
    for point=1:parameters.trialLength
        for channel=1:parameters.totalChannels
            trialBaseline(channel)=trialBaseline(channel)+rawData(trial,point,channel);
            % sum up for each channel */
        end
    end
    for channel=1:parameters.totalChannels
        trialBaseline(channel)=trialBaseline(channel)/parameters.trialLength;
        % get baseline for each channel by taking average */
        % activity during the entire trial in each channel */
    end
    
    for channel=1:parameters.totalChannels
        binBaseline(channel,trialInfo(trial).classification)=binBaseline(channel,trialInfo(trial).classification)+trialBaseline(channel);
        % also sum up so we can create averages later by adding */
        % up now instead of later */
    end
end
for point=1:parameters.trialLength % jt */
    for channel=1:parameters.totalChannels
        if strcmp(remove,'REMOVE')
            rawData(trial,point,channel)=rawData(trial,point,channel)-trialBaseline(channel);
            % substract average if we are removing baseline */
        else
            rawData(trial,point,channel)=rawData(trial,point,channel)+trialBaseline(channel);
            % add averages back in */
        end
    end
end
            
%***********************************************************************/
%  sumTrial()  : sum up a trial for correlation figures                */
%                                                                      */
%  The variance of each channel and the covariance between VEOG, HEOG  */
%  and each channel is computed.  This is done separately for blink    */
%  points and saccade points.                                          */
%  The average of each channel for each trial has already been removed */
%  thus this is based on the data - its average                        */
%                                                                      */
%  Dependent on:                                                       */
%    rawData, markBlink, .sum, .variance, .veogCovar, .heogCovar        */
%    for blinks and saccade                                            */
%                                                                      */
%***********************************************************************/
function [blink, saccade]= sumTrial(trial, blink, saccade, markBlink, rawData, trialInfo, parameters)            
        
for pt=1:parameters.trialLength
    for ch=1:parameters.totalChannels % jt */
        if markBlink(trial,pt)==1
            % if point is marked as a blink */
            if parameters.VEOG_CHANNEL==ch
                blink.pointsInBin(pt,trialInfo(trial).classification)=blink.pointsInBin(pt,trialInfo(trial).classification)+1;
            end
            % increment counter */
            blink.sum(pt,ch,trialInfo(trial).classification)=blink.sum(pt,ch,trialInfo(trial).classification)+rawData(trial,pt,ch);
            % (blink.sum + saccade.sum)/# points in bin = raw average */
            blink.raw.variance(ch)=blink.raw.variance(ch)+rawData(trial,pt,ch)*rawData(trial,pt,ch);
            % NOTE:  Variance labelling below is inaccurate.  In        */
            %        computeVarForTrialWaves() further operations are performed */
            %        which result in the variances and covariances      */
            % standard variance computation, mean x was removed earlier */
            blink.raw.veogCovar(ch)=blink.raw.veogCovar(ch)+rawData(trial,pt,ch)*rawData(trial,pt,parameters.VEOG_CHANNEL);
            % standard covariance computation, mean x was removed ealier */
            blink.raw.heogCovar(ch)=blink.raw.heogCovar(ch)+rawData(trial,pt,ch)*rawData(trial,pt,parameters.HEOG_CHANNEL);
            % standard covariance computation, mean x was removed earlier */
        else
            % point is a saccade point if not a blink point */
            if(parameters.VEOG_CHANNEL==ch)
                saccade.pointsInBin(pt,trialInfo(trial).classification)=saccade.pointsInBin(pt,trialInfo(trial).classification)+1;
                saccade.sum(pt,ch,trialInfo(trial).classification)=saccade.sum(pt,ch,trialInfo(trial).classification)+rawData(trial,pt,ch);
            end
            % (blink.sum + saccade.sum)/# points in bin = raw average */
            saccade.raw.variance(ch)=saccade.raw.variance(ch)+rawData(trial,pt,ch)*rawData(trial,pt,ch);
            % NOTE:  Variance labelling below is inaccurate.  In        */
            %        computeVarForTrialWaves() further operations are performed */
            %        which result in the variances and covariances      */
            %  standard variance computation, mean x was removed earlier */
            saccade.raw.veogCovar(ch)=saccade.raw.veogCovar(ch)+rawData(trial,pt,ch)*rawData(trial,pt,parameters.VEOG_CHANNEL);
            % standard covariance computation, mean x was removed earlier */
            saccade.raw.heogCovar(ch)=saccade.raw.heogCovar(ch)+rawData(trial,pt,ch)*rawData(trial,pt,parameters.HEOG_CHANNEL);
            % standard covariance computation, mean x was removed earlier */
        end
    end % end for(pt) */
end % end for(ch) */



