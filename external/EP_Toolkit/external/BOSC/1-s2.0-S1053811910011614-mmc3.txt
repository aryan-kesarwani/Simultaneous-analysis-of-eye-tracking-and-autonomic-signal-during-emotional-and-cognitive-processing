    This file is part of the Better OSCillation detection (BOSC) library.

    The BOSC library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The BOSC library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

    Copyright 2010 Jeremy B. Caplan, Adam M. Hughes, Tara A. Whitten
    and Clayton T. Dickson.


The BOSC Method: Brief Tutorial
===============================

Core functions
--------------
BOSC_tf.m - calculates time-frequency spectrograms based on the continuous Morlet wavelet transform

BOSC_bgfit.m - estimates the background spectrum

BOSC_thresholds.m - calculates threshold values as a function of frequency, based on the estimate of the background (BOSC_bgfit). The power threshold is based on a percentile cutoff of the theoretical chi-square distribution of wavelet power values at each frequency. The duration threshold is converted into numbers of samples, scaling correctly for each frequency.

BOSC_detect.m - detects oscillatory episodes in a target signal of interest, based on the thresholds calculated by BOSC_thresholds.m


Important note: avoiding edge artifacts
---------------------------------------

   As is important for data-analyses methods based on wavelet or other Fourier-domain method, data-windowing is an important element of time-frequency analysis. An appropriate data window is necessary to avoid edge artifacts due to not having sufficient sample points 
representing the signal when superimposing particular wavelet used. Because the BOSC method has, in addition to the wavelet window, the duration threshold, it requires a data window that takes into account both  durations, both of which are frequency-dependent. A good rule of thumb is to select a data window that is (at a minimum) defined by the following:

  (Duration threshold + wavelet-window duration) X period of the slowest frequency sampled

Using the variables defined below, we can express this in numbers of samples:

	Fsample*(numcyclesthresh*width/min(F))

Running the Oscillation Analyses
--------------------------------

% INITIAL STEPS: Define the following variables (using example values):

width=6; %  - the wavenumber of the Morlet wavelet (Grossman & Morlet, 1985)

Fsample=250; % Sampling rate of signal, in Hz.

F=(2^(1/4)).^(0:23); % Frequency sampling (resolution) for spectral analysis.  
% Logarithmic sampling of frequencies is recommended for computational purposes when using wavelets. The range can be modified but it must be within the bandpass of the amplifier and/or filter(s).

percentilethresh=.95; % Confidence level (range 0 to 1; percentile of the CDF divided by 100) of the estimated chi-square distribution of spectral power used for the power threshold. Increase to make more conservative.

numcyclesthresh=3; % Duration threshold expressed in numbers of oscillation cycles. Increase to make more conservative.


% STEP ONE: Perform wavelet transform on the signal being used for background power estimates
% NB: Background signal (bgsignal) is the time series (either row or column vector) used to estimate the background spectrum.

[B,T]=BOSC_tf(bgsignal,F,Fsample,width); % Compute the time-frequency (wavelet) spectrogram


% STEP TWO: Fit the background spectrum with a linear regression in log space to estimate mean power X frequency representation.
% Note that we use a shoulder based on the lowest frequency (longest period)

bgshoulder=ceil(width*Fsample/min(F));

[pv,meanpower]=BOSC_bgfit(F,B(:,(bgshoulder+1):(end-bgshoulder)));

% NB: meanpower is the estimated background power spectrum (average power as a function of frequency)
% NB: pv contains the slope and y-intercept of the regression line


% STEP THREE: Calculate the threshold values to use for detection

[powthresh,durthresh]=BOSC_thresholds(Fsample,percentilethresh,numcyclesthresh,F,meanpower);

% *** Hint: At this stage, it is a good idea to cross-check the background power spectrum fit (see PLOT #1: Power spectrum and background spectrum fit)


% STEP FOUR: Set the target signal in which oscillations will be detected.
% The variable called "eegsignal" should contain more EEG signal than you want to analyze. starttime and endtime should be set to mark the bounds of a single trial of interest within eegsignal. Make sure there is enough additional signal before and after starttime and endtime (within eegsignal) to be able to include the shoulder to avoid edge artifacts (see comment above)
% NB: A safe way to calculate the shoulder:

shoulder=ceil((width+numcyclesthresh)*Fsample./F);

eegsignal=  ; % eegsignal is assumed to include the time segment of interest, plus additional signal on either side.
starttime=  ; % starttime contains the sample at which the trial of interest starts
endtime=    ; % endtime contains the sample at which the trial of interest ends

for f=1:length(F) % Loop through all frequencies (if all frequencies are desired; otherwise, specify a subset of frequencies)
  targetsignal=eegsignal((starttime-shoulder(f)):(endtime+shoulder(f)));
  % compute the time-frequency (wavelet) spectrogram, aka "scalogram"
  [Btarget,Ttarget]=BOSC_tf(targetsignal,F(f),Fsample,width); 

  % detect oscillations at frequency F(f)
  detected=BOSC_detect(Btarget,powthresh(f),durthresh(f),Fsample);
  detected=detected((shoulder(f)+1):(end-shoulder(f))); % *** Important: strip off the shoulders
  DETECTED(f,:)=detected; % accumulate all frequencies
end % frequency loop

Pepisode=mean(DETECTED,2); % Pepisode as a function of frequency. This is a useful summary measure.


How to make useful plots
------------------------

% PLOT #1: Power spectrum and background spectrum fit

xf=1:4:length(F); % to get log-sampled frequency tick marks
plot(1:length(F),mean(log10(B(:,(bgshoulder+1):(end-bgshoulder))),2),'ko-',1:length(F),log10(meanpower),'r');
set(gca,'XTick',xf,'XTickLabel',F(1:4:end));
ylabel('Log(Power) [dB]');
xlabel('Frequency [Hz]');


% PLOT #2: Pepisode(f) for a detected segment

bar(1:length(F),Pepisode);
set(gca,'XTick',xf,'XTickLabel',F(1:4:end));
ylabel('P_{episode}'); xlabel('Frequency [Hz]');


% PLOT #3: Raw-trace with superimposed BOSC-detected oscillations at a frequency of interest

finterest= % Define a frequency of interest by indexing into F (i.e., not in Hz)
targetsignal=eegsignal(starttime:endtime);
% create copy of target data and find the detected values that = 1
osc=targetsignal; osc(find(DETECTED(finterest,:)~=1))=NaN; 			
t=(starttime:endtime)/Fsample;
h=plot(t,targetsignal,'b',t,osc,'g'); 
% Plot and label figure
ylabel('Voltage [\mu V]'); xlabel('Time [s]');
set(h(2),'LineWidth',2, 'Color', 'r'); % Colour detected oscillations (at frequency F(finterest)) red and use a thicker line
set(h(1),'LineWidth',1,'Color', 'k');  % Draw remaining signal in a thinner, black line
title(sprintf('P_{episode}(%.2f Hz)=%.2f',F(finterest),Pepisode(finterest))); % Optionally, put the frequency of interest and Pepisode in the title
axis tight;


% PLOT #4: Time-frequency plot of BOSC-detected oscillations

% This plots a black background and white denotes detected oscillations
% To reverse simply invert the values in the matrix: DETECTED=1-DETECTED;
imagesc(t,1:size(DETECTED,1),DETECTED);
colormap(gray);
set(gca,'YTick',xf,'YTickLabel',F(xf));
ylabel('Frequency [Hz]')
xlabel ('Time [s]')
axis xy



References

Caplan, J. B., Madsen, J. R., Raghavachari, S. & Kahana, M. J. (2001) Distinct patterns of brain oscillations underlie two basic parameters of human maze learning. Journal of Neurophysiology, 86, 368-380.

Grossmann, A. and Morlet, J. (1985) Decomposition of functions into wavelets of constant shape, and related transforms. In: Mathematics and physics,
Vol. 1 (Streit L., ed), pp 135–165. Singapore: World Scientific.

van Vugt, M. K., Sederberg, P. B. & Kahana, M. J. (2007) Comparison of spectral analysis methods for characterizing brain oscillations. Journal of Neuroscience Methods, 162, 49-63.

Whitten, T. A., Hughes, A. M., Dickson, C. T. & Caplan, J. B. (submitted) BOSC, a Better OSCillation detection method, robustly extracts EEG rhythms across brain state changes: The human alpha rhythm as a test case.
