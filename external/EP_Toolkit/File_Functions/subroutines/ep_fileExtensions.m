function [fileSuffix,formatName]=ep_fileExtensions(formatCode)
% ep_fileFormats - [fileSuffix,formatName,formatCode]=ep_fileFormats(dataType,formatType) -
% Provides file extension information.
%
%Input:
%   formatCode : Which of the file formats on the list.
%
%Output:
%  fileSuffix : The file name suffix for this type of file format.
%  formatName : Plain English name for the file format.

%History
%  by Joseph Dien (8/4/12)
%  jdien07@mac.com
%
% modified 9/16/13 JD
% Added support for Biosemi bdf files
%
% modified 3/11/14 JD
% Added support for Neuromag fiff files
%
% modified 5/20/14 JD
% Added support for BrainVision EEG files
%
% modified 3/14/20 JD
% Added support for Psychophysiology Toolbox .mat files.
%
% modified 3/5/21 JD
% Added support for multiple suffixes.
%
% modified 7/22/22 JD
% Added support for Matlab .mat files.
%
% modified 10/12/22 JD
% Added support for Neuroelectrics .easy files.
% Added BrainVision DAT and SEG files.
%
% modified 1/1/25 JD
% Added support for BESA .raw files.
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

switch formatCode
    case 'egi_egis'
            fileSuffix='.egis';
            formatName='EGI EGIS session';
    case 'egi_egia'
            fileSuffix='.egis';
            formatName='EGI EGIS average';
    case 'egi_sbin'
        fileSuffix='.sbin';
        formatName='EGI Simple Binary';
    case 'eeglab_set'
        fileSuffix={'.study';'.set'};
        formatName='EEGlab set';
    case 'edf'
        fileSuffix='.edf';
        formatName='Biosig EDF';
    case 'ns_cnt'
        fileSuffix='.cnt';
        formatName='Neuroscan CNT';
    case 'ns_eeg'
        fileSuffix='.eeg';
        formatName='Neuroscan EEG';
    case 'ns_avg'
        fileSuffix='.avg';
        formatName='Neuroscan AVG';
    case 'text'
        fileSuffix={'.txt';'.csv'};
        formatName='text';
    case 'ep_mat'
        fileSuffix='.ept';
        formatName='EP .ept';
    case 'ns_mat'
        fileSuffix='.nsf';
        formatName='NetStation Matlab';
    case 'egi_mff'
        fileSuffix='.mff';
        formatName='NetStation MFF';
    case 'eeglab_erp'
        fileSuffix='.erp';
        formatName='ERPlab erp';
    case 'fieldtrip'
        fileSuffix='.mat';
        formatName='fieldtrip';
    case 'biosemi_bdf'
        fileSuffix='.bdf';
        formatName='Biosemi bdf';
    case 'neuromag_fif'
        fileSuffix='.fif';
        formatName='Neuromag FIFF';
    case 'brainvision_eeg'
        fileSuffix='.eeg';
        formatName='BrainVision EEG';
    case 'brainvision_dat'
        fileSuffix='.dat';
        formatName='BrainVision DAT';
    case 'brainvision_seg'
        fileSuffix='.seg';
        formatName='BrainVision SEG';
    case 'PTB_mat'
        fileSuffix='.mat';
        formatName='PTB mat';
    case 'stacked'
        fileSuffix='.txt';
        formatName='Stacked Text';
    case 'matlab_mat'
        fileSuffix='.mat';
        formatName='Matlab Matrix';
    case 'Neuroelectrics_easy'
        fileSuffix='.easy';
        formatName='Neuroelectrics easy';
    case 'besa_raw'
        fileSuffix='.raw';
        formatName='BESA raw';
    otherwise
        fileSuffix=cell(0);
        formatName='';
        disp(['Oops... ' formatCode ' file format not recognized.']);
end
if ~iscell(fileSuffix)
    tempVar=fileSuffix;
    fileSuffix=cell(0);
    fileSuffix{1}=tempVar;
end


