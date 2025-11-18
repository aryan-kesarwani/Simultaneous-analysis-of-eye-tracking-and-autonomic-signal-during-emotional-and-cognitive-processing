function [fileSuffix,formatName,formatCode]=ep_fileFormats(dataType,formatType)
% ep_fileFormats - [fileSuffix,formatName,formatCode]=ep_fileFormats(dataType,formatType) -
% Provides file format information.
%
%Input:
%	dataType  : The type of the data: 'continuous', 'single_trial', or 'average'
%   formatType: The file format type.
%
%Output:
%  fileSuffix : The file name suffix for this type of file format.
%  formatName : Plain English name for the file format.
%  formatCode : Which of the file formats on the list.

%History
%  by Joseph Dien (5/23/10)
%  jdien07@mac.com
%
% modified 6/16/10 JD
% Changed the suffix of EP files to ".ept".
%
% modified 2/26/11 JD
% Added support for EGI Matlab files
%
% modified 9/16/13 JD
% Added support for Biosemi bdf files
%
% modified 3/11/14 JD
% Added support for Neuromag fiff files
%
% modified 5/12/14 JD
% Changed mff to v2.
%
% modified 5/20/14 JD
% Added support for BrainVision EEG files
%
% bugfix 10/28/19 JD
% Fixed handling of BrainVision SEG and DAT files.
%
% modified 3/14/20 JD
% Added support for Psychophysiology Toolbox .mat files.
%
% modified 7/22/22 JD
% Added support for Matlab .mat files.
%
% modified 10/11/22 JD
% Added support for Neuroelectrics .easy files.
%
% modified 1/1/25 JD
% Added support for BESA .raw files.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileSuffix='';
formatName='';
formatCode='';

if nargin==0
    fileSuffix={'EP (.ept)';'EGI EGIS (.egis)';'EGI Simple Binary (.sbin)';'EEGlab (.set/.study)';'Biosig EDF (.edf)';'Neuroscan (.cnt/.eeg/.avg)';'text (.txt)';'EGI Matlab (.nsf)';'EGI MFF (.mff)';'ERPlab (.erp)';'Biosemi (.bdf)';'Neuromag (.fif)';'PTB (.mat)';'fieldtrip (.mat)';'BrainVision (.eeg/.dat/.seg)';'Stacked (.txt)';'Matlab (.mat)';'Neuroelectrics (.easy)'};
else
    switch formatType
        case 'EGI EGIS (.egis)'
            if strcmp(dataType,'single_trial')
                fileSuffix={'*.egis;*.raw;*.ses'};
                formatName='EGI EGIS session';
                formatCode='egi_egis';
            else
                fileSuffix={'*.egis;*.ave;*.gav;*.raw'};
                formatName='EGI EGIS average';
                formatCode='egi_egia';
            end
        case 'EGI Simple Binary (.sbin)'
            fileSuffix={'*.sbin';'*.raw'};
            formatName='EGI Simple Binary';
            formatCode='egi_sbin';
        case 'EEGlab (.set/.study)'
            fileSuffix={'*.set;*.study'};
            formatName='EEGlab set';
            formatCode='eeglab_set';
        case 'Biosig EDF (.edf)'
            fileSuffix={'*.edf'};
            formatName='Biosig EDF';
            formatCode='edf';
        case 'Neuroscan (.cnt/.eeg/.avg)'
            if strcmp(dataType,'continuous')
                fileSuffix={'*.cnt'};
                formatName='Neuroscan CNT';
                formatCode='ns_cnt';
            elseif strcmp(dataType,'single_trial')
                fileSuffix={'*.eeg'};
                formatName='Neuroscan EEG';
                formatCode='ns_eeg';
            else
                fileSuffix={'*.avg'};
                formatName='Neuroscan AVG';
                formatCode='ns_avg';
            end
        case 'text (.txt)'
            fileSuffix={'*.txt'};
            formatName='text';
            formatCode='text';
        case 'EP (.ept)'
            fileSuffix={'*.ept'};
            formatName='EP .ept';
            formatCode='ep_mat';
        case 'EGI Matlab (.nsf)'
            fileSuffix={'*.nsf';'*.mat'};
            formatName='NetStation Matlab';
            formatCode='ns_mat';
        case 'EGI MFF (.mff)'
            fileSuffix='*.mff';
            formatName='NetStation MFF';
            formatCode='egi_mff';
        case 'ERPlab (.erp)'
            fileSuffix={'*.erp'};
            formatName='ERPlab erp';
            formatCode='eeglab_erp';
        case 'fieldtrip (.mat)'
            fileSuffix={'*.mat'};
            formatName='fieldtrip';
            formatCode='fieldtrip';
        case 'Biosemi (.bdf)'
            fileSuffix={'*.bdf'};
            formatName='Biosemi bdf';
            formatCode='biosemi_bdf';
        case 'Neuromag (.fif)'
            fileSuffix={'*.fif'};
            formatName='Neuromag FIFF';
            formatCode='neuromag_fif';
        case 'BrainVision (.eeg/.dat/.seg)'
            if strcmp(dataType,'continuous')
                fileSuffix={'*.eeg'};
                formatName='BrainVision EEG';
                formatCode='brainvision_eeg';
            elseif strcmp(dataType,'single_trial')
                fileSuffix={'*.seg'};
                formatName='BrainVision SEG';
                formatCode='brainvision_seg';
            else
                fileSuffix={'*.dat'};
                formatName='BrainVision DAT';
                formatCode='brainvision_dat';
            end
        case 'PTB (.mat)'
            fileSuffix={'*.mat'};
            formatName='Psychophysiology Toolbox MAT';
            formatCode='PTB_mat';
        case 'Stacked (.txt)'
            fileSuffix={'*.txt'};
            formatName='Stacked Text';
            formatCode='stacked';
        case 'Matlab (.mat)'
            fileSuffix={'*.mat'};
            formatName='Matlab Matrix';
            formatCode='matlab_mat';
        case 'Neuroelectrics (.easy)'
            fileSuffix={'*.easy'};
            formatName='Neuroelectrics easy';
            formatCode='Neuroelectrics_easy';
        case 'BESA (.raw)'
            fileSuffix={'*.raw'};
            formatName='BESA raw';
            formatCode='besa_raw';
        otherwise
            fileSuffix='';
            formatName='';
            formatCode='';
            disp(['Oops... ' formatType ' file format not recognized.']);
    end
end



