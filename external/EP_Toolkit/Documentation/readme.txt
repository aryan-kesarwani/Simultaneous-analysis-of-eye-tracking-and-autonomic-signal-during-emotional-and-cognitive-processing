ERP PCA Toolkit


Copyright (C) 1999-2025  Joseph Dien

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
 
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  Programs are not for medical use or any other application where use or misuse could cause injury.  See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

The rotations other than Varimax, Promax, and Infomax were based on code adapted from Matlab code made publicly available by Coen A. Bernaards and Robert I. Jennrich with their permission.
Website: http://www.stat.ucla.edu/research

The files ave_hdr_offsets_v.m, ave_hdr_offsets_v.m, ses_hdr_offsets_v.m, rd_onetr_allch.m, get_cell_offsets.m, get_fid, and wt_ses_hdr_v2.m are from EGI's EGI Toolbox with permission.  The files rd_PCAegis_hdr_v.m and wt_PCAave_hdr_v.m were modified from this same source.

The files el_topoplot.m, el_traditionaldipfit.m, el_sobi.m, and el_runica.m in the EEGlab folder in the External folder are renamed topoplot.m, traditionaldipfit.m, sobi.m, and runica.m files from EEGlab.  Changes were to the names of the functions, the commenting out of one line and addition of ep_tictoc call in runica, a change from "icadefs" to "evalc('icadefs') in el_topoplot.m, and commenting out of mean correction and whitening and code to accommodate NaN values, and are provided under the GPL license.  This paragraph constitutes the required declaration that it has been modified from the original.  The topoplot modification was made so that calls could be made to the function without confusion occurring as to whether EEGlab's topoplot was intended or FieldTrip's topoplot was intended and to silence unneeded text output from calls to icadefs.  The traditionaldipfit modification was made so that the function would be available for use since it is not visible to outside programs under EEGlab’s standard installation procedure.  The commented out line in runica was performed so that output from repeated applications to a given dataset will fully replicate.  The changes to sobi.m were made to allow for whitening and mean correction to be made in the calling function and the modified code were to allow for NaN values during calculation of the correlation matrices and were marked with a %JD marker. Additionally, ep_headplot.m is a modified version of headplot.m (submitted to the EEGlab group as a revision) which allows for the electrode labels to be individually colored.  The changed lines are marked with a %JD marker.  The file avg152t1.mat is from the diplot plugin and has been unchanged.

The file fmrib_fastr.m in the fmrib1.21 directory in the externals folder has two changes made to it to fix crashing errors.  The changes are documented in the code and included with the EP Toolkit, as permitted under GPL License.

The file amri_eeg_gac.m in the amri_eegfmri_toolbox directory in the externals folder has one set of changes made to it to fix crashing errors.  The changes are documented in the code and included with the EP Toolkit, as permitted under GPL License.

The file ft_sphericalSplineInterpolate.m in the fieldtrip folder in the External folder is unchanged from the fieldtrip-20140807 distribution and is included, as permitted under GPL license, as it is otherwise not directly accessible to EP Toolkit functions.

The file ft_sphericalSplineInterpolate.m in the fieldtrip folder in the External folder is a renamed sphericalSplineInterpolate.m file from the fieldtrip-20140807 distribution.  Only the names of the function was changed and is provided under the GPL license.  This paragraph constitutes the required declaration that it has been modified from the original.  The modification was made so that calls could be made to the function without confusion with the version accompanying fieldtrip.  This was done so that the function would be available for use since it is not visible to outside programs under fieldtrip’s standard installation procedure.

The files in the CCA directory in the externals folder are being included by permission of Maartens De Vos.  The compute_cca.m file has been modified to allow for options parameters to be passed on by the function call.

The files in the invcwt_v1 directory in the externals folder are being redistributed per permissions described in its license.txt file.

The files in the RIDE directory in the externals folder are being redistributed as under GPL License.  The examples directory has been deleted for space reasons.

The files in the dprime directory in the externals folder are being redistributed per permissions described in its license.txt file.

The files in the robertpetermatthew-f_ICC-6bc742d directory in the externals folder are being redistributed per permissions described in its LICENSE file.

The Standard-10-5-Cap385-VEOG.ced file is based on the Standard-10-5-Cap385.spf included with EEGlab and credited to Robert Oostenveld.  This present file has two lower VEOG channels added, using manually generated coordinates.

The files in the BOSC directory in the externals folder are being redistributed as under GPL License.

The files in the eBOSC-master directory in the externals folder are being redistributed as under GPL License.

The files in the fulepo-NIPALS-ac91fbc directory in the externals folder are being redistributed per permissions described in its LICENSE file.

The files in the MFFMatlabIO3.5 directory in the externals folder are being redistributed as under GPL License.

The files in the ReSync1.0 directory in the externals folder are being included by permission of Guang Ouyang.

The files in the FastICA_25 directory in the externals folder are being redistributed as under GPL License, as stated on http://research.ics.aalto.fi/ica/fastica/about.shtml: "The software package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version."

The files in the AMICA1.5.1 directory in the externals folder are being redistributed as under GPL License.

The files in the PrepPipeline0.55.3 directory in the externals folder are being redistributed as under GPL License.

The robust effect size functions in ep_WJGLMml_mat were translated from Rand R. Wilcox's WRS R library with his permission and under the USC-RL v1.0 license that allows for academic and non-commercial usage.

