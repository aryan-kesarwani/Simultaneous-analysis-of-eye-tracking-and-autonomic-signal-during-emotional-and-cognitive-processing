function [experimentFieldNames, experimentData, theHeaders, trialData] = ep_loadSpecFile(specFormat, specFileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [experimentFieldNames, experimentData, theHeaders, trialData] = ep_loadSpecFile(specFormat, specFileName)
% Loads in trial specs from either a regular tabulated text file with column headers or E-Prime .txt file.
%
%Inputs
%  specFileName           : The text file with the trial spec info in it. 
%  specFormat             : The type of trial spec file (EPM=E-Prime, TSP=tab-delimited text file)
%Outputs
%  experimentFieldNames   : The experiment information field names
%  experimentData         : The experiment information
%  fieldNames             : The trial information field names
%  trialData              : The trial data
%
% History:
%
% by Joseph Dien (10/9/18)
% jdien07@mac.com
%
% modified 11/29/19 JD
% Enabled reading of text files with a greater range of variations in character encoding, end-of-line markers, and field separation markers.
%
% bugfix 12/6/21 JD
% Fixed crash when reading tab-text spec files.
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

global EPtictoc

experimentFieldNames=cell(0);
experimentData=cell(0);
theHeaders=cell(0);
trialData=cell(0);

if isempty(specFileName)
    [FileName,PathName,FilterIndex] = uigetfile(['*.txt'],'Load Trial Specs Text File');
    if isnumeric(FileName)
        if FileName == 0
            msg{1}='No file name specified.';
            [msg]=ep_errorMsg(msg);
            return
        end;
    end;
else
    [PathName, FileName, ext] = fileparts(specFileName);
    FileName=[FileName ext];
end

switch specFormat
    case 'EPM' % E-Prime text file
        [experimentFieldNames, experimentData, theHeaders, trialData]=ep_readEprime([PathName filesep FileName]);
    case 'TSP' %tab-delimited text file
        [theHeader, trialData, theDelim] = ep_textScan([PathName filesep FileName],2);
        ep_tictoc;if EPtictoc.stop;EPtictoc.stop=0;ep('start');return;end
        if ~isempty(trialData)
            numSpecs=size(trialData,2);
            if isempty(theHeader)
                msg{1}='No headers.';
                [msg]=ep_errorMsg(msg);
                return
            end
            theHeaders=theHeader{1};
            if length(theHeaders) ~= numSpecs
                msg{1}='Number of column headers is different from the number of data columns.';
                [msg]=ep_errorMsg(msg);
                return
            end
        end
end

