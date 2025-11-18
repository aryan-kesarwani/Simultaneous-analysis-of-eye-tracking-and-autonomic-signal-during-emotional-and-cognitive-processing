function [theHeader, theData, theDelim] = ep_textScan(fileName,firstRow,lastRow,firstCol,lastCol,theDelim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [theHeader, theData, theDelim] = ep_textScan(fileName,firstRow,lastRow,firstCol,lastCol,theDelim)
%
%	Replaces textscan function by also supporting header lines, trimming of the data structure, Unicode text files, different end-of-line encodings, and endian machine coding.
%
%Inputs
%	fileName:   filename and sourcepath.
%   firstRow:   first row of the data.  Any rows before this will be returned in theHeader.
%   lastRow:    last row of the data.  Any rows after this will be discarded.  Zero means up to the very last row.
%   firstCol:   first column of the data.  Any columns before this will be discarded.
%   lastCol:    last column of the data.  Any columns after this will be discarded.  Zero means up to the very last column.
%   theDelim:   the text delimiter.  Optional.  Will guess if not provided.
%
%Outputs
%   theHeader: cell array with cell arrays of the header labels, if any.  Else empty.
%   theData:   the data, organized as a 2D cell array.
%   theDelim:     the delimiter
%
% Re: Unicode byte order marker:
% https://en.wikipedia.org/wiki/Byte_order_mark#UTF-8 11/29/2019
%    The UTF-8 representation of the BOM is the (hexadecimal) byte sequence 0xEF,0xBB,0xBF.
%    In UTF-16, a BOM (U+FEFF) may be placed as the first character of a file or character stream to indicate the endianness (byte order) of all the 16-bit code units of the file or stream. 
%    If an attempt is made to read this stream with the wrong endianness, the bytes will be swapped, thus delivering the character U+FFFE, which is defined by Unicode as a "non character" that should never appear in the text.
% 
%        If the 16-bit units are represented in big-endian byte order, the BOM will appear in the sequence of bytes as 0xFE 0xFF
%        If the 16-bit units use little-endian order, the BOM will appear in the sequence of bytes as 0xFF 0xFE
%
%  UTF-8: 239 187 191
%  UTF-16BE 254 255
%  UTF-16LE 255 254
%
%  Empty data rows, as in excess carriage returns at the end of the file, are ignored.
%  Also, extra delimiters at the end of a row are ignored.

% History:
%
% by Joseph Dien (1/9/20)
% jdien07@mac.com
%
% bugfix 1/20/20 JD
% If text file has an extra delimiter at the start of each row, ignore it.
%
% bugfix 2/14/20 JD
% Fixed crash if firstRow not specified.
%
% bugfix 4/17/20 JD
% Fixed skipping first few characters for UTF-8 text files due to Matlab already dropping the byte order marker when using fgetl.
%
% bugfix 7/11/20 JD
% Fixed not determining number of columns correctly when empty fields at the end of the first data row were not padded out with field delimiters.
% Fixed treating tab and comma delimiters as a single delimiter when repeated due to empty fields.
% Fixed error when there is a final end-of-line character after the final row's end-of-line character.
% Fixed dropping final field if empty by checking against number of fields indicated by the header length.
% Fixed not reading text files correctly when some fields are empty.
% Added option to specify the delimiter.
%
% bugfix 7/29/21 JD
% Fixed not handling blank fields (consecutive delimiters) when reading comma-delimited text files.
% Fixed not dropping the full byte mark code at the start of UTF-8 unicode files.
%
% bugfix 4/29/22 JD
% Fixed not handling space padded numbers.
%
% modified 6/20/22 JD
% Now drops empty rows in data rows, such as when there are excess carriage returns at the end of the file.
%
% bugfix 7/27/22 JD
% Fixed incorrect determination of delimiter when there is a header but no data.
%
% bugfix 5/19/24 JD
% For comma-delimited files, drop single quotes bracketing strings.
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

global EPtictoc

theData=cell(0);
theHeader=cell(0);
if ~exist('theDelim','var')
    theDelim='';
end

if ~exist('firstRow','var')
    firstRow=1;
end
if ~exist('lastRow','var')
    lastRow=0;
end
if ~exist('firstCol','var')
    firstCol=0;
end
if ~exist('lastCol','var')
    lastCol=0;
end

if firstCol==0
    firstCol=1;
end

%determine if there is a unicode byte order marker
fid=fopen(fileName);
if fid==-1
    msg{1}='Problem opening the file.';
    [msg]=ep_errorMsg(msg);
    return
end
tempVar=fread(fid);
if (double(tempVar(1))==239) && (double(tempVar(2))==187) && (double(tempVar(3))==191)
    textType='UTF-8';
elseif (double(tempVar(1))==254) && (double(tempVar(2))==255)
    textType='UTF-16BE';
elseif (double(tempVar(1))==255) && (double(tempVar(2))==254)
    textType='UTF-16LE';
else
    textType='';
end

%determine end-of-line encoding
numLF=length(find(tempVar==10));
numCR=length(find(tempVar==13));
if numLF==numCR
    lineType='windows';
    numRows=numLF;
elseif numLF>numCR
    lineType='unix';
    numRows=numLF;
elseif numLF<numCR
    lineType='OS9';
    numRows=numCR;
else
    msg{1}='Unable to determine end-of-line type.';
    [msg]=ep_errorMsg(msg);
    return
end
if ~any(ismember([10 13],tempVar(end)))
    numRows=numRows+1; %if final row is not terminated with an end-of-line character
end
if isequal(tempVar(end-1:end),[10; 10]) || isequal(tempVar(end-1:end),[13; 13])
    numRows=numRows-1; %if there is an extra end-of-line character after the final row
end
fclose(fid);

%determine the data delimiter
switch textType
    case 'UTF-8'
        fid=fopen(fileName,'r','n','UTF-8');
    case 'UTF-16BE'
        fid=fopen(fileName,'r','b','UTF-8');
    case 'UTF-16LE'
        fid=fopen(fileName,'r','l','UTF-8');
    otherwise
        fid=fopen(fileName,'r','n');
end
if fid==-1
    msg{1}='Problem opening the file.';
    [msg]=ep_errorMsg(msg);
    return
end

if firstRow>1
    for iRow=1:firstRow
        tempVar=EPfgetl(fid,lineType,textType);
        if iRow==(firstRow-1)
            headerVar=tempVar;
        end
    end
else
    headerVar=cell(0);
    tempVar=fgetl(fid);
    %drop the byte order marker
    switch textType
        case 'UTF-8'
            %drops the bytemark when opened with the option.
        case 'UTF-16BE'
            %drops the bytemark when opened with the option.
        case 'UTF-16LE'
            %drops the bytemark when opened with the option.
        otherwise
    end
end

if isempty(theDelim)
    if isempty(tempVar) %case where there is a header but no data
        tempVar2=headerVar;
    else
        tempVar2=tempVar;
    end
    if any(strfind(tempVar2,char(9))) %generally a file with any tabs at all is tab-delimited
        theDelim=char(9);
    else
        whiteSpaceCount=0; %number of consecutive series of spaces so they are treated as a single space each
        whiteSpaceFlag=0;
        for iChar=1:length(tempVar2)
            if strcmp(tempVar2(iChar),' ')
                if ~whiteSpaceFlag %if the last character was not a space, then increment space counter
                    whiteSpaceFlag=1;
                    whiteSpaceCount=whiteSpaceCount+1;
                else
                    whiteSpaceFlag=0;
                end
            end
        end
        if length(strfind(tempVar2,',')) >= whiteSpaceCount
            theDelim=','; %if there are more or equal commas than spaces, assume it is a comma-delimited file.
        else
            theDelim=' '; %if there are more spaces than commas, assume it is a space-delimited file.
        end
    end
end

numcols=0;
if firstRow>1 %if there is a header, then it provides the most reliable indicator of the number of columns
        if strcmp(theDelim,' ')
            rowcols=length(regexp(headerVar,[theDelim '+']))+1; %determine number of data columns based on number of delimiters (treating repeated spaces as a single delimiter)
        else
            rowcols=length(regexp(headerVar,theDelim))+1; %determine number of data columns based on number of delimiters (treating repeats as a single delimiter)
        end
        if regexp(headerVar,[theDelim '$'])
            rowcols=rowcols-1; %if there is an extra delimiter at the end of the line, drop it.  This is especially a problem with EEGlab ced files, which have an excess delimiter at the end of every row for some reason.
        end
%         if regexp(headerVar,['^' theDelim])
%             rowcols=rowcols-1; %if there is an extra delimiter at the start of the line, drop it.  The strtok command will ignore the first delimiter.
%         end
        numcols=rowcols;
else
    endFlag=0;
    while ~endFlag
        if strcmp(theDelim,' ')
            rowcols=length(regexp(tempVar,[theDelim '+']))+1; %determine number of data columns based on number of delimiters (treating repeated spaces as a single delimiter)
        else
            rowcols=length(regexp(tempVar,theDelim))+1; %determine number of data columns based on number of delimiters (treating repeats as a single delimiter)
        end
        if regexp(tempVar,[theDelim '$'])
            rowcols=rowcols-1; %if there is an extra delimiter at the end of the line, drop it.  For text files with no header, this is more likely to denote an excess delimiter than an empty valid field.
        end
%         if regexp(tempVar,['^' theDelim])
%             rowcols=rowcols-1; %if there is an extra delimiter at the start of the line, drop it.  The strtok command will ignore the first delimiter.
%         end
        numcols=max(numcols,rowcols);
        tempVar=fgetl(fid);
        endFlag=feof(fid);
    end
end

%read the data
if lastCol ==0
    lastCol=numcols;
end

%read the header lines
frewind(fid);
%drop the byte order marker
switch textType
    case 'UTF-8'
        theByte=fread(fid,3,'uint8');
    case 'UTF-16BE'
        theByte=fread(fid,2,'uint8');
    case 'UTF-16LE'
        theByte=fread(fid,2,'uint8');
    otherwise
end

%read the header
if firstRow > 1
    theHeader=cell(firstRow-1,1);
    for iRow=1:firstRow-1
        theHeader{iRow}=cell(0);
        tempVar=EPfgetl(fid,lineType,textType);
        while regexp(tempVar,[theDelim '$'])
            tempVar=tempVar(1:end-1); %remove trailing delimiters
        end
        if strcmp(theDelim,char(44)) %if comma-delimited, then need to also recognize quotes for fields with commas
            fieldCounter=0;
            quoteMode=0;
            fieldFlag=0;
            for iChar=1:length(tempVar)
                theChar=tempVar(iChar);
                if ~fieldFlag
                    if double(theChar)==34 %quotes
                        quoteMode=1; %ignore commas till the next quote indicates the end of the field
                    end
                    fieldFlag=1; %now entering a field
                    fieldCounter=fieldCounter+1;
                    theHeader{iRow}{fieldCounter}=theChar;
                else
                    if quoteMode && (double(theChar)==34) %quotes
                        quoteMode=0;
                    elseif double(theChar)==44 %comma
                        fieldFlag=0;
                    else
                        theHeader{iRow}{fieldCounter}(end+1)=theChar;
                    end
                end
            end
            if fieldCounter==numcols
                theHeader{iRow}=theHeader{iRow}(firstCol:lastCol); %if number of header columns equals number of data columns, apply column trimming.
            end
        else
            numHeaderCols=length(strfind(tempVar,theDelim))+1;
            theHeader{iRow}=cell(numHeaderCols,1);
            if numHeaderCols ==1
                theHeader{iRow}{1} = tempVar;
            else
                for iCol = 1:numHeaderCols
                    [theHeader{iRow}{iCol}, tempVar] = strtok(tempVar,theDelim);
                end
                if numHeaderCols==numcols
                    theHeader{iRow}=theHeader{iRow}(firstCol:lastCol); %if number of header columns equals number of data columns, apply column trimming.
                end
            end
        end
    end
    
    if numcols < length(theHeader{1})
        numcols=length(theHeader{1});
    end
end

%read the data
theData=cell(numRows-firstRow+1,numcols);
emptyRow=zeros(size(theData,1),1);
for iRow=1:numRows-firstRow+1
    ep_tictoc;if EPtictoc.stop;return;end
    tempVar=EPfgetl(fid,lineType,textType);
    if isempty(tempVar)
        emptyRow(iRow)=1;
    else
        if strcmp(theDelim,char(44)) %if comma-delimited, then need to also recognize quotes for fields with commas
            fieldCounter=0;
            quoteMode=0;
            fieldFlag=0;
            for iChar=1:length(tempVar)
                theChar=tempVar(iChar);
                if ~fieldFlag
                    if double(theChar)==34 %quotes
                        quoteMode=1; %ignore commas till the next quote indicates the end of the field
                    end
                    if double(theChar)==44 %comma so double delimiters
                        fieldFlag=0; %still just finishing a field
                        fieldCounter=fieldCounter+1;
                        theData{iRow,fieldCounter}=''; %empty field
                    else
                        if (fieldCounter > 0) && (length(theData{iRow,fieldCounter})>1) && (double(theData{iRow,fieldCounter}(1))==39) && (double(theData{iRow,fieldCounter}(end))==39)
                            theData{iRow,fieldCounter}=theData{iRow,fieldCounter}(2:end-1); %drop single-quotes bracketing a string
                        end
                        fieldFlag=1; %now entering a field
                        fieldCounter=fieldCounter+1;
                        theData{iRow,fieldCounter}=theChar;
                    end
                else
                    if quoteMode && (double(theChar)==34) %quotes
                        quoteMode=0;
                    elseif double(theChar)==44 %comma
                        fieldFlag=0;
                    else
                        theData{iRow,fieldCounter}(end+1)=theChar;
                    end
                end
            end
            if double(tempVar(end))==44
                fieldCounter=fieldCounter+1; %if ended on a comma delimliter, add one more field to the count
            end
            if fieldCounter < numcols %extra columns like for subject names is fine
                msg{1}=['Line number ' num2str(iRow) ' did not have the expected number of columns of data.'];
                [msg]=ep_errorMsg(msg);
                return
            end
        else
            numDataCols=length(strfind(tempVar,theDelim))+1;
            if ~isempty(regexp(tempVar,[theDelim '$'])) && (~isempty(theHeader) && (numDataCols ~=length(theHeader{1})))
                numDataCols=numDataCols-1; %if there is an extra delimiter at the end of the line, drop it.
            end
            %         if numDataCols < numcols %extra columns like for subject names is fine
            %             msg{1}=['Line number ' num2str(iRow) ' did not have the expected number of columns of data.'];
            %             [msg]=ep_errorMsg(msg);
            %             return
            %         end
            rowData=split(tempVar,theDelim);
            rowData=rowData(1:numDataCols);
            theData(iRow,1:numDataCols) = rowData;
        end
    end
end
%theData=textscan(fid, repmat('%s',1,numcols),'Delimiter',theDelim, 'MultipleDelimsAsOne', 1,'EndOfLine','\r\n');
fclose(fid);

%drop empty rows
if any(emptyRow)
    disp(['Dropping ' num2str(length(find(emptyRow))) ' empty rows.']);
    theData(find(emptyRow),:)=[];
    numRows=size(theData,1);
end

%drop excess rows
if lastRow ==0
    lastRow=numRows;
end
finalRow=lastRow-firstRow+1;
if finalRow > size(theData,1)
    msg{1}=['Specified final row ' num2str(finalRow) ' was larger than the number of rows of data.'];
    [msg]=ep_errorMsg(msg);
    return
elseif finalRow < size(theData,1)
    theData=theData(1:finalRow,:);
end

%drop excess columns
theData=theData(:,firstCol:lastCol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function theLine=EPfgetl(fid,lineType,textType)
%get one line from the file, taking into account the end-of-line and character encoding schemes.
switch textType
    case 'UTF-16BE'
        outBytes=[];
        eol=0;
        while ~feof(fid) && ~eol
            theByte=fread(fid,1,'uint8');
            if any(ismember([10 13],theByte))
                eol=1;
            else
                outBytes(end+1)=theByte;
            end
        end
        outBytes=outBytes(1:end-1); %discard null byte from the terminator
        if strcmp(lineType,'windows')
            tempVar=fread(fid,2); %discard the carriage return and its null byte
        end
        theLine = native2unicode(outBytes, 'UTF-16BE');
    case 'UTF-16LE'
        outBytes=[];
        eol=0;
        while ~feof(fid) && ~eol
            theByte=fread(fid,1,'uint8');
            if any(ismember([10 13],theByte))
                eol=1;
            else
                outBytes(end+1)=theByte;
            end
        end
        if strcmp(lineType,'windows')
            tempVar=fread(fid,3); %discard remaining null byte from the terminator plus the carriage return and its null byte
        else
            tempVar=fread(fid,1); %discard remaining null byte from the terminator
        end
        theLine = native2unicode(outBytes, 'UTF-16LE');
    otherwise
        theLine=fgetl(fid);
end

