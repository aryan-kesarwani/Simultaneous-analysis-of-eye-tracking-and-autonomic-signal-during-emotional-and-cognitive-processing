function errorFlag=ep_printADFResults(testResults, alpha, followup, outfid, ANOVAoptions)
%errorFlag=ep_printADFResults(testResults, alpha, followup, outfid, ANOVAoptions);
%Prints out statistical results with appropriate formatting.
%
%Inputs
%  testResults : Structured variable with the results of the statistical test
% 	.stat : the value of the test statistic
% 	.numDF : numerator degrees of freedom
% 	.denDF : denominator degrees of freedom  For robust statistic, includes ADF correction.
% 	.pValue : p-value
% 	.epsilon : epsilon correction value
% 	.MSE : mean squared error.  For robust statistic, includes ADF correction.
% 	.statsTest : description of test type.
% 	.statsType 'conventional' or 'robust'
% 	.effectSize : size of the effect size using the appropriate statistic.
%   .effectSizeTest :  type of effect size measure
%   .lowerConfid : lower bound of confidence interval (half alpha).
%   .upperConfid : upper bound of confidence interval (half alpha). 
%   .GPOWER : the effect size statistic needed by GPOWER to perform its calculations (the SPSS version).
%   .GPOWERtest : type of effect size measure for G*POWER analyses.
%   .meanR : mean correlation between within-group levels for Superpower power analyses.
%  alpha :   The alpha thresholds for significance
%    .uncorrected:  Alpha for a priori interest (.05)
%    .corrected:  Alpha for posthoc multiple comparison corrected.
%    .pvar:  twice standard deviation of p-values over reps
%  followup  : whether the output is part of follow-up statistics and hence in a smaller font.
%  outfid:  The fid of the output file.  Will print to screen if zero.
%  ANOVAoptions : preferences for the ANOVA function.
%    .epsilon  : type of epsilon correction for conventional within-group effects ('Greenhouse-Geisser','Huynd-Feldt','Lower Bound').
%Outputs
%   errorFlag : set to 1 if an error has occurred.
%  Prints out statistical results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%History
%  by Joseph Dien (11/8/08)
%  jdien07@mac.com
%
% modified 2/12/09 JD
% Correctly handles NaN results.
%
% modified 3/23/09 JD
% Highlights significant results.  Added option to print to file.
%
% modified 4/21/09 JD
% Green now corresponds to significant 1-tailed test (two times the uncorrected alpha).
%
% modified 7/26/09 JD
% Added indentation for follow-up ANOVAs.
%
% modified 9/20/15 JD
% Added effect sizes.
%
% bugfix 11/30/15 JD
% No longer produces effect sizes for contrasts of more than 1 df.
%
% modified 12/11/15 JD
% Checks if p-value variability exceeds two standard deviations over the threshold.
% p-values exactly equalling the threshold are accepted as significant.
% Prints out number of singular and nearly singular matrices during computation of statistics.
%
% modified 12/27/15 JD
% Changed message to indicate that effect size only available for
% between-group contrasts.
%
% bugfix 6/20/19 JD
% No longer marks NaN results in red.
%
% modified 7/6/19 JD
% Eliminated singularity counts in RESULTS.
%
% modified 7/26/19 JD
% Added MSE to the main results line to make it easier to cut-and-paste into manuscripts.
%
% modified 1/15/21 JD
% Collected test results into a more self-descriptive structured variable and added support for conventional ANOVAs.
%
% bugfix & modified 4/27/21 JD
% Prints last two non-zero digits when MSE is less than 1.
% Fixed p-values sometimes only having one significant digit.
%
% modified 9/27/21 JD
% Changed code to hopefully provide workaround for weird OS X bug with
% outputting with fscanf statement.
%
% bugfix 9/2/22 JD
% Deleted erroneous epsilon output for robust statistics of zero for compiled code and non-zero for Matlab code.
%
% bugfix 7/5/24 JD
% Fixed p-value sometimes displayed with more than two significant digits.
%
% bugfix 8/9/24 JD
% Again fixed p-value sometimes displayed with more than two significant digits.
%
% modified 8/17/24 JD
% Added mean R to output for use with power calculations.
% Added bold for non-significant effects with confidence interval boundaries less than absolute .2.
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

errorFlag=0;

if ~exist('outfid','var')
    outfid =0;
end

%Determine two digits of p-value
P=sprintf('%10.8f', testResults.pValue);
Pp=findstr('.',P);
Pz=findstr(P(Pp+1:length(P)),'0'); %list of zeros after the decimal
if isempty(Pz) || Pz(1)>1
    iSigDigits=1;
else
    for iSigDigits=1:length(Pz)
        if iSigDigits ~= Pz(iSigDigits) %find first non-consecutive zero after decimal
            break
        end
    end
    if (iSigDigits==length(Pz)) && (iSigDigits == Pz(iSigDigits)) %if all the zeroes were consecutive
        iSigDigits=iSigDigits+1;
    end
end

if isnan(testResults.pValue) || isnan(testResults.stat) || isinf(testResults.denDF)
    pSign='=';
    pVal=' could not be calculated';
elseif testResults.pValue >= .00000001
    pSign='=';
    pVal=sprintf(['%1.' num2str(iSigDigits+1) 'f'], testResults.pValue);
else
    pSign='<';
    pVal='0.00000001';
end

%Determine two digits of MSE
if testResults.MSE < 1
    M=sprintf('%10.8f', testResults.MSE);
    Mp=findstr('.',M);
    Mz=findstr(M(Mp+1:length(M)),'0');
    for iSigDigits=1:length(Mz)
        if iSigDigits ~= Mz(iSigDigits)
            break
        end
    end
    if iSigDigits==length(Mz)
        iSigDigits=iSigDigits+1;
    end
    if isempty(iSigDigits)
        iSigDigits=1;
    end
    
    if isnan(testResults.MSE) || isinf(testResults.MSE)
        mSign='=';
        mVal=' could not be calculated';
    elseif testResults.MSE >= .00000001
        mSign='=';
        mVal=sprintf(['%1.' num2str(iSigDigits+1) 'f'], testResults.MSE);
    else
        mSign='<';
        mVal='0.00000001';
    end
else
    mSign='=';
    mVal=sprintf(['%1.2f'], testResults.MSE);
end

if testResults.pValue <= alpha.corrected
    if (testResults.pValue+alpha.pvar) > alpha.corrected
        pVal=[pVal '(not confirmed at the corrected alpha level, consider increasing number of boostrap samples to achieve more stable results)'];
    end
elseif testResults.pValue <= alpha.uncorrected
    if (testResults.pValue+alpha.pvar) > alpha.uncorrected
        pVal=[pVal '(not confirmed at the uncorrected alpha level, consider increasing number of boostrap samples to achieve more stable results)'];
    end
end

effsz='';
gpower='';
if ~isnan(testResults.effectSize)
    if strcmp(testResults.statsType,'robust')
        if isnan(testResults.lowerConfid)
            effsz=sprintf(', %s=%1.2f', testResults.effectSizeTest, testResults.effectSize);
        else
            effsz=sprintf(', %s=%1.2f (CI<sub>95</sub>=%1.2f to %1.2f)', testResults.effectSizeTest, testResults.effectSize, testResults.lowerConfid, testResults.upperConfid);
        end
    else
        switch testResults.statsTest
            case 'paired samples t-test'
                switch testResults.effectSizeTest
                    case 'Hedges grm'
                        effsz=sprintf(', g<sub>rm</sub>=%1.2f (CI<sub>95</sub>=%1.2f to %1.2f)', testResults.effectSize, testResults.lowerConfid, testResults.upperConfid);
                    case 'Hedges gav'
                        effsz=sprintf(', g<sub>av</sub>=%1.2f (CI<sub>95</sub>=%1.2f to %1.2f)', testResults.effectSize, testResults.lowerConfid, testResults.upperConfid);
                    otherwise
                        error('oops - programmer error')
                end
                gpower=sprintf('  G*POWER: d<sub>z</sub>=%1.2f', testResults.GPOWER);
            case 'independent samples t-test'
                effsz=sprintf(', g<sub>s</sub>==%1.2f (CI<sub>95</sub>=%1.2f to %1.2f)', testResults.effectSize, testResults.lowerConfid, testResults.upperConfid);
                gpower=sprintf('  G*POWER: d<sub>s</sub>=%1.2f.', testResults.GPOWER);
            case 'N-way'
                effsz=sprintf(', ɳ<sub>p</sub><sup>2</sup>=%1.2f.', testResults.effectSize);
%                 if length(testResults.GPOWER)  > 1
%                     gpower=sprintf('G*POWER: f=%1.2f, Corr among rep measures=%1.2f', testResults.GPOWER(1), testResults.GPOWER(2));
%                 else
                gpower=sprintf('  G*POWER: f=%1.2f.', testResults.GPOWER);
%                 end
            otherwise
                error('oops - programmer error')
        end
    end
end

meanR='';
if ~isempty(testResults.meanR)
    meanR=sprintf('  Mean R=%1.2f.', testResults.meanR);
end

if strcmp(testResults.statsType,'robust')
    if outfid
        theStat='T<sub>WJt</sub>/c';
    else
        theStat='TWJt/c';
    end
%     theEpsilon=sprintf(', e=%1.2f', testResults.epsilon);
    theEpsilon='';
    df1='%1.1f';
    df2='%1.1f';
else
    if (testResults.numDF==1) && followup
        theStat='t';
        theEpsilon='';
        df1='';
        df2='%1.0f';
    else
        theStat='F';
        if isfield(testResults,'epsilon')
            switch ANOVAoptions.epsilon
                case 'Greenhouse-Geisser'
                    theEpsilon=sprintf(', ê=%1.2f', testResults.epsilon);
                case 'Huynd-Feldt'
                    theEpsilon=sprintf(', ẽ=%1.2f', testResults.epsilon);
                case 'Lower Bound'
                    theEpsilon=sprintf(', e=%1.2f', testResults.epsilon);
                otherwise
                    error('Oops programmer error')
            end
        else
            theEpsilon='';
        end
        df1='%1.0f';
        df2='%1.0f';
    end
end

% theEpsilon='';
% if strcmp(ANOVAoptions.adf,'ADF')
%     df1='%1.1f';
%     df2='%1.1f';
%     if ~statsType
%         testResults.numDF=testResults.numDF*testResults.effectSize;
%         testResults.denDF=testResults.denDF*testResults.effectSize;
%     end
% else
%     df1='%1.0f';
%     df2='%1.0f';
%     if statsType
%         theEpsilon=sprintf(', e=%1.2f', testResults.effectSize);
%         testResults.denDF=testResults.denDF/testResults.effectSize;
%     else
%         switch epsilonType
%             case 'Greenhouse-Geisser'
%                 theEpsilon=sprintf(', ê=%1.2f', testResults.effectSize);
%             case 'Huynd-Feldt'
%                 theEpsilon=sprintf(', ẽ=%1.2f', testResults.effectSize);
%             case 'Lower Bound'
%                 theEpsilon=sprintf(', e=%1.2f', testResults.effectSize);
%             otherwise
%                 error('Oops programmer error')
%         end
%     end
% end

if strcmp(theStat,'t')
    if outfid
        %print to file
        if testResults.pValue <= alpha.corrected
            theText=sprintf([repmat('&nbsp;',1,followup*2) '<font color="red">%s(' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s</font>'], theStat,testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal, gpower, meanR);
        elseif testResults.pValue <= alpha.uncorrected
            theText=sprintf([repmat('&nbsp;',1,followup*2) '<font color="orange">%s(' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s</font>'], theStat,testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal, gpower, meanR);
        elseif testResults.pValue <= 2*alpha.uncorrected
            theText=sprintf([repmat('&nbsp;',1,followup*2) '<font color="green">%s(' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s</font>'], theStat,testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal, gpower, meanR);
        elseif ~isnan(testResults.lowerConfid) && ((testResults.lowerConfid >-.2) && (testResults.upperConfid <.2))
            theText=sprintf([repmat('&nbsp;',1,followup*2) '<b>%s(' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s</b>'], theStat,testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal, gpower, meanR);
        else
            theText=sprintf([repmat('&nbsp;',1,followup*2) '%s(' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s'], theStat,testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal, gpower, meanR);
        end
        fprintf(outfid,[theText '</BR>\n']);
    else
        %print to screen
        if testResults.pValue <= alpha.corrected
            fprintf([repmat('&nbsp;',1,followup*2) '**%s(' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s**\n'], theStat,testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal, gpower, meanR);
        elseif testResults.pValue <= alpha.uncorrected
            fprintf([repmat('&nbsp;',1,followup*2) '*%s(' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s*\n'], theStat,testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal, gpower, meanR);
        else
            fprintf([repmat('&nbsp;',1,followup*2) '%s(' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s\n'], theStat,testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal, gpower, meanR);
        end
    end
else
    if outfid
        %print to file
        if testResults.pValue <= alpha.corrected
            theText=sprintf([repmat('&nbsp;',1,followup*2) '<font color="red">%s(' df1 ',' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s</font>'], theStat,testResults.numDF, testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal, gpower, meanR);
        elseif testResults.pValue <= alpha.uncorrected
            theText=sprintf([repmat('&nbsp;',1,followup*2) '<font color="orange">%s(' df1 ',' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s</font>'], theStat,testResults.numDF, testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal, gpower, meanR);
        elseif testResults.pValue <= 2*alpha.uncorrected
            theText=sprintf([repmat('&nbsp;',1,followup*2) '<font color="green">%s(' df1 ',' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s</font>'], theStat,testResults.numDF, testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal, gpower, meanR);
        elseif ~isnan(testResults.lowerConfid) && ((testResults.lowerConfid >-.2) && (testResults.upperConfid <.2))
            theText=sprintf([repmat('&nbsp;',1,followup*2) '<b>%s(' df1 ',' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s</b>'], theStat,testResults.numDF, testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal, gpower, meanR);
        else
            theText=sprintf([repmat('&nbsp;',1,followup*2) '%s(' df1 ',' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s'], theStat,testResults.numDF, testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal, gpower, meanR);
        end
        fprintf(outfid,[theText '</BR>\n']);
    else
        %print to screen
        if testResults.pValue <= alpha.corrected
            fprintf([repmat('&nbsp;',1,followup*2) '**%s(' df1 ',' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s**'], theStat,testResults.numDF, testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal);
        elseif testResults.pValue <= alpha.uncorrected
            fprintf([repmat('&nbsp;',1,followup*2) '*%s(' df1 ',' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s*'], theStat,testResults.numDF, testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal);
        else
            fprintf([repmat('&nbsp;',1,followup*2) '%s(' df1 ',' df2 ')=%1.2f, p%s%s%s%s, MSe%s%s, %s%s'], theStat,testResults.numDF, testResults.denDF, testResults.stat, pSign, pVal, theEpsilon, effsz, mSign, mVal);
        end
        fprintf(outfid,'</BR>');
    end
end
