function [MUHAT, STDIZER, RESULTS, MSE]=ep_WJGLMml_mat(Y, NX, C, U, OPT1, PER, OPT2, NUMSIM, SEED, MISSING, OPT3, ALPHA, SCALE, LOC1, LOC2)
%function [MUHAT, STDIZER, RESULTS, MSE]=ep_WJGLMml_mat(Y, NX, C, U, OPT1, PER, OPT2, NUMSIM, SEED, MISSING, OPT3, ALPHA, SCALE, LOC1, LOC2);
%Function for calculating robust statistics using Welch-James ADF, boostrapping, and trimmed
%means plus winsorized variances and covariances.
%assumes there is only one random factor of subjects.
%
%Based on SAS/IML code made available by Lisa Lix at:
%http://www.usaskhealthdatalab.ca/sas-programs/
%Robust effect size functions based on R code made available by Rand R. Wilcox at:
%https://osf.io/xhe8u/
%
%Keselman, H. J., Wilcox, R. R., & Lix, L. M. (2003). A generally robust approach
%to hypothesis testing in independent and correlated groups designs
%Psychophysiology, 40, 586-596.
%
%Keselman, H. J., Algina, J., Lix, L. M., Wilcox, R. R., & Deering, K. N. (2008). 
%A generally robust approach for testing hypotheses and setting confidence intervals for effect sizes. 
%Psychological Methods, 13(2), 110.
%
%Wilcox, R. R. (2001). Fundamentals of modern statistical methods. New York: Springer-Verlag.
%
%Wilcox, R. R. (2022).  Introduction to Robust Estimation and Hypothesis Testing, 5th Edition.  London: Academic Press.
%Algina, J., Keselman, H. J., & Penfield, R. D. (2005a). An Alternative to Cohen’s Standardized Mean Difference Effect Size: A Robust Parameter and Confidence Interval in the Two Independent Groups Case. Psychological Methods, 10(3), 317–328. https://doi.org/10.1037/1082-989X.10.3.317
%Algina, J., Keselman, H. J., & Penfield, R. D. (2005b). Effect Sizes and their Intervals: The Two-Level Repeated Measures Case. Educational and Psychological Measurement, 65(2), 241–258. https://doi.org/10.1177/0013164404268675
%
%The MULTP/cterm correction to the location-based effect size measures (d types, not the variance-based xi or D measures) to ensure the scaling is comparable to Cohen's d and thus can be interpreted similarly.
%The lambda correction to the SD and SE is to take into account the fact that the trimming procedure causes the variables to become dependent (correlated) with each other.
%The Cousineau-Moray correction (in ep_ADF) is to ensure that the SE measures reflect the way repeated measure ANOVAs remove the correlated variance between the variables.

%Inputs
%  Y	    : Input data (rows=subjects, columns = cells).  The first NX rows will be assigned to the first group and so forth.
%  NX   	: Row vector with number of subjects in each group.  If empty set then will assume a single group.
%  C    	: Contrast row vector for between factors.  Set to 1 in the case where there is only one group.
%  U    	: Contrast column vector for within factors.  Numbers should sum to zero.  Set to empty set [] for analyses with no within-factors.
%  OPT1     : Activate trimming option (rounded down).  0 = no and 1 = yes.
%  PER  	: Percentage to trim the means.  .05 is the number recommended for ERP data by Dien.
%  OPT2 	: Activate Welch-James and ADF and bootstrapping statistic.  0 = no and 1 = yes.
%  NUMSIM	: Number of simulations used to generate bootstrapping statistic.  p-values will be unstable if too low.  50000 informally recommended.
%  SEED 	: Seed for random number generation.  0 specifies random SEED. 1000 arbitrarily suggested as SEED to ensure RESULTS are replicable.
%  MISSING  : Number to be treated as a missing value.  Observations with missing values are dropped from the analysis.
%  OPT3     : Provide effect sizes. 0 = no and 1 = yes.
%  ALPHA    : Alpha significance level.
%	.corrected : Alpha level corrected for multiple comparisons (not used)
%	.uncorrected : Alpha level not corrected for multiple comparisons
%  SCALE    : "SCALE is a scalar indicator to control the use of a scaling factor for the effect size estimator 
%             (i.e., .642 for 20% symmetric trimming) when robust estimators are adopted. 
%             It takes a value of 0 or 1; a zero indicates that no scaling factor will be used, 
%             while a 1 indicates that a scaling factor will be adopted. The default is SCALE=1.
%  LOC1&LOC2: "If the user specifies LOC1 = 0 and LOC2 = 0, the square root of the average of the variances over the cells 
%             involved in the contrast is used as the standardizer. If the
%             user specifies LOC1 = 99 and LOC2 = 99, no standardizer is selected."

%Outputs
%  MUHAT	: Vector of trimmed means used in calculations. 
%             For between group analyses, within group conditions vary fastest and between group conditions vary slowest.
%  STDIZER	: Winsorized sample variance-covariance matrix
%  RESULTS	: (1) = test statistic
%  RESULTS	: (2) = Numerator DF
%  RESULTS	: (3) = Denominator DF
%  RESULTS	: (4) = Significance
%  RESULTS	: (5) = Effect size using explanatory measure for between groups and projection distance for within groups and dR for 1df tests.
%  RESULTS	: (6) = Lower confidence limit (only for between contrasts with 1 df).
%  RESULTS	: (7) = Upper confidence limit (only for between contrasts with 1 df).
%  RESULTS	: (8) = Scaling factor for effect size estimator (SD units) if SCALE option is chosen (only for between contrasts with 1 df).
%  MSE      : Mean Squared Error.

% modified 10/7/07 JD
% Missing number feature added.
%
% modified 2/11/09 JD
% Warning messages only output once.
%
% modified 6/20/11 JD
% Handles new warning code for ill-conditioned matrices.
%
% modified 3/18/14 JD
% Changed uses of "temp" as a variable name to "tempVar" due to other Matlab programmers often using it as a function
% name, resulting in collisions.
%
% modified 9/14/15 JD
% Revised based on updated SAS/IML code by Lisa Lix to provide effect sizes.
%
% bugfix 10/27/15 JD
% Restored warning messages that were erroneously deleted for problematic analyses and sets p-value output to NaN (not a number).
% Handles case where the entire bootstrap sample is drawn from the same observation and treats as an F of inf.
%
% bugfix 11/30/15 JD
% No longer produces effect sizes for contrasts of more than 1 df.
% Handles case where almost the entire bootstrap sample is drawn from the same observation and was generating an NaN value. 
%
% bugfix 12/6/15 JD
% Fixed NaN effect sizes when PER is set to zero.
%
% modified 12/10/15 JD
% Modified so that singular and near-singular runs are dropped from the
% bootstrapping run and a warning is issued.
% If effect size d* is made positive, confidence interval signs also flipped.
% Added information on number of singular and nearly singular matrices during calculation of statistics to RESULTS.
%
% modified 12/27/15 JD
% Per new information from Lisa Lix, effect size output disabled for within-group
% contrasts.
%
% bugfix 6/16/17 JD
% Fixed apparently providing effect size for within contasts of more than two levels when should have been disabled.
%
% modified 8/8/18 JD
% Outputs MSE.
%
% modified 9/9/18 JD
% If SEED set to empty, then set to random number based on time rather than to zero.
%
% modified 7/5/19 JD
% Changed inv calls to pinv to better handle singular and near-singular data.
% Eliminated singularity counts in RESULTS.
%
% bugfix 9/2/22 JD
% Dropped some incomplete trial code for generating robust epsilon parameter inadvertently left active. 
%
% bugfix 9/2/24 JD
% Fixed outputing Winsorized squared standard error matrix instead of Winsorized variance-covariance matrix.
% Fixed confidence interval lower and upper bounds switched when effect size flipped to positive.
% Added d* effect sizes and confidence intervals for pairwise within-group effects as well adding effect sizes for between and within group effects with more than two levels.
%
% bugfix 5/23/25 JD
% Fixed confidence intervals for pairwise within contrasts when nboot=1 by setting to 1000 instead.
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

%****DEFINE MODULES TO PERFORM ALL CALCULATIONS****;
%****COMPUTE WELCH-JAMES STATISTIC****;

% if ~exist('epsilonType','var')
%     epsilonType='GG';
% end

MUHAT=[];
STDIZER=[];
RESULTS=[];
MSE=[];

errorflag.TooFewSubjects=0;

if (SEED == 0) || isempty(SEED)
    SEED=sum(100*clock);
end

sprev = rng(SEED,'twister');

numsim_b=NUMSIM;
numsim_es=NUMSIM;
numsim_bc=NUMSIM;

alphaThresh=ALPHA.uncorrected;

if isempty(NX)
    NX = size(Y,1);
end
if isempty(C)
    C = [1];
end

if isempty(U)
    U = [1];
end

if size(C,2) ~= length(NX)
    msg{1}=['Number of between group cells (' num2str(length(NX)) ') must equal number of terms in contrast C (' num2str(size(C,2)) ').'];
    [msg]=ep_errorMsg(msg);
    return
end

if size(U,1) ~= size(Y,2)
    msg{1}=['Number of within group cells (' num2str(size(Y,2)) ') must equal number of terms in contrast U (' num2str(size(U,1)) ').'];
    [msg]=ep_errorMsg(msg);
    return
end

if OPT1 ~=0 && OPT1 ~= 1
    msg{1}='OPT1 must equal zero or one.';
    [msg]=ep_errorMsg(msg);
    return
end

if OPT3 ~=0 && OPT3 ~= 1
    msg{1}='OPT3 must equal zero or one.';
    [msg]=ep_errorMsg(msg);
    return
end

% if (size(U,2) > 0) || (size(C,1) > 1)
%     OPT3=0; %can only calculate effect sizes for one degree between contrasts
% end

if ~isempty(MISSING) %drop observations with missing data points
    good = [];
    newNX= zeros(1,length(NX));
    count=0;
    for group = 1:length(NX)
        for obs = 1:NX(group)
            count=count+1;
            if isempty(find(Y(count,:)==MISSING))
                good=[good count];
                newNX(group)=newNX(group)+1;
            end
        end
    end
    if ~isempty(find(newNX==0))
        msg{1}='A group has zero members after dropping missing data points';
        [msg]=ep_errorMsg(msg);
        return
    end
    NX=newNX;
    Y=Y(good,:);
end

%**define module to check initial specifications****;
if isempty(U)
    U=eye(size(Y,2));
end
if size(U,2)>size(U,1)
    msg{1}='Possible Error: Number Of Columns Of U Exceeds Number Of Rows';
    [msg]=ep_errorMsg(msg);
    return
end

if OPT1==1
    if isempty(PER)
        PER=.20;
    end
    if PER > .49
        msg{1}='Error: Percentage Of Trimming Exceeds Upper Limit';
        [msg]=ep_errorMsg(msg);
        return
    end
end
    
if OPT2==1
    if isempty(numsim_b)
        numsim_b=999;
    end
end

if OPT3==1
    if isempty(numsim_es)
        numsim_es=999;
    end
    if isempty(LOC1)
        LOC1=1;
    end
    if isempty(LOC2)
        LOC2=1;
    end
end
        
if isempty(alphaThresh)
    alphaThresh=.05;
end
if isempty(SCALE)
    SCALE=1;
end
if isempty(numsim_bc)
    numsim_bc=699;
end

if sum(NX) ~= size(Y,1)
    msg{1}='Error: Total number in NX does not match Y data matrix.';
    [msg]=ep_errorMsg(msg);
    return
end

for I=1:size(NX,2)
    X1=ones(NX(I),1)*I;
    if I==1
        X=X1;
    else
        X=[X; X1];
    end
end
X=design(X); %full design matrix
NTOT=size(Y,1); %number of subjects
WOBS=size(Y,2); %total number of within cells
BOBS=size(X,2); %number of between groups
WOBS1=WOBS-1; %one less than total number of within cells
R=kron(C,U'); %contrast vector combining within and between contrasts
RESULTS=zeros(8,1);

%****compute Welch-James statistic****;
[MUHAT, BHAT, BHATW, YT, DF, errorflag] = mnmod(Y, OPT1, BOBS, WOBS, NTOT, NX, PER, X, errorflag);
if errorflag.TooFewSubjects
    return;
end
[SIGMA,STDIZER] = sigmod(YT,X,BHATW,DF,BOBS,WOBS,WOBS1,NX);
[FSTAT,DF1,DF2, MSE, errorflag] = testmod(SIGMA,MUHAT,R,DF, BOBS, WOBS, WOBS1, errorflag);

if OPT2==1
    %implements percentile bootstrap (Keselman, Algina, Lix, Wilcox, & Deering, 2008, p. 119, Footnote 7)
    FMAT=zeros(numsim_b,1);
    for simloop=1:numsim_b
        [YB1] = bootdat(Y, BHAT, BOBS, NX);
%         if ~any(std(YB1) > 10^(-12)) %Matlab is showing errors at about the 10^(-14) level
%             FSTATB=inf; %all of the bootstrap observations were identical
%         else
            [YB]=bootcen(YB1,BHAT,BOBS,NX);
            [FSTATB, errorflag] = bootstat(YB, OPT1, R, BOBS, WOBS, WOBS1, NTOT, NX, PER, X, SEED, errorflag);
%         end
        FMAT(simloop)=FSTATB;
    end
    FMAT=sort(FMAT);
    avec=(FSTAT<=FMAT);
    numNaN=sum(isnan(FMAT));
    pval=sum(avec)/(numsim_b-numNaN);
%     theEpsilon=DF2/(DF1*(NTOT-BOBS)); %compute equivalent of the epsilon correction factor for conventional ANOVAs.
else
%     switch epsilonType
%         case 'GG'
%             theEpsilon=epsilonGG(SIGMA,U);
%         case 'HF'
%             
%         case 'LB'
%             
%         otherwise
%             error('Oops - programmer error.')
%     end
%     DF2=DF1*(NTOT-BOBS);
%     
%     pval=1-probf(FSTAT,DF1*theEpsilon,DF2*theEpsilon);

    return; %not ready yet
end

%***calculate significance level for welch-james statistic****;
RESULTS(1)=FSTAT;
RESULTS(2)=DF1;
RESULTS(3)=DF2;
RESULTS(4)=pval;
RESULTS(5) = NaN;
RESULTS(6) = NaN;
RESULTS(7) = NaN;
RESULTS(8) = NaN;
% RESULTS(9) = theEpsilon;

if OPT3==1
    nullValue=0;
    if (BOBS==2) && (WOBS==1)
        %dR Algina,Keselman,Penfield (2005a)
        [AkpEffect,ci,~]=akpEffectCi(Y(1:NX(1),1),Y(NX(1)+1:end,1),alphaThresh,PER,numsim_es,SEED,nullValue);
            RESULTS(5) = AkpEffect;
            RESULTS(6) = ci(1);
            RESULTS(7) = ci(2);
    elseif (BOBS==1) && (WOBS==2)
        %dR Algina,Keselman,Penfield (2005b)
        [EffectSize,ci,~] = DakpEffectCi(Y(:,1),Y(:,2),nullValue,alphaThresh,PER,numsim_es,SEED);
            RESULTS(5) = EffectSize;
            RESULTS(6) = ci(1);
            RESULTS(7) = ci(2);
    elseif (BOBS>1) && (WOBS==1)
        %explanatory measure of effect size (Wilcox & Tian, 2011)
        YR=IML2R(Y,NX);
        [EffectSize,ci]=t1wayEXESci(YR,MUHAT,NX,alphaThresh,PER,numsim_es,SEED,numsim_es,true,pval);
        RESULTS(5)=EffectSize;
        RESULTS(6) = ci(1);
        RESULTS(7) = ci(2);
    elseif (BOBS==1) && (WOBS>1)
        %projection distance measure of effect size (Wilcox, 2022, p.496)
        RESULTS(5)=rmESdifPro(Y,PER);
    elseif (BOBS==1) && (WOBS==1)
        %effect size from updated SAS/IML code by Lisa Lix and colleagues, for use with generating error bars
        esmat=zeros(numsim_es,1);
        for simloop=1:numsim_es
            [YB1] = bootdat(Y, BHAT, BOBS, NX);
            [MULTP,EFFSZB,errorflag]=bootes(YB1,X,LOC1,LOC2,SCALE,OPT1,PER,R,BOBS,WOBS,WOBS1,NTOT,NX,errorflag);
            esmat(simloop)=EFFSZB;
        end
        esmat=sort(esmat);
        numNaN=sum(isnan(esmat));
        numsim_esGood=numsim_es-numNaN;
        if numsim_esGood>0
            ind1=intIML(numsim_esGood*(alphaThresh/2))+1;
            ind2=numsim_esGood-intIML((numsim_esGood*(alphaThresh/2)));
            lcl=esmat(ind1);
            ucl=esmat(ind2);
            [MULTP,EFFSZB,errorflag]=bootes(Y,X,LOC1,LOC2,SCALE,OPT1,PER,R,BOBS,WOBS,WOBS1,NTOT,NX,errorflag);
            RESULTS(5) = abs(EFFSZB); %convention for Cohen's d is to provide absolute value JD
            RESULTS(6) = lcl;
            RESULTS(7) = ucl;
            if sign(EFFSZB) == -1
                RESULTS(6)=-ucl;
                RESULTS(7)=-lcl;
            end
            RESULTS(8) = MULTP;
        end
    end
end

function a = design(X)
a = zeros(size(X,1),length(unique(X)));
for i=1:length(X)
    a(i,X(i)) = 1;
end

function a = intIML(X)
%perform equivalent of SAS/IML int function
a=X-rem(X,1);

%****define module to compute least squares or trimmed means****;
function [MUHAT, BHAT, BHATW, YT, DF, errorflag] = mnmod(Y, OPT1, BOBS, WOBS, NTOT, NX, PER, X, errorflag)
    
if OPT1==0
    BHAT=pinv(X'*X)*X'*Y;
    BHATW=BHAT;
    YT=Y;
    DF=NX-1;
end

if OPT1==1
    BHAT=zeros(BOBS,WOBS);
    BHATW=BHAT;
    YT=zeros(NTOT,WOBS);
    DF=zeros(1,BOBS);
    F=1;
    M=0;
    for J=1:size(NX,2) %loop through between groups
        SAMP=NX(J);%size of between group
        L=M+SAMP;
        G=intIML(PER.*SAMP);
        DF(J)=SAMP-2*G-1;
        for K=1:size(Y,2) %loop through within groups
            tempVar=Y(F:L,K);
            NV=tempVar;
            [TEMP2 index] = sort(NV);
            TRIMY=TEMP2(G+1:SAMP-G,:); %trimmed cell
            TRIMMN=sum(sum(TRIMY))/(DF(J)+1); %trimmed mean
            BHAT(J,K)=TRIMMN; %matrix of trimmed means
            MINT=min(min(TRIMY));
            MAXT=max(max(TRIMY));
            for P=1:size(NV,1)
                if NV(P)<=MINT
                    NV(P)=MINT;
                end
                if NV(P)>=MAXT
                    NV(P)=MAXT;
                end
            end
            YT(F:L,K)=NV; %winsorized sample
            WINMN=sum(sum(NV))/SAMP;
            BHATW(J,K)=WINMN;
        end
        M=L;
        F=F+NX(J);
    end
end

if any(DF == 0)
    errorflag.TooFewSubjects=1;
    MUHAT=0;
    BHAT=0;
    BHATW=0;
    YT=0;
    msg{1}='Error, too few subjects.  Degrees of freedom is zero.';
    [msg]=ep_errorMsg(msg);
    return
end

MUHAT=reshape(BHAT',[],BOBS*WOBS)';

%***DEFINE MODULE TO COMPUTE SIGMA MATRIX****;
function [SIGMA,STDIZER] = sigmod(YT,X,BHATW,DF,BOBS,WOBS,WOBS1,NX)
SIGMA=zeros(WOBS.*BOBS,WOBS.*BOBS);
STDIZER=SIGMA;
for I=1:BOBS
    %Squared Winsorized standard error of the sample mean (3.8).  Wilcox (2022) p. 63, except with Lambda=g/n swapped in..
    SIGB=(diag(X(:,I))*YT-X(:,I)*BHATW(I,:))'*(diag(X(:,I))*YT-X(:,I)*BHATW(I,:))/((DF(I)+1)*DF(I));
    F=I*WOBS-WOBS1;
    L=I*WOBS;
    SIGMA(F:L,F:L)=SIGB;
    STDIZER(F:L,F:L)=SIGB*((DF(I)+1)*DF(I))/(NX(I)-1); %sample Winsorized variance (3.9).  Wilcox (2022) p. 63.
end

%****DEFINE MODULE TO COMPUTE TEST STATISTIC****;
function [FSTAT, DF1, DF2, MSE, errorflag] = testmod(SIGMA, MUHAT, R, DF, BOBS, WOBS, WOBS1, errorflag)
T=(R*MUHAT)'*pinv(R*SIGMA*R')*(R*MUHAT); %T stat squared and without df correction Twj statistic. Johansen (1980)
A=0;
IMAT=eye(WOBS);
for I=1:BOBS
    QMAT=zeros(BOBS*WOBS,BOBS*WOBS);
    F=I*WOBS-WOBS1;
    L=I*WOBS;
    QMAT(F:L,F:L)=IMAT;
    PROD=(SIGMA*R')*pinv(R*SIGMA*R')*R*QMAT;
    A=A+(trace(PROD*PROD)+trace(PROD)^2)/DF(I);  %ADF statistic Lix & Keselman (1995)
end
A=A/2; %for adjusted degrees of freedom
DF1=size(R,1);
DF2=DF1*(DF1+2)/(3*A);
CVAL=DF1+2*A-6*A/(DF1+2);
FSTAT=T/CVAL;
SST=0;
for iContrast=1:DF1
    SST=SST+((R(iContrast,:)*MUHAT)'*(R(iContrast,:)*MUHAT))/sum(R(iContrast,:).^2./kron(DF+1,ones(1,WOBS)));
end
MST=SST/DF1;
MSE=MST/FSTAT;

%****DEFINE MODULES TO PERFORM BOOTSTRAP****;
%***DEFINE MODULE TO GENERATE BOOTSTRAP DATA AND CENTRE DATA****;
%resamples with replacement the observations for each between-group, resulting in a resampled dataset with the same overall dimensions.
function [YB] = bootdat(Y, BHAT, BOBS, NX)
F=1;
M=0;
tempVar=[];
YB=zeros(size(Y));
for J=1:BOBS
    L=M+NX(J);
    tempVar=Y(F:L,:);
    BVAL=tempVar;
    for P=1:size(tempVar,1)
        RVAL=0;
        while RVAL ==0
            RVAL=rand(1);
        end
        BVAL(P,:)=tempVar(ceil(size(tempVar,1)*RVAL),:);
    end
    YB(F:L,:)=BVAL;
    M=L;
    F=F+NX(J);
end

%****CENTRE THE BOOTSTRAP DATA****;
%center the bootstrap sample based on the means from the full sample
function [YB]=bootcen(YB1,BHAT,BOBS,NX)
M=0;
F=1;
YB=YB1;
for I=1:BOBS
    L=M+NX(I);
    MVAL=BHAT(I,:);
    for Q=F:1:L
        YB(Q,:)=YB1(Q,:)-MVAL;
    end
    M=L;
    F=F+NX(I);
end

%****DEFINE MODULE TO COMPUTE BOOTSTRAP STATISTIC****;
function [FSTATB, errorflag] = bootstat(YB,OPT1,R, BOBS, WOBS, WOBS1, NTOT, NX, PER, X, SEED, errorflag)
[MUHATB, BHATB, BHATBW, YTB, DFB, errorflag] = mnmod(YB, OPT1, BOBS, WOBS, NTOT, NX, PER, X, errorflag);
[SIGMAB,STDIZERB] = sigmod(YTB,X,BHATBW,DFB,BOBS,WOBS,WOBS1,NX);
[FSTATB,DF1B,DF2B, MSE, errorflag] = testmod(SIGMAB,MUHATB,R,DFB, BOBS, WOBS, WOBS1, errorflag);

%****define module to compute bootstrap effect size****;
function [MULTP,EFFSZB,errorflag]=bootes(YB,X,LOC1,LOC2,SCALE,OPT1,PER,R,BOBS,WOBS,WOBS1,NTOT,NX,errorflag)
[MUHATB, BHATB, BHATBW, YTB, DFB, errorflag] = mnmod(YB, OPT1, BOBS, WOBS, NTOT, NX, PER, X, errorflag);
[SIGMAB,STDIZERB] = sigmod(YTB,X,BHATBW,DFB,BOBS,WOBS,WOBS1,NX);
[MULTP,EFFSZB]=wjeffsz(BOBS,WOBS,WOBS1,NTOT,NX,LOC1,LOC2,SCALE,OPT1,PER,R,MUHATB,STDIZERB);

%****compute measure of effect size and bootstrap confidence interval****;
function [MULTP,EFFSZ]=wjeffsz(BOBS,WOBS,WOBS1,NTOT,NX,LOC1,LOC2,SCALE,OPT1,PER,R,MUHAT,STDIZER)
if OPT1==0 
    MULTP=1;
end
if OPT1==1
    if SCALE==0
        MULTP=1;
    end
    if SCALE==1
        if PER ~= 0
            cut=sqrt(2)*erfinv(2*PER-1); %probit per wikipedia https://en.wikipedia.org/wiki/Probit 12/6/2015
            %**fun defines the probability density function of the standard normal distribution**;
            fun = @(Z) (Z.^2).*(1/(sqrt(2*3.141592653589793)).*(exp(-(1/2).*Z.^2)));
            int=integral(fun,cut,-cut);
            b=sum(int);
            winvar=b+PER*(2*(cut^2));
            MULTP=sqrt(winvar);
        else
            MULTP=1;
        end
    end
end
num=R*MUHAT; %will always be a scalar since only single df can be calculated at present.
if LOC1==99
    stdz=1;
end
if LOC1==0
    r2=R.^2;
    rvec=reshape(r2',[],BOBS*WOBS);
    stdz1=diag(sqrt(diag(STDIZER)));
    stdz=(rvec*stdz1*rvec')/sum(rvec); %average of square root of variances
end
if LOC1>0
    if LOC1 <99
        loc=LOC1*WOBS-(WOBS-LOC2);
        stdz1 = STDIZER(loc,loc);
        if stdz1> 0
            stdz=sqrt(stdz1);
        end
        if stdz1==0
            stdz=.00001;
        end
    end
end

if (BOBS>1) && (WOBS > 1)
    %Effect sizes not available for between-within contrasts, pending further research.
    EFFSZ=NaN; 
else
    EFFSZ=MULTP*(num/stdz);
end

function x = probf(f,d1,d2)
%Based on code at https://matrixlab-examples.com/f-distribution, implements SAS probf function.
x = 1;
% Computes using inverse for small F-values
if f < 1
    s = d2;
    t = d1;
    z = 1/f;
else
    s = d1;
    t = d2;
    z = f;
end
j = 2/(9*s);
k = 2/(9*t); 

% Uses approximation formulas
y = abs((1 - k)*z^(1/3) - 1 + j)/sqrt(k*z^(2/3) + j);
if t < 4
    y = y*(1 + 0.08*y^4/t^3);
end 

a1 = 0.196854;
a2 = 0.115194;
a3 = 0.000344;
a4 = 0.019527;
x = 0.5/(1 + y*(a1 + y*(a2 + y*(a3 + y*a4))))^4;
x = floor(x*10000 + 0.5)/10000; 

% Adjusts if inverse was computed
if f < 1
    x = 1 - x;
end 

x = 1 - x; %SAS probf provides percentile not tail end value.

% function theEpsilon=epsilonGG(SIGMA,U)
% %calculate Greenhouse-Geisser epsilon-hat for conventional ANOVAs
% U=orth(U);
% theEpsilon=(trace(U'*SIGMA*U))^2/(size(U,2)*trace((U'*SIGMA*U)^2));
% 

function YR=IML2R(Y,NX)
%converts the data format expected by the SAS/IML code to a cell array equivalent of the R code's list format.
YR=cell(length(NX),1);
Ycounter=0;
for iGroup=1:length(NX)
    YR{iGroup}=Y(Ycounter+1:Ycounter+NX(iGroup),:);
    Ycounter=Ycounter+NX(iGroup);
end

function WinvarN = winvarN(x, tr)
% rescale the Winsorized variance so that it equals one for the standard
% normal distribution
% Translated from R library WRS (v43) winvarN

if ~exist('tr','var')
    tr=.2;
end
x=x(~isnan(x));

cterm=[];
if tr==0
    cterm=1;
end
if tr==0.1
    cterm=0.6786546;
end
if tr==0.2
    cterm=0.4120867;
end
if isempty(cterm)
    cut=sqrt(2)*erfinv(2*tr-1); %probit per wikipedia https://en.wikipedia.org/wiki/Probit 12/6/2015
    %**fun defines the probability density function of the standard normal distribution**;
    fun = @(Z) (Z.^2).*(1/(sqrt(2*3.141592653589793)).*(exp(-(1/2).*Z.^2)));
    int=integral(fun,cut,-cut);
    b=sum(int);
    cterm=b+tr*(2*(cut^2)); %sqrt of this term is MULTP
end
WinvarN=winvar(x,tr)/cterm;

function Winvar=winvar(x,tr,narm)
%  Compute the gamma Winsorized variance for the data in the vector x.
%  tr is the amount of Winsorization which defaults to .2.
% Translated from R library WRS (v43) winvar

if ~exist('tr','var')
    tr=.2;
end
if ~exist('narm','var')
    narm=false;
end

remx=x;
if narm
    x=x(~isnan(x));
end
y=sort(x);
n=length(x);
ibot=floor(tr*n)+1;
itop=length(x)-ibot+1;
xbot=y(ibot);
xtop=y(itop);
y(y<=xbot)=xbot;
y(y>=xtop)=xtop;
Winvar=var(y);
if ~narm
    if any(isnan(remx))
        Winvar=NaN;
    end
end

function Wincor=WINCOR(x,tr)
% For convenience, compute Winsorized correlation matrix only.
% Translated from R library WRS (v43) WINCOR

if ~exist('tr','var')
    tr=.2;
end
[Wincor, ~, ~, ~]=winall(x,tr);

function [n, corOut, covOut, pvalue]=wincor(x,y,tr)
%   Compute the Winsorized correlation between x and y.
%   tr is the amount of Winsorization
%   This function also returns the Winsorized covariance
%    Pairwise deletion of missing values is performed.
%   x is a vector, or it can be a matrix with two columns when y=NULL
% Translated from R library WRS (v43) wincor

if ~exist('tr','var')
    tr=.2;
end

if ~exist('y','var')
    y=[];
end

if ~isempty(y)
    m=[x,y];
else
    m=x;
end

m=m(~any(isnan(x),2),:);

nval=size(m,1);
ncol=size(m,2);
if ncol==2
    [a_cor, a_cov, a_pvalue]=wincorSub(m(:,1),m(:,2),tr);
    wcor=a_cor;
    wcov=a_cov;
    sig=a_pvalue;
end

if ncol>2
    wcor=ones(ncol,ncol);
    wcov=zeros(ncol,ncol);
    siglevel=nan(ncol,ncol);
    for i=1:ncol
        ip=i;
        for j=ip:ncol
            [val_cor, val_cov, val_pvalue]=wincorSub(m(:,i),m(:,j),tr);
            wcor(i,j)=val_cor;
            wcor(j,i)=wcor(i,j);
            if i==j
                wcor(i,j)=1;
            end
            wcov(i,j)=val_cov;
            wcov(j,i)=wcov(i,j);
            if i~=j
                siglevel(i,j)=val_pvalue;
                siglevel(j,i)=siglevel(i,j);
            end
        end
    end
    sig=siglevel;
end

n=nval;
corOut=wcor;
covOut=wcov;
pvalue=sig;

function [corOut, covOut, center, pvalues] = winall(m,tr)
%    Compute the Winsorized correlation and covariance matrix for the
%    data in the n by p matrix m.
%
%    This function also returns the two-sided significance level
%    Translated from R library WRS (v43) winall

if ~exist('tr','var')
    tr=.2;
end

ncol=size(m,2);
nrow=size(m,1);
wcor=ones(ncol,ncol);
wcov=zeros(ncol,ncol);
siglevel=nan(ncol,ncol);
for i=1:ncol
    ip=i;
    for j=ip:ncol
        [val_cor, val_cov, ~, val_pvalues]=wincor(m(:,i),m(:,j),tr);
        wcor(i,j)=val_cor;
        wcor(j,i)=wcor(i,j);
        if i==j 
            wcor(i,j)=1;winval
        end
        wcov(i,j)=val_cov;
        wcov(j,i)=wcov(i,j);
        if i~=j 
            siglevel(i,j)=val_pvalues;
            siglevel(j,i)=siglevel(i,j);
        end
    end
end
G=intIML(tr*nrow);
cent=mean(m(G+1:nrow-G,:),2,"omitmissing");
corOut=wcor;
covOut=wcov;
center=cent;
pvalues=siglevel;


function [corOut, covOut, pvalue]=wincorSub(x,y,tr)
%  Translated from R library WRS (v43) wincor.sub

sig=NaN;
g=floor(tr*size(x,1));
xvec=winval(x,tr);
yvec=winval(y,tr);
wcor=corrcoef(xvec,yvec);
wcor=wcor(1,1);
wcov=cov(xvec,yvec);
if sum(x==y) ~= size(x,1)
    test=wcor*sqrt((length(x)-2)/(1.-wcor^2));
    sig=2*(1-tcdf(abs(test),length(x)-2*g-2));
end
corOut=wcor;
covOut=wcov;
pvalue=sig;

function Winval=winval(x,tr)
%  Winsorize the data in the vector x.
%  tr is the amount of Winsorization which defaults to .2.
%
%  This function is used by several other functions that come with this book.
%  Translated from R library WRS (v43) winval

if ~exist('tr','var')
    tr=.2;
end

y=sort(x);
n=length(x);
ibot=floor(tr*n)+1;
itop=length(x)-ibot+1;
xbot=y(ibot);
xtop=y(itop);
x(x<=xbot)=xbot;
x(x>=xtop)=xtop;
Winval=x;

function [x,grpn]=selby(m,grpc,coln)
%
%  A commmon situation is to have data stored in an n by p matrix where
%  one or more of the columns are  group identification numbers.
%  This function groups  all values in column coln according to the
%  group numbers in column grpc and stores the  results in list mode.
%
%  More than one column of data can sorted
%
% grpc indicates the column of the matrix containing group id number
%  Translated from R library WRS (v43) selby

if isempty(m)
    error("Data must be stored in a matrix");
end
if isnan(grpc(1))
    error("The argument grpc is not specified");
end
if isnan(coln(1))
    error("The argument coln is not specified");
end
if length(grpc)~=1
    error("The argument grpc must have length 1");
end

grpn=sort(unique(m(:,grpc)));
x=cell(grpn,1);
it=0;
for ig=1:length(grpn)
    for ic=1:length(coln)
        it=it+1;
        flag=(m(:,grpc)==grpn(ig));
        x{it}=m(flag,coln(ic));
    end
end

function dis=rmESdifPro(x,tr)
%
%  Global measure of effect size,
%   based on difference scores,
%  relative to the  null distribution
%
%  Translated from R library WRS (v43) rmES.dif.pro

if ~exist('tr','var')
    tr=.2;
end

n=size(x,1);
n1=n+1;
J=size(x,2);
ALL=(J^2-J)/2;
M=nan(n,ALL);
ic=0;
for j=1:J
    for k=1:J
        if j<k
            ic=ic+1;
            M(:,ic)=x(:,j)-x(:,k);
        end
    end
end
ES=tmean(M,tr);
Z=zeros(1,ALL);
dis=pdis(M,ES,[],[],[],Z,[]);

function Pdis=pdis(m,pts,MM,cop,dop,center,naRm)
%
% Compute projection distances for points in pts relative to points in m
%  That is, the projection distance from the center of m
%
%
%  MM=F  Projected distance scaled
%  using interquatile range.
%  MM=T  Scale projected distances using MAD.
%
%  There are five options for computing the center of the
%  cloud of points when computing projections:
%  cop=1 uses Donoho-Gasko median
%  cop=2 uses MCD center
%  cop=3 uses median of the marginal distributions.
%  cop=4 uses MVE center
%  cop=5 uses skipped mean
%
%  Translated from R library WRS (v43) pdis

if ~exist('pts','var') || isempty(pts)
    pts=m;
end

if ~exist('MM','var') || isempty(MM)
    MM=false;
end

if ~exist('cop','var') || isempty(cop)
    cop=3;
end

if ~exist('dop','var') || isempty(dop)
    dop=1;
end

if ~exist('center','var') || isempty(center)
    center=nan;
end

if ~exist('naRm','var') || isempty(naRm)
    naRm=true;
end

% if naRm
%     m=m(~isnan(m)); % Remove missing values
%     pts=pts(~isnan(pts));
% end

if iscell(m)
    m=cell2mat(m);
end
nm=size(m,1);

if iscell(pts)
    pts=cell2mat(pts);
end

if size(m,2)>1
    if size(pts,2)==1
        pts=pts';
    end
end

npts=size(pts,1);
mp=[m; pts];
np1=size(m,1)+1;

if size(m,2)==1
    if isnan(center(1))
        center=median(m);
    end
    dis=abs(pts-center);
    disall=abs(m-center);
    temp=idealf(disall);
    if ~MM
        Pdis=dis/(temp.qu-temp.ql);
    end
    if(MM)
        Pdis=dis/mad(disall);
    end
else
    if(isnan(center(1)))
        switch cop
            case 1
                % center=dmean(m,tr=.5,dop=dop);
            case 2
                % center=cov.mcd(m)$center;
            case 3
                center=median(m,2);
            case 4
                % center=cov.mve(m)$center;
            case 5
                % center=smean(m);
            otherwise
        end
    end
    dmat=nan(size(mp,1),size(mp,1));
    for i=1:size(mp,1)
        B=mp(i,:)-center;
        dis=nan;
        BB=B.^2;
        bot=sum(BB);
        if bot~=0
            for j=1:size(mp,1)
                A=mp(j,:)-center;
                temp=sum(A.*B).*B./bot;
                dis(1,j)=sqrt(sum(temp.^2));
            end
            disM=dis(1:nm);
            if ~MM
                temp=idealf(disM);
                dmat(:,i)=dis/(temp.qu-temp.ql);
            end
            if(MM)
                dmat(:,i)=dis/mad(disM,1);
            end
        end
    end
    Pdis=max(dmat,[],2,"omitnan");
    Pdis=Pdis(np1:size(mp,1));
end

function [temp]=idealf(x,naRm)
%
% Compute the ideal fourths for data in x
%
%  Translated from R library WRS (v43) idealf

if ~exist('naRm','var')
    naRm=false;
end

if naRm
    x=x(~isnan(x));
end
j=floor(length(x)/4 + 5/12);
y=sort(x);
g=(length(x)/4)-j+(5/12);
temp.ql=(1-g)*y(j)+g*y(j+1);
k=length(x)-j+1;
temp.qu=(1-g)*y(k)+g*y(k-1);

function dval=DakpEffect(x,y,nullValue,tr)
%
% Computes the robust effect size for one-sample case using
% a simple modification of
% Algina, Keselman, Penfield Pcyh Methods, 2005, 317-328
%
%  When comparing two dependent groups, data for the second group can be stored in
%  the second argument y. The function then computes the difference scores x-y
%  Translated from R library WRS (v43) D.akp.effect

if ~exist('y','var')
    y=[];
end
if ~exist('nullValue','var')
    nullValue=0;
end
if ~exist('tr','var')
    tr=.2;
end

if ~isempty(y)
    x=x-y;
end
x=x(~isnan(x));
s1sq=winvar(x,tr);
cterm=1;
if tr>0
    dnormvar= @(x) x.^2.*normpdf(x); %  Translated from R library WRS (v43) dnormvar
    cterm=integral(dnormvar,norminv(tr),norminv(1-tr))+2*(norminv(tr)^2)*tr;
end
cterm=sqrt(cterm);
dval=cterm*(tmean(x,tr)-nullValue)/sqrt(s1sq);
% change by JD, after consultation with Wilcox (personal communication, 2024)
if (cterm*(tmean(x,tr)-nullValue))==0
    dval=0;
elseif sqrt(s1sq)==0
    dval=Inf; %to indicate that the EF measure has broken down and is large but otherwise uninterpretable.
end
%*************


function [EffectSize,ci,pValue] = DakpEffectCi(x,y,nullValue,alpha,tr,nboot,SEED)
%
% Computes the robust effect size for one-sample case using
% a simple modification of
% Algina, Keselman, Penfield Pcyh Methods, 2005, 317-328
%
%  When comparing two dependent groups, data for the second group can be stored in
%  the second argument y. The function then computes the difference scores x-y
%  Translated from R library WRS (v43) D.akp.effect.ci
if ~exist('y','var')
    y=[];
end
if ~exist('nullValue','var')
    nullValue=0;
end
if ~exist('alpha','var')
    alpha=.05;
end
if ~exist('tr','var')
    tr=.2;
end
if ~exist('nboot','var') || (nboot==1) %JD
    nboot=1000;
end
if ~exist('SEED','var')
    SEED=true;
end

% if(SEED)set.seed(2)

if ~isempty(y)
    x=x-y;
end
x=x(~isnan(x));
a=DakpEffect(x,[],nullValue,tr);
v=nan;
for i = 1:nboot
    X=x(randi(length(x),length(x),1));
    v(i)=DakpEffect(X,[],nullValue,tr);
end
v=sort(v);
ilow=round((alpha/2) * nboot);
ihi=nboot-ilow;
ilow=ilow+1;
ci=v(ilow);
ci(2)=v(ihi);
pv=mean(v<0)+.5*mean(v==0);
pv=2*min(pv,1-pv);
EffectSize=a;
pValue=pv;

function dval=akpEffect(x,y,EQVAR,tr)
%
% Computes the robust effect size suggested by
% Algina, Keselman, Penfield Psych Methods, 2005, 317-328
%  Translated from R library WRS (v43) akp.effect
% Between-group dR measure.

if ~exist('EQVAR','var') || isempty(EQVAR)
    EQVAR=true;
end
if ~exist('tr','var')
    tr=.2;
end
x=x(~isnan(x));
y=y(~isnan(y));

n1=length(x);
n2=length(y);
s1sq=winvar(x,tr);
s2sq=winvar(y,tr);
spsq=(n1-1)*s1sq+(n2-1)*s2sq;
sp=sqrt(spsq/(n1+n2-2));
cterm=1;
if tr>0
    dnormvar= @(x) x.^2.*normpdf(x); %  Translated from R library WRS (v43) dnormvar
    cterm=integral(dnormvar,norminv(tr),norminv(1-tr))+2*(norminv(tr)^2)*tr;
end
cterm=sqrt(cterm);
if EQVAR
    dval=cterm*(tmean(x,tr)-tmean(y,tr))/sp;
    % change by JD, after consultation with Wilcox (personal communication, 2024)
    if (cterm*(tmean(x,tr)-tmean(y,tr)))==0
        dval=0;
    elseif sp==0
        dval=Inf; %to indicate that the EF measure has broken down and is large but otherwise uninterpretable.
    end
    %*************
else
    dval(1)=cterm*(tmean(x,tr)-tmean(y,tr))/sqrt(s1sq);
    dval(2)=cterm*(tmean(x,tr)-tmean(y,tr))/sqrt(s2sq);
    % change by JD, after consultation with Wilcox (personal communication, 2024)
    if (cterm*(tmean(x,tr)-tmean(y,tr)))==0
        dval(1)=0;
        dval(2)=0;
    elseif sqrt(s1sq)==0
        dval(1)=Inf; %to indicate that the EF measure has broken down and is large but otherwise uninterpretable.
        dval(2)=Inf;
    end
    %*************
end

function [AkpEffect,ci,pValue]=akpEffectCi(x,y,alpha,tr,nboot,SEED,nullValue)
%
% Computes the robust effect size for two-sample case using
% Algina, Keselman, Penfield Pcyh Methods, 2005, 317-328
%
%  Translated from R library WRS (v43) akp.effect.ci

if ~exist('y','var')
    y=[];
end
if ~exist('alpha','var')
    alpha=.05;
end
if ~exist('tr','var')
    tr=.2;
end
if ~exist('nboot','var')
    nboot=1000;
end
if ~exist('SEED','var')
    SEED=true;
end
if ~exist('nullValue','var')
    nullValue=0;
end

%if(SEED)set.seed(2)

x=x(~isnan(x));
n1=length(x);
n2=length(y);
beF=nan;
for i=1:nboot
    X=x(randi(n1,n1,1));
    Y=y(randi(n2,n2,1));
    beF(i)=akpEffect(X,Y,true,tr);
end
L=floor(alpha*nboot/2);
U=nboot-L;
beF=sort(beF);
ci=beF(L+1);
ci(2)=beF(U);
est=akpEffect(x,y,true,tr);
pv=mean(beF<nullValue)+mean(beF==nullValue);
pv=2*min(pv,1-pv);
AkpEffect=est;
pValue=pv;


function [EffectSize,ci]=t1wayEXESci(x,xbar,nval,alpha,tr,nboot,SEED,ITER,adj,pValue)
%
% Confidence interval for explanatory measure of effect size
%
%  ITER:  yuenv2, for unequal sample sizes.  iterates to get estimate
%
%  Translated from R library WRS (v43) t1way.EXES.ci

if ~exist('alpha','var')
    alpha=.05;
end
if ~exist('tr','var')
    tr=0;
end
if ~exist('nboot','var')
    nboot=500;
end
if ~exist('SEED','var')
    SEED=true;
end
if ~exist('ITER','var')
    ITER=5;
end
if ~exist('adj','var')
    adj=true;
end

J=length(nval);
% if(SEED)set.seed(2)
[~,chk_EffectSize]=t1wayv2EP(x,tr,nval,nboot,'median');
v=cell(0);
val=nan;
for iGroup=1:J
    x{iGroup}=x{iGroup}(~isnan(x{iGroup}));
end
for i=1:nboot
    for j=1:J
        n1=length(x{j});
        v{j}=x{j}(randi(n1,n1,1));
    end
    [~,val(i)]=t1wayv2EP(v,tr,nval,nboot,'median');
end
ilow=round((alpha/2) * nboot);                                                     
ihi=nboot - ilow;                                                                  
ilow=ilow+1;                                                                       
val=sort(val);                                                                      
ci=val(ilow);                                                                       
ci(2)=val(ihi);                                                                     
if pValue>alpha
    ci(1)=0;
end
if adj
    fix=[1, 1.268757, 1.467181, 1.628221, 1.763191, 1.856621, 1.993326];
end
% if J>8
%     disp('No adjustment available when J>8');
% end
J1=J-1;
if j<=8
    chk_EffectSize=fix(J1)*chk_EffectSize;
    ci=fix(J1)*ci;
else
    chk_EffectSize=[];
end
EffectSize=chk_EffectSize;


function [ExplanatoryPower,EffectSize]=t1wayv2EP(x,tr,nval,nboot,locFun)
% Translated from a portion of Wilcox's R library WRS (v43) t1wayv2
% Explanatory measure of effect size for between-group effects.

if ~exist('tr','var')
    tr=.2;
end
if ~exist('nboot','var')
    nboot=100;
end
if ~exist('locFun','var')
    locFun='median';
end

pts=[];
nval=0;
J=length(x); %for the purpose of this Matlab translation, it will be assumed that all cells of the data are to be analyzed.
for j=1:J
    val=x{j};
    val=val(~isnan(val));
    nval(j)=length(val);
    pts=[pts; val]; %column vectors are more efficient in Matlab than row vectors.
    xbar(j)=tmean(val,tr);
end

if ~any(diff(nval)) %if the groups are all the same size (between group only case)
    top=var(xbar);
    bot=winvarN(pts,tr);
    ePow=top/bot;
    % change by JD, after consultation with Wilcox (personal communication, 2024)
    if top==0
        ePow=0;
    elseif bot==0
        ePow=Inf; %to indicate that the EF measure has broken down and is large but otherwise uninterpretable.
    end
    %*************
else
    vals=0;
    N=min(nval);
    xdat=zeros(J,1);
    for iRep = 1:nboot
        for iGroup = 1:J
            xdat(iGroup)=randperm(x{iGroup},N); %downsample all the groups to the same size, taking the median of the resulting group means.
        end
        vals(iRep)=t1wayEffectEP(x,xbar,tr,nval);
    end
    switch locFun
        case 'median'
            ePow=median(vals,"omitmissing");
        case 'mean'
            ePow=mean(vals,"omitmissing");
        otherwise
            disp(['The option ' locFun ' for the locFun parameter of the t1wayv2 call not recognized.'])
            return
    end
end
ExplanatoryPower=ePow;
EffectSize=sqrt(ePow);

function VarExplained=t1wayEffectEP(Y,xbar,tr,nval)
% Translated from a portion of Wilcox's R library WRS (v43) t1way.effect

if ~exist('tr','var')
    tr=.2;
end

top=var(xbar);
bot=winvarN(Y,tr);
if  bot>0
    ePow=top/bot;
end
if bot==0
    ePow=1;
end
if ePow>=1
    v1=[];
    v2=[];
    xCounter=0;
    for j=1:J
        v1=[v1; repmat(xbar(j),nval(j),1)];
        v2=[v2; Y(xCounter+1:xCounter+nval(j))];
        xCounter=xCounter+nval(j);
    end
    [~, cor, ~, ~]=wincor(v1,v2,tr);
    ePow=cor^2;
end

% change by JD, after consultation with Wilcox (personal communication, 2024)
if top==0
    ePow=0;
elseif bot==0
    ePow=Inf; %to indicate that the EF measure has broken down and is large but otherwise uninterpretable.
end
%*************

VarExplained=ePow;

function Y=tmean(X,tr,dim)
% R's mean (WRS has tmean built on mean) uses floor(n*tr), consistant with the SAS/IML code, whereas Matlab's trimmean uses round(n*tr)
% also, tmean uses decimal (e.g., .05) of each side whereas trimmean uses total trimming in non-decimal (e.g., 10).
if ~exist('dim','var')
    dim=1;
end
n=size(X,dim);
nTrim=floor(n*tr);
switch dim
    case 1
        Xt=X(nTrim+1:end-nTrim,:);
    case 2
        Xt=X(:,nTrim+1:end-nTrim);
    otherwise
        error('programmer error')
end
Y=mean(Xt,dim);