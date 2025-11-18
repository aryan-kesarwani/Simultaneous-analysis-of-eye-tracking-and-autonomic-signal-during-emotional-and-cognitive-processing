function [outCovMat]=ep_covMat(inCovMat, combineList, theOperation, theWeights, theMeans, sampleSizes)
% [outCovMat]=ep_covMat(inCovMat, combineList, theOperation, theWeights)
% Performs operations on a covariance matrix.
%
%Input:
%    inCovMat       : input 3D covariance matrix (chan,chan,sample).
%    combineList    : The index of levels to be combined.  If empty, then all of them.
%    theOperation	: The operation: 'combineChan' to add a new channel.
%                                    'poolObs' to add covariance matrices from different disjoint samples partialling out mean differences, as in covAVE.
%                                    'combineObs' to add covariance matrices from different disjoint samples taking mean differences into account.
%                                    'rereference' to subtract combined channels from all the channels.
%    theWeights     : List of weights for weighted combinations.
%    theMeans       : List of means for each sample (chan,sample).  Only for combineObs.
%    sampleSizes    : List of n for each sample.  Only for combineObs.  If empty, then equal weighting.
%
%Outputs:
%    outCovMat      : output 3D covariance matrix (chan,chan,sample).  With additional channel for combineChan and collapsed into a single sample for combineObs.
%
% B. O'Neill (2014) Some Useful Moment Results in Sampling Problems, The American Statistician, 68:4, 282-296, DOI: 10.1080/00031305.2014.966589
% Charter, R. A., & Alexander, R. A. (1993). A note on combining correlations. Bulletin of the Psychonomic Society, 31(2), 123?124. https://doi.org/10.3758/BF03334158
%
% Issues to keep in mind: 1) sample versus population.  2) within versus between combinations.  3) pooled versus combined.  4) weighting scheme.

%History
%  by Joseph Dien (12/13/19)
%  jdien07@mac.com
%
% bugfix 2/20/21 JD
% Fixed not calculating covariances correctly when present, for poolObs option.
% Now excludes NaN values from calculations.
%
% bugfix 3/31/21 JD
% Fixed crash when performing a rereference or poolObs with NaN values involved.
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outCovMat=[];
if size(inCovMat,7)>1
    inCovMat=permute(inCovMat,[1 7 2 3 4 5 6]);
end
inCovMat=reshape(inCovMat,size(inCovMat,1),size(inCovMat,2),[]);
numChans=size(inCovMat,1);
numChans2=size(inCovMat,2);
numGroups=size(inCovMat,3);
if isempty(combineList)
    if strcmp(theOperation,{'combineChan','rereference'})
        combineList=[1:numChans];
    else
        combineList=[1:numGroups];
    end
end

if (length(combineList)<2) && ~strcmp(theOperation,'rereference')
    disp('oops - programmer error in ep_covMat.  Too few to combine.')
    return
end

if ~exist('theWeights','var')
    theWeights=ones(length(combineList),1);
    if any(strcmp(theOperation,{'combineChan','rereference'}))
        theWeights=theWeights/length(combineList);
    end
elseif ~strcmp(theOperation,'combineObs')
    theWeights=(theWeights/sum(abs(theWeights)))*length(theWeights); %normalize the weights so their absolute values sum to one, then multiply by length so they correspond to unweighted case.
end
if ~exist('sampleSizes','var') || isempty(sampleSizes)
    sampleSizes=ones(length(combineList),1);
else
    sampleSizes(~sampleSizes)=1;
end
if length(theWeights)~=length(combineList)
    disp('oops - programmer error in ep_covMat.  Weights and list of combination levels do not match.')
    return
end
if strcmp(theOperation,'combineObs')
    if ~exist('theMeans','var') || (size(theMeans,2)~=length(theWeights)) || ~exist('sampleSizes','var') || (length(sampleSizes)~=length(theWeights))
        disp('oops - programmer error in ep_covMat.  Wrong number of weights or sample sizes.')
        return
    end
end
if (numChans~=numChans2) && (numChans2~=1)
    disp('oops - programmer error in ep_covMat.  The two channel dimensions are not equal size or a vector.')
    return
end

combineList=combineList(:);
theWeights=theWeights(:);
sampleSizes=sampleSizes(:);

switch theOperation
    case 'combineChan'
        if (numChans~=numChans2)
            disp('oops - programmer error in ep_covMat.  The two channel dimensions are not equal size.')
            return
        end
        outCovMat=NaN(numChans+1,numChans2+1,numGroups);
        for iObs=1:numGroups
            outCovMat(1:numChans,1:numChans2,iObs)=inCovMat(:,:,iObs);
            %The variance of the sum of a set of correlated variables is the sum of their covariances.
            %https://en.wikipedia.org/wiki/Variance#Sum_of_uncorrelated_variables_(Bienaym%C3%A9_formula) 11/21/19
            theMatrix=inCovMat(combineList,combineList,iObs);
            theNonNaNweights=theWeights(~all(isnan(theMatrix),1))/sum(theWeights(~all(isnan(theMatrix),1))); %the weights excluding NaN values. 
            outCovMat(end,end,iObs)=sum(sum(diag(theNonNaNweights)*theMatrix*diag(theNonNaNweights),'omitnan'),'omitnan');
            for iChan=1:numChans
                %Guilford, J. P. (1965). Fundamental statistics in psychology and education. p.427
                %this equation (16.25) translates to simple addition of the weighted covariances.
                theVector=inCovMat(iChan,combineList,iObs);
                theNonNaNweights=theWeights(~isnan(theVector))/sum(theWeights(~isnan(theVector))); %the weights excluding NaN values.
                theCovNewChan=sum(theVector(~isnan(theVector))'.*theNonNaNweights,'omitnan'); %covariance between the new channel and the existing channel.
                outCovMat(iChan,end,iObs)=theCovNewChan;
                outCovMat(end,iChan,iObs)=theCovNewChan;
            end
        end
    case 'combineObs'
        if any(isnan(inCovMat),'all')
            error('Have not yet spent the time to adapt this function to handle NaN values since this function is not actually being used yet.')
        end
        
        %Differential weighting, including negatives, are multiplied against the sample sizes and then normalized such that the sum of the absolute values equals the original total sample size.
        %This is necessary so that the Bessel correction (i.e., "-1") for samples still has the correct scale visa-vis the weights.
        %Negative weights are implemented as a sign reversal for the associated correlation matrix.
        WsampleSizes=((theWeights.*sampleSizes)/sum(abs(theWeights.*sampleSizes)))*sum(sampleSizes);
        
        newVar=zeros(numChans,1);
        %first compute the new combined variances.
        for iChan=1:numChans
            if numChans2==1
                theChan2=1;
            else
                theChan2=iChan;
            end
            
            %O'Neill (2014) p. 283 Result 1
            %as implemented by https://stats.stackexchange.com/questions/384941/combining-two-means-and-sds-of-one-group 12/9/19
            %recursively adds each group to the running pooled statistics.
            n1=abs(WsampleSizes(1));
            var1=inCovMat(1,1);
            mean1=theMeans(iChan,1);
            for iObs=2:numGroups
                n2=abs(WsampleSizes(iObs));
                var2=inCovMat(iChan,theChan2,iObs);
                mean2=theMeans(iChan,iObs);
                var1=(1/(n1+n2-1))*((n1-1)*var1+(n2-1)*var2+((n1*n2)/(n1+n2))*(mean1-mean2)^2); %compute pooled variance
                mean1=(1/n1+n2)*(n1*mean1+n2*mean2); %compute pooled mean
                n1=n1+n2; %compute pooled sample size
            end
            newVar(iChan)=var1;
        end
        if numChans2==1
            outCovMat=newVar;
        else
            %then convert the covariance matrices into correlation matrices
            R=eye(numChans);
            VAR=zeros(numChans,numGroups);
            n=abs(WsampleSizes);
            for i=1:numGroups
                VAR(:,i)=diag(inCovMat(:,:,i));
                R(:,:,i)=diag(sqrt(diag(inCovMat(:,:,i))).^(-1))*inCovMat(:,:,i)*diag(sqrt(diag(inCovMat(:,:,i))).^(-1))*sign(WsampleSizes(i)); %invert correlation signs for negative sample weights;
            end
            %compute the combined sample correlations per Charter & Alexander (1993)
            combinedR=eye(numChans);
            N=sum(n);
            for X=1:numChans
                EX=zeros(numGroups,1);
                EX2=zeros(numGroups,1);
                for i=1:numGroups
                    XM=theMeans(X,i);
                    XS2=VAR(X,i)*((n(i)-1)/n(i)); %adjust for sample rather than population
                    XS=sqrt(XS2);
                    EX(i)=n(i)*XM;
                    EX2(i)=n(i)*(XM^2+XS2);
                end
                for Y=X+1:numChans
                    EY=zeros(numGroups,1);
                    EY2=zeros(numGroups,1);
                    EXY=zeros(numGroups,1);
                    for i=1:numGroups
                        YM=theMeans(Y,i);
                        YS2=VAR(Y,i)*((n(i)-1)/n(i)); %adjust for sample rather than population
                        YS=sqrt(YS2);
                        r=R(X,Y,i);
                        EY(i)=n(i)*YM;
                        EY2(i)=n(i)*(YM^2+YS2);
                        EXY(i)=n(i)*(r*XS*YS+XM*YM);
                    end
                    theCor=(N*sum(EXY)-sum(EX)*sum(EY))/((N*sum(EX2)-sum(EX)^2)*(N*sum(EY2)-sum(EY)^2))^(.5);
                    combinedR(X,Y)=theCor;
                    combinedR(Y,X)=theCor;
                end
            end
            %then convert the correlations back into covariances
            outCovMat=sqrt(diag(newVar))*combinedR*sqrt(diag(newVar));
        end
        
    case 'poolObs'
        %Here the goal is to combine the covariance matrix wherein the deviations are from each sample's own mean rather than the combined mean.
        %The math is relatively simple, just weighted combination of the variances and covariances.  The only tricky thing is the weighting scheme.
        %https://en.wikipedia.org/wiki/Sample_mean_and_covariance#Weighted_samples
        
        %The covariances are weighted such that sum(weights)=1 and are >0.
        %Negative weights are implemented as a sign reversal for the associated correlation matrix.
        normWeights=abs((theWeights.*sampleSizes)/sum(abs(theWeights.*sampleSizes)));
        
        outCovMat=NaN(numChans,numChans2);
        
        %compute the variances
        for iChan=1:numChans
            if numChans2==1
                theChan2=1;
            else
                theChan2=iChan;
            end
            theVector=inCovMat(iChan,theChan2,combineList);
            theNonNaNweights=normWeights(~isnan(theVector))/sum(normWeights(~isnan(theVector))); %the weights excluding NaN values.
            outCovMat(iChan,theChan2)=sum(squeeze(theVector(~isnan(theVector))).*theNonNaNweights,'omitnan');
        end
        if numChans2~=1
            %compute the covariances
            for iChan=1:numChans
                for iChan2=iChan+1:numChans
                    theCov=sum(squeeze(inCovMat(iChan,iChan2,combineList)).*normWeights.*sign(theWeights),'omitnan');
                    outCovMat(iChan,iChan2)=theCov;
                    outCovMat(iChan2,iChan)=theCov;
                end
            end
        end
        
    case 'rereference'
        %this operation consists of first computing the virtual reference site and then subtracting it from each channel.
        %does not change the number of channels.  Implicit channels should be added beforehand.
        %if no weights were specified, the default is for normalized equal weighting for all the combineList (reference) channels.
        if (numChans~=numChans2)
            disp('oops - programmer error in ep_covMat.  The two channel dimensions are not equal size.')
            return
        end
        outCovMat=NaN(size(inCovMat));
        for iObs=1:numGroups
            %The variance of the sum of a set of correlated variables is the sum of their covariances.
            %https://en.wikipedia.org/wiki/Variance#Sum_of_uncorrelated_variables_(Bienaym%C3%A9_formula) 11/21/19
            theMatrix=inCovMat(combineList,combineList,iObs);
            theNonNaNweights=theWeights(~all(isnan(theMatrix),1))/sum(theWeights(~all(isnan(theMatrix),1))); %the weights excluding NaN values.            
            theRefVar=sum(sum(diag(theNonNaNweights)*theMatrix*diag(theNonNaNweights),'omitnan'),'omitnan');
            for iChan=1:numChans
                %Guilford, J. P. (1965). Fundamental statistics in psychology and education. p.427
                %this equation (16.25) translates to simple addition of the weighted covariances.
                theVector=inCovMat(iChan,combineList,iObs);
                theNonNaNweights=theWeights(~isnan(theVector))/sum(theWeights(~isnan(theVector))); %the weights excluding NaN values.
                theCovRefChan=sum(theVector(~isnan(theVector))'.*theNonNaNweights,'omitnan'); %covariance between reference and the channel to be rereferenced.
                %the variance of the rereferenced channel:var(X-R)=var(X)+var(R)-2*cov(X,R)
                %https://en.wikipedia.org/wiki/Variance#Sum_of_uncorrelated_variables_(Bienaym%C3%A9_formula) 12/5/19
                outCovMat(iChan,iChan,iObs)=inCovMat(iChan,iChan,iObs)+theRefVar-(2*theCovRefChan);
                %fill in the rest of the covariance matrix
                for iChan2=iChan+1:numChans
                    %Guilford, J. P. (1965). Fundamental statistics in psychology and education. p.427
                    %this equation (16.25) translates to simple addition of the weighted covariances.
                    theVector=inCovMat(iChan2,combineList,iObs);
                    theNonNaNweights=theWeights(~isnan(theVector))/sum(theWeights(~isnan(theVector))); %the weights excluding NaN values.
                    theCovRefChan2=sum(theVector(~isnan(theVector))'.*theNonNaNweights,'omitnan'); %covariance between reference and the second channel to be rereferenced.
                    %cov(X-R,Y-R)=cov(X,Y)-cov(X,R)-cov(Y,R)+VAR(R) by covariance math
                    theCov=inCovMat(iChan,iChan2,iObs)-theCovRefChan-theCovRefChan2+theRefVar;
                    outCovMat(iChan,iChan2,iObs)=theCov;
                    outCovMat(iChan2,iChan,iObs)=theCov;
                end
            end
        end
    otherwise
        disp('oops - programmer error in ep_covMat.  theOperation not recognized.')
        return
end
        
        
