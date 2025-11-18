/*==========================================================
 * ep_WJGLMml.cc - [MUHAT, STDIZER, RESULTS, MSE]=ep_WJGLMml(Y, NX, C, U, OPT1, PER, OPT2, NUMSIM, SEED, MISSING, OPT3, ALPHA, SCALE, LOC1, LOC2);
 %function [MUHAT, STDIZER, RESULTS, MSE]=ep_WJGLMml(Y, NX, C, U, OPT1, PER, OPT2, NUMSIM, SEED, MISSING, OPT3, ALPHA, SCALE, LOC1, LOC2);
 %Function for calculating robust statistics using Welch-James ADF, boostrapping, and trimmed
 %means plus winsorized variances and covariances.
 %This is a MEX-file for MATLAB.
 % compiled on OS X: mex ep_WJGLMml.cc -I/usr/local/include/eigen3/
 % compiled on Windows: mex ep_WJGLMml.cc -IY:\parallels\windows\eigen3\
 % compiled on Linux: mex ep_WJGLMml.cc -I/usr/local/include/eigen-3.4.0/
 % with boost, compiled on OS X with mex ep_WJGLMml.cc -I/usr/local/include/boost168 -I/usr/local/include/eigen3/
 %
 %Based on SAS/IML code made available by Lisa Lix at:
 %http://www.usaskhealthdatalab.ca/sas-programs/
 %
 %Converted to Matlab 9/14/15 by Joseph Dien (jdien07@mac.com)
 %Changed inverse calls to pinv call 6/20/19 JD
 %
 %Wilcox, R. R. (2001). Fundamentals of modern statistical methods. New York: Springer-Verlag.
 %
 %Keselman, H. J., Wilcox, R. R., & Lix, L. M. (2003). A generally robust approach
 %to hypothesis testing in independent and correlated groups designs
 %Psychophysiology, 40, 586-596.
 %
 %Keselman, H. J., Algina, J., Lix, L. M., Wilcox, R. R., & Deering, K. N. (2008).
 %A generally robust approach for testing hypotheses and setting confidence intervals for effect sizes.
 %Psychological Methods, 13(2), 110.
 
 %Inputs
 %  Y        : Input data (rows=subjects, columns = cells).  The first NX rows will be assigned to the first group and so forth.
 %  NX       : Number of subjects in each group as row vector.  If empty set then will assume a single group.
 %  C        : Contrast row vector for between factors (contrasts,between groups).  Set to 1 in the case where there is only one group.
 %  U        : Contrast column vector for within factors (variables, contrasts).  Numbers should sum to zero.  Set to empty set [] for analyses with no within-factors.
 %  OPT1     : Activate trimming option (rounded down).  0 = no and 1 = yes.
 %  PER      : Percentage to trim the means.  .05 is the number recommended for ERP data by Dien.
 %  OPT2     : Activate Welch-James and ADF and bootstrapping statistic.  0 = no and 1 = yes.
 %  NUMSIM    : Number of simulations used to generate bootstrapping statistic.  p-values will be unstable if too low.  50000 informally recommended.
 %  SEED     : Seed for random number generation.  0 specifies random SEED. 1000 arbitrarily suggested as SEED to ensure RESULTS are replicable.
 %  MISSING  : Number to be treated as a missing value.  Observations with missing values are dropped from the analysis.
 %  OPT3     : Provide effect sizes. 0 = no and 1 = yes.  Disabled.
 %  ALPHA    : Alpha significance level.
 %    .corrected : Alpha level corrected for multiple comparisons (not used)
 %    .uncorrected : Alpha level not corrected for multiple comparisons
 %  SCALE    : "SCALE is a scalar indicator to control the use of a scaling factor for the effect size estimator
 %             (i.e., .642 for 20% symmetric trimming) when robust estimators are adopted.
 %             It takes a value of 0 or 1; a zero indicates that no scaling factor will be used,
 %             while a 1 indicates that a scaling factor will be adopted. The default is SCALE=1.
 %  LOC1&LOC2: "If the user specifies LOC1 = 0 and LOC2 = 0, the square root of the average of the variances over the cells
 %             involved in the contrast is used as the standardizer. If the
 %             user specifies LOC1 = 99 and LOC2 = 99, no standardizer is selected."
 
 %Outputs
 %  MUHAT    : Vector of trimmed means used in calculations.
 %             For between group analyses, within group conditions vary fastest and between group conditions vary slowest.
%  STDIZER	: Winsorized sample variance-covariance matrix
 %  RESULTS    : (1) = test statistic
 %  RESULTS    : (2) = Numerator DF
 %  RESULTS    : (3) = Denominator DF
 %  RESULTS    : (4) = Significance
 %  RESULTS    : (5) = Effect size using delta-hat-star standardizer (only for contrasts with 1 df).
 %  RESULTS    : (6) = Lower confidence limit
 %  RESULTS    : (7) = Upper confidence limit
 %  RESULTS    : (8) = Scaling factor for effect size estimator if SCALE option is chosen.
 %  MSE      : Mean Squared Error.
 
 % bugfix & modified 7/5/19 JD
 % Changed inv to pinv to better handle singularity.
 % Fixed case where last observation is missing data would cause one of the reps to register as significant.
 % Fixed incorrect calculation of F using non-bootstrap option.
 % Fixed case where bootstrap sample has all the same values being treated as a significant sample rather than insignificant when calculating the bootstrap distribution.
 % Eliminated singularity counts in RESULTS.
 %
 % bugfix & modified 7/4/24 JD
 % Updated code to compile with Eigen 3.40.
 % Fixed outputting squared standard error matrix rather than variance-covariance matrix.
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
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#include "mex.h"
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <unsupported/Eigen/KroneckerProduct>
#include <cstdlib>
#include <random>
#include <cmath>
//#include <boost/math/special_functions/erf.hpp>
//#include <boost/math/distributions/normal.hpp>
#include <iostream>
#include <limits>
#include <time.h>

using namespace Eigen;

typedef Map<MatrixXd> MexMat;
typedef Map<VectorXd> MexVec;

template<typename MatrixType, int QRPreconditioner,
bool IsComplex = NumTraits<typename MatrixType::Scalar>::IsComplex>
struct svd_precondition_2x2_block_to_be_real {};

struct errorType {
    double TooFewSubjects;
    
    errorType()
    {
        TooFewSubjects=0;
    }
};

struct mnmodStruct {
    MatrixXd MUHAT;
    MatrixXd BHAT;
    MatrixXd BHATW;
    MatrixXd YT;
    VectorXd DF;
    
    mnmodStruct()
    {
        MUHAT.setZero(1,1);
        BHAT.setZero(1,1);
        BHATW.setZero(1,1);
        YT.setZero(1,1);
        DF.setZero(1);
    }
    
    mnmodStruct(MatrixXd a, MatrixXd b, MatrixXd c, MatrixXd d, VectorXd e)
    {
        MUHAT=a;
        BHAT=b;
        BHATW=c;
        YT=d;
        DF=e;
    }
};

struct sigmodStruct {
    MatrixXd SIGMA;
    MatrixXd STDIZER;
    
    sigmodStruct()
    {
        SIGMA.setZero(1,1);
        STDIZER.setZero(1,1);
    }
    
    sigmodStruct(MatrixXd a, MatrixXd b)
    {
        SIGMA=a;
        STDIZER=b;
    }
};

struct testmodStruct {
    double FSTAT;
    double DF1;
    double DF2;
    double MSE;
    
    testmodStruct()
    {
        FSTAT=0;
        DF1=0;
        DF2=0;
        MSE=0;
    }
    
    testmodStruct(double a, double b, double c, double d)
    {
        FSTAT=a;
        DF1=b;
        DF2=c;
        MSE=d;
    }
};

struct wjeffszStruct {
    double MULTP;
    double EFFSZ;
    
    wjeffszStruct()
    {
        MULTP=0;
        EFFSZ=0;
    }
    
    wjeffszStruct(double a, double b)
    {
        MULTP=a;
        EFFSZ=b;
    }
};

struct bootesStruct {
    double MULTP;
    double EFFSZ;
    
    bootesStruct()
    {
        MULTP=0;
        EFFSZ=0;
    }
    
    bootesStruct(double a, double b)
    {
        MULTP=a;
        EFFSZ=b;
    }
};

template<class Target, class Source>
Target NarrowCast(Source v)
{
    auto r = static_cast<Target>(v);
    if (static_cast<Source>(r) != v)
        mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:narrowCastFail","Narrow cast failed.");
    return r;
}

//this code snippet was provided by the following response on stackoverflow and is therefore covered by the Creative Commons license rather than GPL
//https://stackoverflow.com/questions/56877397/mex-file-implementing-eigen-library-pseudo-inverse-function-crashes
Eigen::MatrixXd PseudoInverse(Eigen::MatrixXd matrix) {
    Eigen::JacobiSVD< Eigen::MatrixXd > svd( matrix, Eigen::ComputeThinU | Eigen::ComputeThinV );
    float tolerance = 1.0e-6f * float(std::max(matrix.rows(), matrix.cols())) * svd.singularValues().array().abs()(0);
    return svd.matrixV()
    * (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal()
    * svd.matrixU().adjoint();
}
//end code snippet

//Declare functions
MatrixXd bootdat(MatrixXd Y, MatrixXd  BHAT, mwSize BOBS, VectorXd NX, std::mt19937::result_type& SEED);
MatrixXd design(MatrixXd X, VectorXd NX);
mnmodStruct mnmod(MatrixXd Y, mwSize OPT1, mwSize BOBS, mwSize WOBS, mwSize NTOT, VectorXd NX, double PER, MatrixXd X, errorType* errorflag);
sigmodStruct sigmod(MatrixXd YT, MatrixXd X, MatrixXd BHATW, VectorXd DF, mwSize BOBS, mwSize WOBS, mwSize WOBS1, VectorXd NX);
testmodStruct testmod(MatrixXd SIGMA, MatrixXd MUHAT, MatrixXd R, VectorXd DF, mwSize BOBS, mwSize WOBS, mwSize WOBS1, errorType* errorflag);
wjeffszStruct wjeffsz(mwSize BOBS, mwSize WOBS, mwSize WOBS1, mwSize NTOT, VectorXd NX, mwSize LOC1, mwSize LOC2, mwSize SCALE, mwSize OPT1, double PER, MatrixXd R, MatrixXd MUHAT, MatrixXd STDIZER);
double probf(double f,double d1, double d2);
MatrixXd bootcen(MatrixXd YB1, MatrixXd BHAT, mwSize BOBS, VectorXd NX);
double bootstat(MatrixXd YB, mwSize OPT1, MatrixXd R, mwSize BOBS, mwSize WOBS, mwSize WOBS1, mwSize NTOT, VectorXd NX, double PER, MatrixXd X, std::mt19937::result_type& SEED, errorType* errorflag);
bootesStruct bootes(MatrixXd YB, MatrixXd X, mwSize LOC1, mwSize LOC2, mwSize SCALE, mwSize OPT1, double PER, MatrixXd R, mwSize BOBS, mwSize WOBS, mwSize WOBS1, mwSize NTOT, VectorXd NX, errorType* errorflag);
mwSize random_in_range(mwSize min, mwSize max, std::mt19937::result_type& SEED);
MatrixXd PseudoInverse(MatrixXd matrix);

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    try {
        //std::cout << "Here is the matrix tempMat1:\n" << tempMat1 << std::endl;
        //mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:badInputs","test.");*/
        
        MatrixXd Y;                    /*input data matrix*/
        VectorXd NX;                    /* Nx1 input matrix of between group sizes */
        MatrixXd C;                     /* row input matrix of between contrast */
        MatrixXd U;                     /* column input matrix of within contrast */
        mwSize OPT1;                   /* input scalar */
        double PER;                    /* input scalar */
        mwSize OPT2;                   /* input scalar */
        mwSize NUMSIM;                 /* input scalar */
        std::mt19937::result_type SEED;/* input scalar seed for random number generation */
        double MISSING;                /* input scalar */
        mwSize OPT3;                   /* input scalar */
        double SCALE;                  /* input scalar */
        mwSize LOC1;                   /* input scalar */
        mwSize LOC2;                   /* input scalar */
        double alphaThresh;
        size_t Yr;                   /* rows of Y */
        size_t Yc;                   /* cols of Y */
        size_t NXc;                   /* rows of NX */
        VectorXd  RESULTS;             /*output matrix*/
        MatrixXd YB;                   /*bootstrapped Y matrix*/
        mwSize numsim_b=0;
        mwSize numsim_es=0;
        mwSize numsim_bc=0;
        VectorXd esmat;
        MatrixXd tempMat1;
        MatrixXd tempMat2;
        MatrixXd tempMat3;
        MatrixXd tempMat4;
        errorType errorflagStruct;
        errorType * errorflag;
        errorflag=&errorflagStruct;
        
        /* check for proper number of arguments */
        if(nrhs!=15) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:nrhs","Fifteen inputs required.");
        }
        if(nlhs!=4) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:nlhs","Four outputs required.");
        }
        /* make sure Y is type double */
        if( !mxIsDouble(prhs[0]) ||
           mxIsComplex(prhs[0])) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notDouble","Y must be type double.");
        }
        
        /* make sure NX is type double */
        if( !mxIsDouble(prhs[1]) ||
           mxIsComplex(prhs[1])) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notDouble","NX must be type double.");
        }
        
        /* check that number of cols in NX is 1 */
        if(mxGetM(prhs[1]) >1) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notColVector","NX must be a row vector.");
        }
        
        /* make sure C is type double */
        if( !mxIsDouble(prhs[2]) ||
           mxIsComplex(prhs[2])) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notDouble","C must be type double.");
        }
        
        /* make sure U is type double */
        if( !mxIsDouble(prhs[3]) ||
           mxIsComplex(prhs[3])) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notDouble","U must be type double.");
        }
        
        /* make sure OPT1 is scalar */
        if( !mxIsDouble(prhs[4]) ||
           mxIsComplex(prhs[4]) ||
           mxGetNumberOfElements(prhs[4])>1 ) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notScalar","OPT1 must be a scalar.");
        }
        
        /* make sure PER is scalar */
        if( !mxIsDouble(prhs[5]) ||
           mxIsComplex(prhs[5]) ||
           mxGetNumberOfElements(prhs[5])>1 ) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notScalar","PER must be a scalar.");
        }
        
        /* make sure OPT2 is scalar */
        if( !mxIsDouble(prhs[6]) ||
           mxIsComplex(prhs[6]) ||
           mxGetNumberOfElements(prhs[6])>1 ) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notScalar","OPT2 must be a scalar.");
        }
        
        /* make sure NUMSIM is scalar */
        if( !mxIsDouble(prhs[7]) ||
           mxIsComplex(prhs[7]) ||
           mxGetNumberOfElements(prhs[7])>1 ) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notScalar","NUMSIM must be a scalar.");
        }
        
        /* make sure SEED is scalar */
        if( !mxIsDouble(prhs[8]) ||
           mxIsComplex(prhs[8]) ||
           mxGetNumberOfElements(prhs[8])>1 ) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notScalar","SEED must be a scalar.");
        }
        
        /* make sure MISSING is scalar */
        if( !mxIsDouble(prhs[9]) ||
           mxIsComplex(prhs[9]) ||
           mxGetNumberOfElements(prhs[9])>1 ) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notScalar","MISSING must be a scalar.");
        }
        
        /* make sure OPT3 is scalar */
        if( !mxIsDouble(prhs[10]) ||
           mxIsComplex(prhs[10]) ||
           mxGetNumberOfElements(prhs[10])>1 ) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notScalar","OPT3 must be a scalar.");
        }
        
        /* make sure SCALE is scalar */
        if( !mxIsDouble(prhs[12]) ||
           mxIsComplex(prhs[12]) ||
           mxGetNumberOfElements(prhs[12])>1 ) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notScalar","SCALE must be a scalar.");
        }
        
        /* make sure LOC1 is scalar */
        if( !mxIsDouble(prhs[13]) ||
           mxIsComplex(prhs[13]) ||
           mxGetNumberOfElements(prhs[13])>1 ) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notScalar","LOC1 must be a scalar.");
        }
        
        /* make sure LOC2 is scalar */
        if( !mxIsDouble(prhs[14]) ||
           mxIsComplex(prhs[14]) ||
           mxGetNumberOfElements(prhs[14])>1 ) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notScalar","LOC2 must be a scalar.");
        }
        
        /* get the value of the inputs  */
        Yr = mxGetM(prhs[0]);
        Yc = mxGetN(prhs[0]);
        MexMat mapY(mxGetPr(prhs[0]), Yr, Yc);
        Y=mapY;
        
        MexVec mapNX(mxGetPr(prhs[1]), mxGetN(prhs[1]));
        NX=mapNX;
        NXc=mxGetN(prhs[1]);
        
        MexMat mapC(mxGetPr(prhs[2]), mxGetM(prhs[2]), mxGetN(prhs[2]));
        C=mapC;
        
        MexMat mapU(mxGetPr(prhs[3]), mxGetM(prhs[3]), mxGetN(prhs[3]));
        U=mapU;
        
        OPT1 = mxGetScalar(prhs[4]);
        
        PER = mxGetScalar(prhs[5]);
        
        OPT2 = mxGetScalar(prhs[6]);
        
        NUMSIM = mxGetScalar(prhs[7]);
        
        SEED = mxGetScalar(prhs[8]);
        
        MISSING = mxGetScalar(prhs[9]);
        
        OPT3 = mxGetScalar(prhs[10]);
        
        numsim_b=NUMSIM;
        numsim_es=NUMSIM;
        numsim_bc=NUMSIM;
        
        /*        StructArray const matlabStructArray = prhs[11];
         if (matlabStructArray.getType() != ArrayType::STRUCT) {
         mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:notScalar","ALPHA must be a structure.");
         }
         size_t nfields = matlabStructArray.getNumberOfFields();*/
        /*const char *fieldNames[] = {"corrected", "uncorrected"};
         mxArray* positionArray = mxGetField(prhs[11], 0, fieldNames[1] );
         alphaThresh  = (double) mxGetScalar(positionArray);*/
        
        alphaThresh = mxGetScalar(prhs[11]);
        
        SCALE = mxGetScalar(prhs[12]);
        
        LOC1 = mxGetScalar(prhs[13]);
        
        LOC2 = mxGetScalar(prhs[14]);
        
        /* create the output matrix */
        plhs[0] = mxCreateDoubleMatrix(NXc*Yc,1,mxREAL);
        Eigen::Map<Eigen::MatrixXd> outputMUHAT(mxGetPr(plhs[0]), NXc*Yc, 1); // Map the array
        
        plhs[1] = mxCreateDoubleMatrix(NXc*Yc,NXc*Yc,mxREAL);
        Eigen::Map<Eigen::MatrixXd> outputSIGMA(mxGetPr(plhs[1]), NXc*Yc, NXc*Yc); // Map the array
        
        plhs[2] = mxCreateDoubleMatrix(8,1,mxREAL);
        Eigen::Map<Eigen::MatrixXd> outputRESULTS(mxGetPr(plhs[2]), 8, 1); // Map the array
        
        plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
        Eigen::Map<Eigen::MatrixXd> outputMSE(mxGetPr(plhs[3]), 1, 1); // Map the array
        
        if ((SEED==0) || (mxIsEmpty(prhs[8]))) {
            SEED = time(NULL);
        }
        
        if ((mxIsEmpty(prhs[1])) || (NX(0)==0)) {
            NX.resize(1);
            NX[0]=NarrowCast<double>(mxGetM(prhs[1]));
        }
        
        if ((mxIsEmpty(prhs[2])) || ((C.rows()==1 && C.cols()==1 && C(0,0)==0))) {
            C.resize(1,1);
            C(0,0)=1;
        }
        
        if ((mxIsEmpty(prhs[3])) || ((U.rows()==1 && U.cols()==1 && U(0,0)==0))) {
            U.resize(Y.cols(),1);
            U.setIdentity();
        }
        
        if (NX.size() != C.cols()) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:badInputs","Number of between group cells (%d) must equal number of terms in contrast C (%d).",mxGetM(prhs[1]),mxGetN(prhs[2]));
        }
        
        if (Y.cols() != U.rows()) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:badInputs","Number of within group cells (%d) must equal number of terms in contrast U (%d).",mxGetN(prhs[0]),mxGetM(prhs[3]));
        }
        
        if ((OPT1 != 0) && (OPT1 != 1)) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:badInputs","OPT1 must equal zero or one.");
        }
        
        if ((OPT3 != 0) && (OPT3 != 1)) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:badInputs","OPT3 must equal zero or one.");
        }
        
        if (OPT3 == 1) {
            OPT3=0;
            std::cout << "OPT3 effect size option is disabled as it is too limited to be worth implementing.\n" << std::endl;
        }
        
        if ((U.cols() > 1) || (C.rows() > 1)) {
            OPT3=0; //cannot calculate effect sizes for more than 1 degree of freedom contrasts
        }
        
        mwSize iNX=0;
        mwSize iGroup=0;
        mwSize count=0;
        mwSize iCol=0;
        mwSize iRow=0;
        mwSize badFlag=0;
        mwSize simloop=0;
        
        if (!mxIsEmpty(prhs[9])) { //drop observations with missing data points
            VectorXd good = VectorXd::Zero(Y.rows());
            VectorXd newNX = VectorXd::Zero(NX.size());
            mwSize iObs=0;
            mwSize goodCount=0;
            mwSize newRow=0;
            
            for (iGroup=0; iGroup<NX.size(); iGroup++) {
                for (iObs=0; iObs<NX(iGroup); iObs++) {
                    badFlag=0;
                    for (iCol=0; iCol<Y.cols(); iCol++) {
                        if (Y(count,iCol)==MISSING) {
                            badFlag=1;
                        }
                    }
                    if (badFlag ==0) {
                        good(count)=1;
                        goodCount++;
                        newNX(iGroup)++;
                    }
                    count++;
                }
                if (newNX(iGroup)==0) {
                    mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:badInputs","A group has zero members after dropping missing data points.");
                }
            }
            
            if (goodCount<Y.rows()) {
                NX=newNX;
                MatrixXd newY = MatrixXd::Zero(goodCount,Y.cols());
                newRow=0;
                for (iRow=0; iRow<Y.rows(); iRow++) {
                    if (good(iRow)==1) {
                        newY.row(newRow)=Y.row(iRow);
                        newRow++;
                    }
                }
                Y=newY;
            }
        }
        
        //**define module to check initial specifications****;
        
        if (U.cols() > U.rows()) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:badInputs","Possible Error: Number Of Columns Of U Exceeds Number Of Rows.");
        }
        
        if (OPT1==1) {
            if (mxIsEmpty(prhs[5])) {
                PER=.20;
            }
            if (PER > .49) {
                mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:badInputs","Error: Percentage Of Trimming Exceeds Upper Limit");
            }
        }
        
        if (OPT2==1) {
            if (mxIsEmpty(prhs[5])) {
                numsim_b=999;
            }
        }
        
        if (OPT3==1) {
            if (mxIsEmpty(prhs[5])) {
                numsim_es=999;
            }
            if (mxIsEmpty(prhs[13])) {
                LOC1=1;
            }
            if (mxIsEmpty(prhs[14])) {
                LOC2=1;
            }
        }
        
        if (mxIsEmpty(prhs[11])) {
            alphaThresh=.05;
        }
        
        if (mxIsEmpty(prhs[12])) {
            SCALE=1;
        }
        
        if (mxIsEmpty(prhs[5])) {
            numsim_bc=699;
        }
        
        int totalNX=0;
        for (iNX=0; iNX<NX.size(); iNX++) {
            totalNX=totalNX+NX(iNX);
        }
        
        if (totalNX != Y.rows()) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:badInputs","Error: Total number in NX does not match Y data matrix.");
        }
        
        MatrixXd X = MatrixXd::Zero(Y.rows(),1);
        count=0;
        for (iNX=0; iNX<NX.size(); iNX++) {
            for (iGroup=0; iGroup<NX(iNX); iGroup++) {
                X(count,0)=iNX;
                count++;
            }
        }
        X=design(X, NX);
        mwSize NTOT=Y.rows(); //number of subjects
        mwSize WOBS=Y.cols(); //total number of within cells
        mwSize BOBS=X.cols(); //number of between groups
        mwSize WOBS1=WOBS-1; //one less than total number of within cells
        MatrixXd R=kroneckerProduct(C,U.transpose()); //contrast vector combining within and between contrasts
        
        double MULTP=0;
        double EFFSZ=0;
        double EFFSZB=0;
        double nanCount=0;
        double sigCount=0;
        
        //****compute Welch-James statistic****;
        mnmodStruct mnmodOut = mnmod(Y, OPT1, BOBS, WOBS, NTOT, NX, PER, X, errorflag);
        
        MatrixXd MUHAT = mnmodOut.MUHAT;
        MatrixXd BHAT = mnmodOut.BHAT;
        MatrixXd BHATW = mnmodOut.BHATW;
        MatrixXd YT = mnmodOut.YT;
        VectorXd DF = mnmodOut.DF;
        
        if (errorflag->TooFewSubjects==1) {
            mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:badInputs","Error: Too few subjects.");
        }
        
        sigmodStruct sigmodOut = sigmod(YT,X,BHATW,DF,BOBS,WOBS,WOBS1,NX);
        MatrixXd SIGMA = sigmodOut.SIGMA;
        MatrixXd STDIZER = sigmodOut.STDIZER;
        
        testmodStruct testmodOut = testmod(SIGMA,MUHAT,R,DF, BOBS, WOBS, WOBS1, errorflag);
        
        double FSTAT = testmodOut.FSTAT;
        double DF1 = testmodOut.DF1;
        double DF2 = testmodOut.DF2;
        double MSE = testmodOut.MSE;
        
        //       std::cout << "Here is the matrix SIGMA:\n" << SIGMA << std::endl;
        //       mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:badInputs","test.");
        
        RESULTS.setZero(10);
        
        if (OPT2==1) {
            MatrixXd YB1;
            double FSTATB=0;
            VectorXd FMAT;
            FMAT.setZero(numsim_b);
            double theVal=0;
            mwSize sameFlag=0;
            for (simloop=0; simloop < numsim_b; simloop++) {
                YB1 = bootdat(Y, BHAT, BOBS, NX, SEED);
                //check to see if the bootstrap randomly picked all the same observations for the simulation
                for (iCol=0; iCol < Y.cols(); iCol++) {
                    sameFlag=1;
                    for (iRow=0; iRow < Y.rows(); iRow++) {
                        if (iRow==0) {
                            theVal=Y(0,iCol);
                        }
                        else
                        {
                            if (theVal != Y(iRow,iCol)) {
                                sameFlag=0;
                            }
                        }
                    }
                    if (sameFlag==1) {
                        badFlag=1; //all of the bootstrap observations were identical
                    }
                }
                YB=bootcen(YB1,BHAT,BOBS,NX);
                FSTATB = bootstat(YB, OPT1, R, BOBS, WOBS, WOBS1, NTOT, NX, PER, X, SEED, errorflag);
                FMAT(simloop)=FSTATB;
            }
            std::sort(FMAT.data(),FMAT.data()+FMAT.size());
            FMAT.reverse();
            Matrix<bool, Dynamic, 1> avec=FMAT.array()>=FSTAT;
            nanCount=0;
            sigCount=0;
            for (simloop=0; simloop < FMAT.size(); simloop++) {
                if (FMAT(simloop)==mxGetNaN()) {
                    nanCount++;
                }
                if (avec(simloop)) {
                    sigCount++;
                }
            }
            if (numsim_b==nanCount) {
                RESULTS(3)=mxGetNaN();
            }
            else
            {
                RESULTS(3)=sigCount/(numsim_b-nanCount);
            }
        }
        
        //***calculate significance level for welch-james statistic****;
        RESULTS(0)=FSTAT;
        RESULTS(1)=DF1;
        RESULTS(2)=DF2;
        if (OPT2==0) {
            RESULTS(3)=1-probf(RESULTS(0),DF1,DF2);
        }
        RESULTS(4) = mxGetNaN();
        RESULTS(5) = mxGetNaN();
        RESULTS(6) = mxGetNaN();
        RESULTS(7) = mxGetNaN();
        
        if (OPT3==1) {
            if ((DF1>1) || (WOBS>1)) {
                //Effect sizes only available for one DF between-group contrasts, pending further research by Dr. Lix.
                RESULTS(4) = mxGetNaN();
                RESULTS(5) = mxGetNaN();
                RESULTS(6) = mxGetNaN();
                RESULTS(7) = mxGetNaN();
            }
            else
            {
                //Effect sizes only available for one DF between-group contrasts, pending further research by Dr. Lix.
                esmat.setZero(numsim_es);
                MatrixXd YB1;
                for (simloop=0; simloop<numsim_es; simloop++) {
                    YB1 = bootdat(Y, BHAT, BOBS, NX, SEED);
                    bootesStruct bootesOut=bootes(YB1,X,LOC1,LOC2,SCALE,OPT1,PER,R,BOBS,WOBS,WOBS1,NTOT,NX,errorflag);
                    MULTP = bootesOut.MULTP;
                    EFFSZB = bootesOut.EFFSZ;
                    esmat(simloop)=EFFSZB;
                }
                std::sort(esmat.data(),esmat.data()+esmat.size());
                esmat.reverse();
                nanCount=0;
                for (simloop=0; simloop < esmat.size(); simloop++) {
                    if (esmat(simloop)==mxGetNaN()) {
                        nanCount++;
                    }
                }
                mwSize numsim_esGood=numsim_es-nanCount;
                if (numsim_esGood>0) {
                    Index ind1=trunc(numsim_esGood*(alphaThresh/2))+1;
                    Index ind2=numsim_esGood-trunc((numsim_esGood*(alphaThresh/2)));
                    double lcl=esmat(ind1);
                    double ucl=esmat(ind2);
                    bootesStruct bootesOut=bootes(Y,X,LOC1,LOC2,SCALE,OPT1,PER,R,BOBS,WOBS,WOBS1,NTOT,NX,errorflag);
                    MULTP = bootesOut.MULTP;
                    EFFSZB = bootesOut.EFFSZ;
                    RESULTS(4) = fabs(EFFSZB); //convention for Cohen's d is to provide absolute value JD
                    RESULTS(5) = lcl;
                    RESULTS(6) = ucl;
                    if (EFFSZB<0) {
                        RESULTS(5)=-RESULTS(5);
                        RESULTS(6)=-RESULTS(6);
                    }
                    RESULTS(7) = MULTP;
                }
            }
        }
        outputMUHAT=MUHAT;
        outputSIGMA=STDIZER;
        outputRESULTS=RESULTS;
        outputMSE(0)=MSE;
        
    } catch (std::exception& ex) {
        /* In case of any exception, issue a MATLAB exception with
         * the text of the C++ exception. This terminates the MEX
         * file without returning any data back to MATLAB.
         */
        mexErrMsgTxt(ex.what());
    }
}

//Design function
MatrixXd design(MatrixXd X, VectorXd NX)
{
    Index iObs=0;
    Index theX=0;
    MatrixXd a = MatrixXd::Zero(X.rows(),NX.size());
    for (iObs=0; iObs<X.rows(); iObs++) {
        theX=X(iObs);
        a(iObs,theX) = 1;
    }
    return a;
}

mnmodStruct mnmod(MatrixXd Y, mwSize OPT1, mwSize BOBS, mwSize WOBS, mwSize NTOT, VectorXd NX, double PER, MatrixXd X, errorType* errorflag)
//****define module to compute least squares or trimmed means****;
{
    MatrixXd MUHAT=MatrixXd::Zero(BOBS*WOBS,1);
    MatrixXd BHAT=MatrixXd::Zero(BOBS, WOBS);
    MatrixXd BHATW=MatrixXd::Zero(BOBS, WOBS);
    MatrixXd YT=MatrixXd::Zero(NTOT, WOBS);
    VectorXd DF=VectorXd::Zero(BOBS);
    MatrixXd tempMat1;
    MatrixXd tempMat2;
    mwSize J=0;
    mwSize K=0;
    
    if (OPT1==0) {
        tempMat1 = X.transpose()*X;
        tempMat2 = PseudoInverse(tempMat1);
        BHAT=tempMat2*X.transpose()*Y;
        BHATW=BHAT;
        YT=Y;
        DF=NX.array()-1;
    }
    
    if (OPT1==1) {
        mwSize F=1;
        mwSize M=0;
        mwSize L=0;
        mwSize P=0;
        mwSize G=0;
        mwSize SAMP=0;
        VectorXd NV;
        VectorXd NV2;
        VectorXd TRIMY;
        double TRIMMN=0;
        double MINT=0;
        double MAXT=0;
        double WINMN=0;
        
        for (J=0; J<NX.size(); J++) { //loop through between groups
            SAMP=NX(J);//size of between group
            L=M+SAMP;
            G=trunc(PER*SAMP);
            DF(J)=SAMP-2*G-1;
            for (K=0; K<Y.cols(); K++) { //loop through within groups
                NV=Y.block(F-1,K,L-F+1,1);
                NV2=NV;
                std::sort(NV2.data(),NV2.data()+NV2.size());
                NV2.reverse();
                TRIMY=NV2.segment(G,SAMP-(2*G)); //trimmed cell
                TRIMMN=TRIMY.sum()/(DF(J)+1); //trimmed mean
                BHAT(J,K)=TRIMMN; //matrix of trimmed means
                MINT=TRIMY.minCoeff();
                MAXT=TRIMY.maxCoeff();
                for (P=0; P<NV.size(); P++) {
                    if (NV(P)<=MINT) {
                        NV(P)=MINT;
                    }
                    if (NV(P)>=MAXT) {
                        NV(P)=MAXT;
                    }
                }
                YT.block(F-1,K,L-F+1,1)=NV; //winsorized sample
                WINMN=NV.sum()/SAMP;
                BHATW(J,K)=WINMN;
            }
            M=L;
            F=F+NX(J);
        }
    }
    
    Matrix<bool, Dynamic, 1> testBad = DF.array()==0;
    if (testBad.any()) {
        errorflag->TooFewSubjects=1;
        MUHAT.resize(1,1);
        MUHAT(0,0)=0;
        BHAT.resize(1,1);
        BHAT(0,0)=0;
        BHATW.resize(1,1);
        BHATW(0,0)=0;
        YT.resize(1,1);
        YT(0,0)=0;
        mexErrMsgIdAndTxt("MyToolbox:ep_WJGLMml:badInputs","Error, too few subjects.  Degrees of freedom is zero.");
    }
    
    for (J=0; J<BOBS; J++) { //loop through between groups
        for (K=0; K<WOBS; K++) { //loop through within groups
            MUHAT((J*WOBS)+K,0)=BHAT(J,K);
        }
    }
    return {MUHAT,BHAT,BHATW,YT,DF};
}

sigmodStruct sigmod(MatrixXd YT, MatrixXd X, MatrixXd BHATW, VectorXd DF, mwSize BOBS, mwSize WOBS, mwSize WOBS1, VectorXd NX)
//***DEFINE MODULE TO COMPUTE SIGMA MATRIX****;
{
    mwSize I=0;
    mwSize F=0;
    mwSize L=0;
    MatrixXd SIGB;
    MatrixXd SIGMA=MatrixXd::Zero(WOBS*BOBS, WOBS*BOBS);
    MatrixXd STDIZER=SIGMA;
    MatrixXd tempMat1;
    MatrixXd tempMat2;
    MatrixXd tempMat3;
    MatrixXd tempMat4;
    
    for (I=0; I<BOBS; I++) {
        tempMat1=X.col(I);
        tempMat2=tempMat1.asDiagonal()*YT-X.col(I)*BHATW.row(I);
        tempMat3=X.col(I);
        tempMat4=tempMat3.asDiagonal();
        SIGB=tempMat2.transpose()*(tempMat4*YT-X.col(I)*BHATW.row(I))/((DF(I)+1)*DF(I));
        F=(I+1)*WOBS-WOBS1-1;
        L=(I+1)*WOBS-1;
        SIGMA.block(F,F,L-F+1,L-F+1)=SIGB;
        STDIZER.block(F,F,L-F+1,L-F+1)=SIGB*((DF(I)+1)*DF(I))/(NX(I)-1);
    }
    
    return {SIGMA,STDIZER};
}


testmodStruct testmod(MatrixXd SIGMA, MatrixXd MUHAT, MatrixXd R, VectorXd DF, mwSize BOBS, mwSize WOBS, mwSize WOBS1, errorType* errorflag)
//****DEFINE MODULE TO COMPUTE TEST STATISTIC****;
{
    mwSize I=0;
    mwSize J=0;
    mwSize F=0;
    mwSize L=0;
    mwSize C=0;
    mwSize iContrast=0;
    MatrixXd IMAT;
    MatrixXd QMAT;
    MatrixXd PROD;
    double SST=0;
    double MST=0;
    double MSE=0;
    double FSTAT=0;
    double T=0;
    double CVAL=0;
    double DF1=0;
    double DF2=0;
    double A=0;
    MatrixXd tempMat1;
    MatrixXd tempMat2;
    MatrixXd tempMat3;
    MatrixXd tempMat4;
    MatrixXd tempMat5;
    MatrixXd tempMat6;
    
    tempMat1=R*MUHAT;
    tempMat2=R*SIGMA*R.transpose();
    tempMat3=tempMat2.inverse();
    tempMat4=tempMat1.transpose()*tempMat3*(R*MUHAT); //T stat squared and without df correction Twj statistic. Johansen (1980)
    T=tempMat4(0,0);
    IMAT=MatrixXd::Identity(WOBS, WOBS);
    for (I=0; I<BOBS; I++) {
        QMAT=MatrixXd::Zero(BOBS*WOBS, BOBS*WOBS);
        F=(I+1)*WOBS-WOBS1-1;
        L=(I+1)*WOBS-1;
        QMAT.block(F,F,L-F+1,L-F+1)=IMAT;
        tempMat1=R*SIGMA*R.transpose();
        tempMat2=PseudoInverse(tempMat1);
        PROD=(SIGMA*R.transpose())*tempMat2*R*QMAT;
        tempMat1=PROD*PROD;
        //ADF statistic Lix & Keselman (1995)
        A=A+(tempMat1.trace()+pow(PROD.trace(),2))/DF(I);
    }
    A=A/2; //for adjusted degrees of freedom
    DF1=R.rows();
    DF2=DF1*(DF1+2)/(3*A);
    CVAL=DF1+2*A-6*A/(DF1+2);
    FSTAT=T/CVAL;
    SST=0;
    
    for (iContrast=0; iContrast<DF1; iContrast++) {
        tempMat1=R.row(iContrast)*MUHAT;
        tempMat2=R.row(iContrast);
        tempMat3=tempMat2.array().pow(2);
        tempMat4=MatrixXd::Ones(1, WOBS*BOBS);
        C=0;
        for (I=0; I<BOBS; I++) {
            for (J=0; J<WOBS; J++) {
                tempMat4(C)=DF(I)+1;
                C++;
            }
        }
        tempMat5=tempMat3.array()/tempMat4.array();
        tempMat6=(tempMat1.transpose()*(R.row(iContrast)*MUHAT));
        SST=SST+tempMat6(0)/tempMat5.sum();
    }
    MST=SST/DF1;
    MSE=MST/FSTAT;
    
    return {FSTAT,DF1,DF2,MSE};
    
}

//****DEFINE MODULES TO PERFORM BOOTSTRAP****;
//***DEFINE MODULE TO GENERATE BOOTSTRAP DATA AND CENTRE DATA****;
//resamples with replacement the observations for each between-group, resulting in a resampled dataset with the same overall dimensions.
MatrixXd bootdat(MatrixXd Y, MatrixXd  BHAT, mwSize BOBS, VectorXd NX, std::mt19937::result_type& SEED)
//resamples with replacement the observations for each between-group, resulting in a resampled dataset with the same overall dimensions.
{
    MatrixXd BVAL;
    mwSize J=0;
    mwSize P=0;
    mwSize tempRows=0;
    
    mwSize F=1;
    mwSize M=0;
    mwSize L=0;
    mwSize RVAL=0;
    MatrixXd tempMat1;
    
    MatrixXd YB(Y.rows(),Y.cols()); //output variable
    
    for (J=0; J<BOBS; J++) {
        L=M+NX[J];
        tempRows=L-F+1;
        tempMat1.resize(tempRows,Y.cols());
        tempMat1.block(0,0,tempRows,Y.cols())=Y.block(F-1,0,tempRows,Y.cols());
        BVAL.resize(tempRows,Y.cols());
        for (P=0; P<(tempRows); P++) {
            RVAL=random_in_range(0, tempRows-1, SEED);
            BVAL.block(P,0,1,Y.cols())=tempMat1.block(RVAL,0,1,Y.cols());
        }
        
        mwSize BVALr=BVAL.rows();
        mwSize BVALc=BVAL.cols();
        YB.block(F-1,0,tempRows,Y.cols())=BVAL;
        M=L;
        F=F+NX[J];
    }
    return YB;
}

//****CENTRE THE BOOTSTRAP DATA****;
//center the bootstrap sample based on the means from the full sample
MatrixXd bootcen(MatrixXd YB1, MatrixXd BHAT, mwSize BOBS, VectorXd NX)
{
    mwSize M=0;
    mwSize F=0;
    int L=-1;
    mwSize I=0;
    mwSize Q=0;
    MatrixXd YB;
    MatrixXd MVAL;
    
    YB=YB1;
    for (I=0; I<BOBS; I++) {
        L=M+NX(I);
        MVAL=BHAT.row(I);
        for (Q=F; Q<L; Q++) {
            YB.row(Q)=YB1.row(Q)-MVAL;
        }
        M=L;
        F=F+NX(I);
    }
    return YB;
}

//****DEFINE MODULE TO COMPUTE BOOTSTRAP STATISTIC****;
double bootstat(MatrixXd YB, mwSize OPT1, MatrixXd R, mwSize BOBS, mwSize WOBS, mwSize WOBS1, mwSize NTOT, VectorXd NX, double PER, MatrixXd X, std::mt19937::result_type& SEED, errorType* errorflag)
{
    mnmodStruct mnmodOut = mnmod(YB, OPT1, BOBS, WOBS, NTOT, NX, PER, X, errorflag);
    MatrixXd MUHATB = mnmodOut.MUHAT;
    MatrixXd BHATB = mnmodOut.BHAT;
    MatrixXd BHATBW = mnmodOut.BHATW;
    MatrixXd YTB = mnmodOut.YT;
    VectorXd DFB = mnmodOut.DF;
    
    sigmodStruct sigmodOut = sigmod(YTB,X,BHATBW,DFB,BOBS,WOBS,WOBS1,NX);
    MatrixXd SIGMAB =sigmodOut.SIGMA;
    MatrixXd STDIZER = sigmodOut.STDIZER;
    
    
    testmodStruct testmodOut = testmod(SIGMAB,MUHATB,R,DFB, BOBS, WOBS, WOBS1, errorflag);
    double FSTATB = testmodOut.FSTAT;
    double DF1 = testmodOut.DF1;
    double DF2 = testmodOut.DF2;
    double MSE = testmodOut.MSE;
    
    return FSTATB;
}

//****define module to compute bootstrap effect size****;
bootesStruct bootes(MatrixXd YB, MatrixXd X, mwSize LOC1, mwSize LOC2, mwSize SCALE, mwSize OPT1, double PER, MatrixXd R, mwSize BOBS, mwSize WOBS, mwSize WOBS1, mwSize NTOT, VectorXd NX, errorType* errorflag)
{
    mnmodStruct mnmodOut = mnmod(YB, OPT1, BOBS, WOBS, NTOT, NX, PER, X, errorflag);
    MatrixXd MUHATB = mnmodOut.MUHAT;
    MatrixXd BHATB = mnmodOut.BHAT;
    MatrixXd BHATBW = mnmodOut.BHATW;
    MatrixXd YTB = mnmodOut.YT;
    VectorXd DFB = mnmodOut.DF;
    
    sigmodStruct sigmodOut = sigmod(YTB,X,BHATBW,DFB,BOBS,WOBS,WOBS1,NX);
    MatrixXd SIGMAB =sigmodOut.SIGMA;
    MatrixXd STDIZERB = sigmodOut.STDIZER;
    
    wjeffszStruct wjeffszOut = wjeffsz(BOBS,WOBS,WOBS1,NTOT,NX,LOC1,LOC2,SCALE,OPT1,PER,R,MUHATB,STDIZERB);
    double MULTP = wjeffszOut.MULTP;
    double EFFSZ = wjeffszOut.EFFSZ;
    
    return {MULTP,EFFSZ};
}


wjeffszStruct wjeffsz(mwSize BOBS, mwSize WOBS, mwSize WOBS1, mwSize NTOT, VectorXd NX, mwSize LOC1, mwSize LOC2, mwSize SCALE, mwSize OPT1, double PER, MatrixXd R, MatrixXd MUHAT, MatrixXd STDIZER)
//****compute measure of effect size and bootstrap confidence interval****;
{
    double MULTP=0;
    double EFFSZ=0;
    mwSize J=0;
    double stdz=0;
    MatrixXd tempMat1;
    
    if (OPT1==0) {
        MULTP=1;
    }
    
    if (OPT1==1) {
        if (SCALE==0) {
            MULTP=1;
        }
        /*       if (SCALE==1) {
         if (PER != 0) {
         double cut=sqrt(2)*boost::math::erf_inv(2*PER-1); //probit per wikipedia https://en.wikipedia.org/wiki/Probit 12/6/2015
         // standard normal distribution object:  should actually have (Z.^2).*(1/(sqrt(2*3.141592653589793)).*(exp(-(1/2).*Z.^2))); term but not worth figuring out how to implement
         boost::math::normal normDist;
         double b=cdf(complement(normDist, cut))-cdf(complement(normDist, -cut));
         double winvar=b+PER*(2*(pow(cut,2)));
         MULTP=sqrt(winvar);
         //std::cout << "Here is the matrix winvar:\n" << winvar << std::endl;//
         //std::cout << "Here is the matrix b:\n" << b << std::endl;//
         //std::cout << "Here is the matrix PER:\n" << PER << std::endl;//
         //std::cout << "Here is the matrix cut:\n" << cut << std::endl;//
         }
         else
         {
         MULTP=1;
         }
         }*/
    }
    
    tempMat1=R*MUHAT;
    double num=tempMat1(1); // will always be a scalar since only single df can be calculated at present.
    if (LOC1==99) {
        stdz=1;
    }
    if (LOC1==0) {
        mwSize K=0;
        MatrixXd r2=R.array().pow(2);
        MatrixXd rvec=MatrixXd::Zero(1,BOBS*WOBS);
        for (J=0; J<BOBS; J++) { //loop through between groups
            for (K=0; K<WOBS; K++) { //loop through within groups
                rvec(0,(J*WOBS)+K)=r2(J,K);
            }
        }
        MatrixXd stdz1=STDIZER.diagonal();
        MatrixXd stdz2=stdz1.cwiseSqrt();
        
        tempMat1=(rvec*stdz2.asDiagonal()*rvec.transpose())/rvec.sum(); //average of square root of variances
        stdz=tempMat1(1);
    }
    if (LOC1>0) {
        if (LOC1 <99) {
            mwSize loc=LOC1*WOBS-(WOBS-LOC2);
            double stdz3 = STDIZER(loc,loc);
            if (stdz3> 0) {
                stdz=sqrt(stdz3);
            }
            if (stdz3==0) {
                stdz=.00001;
            }
        }
    }
    
    //if ((num.size()>1) || (WOBS > 1)) {
    if (WOBS > 1) {
        //Effect sizes only available for between-group contrasts, pending further research by Dr. Lix.
        EFFSZ=mxGetNaN();
    }
    else
    {
        EFFSZ=MULTP*(num/stdz);
    }
    
    return {MULTP, EFFSZ};
}

double probf(double f,double d1, double d2)
//Based on code from matrixlab-examples.com, implements SAS probf function.
{
    double x = 1;
    double s = 0;
    double t = 0;
    double z = 0;
    double j = 0;
    double k = 0;
    double y = 0;
    double a1 = 0.196854;
    double a2 = 0.115194;
    double a3 = 0.000344;
    double a4 = 0.019527;
    //Computes using inverse for small F-values
    if (f < 1) {
        s = d2;
        t = d1;
        z = 1/f;
    }
    else
    {
        s = d1;
        t = d2;
        z = f;
    }
    j = 2/(9*s);
    k = 2/(9*t);
    
    //Uses approximation formulas
    y = fabs((1 - k)*pow(z,(1.0/3.0)) - 1 + j)/sqrt(k*pow(z,(2.0/3.0)) + j);
    if (t < 4) {
        y = y*(1 + 0.08*pow(y,4)/pow(t,3));
    }
    
    x = 0.5/pow((1 + y*(a1 + y*(a2 + y*(a3 + y*a4)))),4);
    x = floor(x*10000 + 0.5)/10000;
    
    //Adjusts if inverse was computed
    if (f < 1) {
        x = 1 - x;
    }
    
    x = 1 - x; //SAS probf provides percentile not tail end value.
    
    return x;
}

mwSize random_in_range(mwSize min, mwSize max, std::mt19937::result_type& SEED) {
    //codereview.stackexchange.com/questions/101525/random-number-generation-seeding-in-c
    //random number generator
    
    thread_local static std::mt19937 mt(SEED);
    thread_local static std::uniform_int_distribution<mwSize> pick;
    
    // assuming param_type is lighter weight to construct
    // than a uniform_int_distribution
    return pick(mt, decltype(pick)::param_type{min, max});
}
