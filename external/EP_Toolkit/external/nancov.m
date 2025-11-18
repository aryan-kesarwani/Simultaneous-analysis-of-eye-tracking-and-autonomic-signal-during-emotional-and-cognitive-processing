%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NANCOV.M
% Program to compute a covariance matrix ignoring NaNs
%
% Usage: C = nancov(A,B)
%
% NANCOV calculates the matrix product A*B ignoring NaNs by replacing them
% with zeros, and then normalizing each element C(i,j) by the number of 
% non-NaN values in the vector product A(i,:)*B(:,j).
%
% A - left hand matrix to be multiplied
% B - right hand matrix to be multiplied
% C - resultant covariance matrix
% Example: A = [1 NaN 1] , B = [1
%                               1
%                               1]
% then nancov(A,B) is 2/2 = 1
%
% This is much better than Josh's stupid way (his words, not mine).
%  (© Kara Lavender) http://pordlabs.ucsd.edu/matlab/nan.htm

function [C]=nancov(A,B);

NmatA=~isnan(A); % Counter matrix
NmatB=~isnan(B);

A(isnan(A))=0; % Replace NaNs in A,B, and counter matrices
B(isnan(B))=0; % with zeros

Npts=NmatA*NmatB;
C=(A*B)./Npts;