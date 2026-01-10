function [RA,RB,ANGLES ]= AlignedSubspaces(A,B,TOL_SVD,ANGLE)
% Adaptationof the function IntersectionSubspaces.m
% Given the column space of two matrices A and B, this function
% returns the  subspace of  col(A) and col(B) that subtends angle below
% a given tolerance
% OUTPUT: Orthogonal matrix for the column space of the intersection
% JAHO, 11-Apr-2019, 13-May-2020 . Source: Golub

if nargin == 0
    load('tmp.mat')
    B = B(:,1:end-1) ;
end

      
DATALOC.RELATIVE_SVD = 1 ;
[QA,SA,~] = SVDT(A,TOL_SVD,DATALOC) ;
[QB,SB,~] = SVDT(B,TOL_SVD,DATALOC) ;
[Yab,CTab,Zab] = SVDT(QA'*QB) ;
ANGLES = real(acos(CTab))*180/pi ;   % It gives the  same with 
nREACT = length(find(ANGLES <= ANGLE));
RA = QA*Yab(:,1:nREACT) ;
RB = QB*Zab(:,1:nREACT) ;

    
%      DATALOC.RELATIVE_SVD = 1 ;
%     [QA,SA,~] = SVDT(RA,TOL_SVD,DATALOC) ;
%     [QB,SB,~] = SVDT(RB,TOL_SVD,DATALOC) ;
%     [Yab,CTab,Zab] = SVDT(QA'*QB) ;
%      [Yba,CTba,Zba] = SVDT(QB'*QA) ;
%     ANGLES = real(acos(CTab))*180/pi ;
%     
%    
% 
% 
% sa = size(QA,2);
% sb = size(QB,2);


