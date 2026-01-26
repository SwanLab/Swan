function [qRB,qDEF,rDEF]  =    AssemblyUnconstrainedMinimization(nDOM,Pcomp, bCOMP,INDrig,INDdef,...
    KdomREDglo,Treac,Hqr,Fbar,cREAC,DATAIN)

if nargin == 0
    load('tmp.mat')
end


PcompRR = Pcomp(INDrig,INDrig) ;
PcompRD = Pcomp(INDrig,INDdef) ;
PcompDD = Pcomp(INDdef,INDdef) ;
bCOMPr = bCOMP(INDrig);
bCOMPd = bCOMP(INDdef) ;


nrows = size(PcompRR,1) +  size(PcompDD,1) + size(Treac,1)  ;
A = sparse(nrows,nrows) ;
b = sparse(nrows,1) ;
% FIRST ROW
% -----------------------
iacum = 0;
indROWS = 1:size(PcompRR,1) ;
indCOLS = 1:(size(PcompRR,2) + size(PcompRD,2)  );
A(indROWS,indCOLS) = [PcompRR,PcompRD] ;
b(indROWS) = bCOMPr ;
% SECOND ROW
% -------------
iacum = indROWS(end);
indROWS =  iacum +[ 1:size(PcompDD,1)] ;
indCOLS = 1:(size(PcompDD,2) + size(PcompRD,1)  );
A(indROWS,indCOLS) = [PcompRD',PcompDD+KdomREDglo] ;
indCOLS = (indCOLS(end)+1):nrows;
A(indROWS,indCOLS) = -Hqr ;
b(indROWS) = bCOMPd + Fbar ;

% THIRD ROW
% -------------
iacum = indROWS(end);
indROWS =  iacum +[ 1:size(Treac,1)] ;
indCOLS = (nrows-size(Treac,2)-size(Hqr',2)+1):nrows;
A(indROWS,indCOLS) = [-Hqr',Treac] ;
b(indROWS) = -cREAC;

 

%%% SOLUTION
disp(['Solving reduced order system (number of eq. =',num2str(length(b)),')'])
tic
x = full(A\b) ;
toc
%
iacum = 0 ;
qRB = x(1:size(PcompRR,1)) ;
iacum =  size(PcompRR,1);
qDEF = x((iacum+1):(iacum+ size(PcompDD,2))) ;
iacum = iacum+size(PcompDD,2);
rDEF = x((iacum+1):(iacum+size(Treac,2))) ;
