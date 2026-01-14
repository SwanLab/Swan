function [qRB,qDEF,rDEF] =DomEqAsseConstrReaction(nDOM,Pcomp, bCOMP,INDrig,INDdef,...
    KdomREDglo,Treac,Hqr,Fbar,cREAC,DATAIN,DATAwmethod)

if nargin  == 0
    load('tmp.mat')
end

PcompRR = Pcomp(INDrig,INDrig) ;
PcompRD = Pcomp(INDrig,INDdef) ;
PcompDD = Pcomp(INDdef,INDdef) ;
bCOMPr = bCOMP(INDrig);
bCOMPd = bCOMP(INDdef) ;
JintDEF = DATAwmethod.JintDEF ; 
Jr_r = DATAwmethod.Jr_r ; 
INDrigR = DATAwmethod.INDrigR ; 
INDdefR = DATAwmethod.INDdefR ; 


nrows = size(PcompRR,1) +  size(PcompDD,1) + size(Treac,1) + size(KdomREDglo,1) + size(JintDEF,1) ;
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
A(indROWS,indCOLS) = [PcompRD',PcompDD] ;
indCOLS = (nrows-size(KdomREDglo,2)+1):nrows;
A(indROWS,indCOLS) = KdomREDglo' ;
b(indROWS) = bCOMPd ;

% THIRD ROW
% -------------
iacum = indROWS(end);
indROWS =  iacum +[ 1:size(Treac,1)] ;
indCOLS = (nrows-size(Treac,2)-size(Hqr',2)-size(JintDEF,1)+1):nrows;
A(indROWS,indCOLS) = [Treac,-Hqr',JintDEF'] ;
b(indROWS) = -cREAC;

% FOURTH ROW
% --------------
iacum = indROWS(end);
indROWS =  iacum +[ 1:size(KdomREDglo,1)] ;
indCOLS = (size(PcompRR,2)+1):(size(PcompRR,2)+size(KdomREDglo,2) + size(Hqr,2)) ;
A(indROWS,indCOLS) = [KdomREDglo,-Hqr] ;
b(indROWS) = Fbar;

% FIFTH ROW
% --------------
iacum = indROWS(end);
indROWS =  iacum +[ 1:size(JintDEF,1)] ;
indCOLS = (size(PcompRR,2)+size(PcompRD,2)+1):(size(PcompRR,2)+size(PcompRD,2)+size(JintDEF,2)) ;
A(indROWS,indCOLS) = [JintDEF] ;
b(indROWS) = -Jr_r;

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