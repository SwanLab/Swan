function [qRB,qDEF,rDEF] = AssemblyStandardMINBWORK(nDOM,Pcomp, bCOMP,INDrig,INDdef,...
    KdomREDglo,Treac,Hqr,Fbar,cREAC,DATAwmethod,rRB)

if nargin ==0 
    load('tmp.mat')
end

E = DATAwmethod.Ework ; 
INDrigR = DATAwmethod.INDrigR ; 
INDdefR = DATAwmethod.INDdefR ; 
Err = E(INDrigR,INDrig); 
Erd = E(INDrigR,INDdef); 
Edr = E(INDdefR,INDrig); 
Edd = E(INDdefR,INDdef); 


PcompRR = Pcomp(INDrig,INDrig) ;
PcompRD = Pcomp(INDrig,INDdef) ;
PcompDD = Pcomp(INDdef,INDdef) ;
bCOMPr = bCOMP(INDrig);
bCOMPd = bCOMP(INDdef) ;


nrows = size(PcompRR,1) +  size(PcompDD,1) + size(Treac,1) + size(KdomREDglo,1) ;
A = sparse(nrows,nrows) ;
b = sparse(nrows,1) ;
% FIRST ROW
% -----------------------
iacum = 0;
indROWS = 1:size(PcompRR,1) ;
indCOLS = 1:(size(PcompRR,2) + size(PcompRD,2) +size(Edr',2) );
A(indROWS,indCOLS) = [PcompRR,PcompRD,Edr'] ;
b(indROWS) = bCOMPr -Err'*rRB ;
% SECOND ROW
% -------------
iacum = indROWS(end);
indROWS =  iacum +[ 1:size(PcompDD,1)] ;
indCOLS = 1:(size(PcompDD,2) + size(PcompRD,1)  +size(Edd',2) );
A(indROWS,indCOLS) = [PcompRD',PcompDD,Edd'] ;
indCOLS = (nrows-size(KdomREDglo,2)+1):nrows;
A(indROWS,indCOLS) = KdomREDglo' ;
b(indROWS) = bCOMPd-Erd'*rRB ;

% THIRD ROW
% -------------
iacum = indROWS(end);
indROWS =  iacum +[ 1:size(Treac,1)] ;
%indCOLS = (nrows-size(Treac,2)-size(Hqr',2)+1):nrows;
A(indROWS,:) = [Edr,Edd,Treac,-Hqr'] ;
b(indROWS) = -cREAC;

% FOURTH ROW
% --------------
iacum = indROWS(end);
indROWS =  iacum +[ 1:size(KdomREDglo,1)] ;
indCOLS = (size(PcompRR,2)+1):(size(PcompRR,2)+size(KdomREDglo,2) + size(Hqr,2)) ;
A(indROWS,indCOLS) = [KdomREDglo,-Hqr] ;
b(indROWS) = Fbar;

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