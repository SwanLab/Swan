function   [qRB,qDEF,rDEF] = AssemblySimpleProjection(nDOM,DATAwmethod,INDrig,INDdef,...
    KdomREDglo,Hqr,Fbar)

if nargin == 0
    load('tmp.mat')
end

bP =  DATAwmethod.b ;
Pd =  DATAwmethod.Pd ;
Pr =  DATAwmethod.Pr ;
Td =  DATAwmethod.Td ;
cTr =  DATAwmethod.cTr ;
INDdefR =  DATAwmethod.INDdefR ;
INDrigR =   DATAwmethod.INDrigR ;



nrows = size(Td,1) +  size(Pr,1) +size(Pr,2)+ size(Hqr,2) + size(KdomREDglo,1) ;
A = sparse(nrows,nrows) ;
b = sparse(nrows,1) ;
% FIRST ROW
% -----------------------
iacum = 0;
indROWS = 1:size(KdomREDglo,1) ;
indCOLS = 1:(size(KdomREDglo,2) + size(Hqr,2)  );
A(indROWS,indCOLS) = [KdomREDglo,-Hqr] ;
b(indROWS) = Fbar ;
indCOLS= (indCOLS(end)+size(Pr,2) +1):(indCOLS(end)+size(Pr,2)+size(Pd,1));
A(indROWS,indCOLS) = Pd' ; 
% SECOND ROW
% -------------
indROWS =  indROWS(end)+[ 1:size(Hqr,2)] ;
indCOLS = 1:(size(Hqr,1) );
A(indROWS,indCOLS) = [-Hqr'] ;
indCOLS = (nrows-size(Td,1)+1):nrows;
A(indROWS,indCOLS) = Td' ;

% THIRD ROW
% -------------
 
indROWS =  indROWS(end) +[ 1:size(Pr,2)] ;
indCOLS = (nrows-size(Pr,1)-size(Td',2)+1):(nrows-size(Td',2));
A(indROWS,indCOLS) = [Pr'] ;
 
% FOURTH ROW
% --------------
 
indROWS =  indROWS(end) +[ 1:size(Pd,1)] ;
indCOLS =  1:size(Pd,2) ;
A(indROWS,indCOLS) = [Pd];
b(indROWS) = bP;
indCOLS =  (indCOLS(end)+size(Td,2) +1):(indCOLS(end)+size(Td,2)+size(Pr,2)); ;
 A(indROWS,indCOLS) = [Pr];
 
 % FIFTH ROW
% --------------
 
indROWS =  indROWS(end) +[ 1:size(Td,1)] ;
indCOLS =  (size(Pd,2)+1):(size(Pd,2)+size(Td,2)) ;
A(indROWS,indCOLS) = [Td];
b(indROWS) = -cTr;

 
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