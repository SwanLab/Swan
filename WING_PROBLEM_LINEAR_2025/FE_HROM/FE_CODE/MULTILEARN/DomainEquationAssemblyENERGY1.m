function   [qRB,qDEF,rDEF] =  DomainEquationAssemblyENERGY1(nDOM,DATAwmethod,INDrig,INDdef,...
    KdomREDglo,Hqr,Fbar,DATAIN) ;

if nargin ==0
    load('tmp.mat')
end

INDrigR = DATAwmethod.INDrigR ;
INDdefR = DATAwmethod.INDdefR ;
PcompBC = DATAwmethod.PcompBC ;
bCOMPbc = DATAwmethod.bCOMPbc ;
Ework = DATAwmethod.Ework ;

nrows =  length(INDrig)+ length(INDdef) + length(INDrigR)+ length(INDdefR) ;
A = sparse(nrows,nrows) ;
b = zeros(nrows,1) ;
% FIRST ROW
% -----------------------
iacum = 0;
indROWS = 1:size(KdomREDglo,1) ;
indCOLS = 1:nrows;
G = Ework'-Hqr ; 
A(indROWS,indCOLS) = [KdomREDglo+PcompBC,G ] ;
b(indROWS) = bCOMPbc + Fbar ;
% SECOND ROW
% -------------
iacum = indROWS(end);
indROWS =   (indROWS(end)+1):nrows ;
indCOLS = 1:( size(G',2));
A(indROWS,indCOLS) = [G'] ;
  

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
