function [A,b] = DomainEquationAssemblyMULTI(nDOM,nRB,nDEF,nREAC,PcompRR,PcompRD,PcompDD,...
    KdomREDglo,bCOMPd,bCOMPr,Treac,Hqr,Fbar,cREAC)
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
% EQUATIONS
%-------------------------
%------------------------
 
nrows = size(PcompRR,1) +  size(PcompDD,1) + size(Tread,1) + size(KdomREDglo,1) ; 
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
indCOLS = (nrows-size(Treac,2)-size(Hqr',2)+1):nrows;
A(indROWS,indCOLS) = [Treac,-Hqr'] ;
b(indROWS) = -cREAC;

% FOURTH ROW
% --------------
iacum = indROWS(end);
indROWS =  iacum +[ 1:size(KdomREDglo,1)] ;
indCOLS = (size(PcompRR,2)+1):(size(PcompRR,2)+size(KdomREDglo,2) + size(Hqr,2)) ;
A(indROWS,indCOLS) = [KdomREDglo,-Hqr] ;
b(indROWS) = Fbar;