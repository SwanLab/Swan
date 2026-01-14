function [qRB,qDEF,rDEF] = DomainEquationAssemblySCHUR(nDOM,nRB,nDEF,nREAC,P, b,INDrig,INDdef,...
    K,T,H,Fbar,c,DATAIN)
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
% See Implementation.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ==0
    load('tmp.mat')
end

% Offline stage. If it is already in memory, there is no need to calculate
% again 
ATIME = tic ; 
disp('---------------------------------')
disp('Computation offline matrices ...')
% 1.
 hTIL = [ P\b; -T\c] ;
% 2. Q matrix 
Q = T\H' ;  % full
 % 3. 
S = H*Q ;  
% 4. 
m  = Q'*c ; 
ATIME1 = toc(ATIME) ; 
disp(['...done in ',num2str(ATIME),' seconds'])

%%%%%%%%%%%%%%
disp('---------------------------------')
disp('Computation online matrices ...')
ATIME = tic ; 
% 1. Ktil
nrows = size(K,1) ; 
ncols = length(INDrig) + +length(INDdef) ; 
Ktil = sparse(nrows,ncols) ; 
Ktil(:,INDdef) = K ; 
% 2. 
Rtil = P\Ktil' ; 
% 3. System of equations
R = Ktil*Rtil + S ; 
a =  -Fbar +Rtil'*b + m ; 
lambda = R\a ; 




% DATAIN = DefaultField(DATAIN,'SchurComplementSolution',0) ;
% 
% if DATAIN.SchurComplementSolution == 0
%     nrows = size(PcompRR,1) +  size(PcompDD,1) + size(Treac,1) + size(KdomREDglo,1) ;
%     A = sparse(nrows,nrows) ;
%     b = sparse(nrows,1) ;
%     % FIRST ROW
%     % -----------------------
%     iacum = 0;
%     indROWS = 1:size(PcompRR,1) ;
%     indCOLS = 1:(size(PcompRR,2) + size(PcompRD,2)  );
%     A(indROWS,indCOLS) = [PcompRR,PcompRD] ;
%     b(indROWS) = bCOMPr ;
%     % SECOND ROW
%     % -------------
%     iacum = indROWS(end);
%     indROWS =  iacum +[ 1:size(PcompDD,1)] ;
%     indCOLS = 1:(size(PcompDD,2) + size(PcompRD,1)  );
%     A(indROWS,indCOLS) = [PcompRD',PcompDD] ;
%     indCOLS = (nrows-size(KdomREDglo,2)+1):nrows;
%     A(indROWS,indCOLS) = KdomREDglo' ;
%     b(indROWS) = bCOMPd ;
%     
%     % THIRD ROW
%     % -------------
%     iacum = indROWS(end);
%     indROWS =  iacum +[ 1:size(Treac,1)] ;
%     indCOLS = (nrows-size(Treac,2)-size(Hqr',2)+1):nrows;
%     A(indROWS,indCOLS) = [Treac,-Hqr'] ;
%     b(indROWS) = -cREAC;
%     
%     % FOURTH ROW
%     % --------------
%     iacum = indROWS(end);
%     indROWS =  iacum +[ 1:size(KdomREDglo,1)] ;
%     indCOLS = (size(PcompRR,2)+1):(size(PcompRR,2)+size(KdomREDglo,2) + size(Hqr,2)) ;
%     A(indROWS,indCOLS) = [KdomREDglo,-Hqr] ;
%     b(indROWS) = Fbar;
%     
%     %%% SOLUTION
% disp(['Solving reduced order system (number of eq. =',num2str(length(b)),')'])
% tic
% x = full(A\b) ;
% toc
% %
% iacum = 0 ;
% qRB = x(1:size(PcompRR,1)) ;
% iacum =  size(PcompRR,1);
% qDEF = x((iacum+1):(iacum+ size(PcompDD,2))) ;
% iacum = iacum+size(PcompDD,2);
% rDEF = x((iacum+1):(iacum+size(Treac,2))) ;
%     
% else
%     
%   [qRB,qDEF,rDEF] =  DomainEquationAssemblySCHUR(nDOM,nRB,nDEF,nREAC,PcompRR,PcompRD,PcompDD,...
%     KdomREDglo,bCOMPd,bCOMPr,Treac,Hqr,Fbar,cREAC,DATAIN) ; 
%     
% end