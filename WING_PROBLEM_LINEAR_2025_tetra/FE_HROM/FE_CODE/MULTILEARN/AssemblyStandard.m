function [qRB,qDEF,rDEF] = AssemblyStandard(nDOM,Pcomp, bCOMP,INDrig,INDdef,...
    KdomREDglo,Treac,Hqr,Fbar,cREAC,DATAIN) 



% DATAIN = DefaultField(DATAIN,'DO_NO_INCLUDE_MIN_REACTIONS',0) ; 
% if DATAIN.DO_NO_INCLUDE_MIN_REACTIONS == 1
%     Treac = sparse(size(Treac,1),size(Treac,2)) ; 
%     cREAC = sparse(size(cREAC,1),size(cREAC,2)) ;
% end

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
    
    %%% SOLUTION
    disp(['Solving reduced order system (number of eq. =',num2str(length(b)),')'])
    tic
    
    DATAIN = DefaultField(DATAIN,'TYPESOLVER',0) ; 
    
    if DATAIN.TYPESOLVER  == 0
        x = full(A\b) ;
    elseif DATAIN.TYPESOLVER  ==1
        disp('Conjugated gradient')
        DATAIN = DefaultField(DATAIN,'USE_PRECOND','NO') ; % = 'CHOLESKY' ; )
                DATAIN = DefaultField(DATAIN,'tolCONJG',1e-6) ; % = 'CHOLESKY' ; )
                DATAIN = DefaultField(DATAIN,'niterCONJG',6000) ; % = 'CHOLESKY' ; )

        switch DATAIN.USE_PRECOND
            case 'NO'
                x=    pcg(A,b,DATAIN.tolCONJG,DATAIN.niterCONJG) ;
            case 'CHOLESKY'
%                 LL = ichol(Kll,struct('michol','on'));
%                 dL=    pcg(Kll,Fl-KlrU,DATA.tolCONJG,DATA.niterCONJG,LL,LL') ;
        end
    end
    
  
    toc
    %
    iacum = 0 ;
    qRB = x(1:size(PcompRR,1)) ;
    iacum =  size(PcompRR,1);
    qDEF = x((iacum+1):(iacum+ size(PcompDD,2))) ;
    iacum = iacum+size(PcompDD,2);
    rDEF = x((iacum+1):(iacum+size(Treac,2))) ;