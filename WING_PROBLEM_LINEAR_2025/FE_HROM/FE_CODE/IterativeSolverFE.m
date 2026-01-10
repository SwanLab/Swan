function dL =  IterativeSolverFE(DATA,Kll,Fl,KlrU)


DATA =DefaultField(DATA,'tolCONJG',[]) ; 
DATA =DefaultField(DATA,'niterCONJG',[]) ; 

if DATA.TYPESOLVER  ==1
    
    
%    error('This solver has not proved efficient. Set   DATA.TYPESOLVER = 0 (direct) ir DATA.TYPESOLVER = 2 (minres) ')
    disp('Conjugated gradient')
    DATA = DefaultField(DATA,'USE_PRECOND','NO') ; % = 'CHOLESKY' ; )
    
    switch DATA.USE_PRECOND
        case 'NO'
            dL=    pcg(Kll,Fl-KlrU,DATA.tolCONJG,DATA.niterCONJG) ;
        case 'CHOLESKY'
            disp('Incomplete cholesky')
            try
                LL = ichol(Kll,struct('michol','on'));
            catch
                alpha = .1;
                alpha = max(sum(abs(Kll),2)./diag(Kll))-2 ;
                LL = ichol(Kll, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
            end
            disp('Solving pcg ....')
            dL=    pcg(Kll,Fl-KlrU,DATA.tolCONJG,DATA.niterCONJG,LL,LL') ;
            
                
    end
    
elseif   DATA.TYPESOLVER  ==2
     disp('MINRES solver')
    Fl = Fl-KlrU; 
    dL = minres(Kll,Fl,DATA.tolCONJG,DATA.niterCONJG) ; 
    
else('Option not implemented')
    
end