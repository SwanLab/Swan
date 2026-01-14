function xiD = CalcDampingRatios(DATA,FREQ) ;

METHOD = 2;

if METHOD ==1
    xiDmin = DATA.xiDmin ;
    freqNmin = FREQ(1) ;
    freqNmax = FREQ(end) ;
    % Least squares
    b = [xiDmin 1]' ;
    A = [0.5/freqNmin 0.5*freqNmin
        0.5/freqNmax  0.5*freqNmax]
    
    SOL =  lsqnonneg(A,b)
    xi = (A*SOL)
    
    betaD = SOL(2) ;
    alphaD = SOL(1) ;
    % Damping ratios
    xiD = 0.5*(alphaD./FREQ + betaD*FREQ) ;
    aaaa = find(xiD >1);
    
    if ~isempty(aaaa)
        xiD(aaaa) = 1 ;
    end
    
elseif METHOD ==2
        xiDmin = DATA.xiDmin ;
 
    % Least squares
    b = xiDmin*ones(size(FREQ)); %[xiDmin 1]' ;
    A =  0.5*[1./FREQ FREQ] ; 
    
    SOL =  lsqnonneg(A,b) ; 
    xi = (A*SOL)
    
    betaD = SOL(2) ;
    alphaD = SOL(1) ;
    % Damping ratios
    xiD = 0.5*(alphaD./FREQ + betaD*FREQ) ;
    aaaa = find(xiD >1);
    
    if ~isempty(aaaa)
        xiD(aaaa) = 1 ;
    end
    
end