function  ResultantReactionsBeams(DATA,DATAOUT)

if nargin == 0
    load('tmp.mat')
end

% Computing resultant of reaction forces
DATA = DefaultField(DATA,'BasisUrb',[]) ;
if isempty(DATA.BasisUrb)
    disp('-----------------------------------------------------------------')
    disp('WARNING: Resultant of reaction forces cannot be computed (DATA.BasisUrb is absent)')
    disp('-----------------------------------------------------------------')
else
    
    if ~isempty(DATA.DOFA)
        disp('-----------------------------------------------------------------')
        disp('REACTIONS END = 1')
        React = DATAOUT.React(DATA.DOFA) ;
        
      %  if ~isempty(DATA.RotMatrixA)
         React = RotateMatrix(DATA.RotMatrixA,React)  ; 
      %  end  
        DATA = DefaultField(DATA,'BasisUrbA',DATA.BasisUrb) ;
        Resultant = DATA.BasisUrbA'*React ;
        
        disp(['Fx = ',num2str(Resultant(1))]) ; 
        disp(['Fy = ',num2str(Resultant(2))]) ;
        
        if length(Resultant) == 3
            disp(['M = ',num2str(Resultant(3))]) ; 
        else
            disp(['Fz = ',num2str(Resultant(3))]) ;
            disp(['Mx = ',num2str(Resultant(4))]) ;
            disp(['My = ',num2str(Resultant(5))]) ;
            disp(['Mz = ',num2str(Resultant(6))]) ;
        end
        disp('-----------------------------------------------------------------')
        ResultantA = Resultant ; 
        
        save(DATA.nameWORKSPACE,'ResultantA','-append')
    end
    
    
     if ~isempty(DATA.DOFB)
         
        disp('REACTIONS END = 2')
        React = DATAOUT.React(DATA.DOFB) ;
           React = RotateMatrix(DATA.RotMatrixB,React)  ; 
        DATA = DefaultField(DATA,'BasisUrbB',DATA.BasisUrb) ;
        Resultant = DATA.BasisUrbB'*React ;
        
         disp(['Fx = ',num2str(Resultant(1))]) ; 
        disp(['Fy = ',num2str(Resultant(2))]) ;
        
        if length(Resultant) == 3
            disp(['M = ',num2str(Resultant(3))]) ; 
        else
            disp(['Fz = ',num2str(Resultant(3))]) ;
            disp(['Mx = ',num2str(Resultant(4))]) ;
            disp(['My = ',num2str(Resultant(5))]) ;
            disp(['Mz = ',num2str(Resultant(6))]) ;
        end
       ResultantB = Resultant ;  
        save(DATA.nameWORKSPACE,'ResultantB','-append')
    end
    disp('-----------------------------------------------------------------')
end