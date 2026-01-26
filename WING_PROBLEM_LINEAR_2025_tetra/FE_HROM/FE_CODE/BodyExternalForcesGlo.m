function  [Fb, Nst,  DATA, wSTs_RHS, posgp_RHS] = BodyExternalForcesGlo(COOR,CN,TypeElement, fNOD,DATA)

if DATA.VECTcode  == 0
    % Standard way (elementwise)
    Fb = ComputeFb(COOR,CN,TypeElement, fNOD); Nst = [] ;
else
    % Vectorized code
    % dbstop('90')
    if   DATA.RECALCULATE_STIFFNESS == 1
        Nst = [] ; wSTs_RHS = [] ; posgp_RHS = [] ;
    else
        disp('Retrieving Nst Matrix...')
        %  dbstop('69')
        Nst = [] ; wSTs_RHS = [] ; posgp_RHS = [] ;
        load(DATA.nameWORKSPACE_Kstiff,'Nst','posgp_RHS','wSTs_RHS');
        disp('Done')
    end
    [Fb, Nst,  DATA, wSTs_RHS, posgp_RHS]= ComputeFbVect(COOR,CN,TypeElement, fNOD,DATA,Nst,wSTs_RHS,posgp_RHS);
    
    if exist(DATA.nameWORKSPACE)==0
        APPEND = '' ;
    else
        APPEND = '-append' ;
    end
    
    if  DATA.STORE_STIFFNESS ==2 
        
        if DATA.RECALCULATE_STIFFNESS == 1  
            disp('Storing Nst Matrix...')
            
            save(DATA.nameWORKSPACE,'posgp_RHS','wSTs_RHS','fNOD','Fb',APPEND);
          %  if       DATA.DO_NOT_STORE_STIFFNESS_AND_B_MATRIX ==0 
                save(DATA.nameWORKSPACE,'Nst','-append');
           % end

            
            disp('Done')
        else
            save(DATA.nameWORKSPACE,'fNOD','Fb',APPEND);
            
        end
    end
end
