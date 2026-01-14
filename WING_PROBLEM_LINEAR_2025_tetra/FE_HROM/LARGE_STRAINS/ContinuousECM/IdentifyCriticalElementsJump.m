function  DATALOC =  IdentifyCriticalElementsJump(DATALOC,iterWEIGHTS,VARCnew)


 if DATALOC.EXCLUDE_ELEMENTS_TRIGGERED_NONCONVERGENCE_ITERATION >0  && iterWEIGHTS >= DATALOC.EXCLUDE_ELEMENTS_TRIGGERED_NONCONVERGENCE_ITERATION
        IISSS =cellfun(@isempty,VARCnew.ListElementsInTransitionINNERloop)  ;
        INDEX_ITERATION_WITH_elementJUMPS = find(IISSS == 0) ;
        if ~isempty(INDEX_ITERATION_WITH_elementJUMPS)
            % First iteration with jump
            SALIRloc = 0 ; 
            iterLOCsalir = 1; 
            maxITERLOCrestricted = min(length(INDEX_ITERATION_WITH_elementJUMPS),3) ;  
            while SALIRloc == 0 && iterLOCsalir<=maxITERLOCrestricted
                Ifirst = INDEX_ITERATION_WITH_elementJUMPS(iterLOCsalir)  ;
                ElementsInTheJump = unique(cell2mat(VARCnew.ListElementsInTransitionINNERloop(Ifirst)')) ;
                %         % Now these elements are part of the critical element
                %         jumps. Accordingly, we expand
                
                NewElementsBad = setdiff( ElementsInTheJump,DATALOC.EXCLUDE_ELEMENT_TRANSITION_LIST) ;
                
                 NewElementsBadNOB = setdiff( NewElementsBad,DATALOC.BenignJumps) ;
                if  isempty(NewElementsBadNOB)
                    iterLOCsalir = iterLOCsalir +1 ; 
                else
                    
                     SALIRloc = 1; 
                    
                    disp('---------------------------------------------------')
                    disp('Updating the list of constrained elements')
                    disp(num2str(NewElementsBadNOB'))
                    disp('---------------------------------------------------')
                    DATALOC.EXCLUDE_ELEMENT_TRANSITION_LIST = [DATALOC.EXCLUDE_ELEMENT_TRANSITION_LIST;NewElementsBadNOB ] ;
                end
            end
        end
        
    end