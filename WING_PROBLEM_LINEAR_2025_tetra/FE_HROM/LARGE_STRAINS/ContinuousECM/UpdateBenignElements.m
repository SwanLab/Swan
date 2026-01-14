function DATALOC = UpdateBenignElements(VARCnew,DATALOC)


 % CONVERGENCE HAS BEEN ACHIEVED
        IISSS =cellfun(@isempty,VARCnew.ListElementsInTransitionINNERloop)  ;
        INDEX_ITERATION_WITH_elementJUMPS = find(IISSS == 0) ;
        if ~isempty(INDEX_ITERATION_WITH_elementJUMPS)
            % First iteration with jump
            
            ElementsInTheJump = unique(cell2mat(VARCnew.ListElementsInTransitionINNERloop(INDEX_ITERATION_WITH_elementJUMPS)')) ;
            % Now these elements are part of the set of "benign" jumps
             NewElements = setdiff( ElementsInTheJump,DATALOC.BenignJumps) ;
            if  ~isempty(NewElements)
            disp('---------------------------------------------------')
            disp('Updating the list of benign elements')
            disp(num2str(NewElements'))
            disp('---------------------------------------------------')
            DATALOC.BenignJumps = [DATALOC.BenignJumps;NewElements ] ; 
            end
        %    DATALOC.BenignJumps = unique([DATALOC.BenignJumps;ElementsInTheJump(:)]) ;
        end
    