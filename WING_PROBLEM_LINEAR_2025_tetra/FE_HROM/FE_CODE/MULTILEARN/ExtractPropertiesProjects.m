function [KdomRED,FORCE_PROJECTS,DATA,Cglo,MaterialType,nMAT,K] = ...
    ExtractPropertiesProjects(PROJECT_LOADS,COORref,BdomRED,Wdom,DATAINM,setPoints,WdomRED)


FORCE_PROJECTS = {} ;
for iproject = 1:length(PROJECT_LOADS)
    FUNinput = [] ;     DATALOC = [] ; 
    eval(PROJECT_LOADS{iproject}) ;
    nMAT = length(MATERIAL.PLY) ;
    DATA.MATERIAL = MATERIAL ; 
    DATALOC.INPUTDATAfile = PROJECT_LOADS{iproject} ;
    DATALOC.NOCALCULATE_DISPLACEMENTS = 1 ;
    % Calling Finite Element elastostatic program (but only for computing external forces)
    DATAOUT = FE_ELASTOSTATIC(FUNinput,DATALOC) ;
    load(DATAOUT.nameWORKSPACE,'Fb','Ftrac','CN','COOR','MaterialType','K') ;    
    %  load(DATAOUT.nameWORKSPACE,'Fb','Ftrac','CN','COOR','MaterialType','K') ;
    [IDX D]= knnsearch(COOR,COORref) ;    
    if any(abs(D) >1e-16 )
      %  dbstop('51')
        error('Non-conforming meshes')
    end
    IDXdofs = Nod2DOF(IDX,size(COOR,2)) ;
    % All properties are set in terms of the numbering of reference mesh
    %  dbstop('55')
    if iproject == 1
        load(DATAOUT.nameWORKSPACE,'Cglo') ; % global elasticity matrix
        KdomRED = {} ;
        for itype = 1:length(BdomRED)
            nstrain = size(BdomRED{1},1)/length(Wdom) ;
            if  DATAINM.CUBATURE.ACTIVE == 1
                setIndices =  small2large(setPoints{itype},nstrain) ;
                % Cglo is multiplied by Wdom(setPoints). Accordingly, we
                % define
                rW = WdomRED{itype}./Wdom(setPoints{itype}) ;
                % And make
                rwDIAG = CompWeightDiag(rW,nstrain)  ;
                KdomRED{itype} = (rwDIAG*BdomRED{itype}(setIndices,:))'*(Cglo(setIndices,setIndices)*BdomRED{itype}(setIndices,:)) ;
                
            else
                KdomRED{itype} =   BdomRED{itype}'*(Cglo*BdomRED{itype}) ;
                %           KdomRED{itype} =  BasisUdef{itype}'*(K*BasisUdef{itype}) ;
            end
        end
        
    end
    FORCE_PROJECTS{iproject} = Fb(IDXdofs) + Ftrac(IDXdofs) ;   % External forces
end