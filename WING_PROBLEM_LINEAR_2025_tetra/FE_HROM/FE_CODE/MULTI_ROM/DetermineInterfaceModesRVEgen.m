function [DATAOUT,BASES,V,TEXTP,BasisUdef,BasisRdef,DATAIN] = DetermineInterfaceModesRVEgen(DATAIN,SROTrows,DATA_REFMESH,DATAOUT,BasisUdef,SinvVal_Udef,...
    f,fI,SinvVal_Rdef,FACES_GROUPS,BASES,BasisRdef,TEXTP)



if DATAIN.INTERFACE_MODES_REACTIONS_WORK.ACTIVE ==1
    %   METHOD IMPLEMENTED MAY-16-2019
    % --- EQUIVALENT IN BEAMS:
    %  ELASTOSTATIC_GEN/FE_CODE/BeamROM/ReactANDintf_WORKREACT.m
    if any(SROTrows)
        error('Function not adapted to curved surfaces')
    end
    error('Warning: PArt of the code abandoned Sept-2019 in favor of DATAIN.DeformationalInterfaceModes_AlignmentMethod')
    [DATAOUT,BASES,V,TEXTP,BasisUdef] =  ReactANDintf_WORKREACT_RVES(DATA_REFMESH,DATAOUT,DATAIN,BasisUdef,BasisRdef,...
        SinvVal_Udef,f,fI,SinvVal_Rdef,FACES_GROUPS,BASES) ;
elseif DATAIN.DeformationalInterfaceModes_AlignmentMethod ==1
    % METHOD; 21-May-2019
    % Extension DATAIN.DeformationalInterfaceModes_AlignmentMethod to pair
    % of faces
    
    [DATAOUT,BASES,V,TEXTP,BasisUdef] =  ReactANDintf_ALIGN_TWOFACES(DATA_REFMESH,DATAOUT,DATAIN,BasisUdef,BasisRdef,...
        SinvVal_Udef,f,fI,SinvVal_Rdef,FACES_GROUPS,BASES,TEXTP) ;
    
    
else
    
        error('Warning: PArt of the code abandoned Sept-2019 in favor of DATAIN.DeformationalInterfaceModes_AlignmentMethod')

    if any(SROTrows)
        error('Function not adapted to curved surfaces')
    end
    
    % Rigid body modes for interface
    Vrb = DATA_REFMESH.RigidBodyModesInterface ;
    DATAOUT.BasisIntRB = Vrb ;
    DATAIN = DefaultField(DATAIN,'RIGID_BODY_MODES_INTERFACE',1) ;
    %%% Displacements interface boundaries
    %DISP_BOUND = BASES.DISPLACEMENTS_INTERFACES_BOUNDARY ;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Selection of reaction modes to meet stability conditions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DATAIN = DefaultField(DATAIN,'TOL_SINGULAR_VALUES_Hqr',0.1) ;
    %
    DATAIN = DefaultField(DATAIN,'EQUAL_NUMBER_REACTIONS_NUMBER_MODES_REACTION_DRIVEN',0) ;
    TEXTP = {} ;
    if  DATAIN.EQUAL_NUMBER_REACTIONS_NUMBER_MODES_REACTION_DRIVEN == 0
        BasisRdef = SelectReactionModesStability(BasisUdef,BasisRdef,f,DATAIN.TOL_SINGULAR_VALUES_Hqr) ;
        
    else
        [BasisUdef,SinvVal_Udef, TEXTP]= ReactionDispModesEqual(BasisUdef,BasisRdef,f,...
            SinvVal_Udef) ;
        DATAOUT.BasisUdef = BasisUdef ;
        
        BASES.DISPLACEMENTS.U = BasisUdef ;
        BASES.DISPLACEMENTS.S = SinvVal_Udef ;
        
        % PLOT AGAING DISPL. MODES
        
        %if DATAIN.PLOT_MODES == 1
        %  ELEMS =  DATAOUT.DOMAINVAR.ListElements{idom} ;
        
        %     INCLUDE_MASS_MATRIX = 0;
        %     FACTOR_MASS_MATRIX = 1e3 ;
        %     if    INCLUDE_MASS_MATRIX == 1
        %
        %         [BasisUdefPlot,SSS,VVV] =  WSVDT(BasisUdef,DATA_REFMESH.M/FACTOR_MASS_MATRIX) ;
        %
        %     else
        BasisUdefPlot = BasisUdef ;
        
        % end
        
        CNref = DATA_REFMESH.CN  ;
        COOR = DATA_REFMESH.COOR  ;
        NAME_MODES = [DATAIN.NAME_WS_MODES(1:end-4),'DISP_ALIGNED' ];
        DATA= [] ;
        LEGENDG = 'DISP.'
        GidPostProcessModes_dom(COOR,CNref,DATA_REFMESH.TypeElement,BasisUdefPlot,DATA_REFMESH.posgp,...
            NAME_MODES,DATA,LEGENDG);
        
        
        %end
        
        
    end
    
    
    
    DATAOUT.BasisRdef = BasisRdef ;
    
    %%%%%
    
    
    
    % ------------------------
    % Geometric mass matrices
    % ------------------------
    nfaces = length(Vrb);
    M1d = DATA_REFMESH.GeometricMassMatrixInterface ;
    M = cell(nfaces,1) ;
    for iface= 1:nfaces
        M{iface} = sparse(size(M1d{iface},1)*ndim,size(M1d{iface},2)*ndim) ;
        for idim = 1:ndim
            M{iface}(idim:ndim:end,idim:ndim:end) = M1d{iface} ;
        end
    end
    
    DATAOUT.MassMinterfaces = M ;
    
    
    
    
    % DETERMINATION OF INTERFACE MODES.
    % ----------------------------------
    
    %APPROACHES BEFORE April-21st, 2019 (before kinematical constraints)
    if DATAIN.KINEMATIC_CONSTRAINTS_MODES.ACTIVE == 0
        [DATAOUT,V,TEXTP] =  InterfaceRVEModesStandard(DATAIN,Vrb,DATAOUT,nfaces,fI,BasisUdef,BasisRdef,M,...
            SinvVal_Udef,SinvVal_Rdef,TEXTP,FACES_GROUPS) ;
    else
        DATAIN.KINEMATIC_CONSTRAINTS_MODES = DefaultField(DATAIN.KINEMATIC_CONSTRAINTS_MODES,'COMBINED_ENERGETIC_METHOD',0) ;
        if DATAIN.KINEMATIC_CONSTRAINTS_MODES.COMBINED_ENERGETIC_METHOD == 0
            [DATAOUT,V,TEXTP] =  InterfaceRVEKinematicConstraint(DATAIN,Vrb,DATAOUT,nfaces,fI,BasisUdef,BasisRdef,M,...
                SinvVal_Udef,SinvVal_Rdef,TEXTP,FACES_GROUPS,BasisUrb,DATA_REFMESH.M ) ;
        else
            % Combined with the energetic method (InterfaceRVEModesStandard)
            DATAIN.DeformationalInterfaceModes_AlignmentMethod =1;
            [DATAOUT,V,TEXTP] =  InterfaceRVEenergetic(DATAIN,Vrb,DATAOUT,nfaces,fI,BasisUdef,BasisRdef,M,...
                SinvVal_Udef,SinvVal_Rdef,TEXTP,FACES_GROUPS,BasisUrb,DATA_REFMESH.M,SinvVal_Udef) ;
            
            
        end
        
    end
    
    
end