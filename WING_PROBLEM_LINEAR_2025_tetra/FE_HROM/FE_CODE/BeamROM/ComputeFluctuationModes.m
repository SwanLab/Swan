function [Ufluct,ksfluc] = ComputeFluctuationModes(BasisUdef,f1,nBASES_BEAM,DATAIN,Vrb,...
    M,NODES_faces12,DATA_REFMESH,DATAsnap)

if nargin == 0
    load('tmp2.mat')
end
ndim = 3;
DATAIN = DefaultField(DATAIN,'Type_Fluctuation_Modes','Separated_by_tests');
%DATA.Type_Fluctuation_Modes = 'Separated_by_tests' ;  % New variable to separately treat the fluc. mode of each beam test


ksfluc = [] ; 
if ~isempty(DATAsnap)
    DATAIN.Type_Fluctuation_Modes = 'ImposedDisplacementTests' ;
end


DATAIN = DefaultField(DATAIN,'ONLY_DIRICHLET_BCS_JOINTS',1) ;

if DATAIN.ONLY_DIRICHLET_BCS_JOINTS == 1
    DATAIN.Type_Fluctuation_Modes = 'ONLY_DIRICHLET_BCS_JOINTS' ;
end


switch DATAIN.Type_Fluctuation_Modes
    case 'ImposedDisplacementTests'
        if iscell(M)
            M = M{1} ;
            Vrb = Vrb{1} ;
        end
        % Reaction forces face A
        nodesA = DATA_REFMESH.NODES_faces12{1} ;
        DOFA = small2large(nodesA,ndim)  ;
        nodesB = DATA_REFMESH.NODES_faces12{2} ;
        DOFB = small2large(nodesB,ndim)  ;
        Dflex = zeros(6,6) ;
        for itest = 1:length(DATAsnap.Reactions)
            FORCEA = Vrb'*DATAsnap.Reactions{itest}(DOFA)  ;
            FORCEB = Vrb'*DATAsnap.Reactions{itest}(DOFB)  ;
            Dflex(:,itest) = FORCEA' ;
        end
        
        % Coefficients
        invDflex = inv(Dflex) ;
        TEST_names = {'axial','shear_y','shear_z','torsion','pbend_y','pbend_z'};
        coeff_TEST = invDflex ; %
        %%%% COMPUTING FLUCTUATIONS
        
        %          TEST_LOCAl = 'shear_z' ;
        %                     INF = 3 ;
        %                     InputForce = FB(INF) ;
        %                     % ---- REmoved the part corresponding to bending
        %                     Ub = Ufluct.('pbend_y') ;
        %                     dBfluc = dBfluc - Ub*(Ub'*M*dBfluc) ;
        %                     % Norm of the vector
        %                     normD = sqrt(dBfluc'*(M*dBfluc)) ;
        %                     % Fluctuation mode (---> shear z)
        %                     Uloc = dBfluc/normD;   % M-Norm = 1
        %                     Ufluct.(TEST_LOCAl) =  Uloc ; % Matrix of fluctuation modes
        
        Ufluct = [] ;
        SNAP_dA = cell2mat(DATAsnap.Displacements) ;
        SNAP_dA = SNAP_dA(DOFA,:) ; % Snapshots of displacements, face A
        VrbBAR = M*Vrb ;
        PSD = inv(VrbBAR'*Vrb);
        ksfluc = [] ;
        VALUE_FORCE = 1;
        for idim = 1:length(TEST_names)
            dFLUC = SNAP_dA*coeff_TEST(:,idim) ; % Total displacements
            coeffs = PSD*(VrbBAR'*dFLUC) ;
            dFLUC = dFLUC - Vrb*coeffs ;
            normD = sqrt(dFLUC'*(M*dFLUC)) ;
            Uloc = dFLUC/normD;   % M-Norm = 1
            TEST_LOCAL = TEST_names{idim} ;
            Ufluct.(TEST_LOCAL) =  Uloc ; % Matrix of fluctuation modes
            ksfluc.(TEST_LOCAL) = VALUE_FORCE/normD ;
        end
        
        
        
        
        
        
    case 'Separated_by_tests'
        % -------------------------------
        load(DATAIN.NameWs_displacementB_beams,'DATA_DISPLACEMENT_faceB') ;
        % Computing the displacements associated to each beam test
        TESTS = fieldnames(DATA_DISPLACEMENT_faceB) ;
        Ufluct = [] ;
        ksfluc = [] ;
        
        if iscell(M)
            M = M{1} ;
            Vrb = Vrb{1} ;
        end
        
        for itest = 1:length(TESTS)
            TEST_LOCAl = TESTS{itest} ;
            dB = DATA_DISPLACEMENT_faceB.(TEST_LOCAl).DISP;  % Displacement face B
            FB = DATA_DISPLACEMENT_faceB.(TEST_LOCAl).INPUT_FORCE;
            Rbar = M*Vrb ;
            PRB = (Rbar'*Vrb) ;
            coeff = PRB\(Rbar'*dB) ; % Coefficients projection onto span Vrb (M-norm)
            dBfluc = dB - Vrb*coeff ; % Fluctuation mode
            
            switch TEST_LOCAl
                case {'sbend_y'}
                    TEST_LOCAl = 'shear_z' ;
                    INF = 3 ;
                    InputForce = FB(INF) ;
                    % ---- REmoved the part corresponding to bending
                    Ub = Ufluct.('pbend_y') ;
                    dBfluc = dBfluc - Ub*(Ub'*M*dBfluc) ;
                    % Norm of the vector
                    normD = sqrt(dBfluc'*(M*dBfluc)) ;
                    % Fluctuation mode (---> shear z)
                    Uloc = dBfluc/normD;   % M-Norm = 1
                    Ufluct.(TEST_LOCAl) =  Uloc ; % Matrix of fluctuation modes
                case {'sbend_z'}
                    TEST_LOCAl = 'shear_y' ;
                    INF = 2 ;
                    InputForce = FB(INF) ;
                    % ---- REmoved the part corresponding to bending
                    Ub = Ufluct.('pbend_z') ;
                    dBfluc = dBfluc - Ub*(Ub'*M*dBfluc) ;
                    % Norm of the vector
                    normD = sqrt(dBfluc'*(M*dBfluc)) ;
                    % Fluctuation mode (---> shear z)
                    Uloc = dBfluc/normD;   % M-Norm = 1
                    Ufluct.(TEST_LOCAl) =  Uloc ; % Matrix of fluctuation modes
                    
                    
                otherwise
                    normD = sqrt(dBfluc'*(M*dBfluc)) ;
                    Uloc = dBfluc/normD;   % M-Norm = 1
                    Ufluct.(TEST_LOCAl) =  Uloc ; % Matrix of fluctuation modes
                    % Determination stiffness constant associated to
                    % fluctuations
                    INF = find(FB~=0) ;
                    InputForce = FB(INF) ;
                    
            end
            ksfluc.(TEST_LOCAl)  = InputForce/normD ;
            
        end
        %%T%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
        
        
        
        % --------------------------------
    otherwise
         if iscell(M)
            M = M{1} ;
            Vrb = Vrb{1} ;
        end
        BasisU_f1 = BasisUdef(f1,1:nBASES_BEAM.DISPLACEMENTS) ;
        % Projection onto the orthogonal complement of Vrb
        coeff= (Vrb'*M*Vrb)\(Vrb'*M*BasisU_f1) ;
        U = BasisU_f1 - Vrb*coeff ;   % Fluctuation mode
       % TOL = 1e-6 ;
     %   [Ufluct,SINGV ]= SVDT(U,TOL) ;  % Orthogonalization
        %
        % Make them M-orthogonal 
        Mchol = chol(M) ;
        Xbar = Mchol*U ;
        TOL = 0 ;
        [Ubar,S,Vbar] = SVDT(Xbar,TOL) ;
        Ufluct = Mchol\Ubar ;
        
        
        DATAIN = DefaultField(DATAIN,'PostProcessFluctuationModes',1) ;
        
        %%%
        % For prescribing BCs
        %         CHECK = 0;
        %         if CHECK == 1
        %             warning('Checking BCs')
        %             coeff= (Ufluct'*M*Ufluct)\(Ufluct'*M) ;
        %             Q = Ufluct*coeff ;
        %             IQ = eye(size(Q))-Q ;
        %             R = Vrb ;
        %
        %             [aaa l] = licols(IQ') ;
        %             r = 1:size(Q,1) ;
        %             r(l) = [] ;
        %             IQ_ll = eye(length(l)) -Q(l,l) ;
        %
        %
        %
        %         end
        
        %%%
        
end
% -----------------------------------------
DATAIN = DefaultField(DATAIN,'PostProcessFluctuationModes',1) ;


if DATAIN.PostProcessFluctuationModes == 1
    
    if  isstruct(Ufluct)
        UfluctINP = [] ;
        TESTS = fieldnames(Ufluct) ;
        for itest = 1:length(TESTS)
            UfluctINP = [UfluctINP,Ufluct.(TESTS{itest})] ;
        end
    else
        UfluctINP = Ufluct ;
    end
    
    refFACE = 1;
    COOR =DATA_REFMESH.COOR(NODES_faces12{refFACE},:) ;
    CNref =  RenumberConnectivities( DATA_REFMESH.CONNECTb{refFACE},1:length(NODES_faces12{refFACE})) ;
    TypeElementB = DATA_REFMESH.TypeElementB ;
    posgp = [] ;
    LEGENDG= ['FLUCT.'] ;
    NAME_MODES = [DATAIN.NAME_WS_MODES(1:end-4),LEGENDG ];
    
    DATAMODES = [] ;
    
    
    GidPostProcessModes_dom(COOR,CNref,TypeElementB,UfluctINP,posgp,NAME_MODES,DATAMODES,LEGENDG);
    
end