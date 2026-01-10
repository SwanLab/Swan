function  [BasisF, SingVal_F,MSG] = InternalFocesECM_fromStressModes(BASES,DATAIN,MSG,BdomRED,Wdom,DATAROM,Bdom)


    % Old, classical method

    BASES.STRESSES = DefaultField(BASES.STRESSES,'SNAPSHOTS',[]) ;
    if isempty(BASES.STRESSES.SNAPSHOTS)
        BasisS = BASES.STRESSES.U ;
        SingVal_stress =  BASES.STRESSES.S/BASES.STRESSES.S(1);  % 28-May-2019
        if length(SingVal_stress) ~= size(BasisS,2)
            SingVal_stress = ones(size(BasisS,2),1) ;
        end
    else
        BasisS = BASES.STRESSES.SNAPSHOTS ;  % Use directly STRESS SNAPSHOTS FOR CONSTRUCTING BasisF
        
        SingVal_stress =[] ; 
        DATAINloc.CUBATURE.IncludeSingularValuesStress = 0 ; 
    end





DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'AlignStressesWithStrains',0) ;

if DATAIN.CUBATURE.AlignStressesWithStrains == 1   && length(SingVal_stress) > size(BdomRED,2)
    error('Option not properly assessed...')
    MSG{end+1} = '---------------------------------------------------------' ; 
    MSG{end+1} = 'Alignment with strains when computing integration points' ;
    MSG{end+1} = ['Number of stress modes before alignment  = ',num2str(length(SingVal_stress))];
    [BasisS,SingVal_stress] = AlignmentSubspaces(BasisS,SingVal_stress,BdomRED);
else
        MSG{end+1} = ['Number of stress modes  = ',num2str(length(SingVal_stress))];
end

if DATAIN.CUBATURE.IncludeSingularValues_displacements ==1
%  error('OPtion not properly assessed ...')
   SingVal_Udef =  DATAROM.SingVal_Udef/DATAROM.SingVal_Udef(1) ; 
   BasisUdef = bsxfun(@times,DATAROM.BasisUdef',SingVal_Udef)' ; 
   BdomRED_singval = Bdom*BasisUdef ; 
else
    BdomRED_singval = BdomRED; 
end


% CgloDOM = CgloDOM*Bdom ;
% Basis for streses

DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'UseStrainsAsStresses',0) ; % Include Singular Value displacement 

if DATAIN.CUBATURE.UseStrainsAsStresses == 1   
    error('Option not properly assessed')
    BasisS = BdomRED ; 
     SingVal_stress = ones(size(BasisS,2),1) ;    
end


DATAINloc.PLOT_RIGHT_SINGULAR_VECTORS  =0  ;

DATAIN.CUBATURE = DefaultField(DATAIN.CUBATURE,'TOL_LOC_InternalForces',1e-6) ;
DATAINloc.TOL_LOC_InternalForces = DATAIN.CUBATURE.TOL_LOC_InternalForces ;
 
 
[BasisF, SingVal_F]= BasisFfromStress(BasisS,SingVal_stress,...
    BdomRED_singval, Wdom,DATAINloc) ;
