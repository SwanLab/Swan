function STRESSES_FINE = ReconstStressesEIFE(dE,EIFE_prop,stressesREF,ROTATIONstresses,Vrot,DATA,SCALE_FACTOR_ALLS)

if nargin == 0
    load('tmp.mat')
end




if DATA.RecoveryStressesFromElasticityMatrix == 1
    
    if  isfield(EIFE_prop.RECONSTRUCTION.PK2STRESS,'Celas_Bmat_PhiDEF')
        STRESSES_FINE_REF = (EIFE_prop.RECONSTRUCTION.PK2STRESS.Celas_Bmat_PhiDEF)*...
            EIFE_prop.OPER.HdefINV_PsiDEFfT*(Vrot*dE)/SCALE_FACTOR_ALLS;
        
    else
        error('DATA.RecoveryStressesFromElasticityMatrix  = 1  requires the field EIFE_prop.RECONSTRUCTION.PK2STRESS.Celas_Bmat_PhiDEF ')
    end
    
    %     IndexDomainLOC =  TRANSF_COORD.IndexParentDomain ;
    %     EIFEoper = PROPMAT(MaterialType(e)).EIFE_prop(IndexDomainLOC) ;
    
    
    
    %  weig = WEIGHTSinteg.INTforces{e} ;
    % ngaus = length(weig) ;
    nstrain = EIFE_prop.MESH.nstrain ;
    
    EIFE_prop.RECONSTRUCTION.STRAINS.coeff ;
    
    
else
    coeff = EIFE_prop.RECONSTRUCTION.PK2STRESS.coeff*stressesREF(:);
    
    STRESSES_FINE_REF = EIFE_prop.RECONSTRUCTION.PK2STRESS.BASIS*coeff ;
    
    
    % ROTATION
    
    
end


if EIFE_prop.MESH.nstrain == 4
    STRESSES_FINE =  reshape(STRESSES_FINE_REF,EIFE_prop.MESH.nstrain,[])  ;
    STRESSES_FINE(1:3,:) = ROTATIONstresses*STRESSES_FINE(1:3,:) ;
    
    % GID CONVENTION FOR PRINTING
    %   STRESSES_FINE = STRESSES_FINE([1 2 4 3],:) ;
    STRESSES_FINE = STRESSES_FINE(:);
    ngaus =size(EIFE_prop.MESH.posgp,2) ;
    nstrain = (EIFE_prop.MESH.nstrain) ;
    
    STRESSES_FINE =  reshape(STRESSES_FINE,nstrain*ngaus,[])  ;
    
    
else
    error('Option not implemented')
end