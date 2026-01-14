function [DISP,qDEF] = ReconstDisplacementsEIFE(dCOARSEloc,EIFE_prop,Vrot,DATA,SCALE_FACTOR_ALLS,ROTATIONmat,ndim)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/101_MULTI2D_2023/01_HOMOG/07_PostProcess_2D.mlx

 if nargin == 0
     load('tmp1.mat')
 end
 
 factorDEF = DATA.FACTOR_MULTIPLYING_DEFORMATIONAL_RB_amplitudes(1) ;
  factorRB = DATA.FACTOR_MULTIPLYING_DEFORMATIONAL_RB_amplitudes(2) ;
  
  if factorDEF*factorRB ~= 1
      disp(['Warning: deformational/rigid body displacements displayed with scaling factors different from one'])
  end

% deformational part (amplitude modes)
qDEF = factorDEF*EIFE_prop.OPER.HdefINV_PsiDEFfT*Vrot*dCOARSEloc ;
% rigid body part (amplitude modes)
qRB = factorRB*EIFE_prop.OPER.Ginv_PhiRBfT_m_HrbHdefINV_PsiDEFfT*Vrot*dCOARSEloc ;
% Reconstruction (parent domain)
dREF = EIFE_prop.RECONSTRUCTION.DEF_DISP.BASIS*qDEF + EIFE_prop.RECONSTRUCTION.RB_DISP.BASIS*qRB ;

if DATA.UNIFORM_SCALING_REFERENCE_ELEMENT == 1
    
    %    dREFrot = SCALE_FACTOR_ALLS*ROTATIONmat*reshape(dREF,ndim,[]) ;
    
    
    DISP = zeros(size(dREF)) ;
    
    for itime = 1:size(dREF,2)
        
        dREFrot = ROTATIONmat*reshape(dREF(:,itime),ndim,[]) ;
        DISP(:,itime) = dREFrot(:) ;
        
    end
    
    
    
else
    error('Option not implemented')
end