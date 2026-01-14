function [STRESSES_FINE,VMSTRESS] = ReconstStressesEIFE_viar(EIFE_prop,stressesREF,ROTATIONmatrix,DATA)

if nargin == 0
    load('tmp.mat')
end



coeff = EIFE_prop.RECONSTRUCTION.PK2STRESS.coeff*stressesREF;

STRESSES_FINE_REF = EIFE_prop.RECONSTRUCTION.PK2STRESS.BASIS*coeff ;
ndim = size(EIFE_prop.MESH.COOR,2) ;

VMSTRESS = VonMisesStressCOMP(STRESSES_FINE_REF,ndim,DATA) ;



if EIFE_prop.MESH.nstrain == 4
    
    STRESSES_FINE =zeros(size(STRESSES_FINE_REF)) ;
    
    for itime = 1:size(STRESSES_FINE_REF,2)
        STRESSES_FINE_loc =  reshape(STRESSES_FINE_REF(:,itime),EIFE_prop.MESH.nstrain,[])  ;
        
        ROTATION_STRESSES =    RotateStress2Dplanestrain(ROTATIONmatrix(1,1),ROTATIONmatrix(2,1)) ;
        STRESSES_FINE_loc(1:3,:) = ROTATION_STRESSES*STRESSES_FINE_loc(1:3,:) ;
        
        % GID CONVENTION FOR PRINTING
        %   STRESSES_FINE = STRESSES_FINE([1 2 4 3],:) ;
        STRESSES_FINE(:,itime) = STRESSES_FINE_loc(:);
        %  ngaus =size(EIFE_prop.MESH.posgp,2) ;
        %  nstrain = (EIFE_prop.MESH.nstrain) ;
        
        %  STRESSES_FINE_loc =  reshape(STRESSES_FINE_loc,nstrain*ngaus,[])  ; 
        
    end
    
    
else
    STRESSES_FINE = STRESSES_FINE_REF ; 
   % error('Option not implemented')
end