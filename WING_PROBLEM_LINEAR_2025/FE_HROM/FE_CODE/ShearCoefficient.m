function stressAVG = ShearCoefficient(DATA,stressAVG,stressGLOv_W,strainGLOv,VOLUME,dZ)




if  DATA.SPECIAL_BC_SHEAR_ON == 1 ;
    stressAVG(7:8) =  DATA.SPECIAL_BC_SHEAR_FACTOR*stressAVG(7:8) ;
end


%
% Recalculaing 7 and  8 components
% dbstop('98')

if DATA.strainINP(8)~= 0 & any(DATA.strainINP([1:7]))==0
    % Energy associated to components 4 and 5
    enerXZ = stressGLOv_W(5,:).*strainGLOv(5,:) ;
    enerXZ = sum(enerXZ);
    gamma_xz = DATA.strainINP(8) ;
    supRVE = VOLUME/dZ ;
    Qx = enerXZ/gamma_xz/supRVE ;
    
    stressAVG(8) = Qx ;
    
    %      shear_avg = sum(strainGLOv_W(:,5)) ;
    % %
    %      Qx = enerXZ/shear_avg/supRVE ;
    
end


if DATA.strainINP(7)~= 0 & any(DATA.strainINP([1:6,8]))==0
    % Energy associated to components 4 and 5
    enerYZ = stressGLOv_W(4,:).*strainGLOv(4,:) ;
    enerYZ = sum(enerYZ) ;
    gamma_yz = DATA.strainINP(7) ;
    supRVE = VOLUME/dZ ;
    Qy = enerYZ/gamma_yz/supRVE ;
    
    stressAVG(7) = Qy ;
    
    %       shear_avg = sum(strainGLOv_W(:,4)) ;
    %
    %     Qy = enerYZ/shear_avg/supRVE ;
    
end
%
%
%
%
