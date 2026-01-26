clc
clear all

% Compute Failure  Envelope Laminate
% ------------------------------------

RECALCULATE_FAILURE = 1;
NameFileMesh ='lam10lay_8000.msh' ; 'layer10_16rve.msh'; 'layer10_4rve_32000.msh';'layer10_16rve.msh' ;'lam10lay_128.msh' ;'lam10lay_32000.msh' ;
NameWS_unitcell ='DATAWS/Celas_ejem1.msh.mat' ; 'DATAWS/Celas_MYFIRSTMESH.msh.mat' ;
NameWS_failureUnitCell = 'DATAWS/FAILURE_ejem1.msh.mat' ;
% Elasticity matrix for each  layer (and angle of rotation)
% The name of the .mat containg such information is to be provided
% -----------------------
%
nplies = 8 ;
ANG_FIB = [0 0 0 0 0 0 0 0]
for iply = 1:nplies
    MATERIAL.PLY(iply).NAMEWS = NameWS_unitcell;  % Store the Celas in this folder
    MATERIAL.PLY(iply).ANGLE = ANG_FIB(iply) ;  % Angle subtended by the fiber and the x-axis
    MATERIAL.PLY(iply).NAMEWS_FAIL = NameWS_failureUnitCell ;
end
DATA.RECALCULATE_STIFFNESS = 1;   % To avoid computing again the stiffness matrix (if =0) when computing
% the average stresses for different input
% strains, yet similar
% material/geometric properties
DATA.CALCULATE_averageSTRESS = 2; %
 DATA.NUMBER_ELEMENTS_FAIL = 20;
DATA.BOUNDARY_CONDlam =   'ZERO';'ZERO_minTB' ;  %PERIODIC';     % PERIODIC/ZERO
DATA.niterCONJG = 3000 ; % Number of iterations  conjugated gradient
% ---- END INPUTS ----------------------------------
nameFailure= ['DATAWS/FailureLam_',NameFileMesh,'.mat'] ;
%% STRAIN TRAJECTORIES
Curvatures = [1 0 0 1 1 0 1 0.5
    0 1 0 1 0 1  1 -0.7
    0 0 1 0 1 1 1  0.3];
compp = [4 5 6] ;
STRAIN_TRAJECTORIES = zeros(8,size(Curvatures,2)) ;
STRAIN_TRAJECTORIES(compp,:) = Curvatures ;
STRAIN_TRAJECTORIES = [STRAIN_TRAJECTORIES -STRAIN_TRAJECTORIES] ;

if  RECALCULATE_FAILURE ==1
    ntraj = size(STRAIN_TRAJECTORIES,2) ;
    stressFAILURE = zeros(size(STRAIN_TRAJECTORIES))  ;
    for i = 1:8
        strainINP = zeros(8,1) ;
        strainINP(i) = 1;
        [stressMACRO DATAOUT]= LaminateSTRESScal(strainINP,NameFileMesh,MATERIAL,DATA) ;
        
        stressFAILURE(:,i) =DATAOUT.stressAVG_FAIL ;
        DATA.RECALCULATE_STIFFNESS =0 ;
    end
    %------
    
    save(nameFailure,'stressFAILURE') ;
else
    % -------------------
    load(nameFailure)
    % -------------------
end
% ----------------------------------
% Failure envelope for stresses
% --------------------------------
Mx = stressFAILURE(4,:)' ;
My = stressFAILURE(5,:)' ;
Mxy = stressFAILURE(6,:)' ;
LEGENDS  = {'M_x','M_y','M_{xy}'} ;
numfig = 4;
[Fenv radii_stress evecs_stress chi2_stress] = PlotEnvelopeFailure(Mx,My,Mxy,LEGENDS,numfig) ;

save(nameFailure,'stressFAILURE','Fenv','radii_stress','-append') ;
%

% % ----------------------------------
% % Failure envelope for strains
% % --------------------------------
% First, we need to obtain the Celas matrix. How ?

%
%
%
% e1 = strainFAILURE(1,:)' ;
% e2 = strainFAILURE(2,:)' ;
% e12 = strainFAILURE(6,:)' ;
% LEGENDS  = {'\epsilon_1','\epsilon_2','\gamma_{12}'} ;
% numfig = 5;
% [Fenv radii evecs chi2] = PlotEnvelopeFailure(e1,e2,e12,LEGENDS,numfig) ;
%
%
% load(nameCelas,'Celas')
% % We define a larger Fenv
% FenvGlo = zeros(6,6) ;
% comp = [1,2,6] ;
% FenvGlo(comp,comp) = Fenv ;
% % Compliance matrix
% Selas = inv(Celas) ;
% % Thus
% FenvGlo_e = Selas'*FenvGlo*Selas ;
% FenvSTRAIN = FenvGlo_e(comp,comp);
%
% %[ center, radii, evecs, v, chi2 ] = ellipsoid_fit( [x,y,z],'') ;
% center = [0,0,0] ;
% [evecs radii ]= eig(FenvSTRAIN);
% radii = sqrt(diag(radii))
%
% x = strainFAILURE(1,:)' ;
% y = strainFAILURE(2,:)' ;
% z = strainFAILURE(6,:)' ;
% figure(7)
% hold on
% xlabel(LEGENDS{1})
% ylabel(LEGENDS{2})
% zlabel(LEGENDS{3})
% plot3( x, y, z, '.r' );
%       nn = 50 ;
%     [xx yy zz] = ellipsoid(center(1),center(2),center(3),radii(1),radii(2),radii(3),nn) ;
%
%     % ROTATION
%     X = xx(:);
%     Y = yy(:);
%     Z = zz(:) ;
%     XMAT = evecs*[X';Y';Z'];
%     X = XMAT(1,:) ;
%     Y = XMAT(2,:) ;
%     Z = XMAT(3,:) ;
%     xx = reshape(X,nn+1,[]  );
%     yy = reshape(Y,nn+1,[]  );
%     zz = reshape(Z,nn+1,[]  );
%     %
%
%     surfl(xx,yy,zz)
%     shading interp
%     colormap(gray);
%
