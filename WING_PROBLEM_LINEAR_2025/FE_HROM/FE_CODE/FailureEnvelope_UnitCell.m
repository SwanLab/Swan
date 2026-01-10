clc
clear all

% Compute elasticity matrix
% -------------------------

RECALCULATE_FAILURE = 1;
NameFileMesh = 'ejem1.msh'; 'iSO,msh'; 'Lab1_06.msh' ; ; %'mesh40k.msh' ;
 
% Transverse isotropic material
MATERIAL.FIBER.E(1) = 230e3 ; %MPa  % Longitudinal elastic modulus
MATERIAL.FIBER.E(2) = 8e3 ; %MPa   % Transverse elastic modulus
MATERIAL.FIBER.nu(1) = 0.25 ;  % Major poisson's ratio (contraction in the transverse directionupon an extension in the fiber direction)
MATERIAL.FIBER.nu(2) = 0.3 ;  % Transverse poisson's ratio
MATERIAL.FIBER.Gshear(1) = 27.3e3 ; %MPa % In-plane shear modulus
MATERIAL.FIBER.Gshear(2)= 3.08e3 ; % MPa % Transverse shear modulus
MATERIAL.FIBER.INDEX = 2;  % Material index (within GID mesh)
MATERIAL.FIBER.strength  = 3930;  %  % Strength in MPa (same in all directions)
MATERIAL.FIBER.CONSTITUTIVE_MODEL  = 'VONMISES';  % s
%%%
E = 3e3 ; nu = 0.3 ;
MATERIAL.MATRIX.E(1) = E  ; % MPa Elastic modulus   % MATRIX
MATERIAL.MATRIX.nu(1) = nu ; % Poisson's ratio
MATERIAL.MATRIX.Gshear(1) = E/2/(1+nu)  ;  % In-plane shear modulus
MATERIAL.MATRIX.INDEX = 1; % Material index (within GID mesh)
MATERIAL.MATRIX.strength  = 64;  %   % Strength in MPa (same in all directions)
MATERIAL.MATRIX.CONSTITUTIVE_MODEL  = 'VONMISES';  %



% % Transverse isotropic material
% E= 230e3 ; nu = 0.25 ;
% MATERIAL.FIBER.E(1) = E ; %MPa  % Longitudinal elastic modulus
% MATERIAL.FIBER.E(2) =E ; %MPa   % Transverse elastic modulus
% MATERIAL.FIBER.nu(1) = nu ;  % Major poisson's ratio (contraction in the transverse directionupon an extension in the fiber direction)
% MATERIAL.FIBER.nu(2) = nu ;  % Transverse poisson's ratio
% MATERIAL.FIBER.Gshear(1) = E/2/(1+nu); %MPa % In-plane shear modulus
% MATERIAL.FIBER.Gshear(2)= E/2/(1+nu)  ; % MPa % Transverse shear modulus
% MATERIAL.FIBER.INDEX = 2;  % Material index (within GID mesh)
% MATERIAL.FIBER.strength  = 3930;  %  % Strength in MPa (same in all directions)
% MATERIAL.FIBER.CONSTITUTIVE_MODEL  = 'VONMISES';  % s
% %%%
% E = MATERIAL.FIBER.E(1) ; nu =MATERIAL.FIBER.nu(1) ;
% MATERIAL.MATRIX.E(1) = E  ; % MPa Elastic modulus   % MATRIX
% MATERIAL.MATRIX.nu(1) = nu ; % Poisson's ratio
% MATERIAL.MATRIX.Gshear(1) = MATERIAL.FIBER.Gshear(1)   ;  % In-plane shear modulus
% MATERIAL.MATRIX.INDEX = 1; % Material index (within GID mesh)
% MATERIAL.MATRIX.strength  =MATERIAL.FIBER.strength;  %   % Strength in MPa (same in all directions)
% MATERIAL.MATRIX.CONSTITUTIVE_MODEL  = 'VONMISES';  %


%
DATA.CALCULATE_STRENGTH = 1;    % Calculate strength
DATA.NUMBER_ELEMENTS_FAIL = 100; %  At least,  DATA.NUMBER_ELEMENTS_FAIL of all elements should
%fail in order to consider failure of the RVE
DATA.RECALCULATE_STIFFNESS =1 ;  % To avoid computing again the stiffness matrix (if =0) when computing
% the average stresses for different input
% strains, yet similar
% material/geometric properties
%% STRAIN TRAJECTORIES
strainX_Y_XY = [1 0 0 1 1 0 1 0.5
                0 1 0 1 0 1  1 -0.7
                0 0 1 0 1 1 1  0.3];
compp = [1 2 6] ;
STRAIN_TRAJECTORIES = zeros(6,size(strainX_Y_XY,2)) ;
STRAIN_TRAJECTORIES(compp,:) = strainX_Y_XY ;
STRAIN_TRAJECTORIES = [STRAIN_TRAJECTORIES -STRAIN_TRAJECTORIES] ;
% ---- END INPUTS ----------------------------------
nameFAILURE = ['DATAWS/FAILURE_',NameFileMesh,'.mat'] ;
nameCelas = ['DATAWS/Celas_',NameFileMesh,'.mat'] ;
if  RECALCULATE_FAILURE ==1
    %     figure(1)
    %     hold on
    %     title('Failure envelope in stress space')
    %     xlabel('\sigma_1')
    %     ylabel('\sigma_2')
    ntraj = size(STRAIN_TRAJECTORIES,2) ;
    strainFAILURE = zeros(size(STRAIN_TRAJECTORIES)) ;
    stressFAILURE = zeros(size(STRAIN_TRAJECTORIES))  ;
    for i = 1:ntraj
        strainINP  = STRAIN_TRAJECTORIES(:,i);
        [stressMACRO DATAOUT]   = CompHomog_CF(strainINP,NameFileMesh,MATERIAL,DATA) ;
        strainFAILURE(:,i) =DATAOUT.strainAVG_FAIL ;
        stressFAILURE(:,i) =DATAOUT.stressAVG_FAIL ;
        DATA.RECALCULATE_STIFFNESS =0 ;
        %            plot(stressFAILURE(1,i),stressFAILURE(2,i),'r*')
        %         text(stressFAILURE(1,i),stressFAILURE(2,i),num2str(i))
    end
    save(nameFAILURE,'strainFAILURE','stressFAILURE') ;
else
    % -------------------
    load(nameFAILURE)
    % -------------------
    
    
end
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% hold on
% title('Failure envelope in stress space  ')
% xlabel('\sigma_1')
% ylabel('\sigma_2')
%  for i=1:size(stressFAILURE,2)
%      if abs(stressFAILURE(6,i))<1e-2
%          plot(stressFAILURE(1,i),stressFAILURE(2,i),'r*')
%          text(stressFAILURE(1,i),stressFAILURE(2,i),num2str(i))
%      end
%  end
%
% figure(2)
% hold on
% title('Failure envelope in stress space (plane sigma1-sigma2)')
% xlabel('\sigma_1')
% ylabel('\sigma_2')
% for i=1:size(stressFAILURE,2)
%     plot(stressFAILURE(1,:),stressFAILURE(6,:),'r*')
%     text(s;
 %   stressFAILURE = zeros(size(STRAIN_TRAJECTORIES))  stressFAILURE(1,:),stressFAILURE(6,:),num2str(i))
% end

% ----------------------------------
% Failure envelope for stresses
% --------------------------------
sigma1 = stressFAILURE(1,:)' ;
sigma2 = stressFAILURE(2,:)' ;
tau12 = stressFAILURE(6,:)' ;
LEGENDS  = {'\sigma_1','\sigma_2','\tau_{12}'} ;
numfig = 4;
[Fenv radii_stress evecs_stress chi2_stress] = PlotEnvelopeFailure(sigma1,sigma2,tau12,LEGENDS,numfig) ;

save(nameFAILURE,'strainFAILURE','stressFAILURE','Fenv','radii_stress','-append') ;
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
