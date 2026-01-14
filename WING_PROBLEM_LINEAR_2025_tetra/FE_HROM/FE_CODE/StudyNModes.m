clc
clear all

% Compute elasticity matrix
% -------------------------

RECALCULATE  = 0;
NameFileMesh = 'MYFIRSTMESH.msh'; %'mesh40k.msh' ;

% Transverse isotropic material
MATERIAL.FIBER.E(1) = 230e3 ; %MPa  % Longitudinal elastic modulus
MATERIAL.FIBER.E(2) = 8e3 ; %MPa   % Transverse elastic modulus
MATERIAL.FIBER.nu(1) = 0.25 ;  % Poisson's ratio (contraction in the transverse directionupon an extension in the fiber direction)
MATERIAL.FIBER.nu(2) = 0.3 ;  % Transverse poisson's ratio
MATERIAL.FIBER.Gshear(1) = 27.3e3 ; %MPa % In-plane shear modulus  G_LT
MATERIAL.FIBER.Gshear(2)= 3.08e3 ; % MPa % Transverse shear modulus G_TT
MATERIAL.FIBER.INDEX = 1;  % Material index
%%%
E = 3e3 ; nu = 0.3 ;
MATERIAL.MATRIX.E(1) = E  ; % MPa Elastic modulus   % MATRIX
MATERIAL.MATRIX.nu(1) = nu ; % Poisson's ratio
MATERIAL.MATRIX.Gshear(1) = E/2/(1+nu)  ;  % In-plane shear modulus
MATERIAL.MATRIX.INDEX = 2; % Material index


DATA.RECALCULATE_STIFFNESS =1 ;  % To avoid computing again the stiffness matrix (if =0) when computing
% the average stresses for different input
% strains, yet similar
% material/geometric properties
% ---- END INPUTS ----------------------------------

SNAPSHOTS_STRESS = [] ;
SNAPSHOTS_STRAIN = [] ;
SNAPSHOTS_DISP= [] ;
wsSNAPSHOTS = ['DATAWS/SNAPS_',NameFileMesh,'.mat'] ;



if  RECALCULATE  ==1
    for i = 1:6
        for j=1:6
            strainINP = zeros(6,1) ;
            strainINP(i) = 1;
            strainINP(j) = 1;
            [~, DATAOUT]= CompHomog_CF(strainINP,NameFileMesh,MATERIAL,DATA) ;
            DATA.RECALCULATE_STIFFNESS =0 ;
            SNAPSHOTS_STRESS =  [SNAPSHOTS_STRESS DATAOUT.stress] ;
            SNAPSHOTS_STRAIN = [SNAPSHOTS_STRAIN DATAOUT.strain] ;
            SNAPSHOTS_DISP = [SNAPSHOTS_DISP DATAOUT.disp] ;
        end
        
    end
    
    BOUNDATA = DATAOUT.BOUNDDATA ;
    
    save(wsSNAPSHOTS,'SNAPSHOTS_STRESS','SNAPSHOTS_STRAIN','SNAPSHOTS_DISP','BOUNDATA') ;
else
    % -------------------
    load(wsSNAPSHOTS)
    % -------------------
end

%%%%

rSTRESS = rank(SNAPSHOTS_STRESS);
rSTRAIN = rank(SNAPSHOTS_STRAIN);

% SVD
[Us,Ss,Vs] = svd(SNAPSHOTS_STRESS,0) ;
[Ue,Se,Ve] = svd(SNAPSHOTS_STRAIN,0) ;

%%
Ss = diag(Ss) ;
figure(1)
hold on
plot(log10(Ss/Ss(1)))
xlabel('N of   stress modes')
ylabel('log10(Singular values)')

Se = diag(Se) ;
figure(2)
hold on
plot(log10(Se/Se(1)))
xlabel('N of   strains modes')
ylabel('log10(Singular values)')

%%% COMPUTING ENERGY
ENERGY_P = 0.5*SNAPSHOTS_STRAIN.*SNAPSHOTS_STRESS ;
ENERGY = zeros(size(SNAPSHOTS_STRAIN,1)/6,size(SNAPSHOTS_STRESS,2)) ;
for i= 1:6
    indS = i:6:size(SNAPSHOTS_STRAIN,1) ;
    ENERGY = ENERGY + ENERGY_P(indS,:) ;
end

[U,S,V] = svd(ENERGY,0) ;
S = diag(log10(S/S(1))) ;
figure(3)
hold on
plot(S)
xlabel('N of energy modes')
ylabel('log10(Singular values)')


%% And what about internal forces ?
% -------------------------------

DOFm = BOUNDATA.DOFm   ;
DOFr = BOUNDATA.DOFr    ;
DOFf = BOUNDATA.DOFf    ;
Gb   = BOUNDATA.Gb      ;

DOFl = [DOFf; DOFm] ;
NAMEWS = ['DATAWS/',NameFileMesh,'_WS.mat'] ;
load(NAMEWS,'Bst');

Bf = [Bst(:,DOFf) Bst(:,DOFm)+Bst(:,DOFr)*Gb] ;

%% Displacement modes
displac = SNAPSHOTS_DISP(DOFl,:) ;
[U,S,V] = svd(displac,0) ;
S = diag(log10(S/S(1))) ;
BasisU = U(:,1:6) ;
BfRI = Bf*BasisU ;
%%%%%

% Internal forces (reduced)
ngaus = size(BfRI,1)/6 ;
nmodesU = 6;
INT_FORCES = zeros(ngaus*nmodesU,size(SNAPSHOTS_DISP,2)) ;
nstrain = 6 ;
for igaus=1:ngaus
    iini = (igaus-1)*nmodesU+1;
    ifin =  igaus*nmodesU ;
    ig = (igaus-1)*nstrain+1 ;
    fg = (igaus*nstrain);
    
    
    INT_FORCES(iini:ifin,:) = BfRI(ig:fg,:)'*SNAPSHOTS_STRESS(ig:fg,:) ;
    
    
    aaaa =mod(igaus,100) ;
    if aaaa==0
        %         if         igaus == 130000;
        %             dbstop('46')
        %         end
        disp(['igaus=',num2str(igaus),' of ',num2str(ngaus)])
    end
end
[U,S,V] = svd(INT_FORCES,0) ;
S = diag(log10(S/S(1))) ;
figure(4)
hold on
plot(S)
xlabel('N of internal force modes')
ylabel('log10(Singular values)')
 
%
BasisFex = U(:,1:6) ; 
nmodesFINTex = 6 ; 
Jmin = zeros(nmodesU*nmodesFINTex,ngaus);

for imode=1:nmodesU
    imodeG = imode:nmodesU:nmodesU*ngaus ;
    JminLOC = BasisFex(imodeG,:)' ;
    
    igausG = imode:nmodesU:nmodesU*nmodesFINTex ;
    Jmin(igausG,:) = JminLOC ;
    
end
 
S= svd(Jmin,0) ;
%S = diag(log(S/S(1))) ;
figure(6)
hold on
plot(log10(S/S(1)))
xlabel('N of Modes of Jmin')
ylabel('log10(Singular values)')
 

%%%%%