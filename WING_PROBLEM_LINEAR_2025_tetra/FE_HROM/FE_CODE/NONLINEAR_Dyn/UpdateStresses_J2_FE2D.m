function [Ctang,  stressNP1, EPnp1,  alphaNP1,  sigmayNP1,CtangM ] = ...
    UpdateStresses_J2_FE2D(eNP1,PMAT,EPn,sigmayN,alphaN,...
    CALC_CTANG,iter,WEIGHTS)

%dbstop('6')
if nargin == 0
    load('tmp1.mat')
   % sigmayN = 0.0004*ones(size(sigmayN)) ;
    %iter=2 ;
end
ngaus = size(alphaN,1); ngausT = size(eNP1,1); nstrain = ngausT/ngaus ;
ngausE =4 ;
COMP.c1 = 1:nstrain:ngausT ;  %  Indices of Component x of stresses and strains
COMP.c2 = 2:nstrain:ngausT ;  % Indices of Component y of stresses and strains
COMP.c3 = 3:nstrain:ngausT ; % Indices of Component xy of stresses and strains
COMP.c4 = 4:nstrain:ngausT ;  % Indices of Component z of stresses and strains
% muA = PMAT.muA ; kappaA = PMAT.kappaA ; kappa = PMAT.kappa ;
% mu = PMAT.mu ;
% Initializations
EPnp1 = EPn ;
alphaNP1 = alphaN ;
sigmayNP1 = sigmayN;
%CtangNP1_B = sparse(size(BBfNW,1),size(BBfNW,2));
%CtangNP1_B= BBfNW ;
CtangM = PMAT.Celas ;
stressNP1 = zeros(size(EPnp1));
% 1. Trial step
% -------------
% Trial elastic strains
eTRelas = eNP1-EPn ;
% Deviatoric part and spherical parts
tr_elasTR = TracePlaneStrain(eTRelas,COMP) ;  % Trace
dev_elasTR = eTRelas -1/3*tr_elasTR ;
% Deviatoric trial stress
muA =reshape([PMAT.mu'; PMAT.mu'; PMAT.mu'; PMAT.mu'],ngausT,1) ;
sTR = StrainMultiplStress(dev_elasTR,2*muA,COMP) ;
% Spheric trial stress
kappaA =reshape([PMAT.kappa'; PMAT.kappa'; PMAT.kappa'; PMAT.kappa'],ngausT,1) ;
pTR = kappaA.*tr_elasTR ;

% 2. Check yield condition
%-------------------------
normSTR = sqrt(sTR(COMP.c1).^2 +sTR(COMP.c2).^2 +2*sTR(COMP.c3).^2    +sTR(COMP.c4).^2 ) ;
fTR = normSTR -sqrt(2/3)*sigmayN ;
indPLAS = find(fTR>0) ;  % Plastic gauss points
indELAS = 1:ngaus ; indELAS(indPLAS) = [] ;  % Elastic gauss points
%----------------------------------------------
% --- POINTS BEHAVING ELASTICALLY
%----------------------------------------------
%indELASm = INDEXcomponents(indELAS,nstrain)  ;
% By default
%stressNP1(indELASm)  = sTR(indELASm) + pTR(indELASm) ;
stressNP1 = sTR + pTR ;

%  Product CtangST*BBnw
% ----------------------
% METHOD_CTANG = 2 ;
% if  CALC_CTANG == 1
%     if METHOD_CTANG ==1
%         CtangNP1_B(indELASm,:)  = CtangB_elas(PMAT,BBfNW,indELASm,nstrain,indELAS) ;
%     end
% end
%----------------------------------------------
% --- POINTS BEHAVING PLASTICALLY
%----------------------------------------------
%dbstop('64')
niter_elast = 1 ; 
if ~isempty(indPLAS) && iter>niter_elast
    indPLASm = INDEXcomponents(indPLAS,nstrain)  ;
    % Plastic multiplier
    pmult = fTR(indPLAS)./(2*PMAT.mu(indPLAS) + 2/3*PMAT.Hmodulus(indPLAS)) ;
    ngausTg =length(indPLASm) ;
    COMPg.c1 = 1:nstrain:ngausTg ;  %  Indices of Component x of stresses and strains
    COMPg.c2 = 2:nstrain:ngausTg ;  % Indices of Component y of stresses and strains
    COMPg.c3 = 3:nstrain:ngausTg ; % Indices of Component xy of stresses and strains
    COMPg.c4 = 4:nstrain:ngausTg ;  % Indices of Component z of stresses and strains
    pmultA =   reshape([pmult'; pmult'; pmult'; pmult'],ngausTg,1) ; % At All components
    % Internal variables
    alphaNP1(indPLAS) = alphaN(indPLAS) + sqrt(2/3)*pmult;
    sigmayNP1(indPLAS) = sigmayN(indPLAS) + PMAT.Hmodulus(indPLAS).*(alphaNP1(indPLAS)- alphaN(indPLAS));
    %
    % Normal along the plastic flow direction
    normSTRg =  normSTR(indPLAS) ;
    normSTRgA =   reshape([normSTRg'; normSTRg'; normSTRg'; normSTRg'],ngausTg,1) ;  % All components
    normalNP = sTR(indPLASm)./normSTRgA ;
    % Deviatoric stress
    sNP = sTR(indPLASm) - 2*muA(indPLASm).*pmultA.*normalNP ;
    % Total stress
    stressNP1(indPLASm)  =   pTR(indPLASm)  + sNP ;
    % Plastic strains
    devdEP =   StressMultiplStrain(normalNP,pmultA,COMPg) ;
    EPnp1(indPLASm) =  EPn(indPLASm) + devdEP ;
    % Product C*BBfNW
    if  CALC_CTANG == 1
        %      if METHOD_CTANG ==1
        %         CtangNP1_B(indPLASm,:)  = CtangB_plas(PMAT,BBfNW,indPLASm,nstrain,indPLAS,normSTRg,pmult,normalNP) ;
        %    else
        HmodulusG = PMAT.Hmodulus(indPLAS);
        HmodulusG =reshape([HmodulusG'; HmodulusG'; HmodulusG'; HmodulusG'],length(indPLASm),1) ;
        CtangPLAS  =  Ctang_plasMethod2(nstrain,normSTRgA,pmultA...
            ,normalNP,muA(indPLASm),kappaA(indPLASm),HmodulusG,PMAT.identVOL(indPLASm,:),PMAT.identV(indPLASm,:)) ;
        CtangM(indPLASm,:) = CtangPLAS ;
        %   end
        
    end
    
% elseif  ~isempty(indPLAS) && iter<=niter_elast
%     stressNP1(indPLAS)  = sTR(indPLAS) + pTR(indPLAS) ;
%     %   if METHOD_CTANG ==1
%     %     CtangNP1_B(indPLAS,:)  = CtangB_elas(PMAT,BBfNW,indPLASm,nstrain,indPLAS) ;
%     %  end
end

if  CALC_CTANG ==1
    %dbstop('106')
    if ~isempty(WEIGHTS)
        for istrain = 1:size(CtangM,2)
            CtangM(:,istrain) =  CtangM(:,istrain).*WEIGHTS ;
        end
        
    end
    
    
    Ctang =     ConvertCmatSparseMatrix(CtangM,PMAT.indI,PMAT.indJ) ;
  %  CtangNP1_B = Ctang*BBfNW ;
end
% 
% 
% TEST_SOLUTION = 1 ; 
% if TEST_SOLUTION ==1 
%     stress = Ctang*eNP1 ;
% end

