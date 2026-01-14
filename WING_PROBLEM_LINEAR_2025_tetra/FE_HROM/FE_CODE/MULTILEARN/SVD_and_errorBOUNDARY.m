function  [U,S,V,h,h2] = SVD_and_errorBOUNDARY(SNAP,nfigure,LEGENDG,NMODES_SHOW,COLOR,DATA,...
    NODESfaces,ndim,DATAONLINE,DOFS_reference,COLUMNS_RVEloc )



%dbstop('4')
if nargin ==  0
    load('tmp.mat')
    %  DATA.SVD_RANDOM = 1 ;
end
%SNAP = cell2mat(dRVE) ;
if iscell(SNAP)
    SNAP = cell2mat(SNAP) ;
end


NODES = [] ;
%FACES_CONTACT = ones(size()) ;
for i = 1:length(NODESfaces)
    %   if FACES_CONTACT(i)==1
    NODES = [NODES; NODESfaces{i}'] ;
    %  end
end
NODES = unique(NODES) ;
DOFf = small2large(NODES,ndim) ;  % Boundary DOFf
DOFi = 1:(size(SNAP,1));  % Interior DOFf
DOFi(DOFf) = [] ;

DATA = DefaultField(DATA,'MINIMIZATON_SVD_WITH_STIFFNESS_MATRIX',0) ;

if  DATA.MINIMIZATON_SVD_WITH_STIFFNESS_MATRIX == 1
    K = RecoverStiffnessMatrixUnitaryProject(DATAONLINE) ;
    %   u = K(DOFi,DOFi)\ones(length(DOFi),1) ;
    Kffbar = K(DOFf,DOFf) - K(DOFf,DOFi)*(K(DOFi,DOFi)\K(DOFi,DOFf)) ;
    [INTTTG DOFremove DOFremove2] =  intersect(DOFf,DOFS_reference);
    DOFf_loc = 1:length(DOFf) ;
    DOFf_loc(DOFremove) = [] ;
    Kffbar =     Kffbar(DOFf_loc,DOFf_loc)  ;
    DOFf(DOFremove)=[] ;
    DOFi = [DOFi  DOFS_reference'] ;
    
    Hbar = chol(Kffbar) ;  % If its no positive definite, then check what happens with DOFS_reference
else
    Hbar = [] ;
end

DATA = DefaultField(DATA,'SVD_RANDOM',0) ;
DATA = DefaultField(DATA,'TOL_LOC',0) ;

%%
if    (DATA.REACTIONmodes_equal_DISPmodes) ~= 0
    error('Option not available for this type of approximation')
    
end

%%

DATASVD.RELATIVE_SVD = 1 ;

DATA = DefaultField(DATA,'TOLERANCE_SVD_PROJECT',[]); 

if   ~isempty(DATA.TOLERANCE_SVD_PROJECT)
     
       if  isempty(Hbar)
             [Uf,S,V,eSVD] = SVD_accuracy_project(SNAP(DOFf,:),DATA,COLUMNS_RVEloc) ; 
        else
             [U,S,V,eSVD] = SVD_accuracy_project(Hbar*SNAP(DOFf,:),DATA,COLUMNS_RVEloc) ; 
        end
else
    if DATA.SVD_RANDOM == 0
        if  isempty(Hbar)
            [Uf,S,V,eSVD] = SVDT(SNAP(DOFf,:),DATA.TOL_LOC,DATASVD) ;
        else
            [U,S,V,eSVD] = SVDT(Hbar*SNAP(DOFf,:),DATA.TOL_LOC,DATASVD) ;
        end
    else
        if  isempty(Hbar)
            [Uf,S,V,eSVD] = RSVDT(SNAP(DOFf,:),DATA.TOL_LOC,[],0,DATASVD) ;
        else
            [Uf,S,V,eSVD] = RSVDT(Hbar*SNAP(DOFf,:),DATA.TOL_LOC,[],0,DATASVD) ;
        end
    end
end




if eSVD == 0
    eSVD =  eps(1);
end

if isempty(Hbar)
    U = zeros(size(SNAP,1),size(Uf,2)) ;
    U(DOFf,:) = Uf ;
    
    Ui = SNAP(DOFi,:)*V ;
    Ui = bsxfun(@times,Ui',1./S)' ;
    U(DOFi,:) = Ui ;
else
    U = SNAP*V ;
    U = bsxfun(@times,U',1./S)' ;
end

SingVsq =  (S.*S) ;
nTOTAL = sqrt(sum(SingVsq) ); % + errorHOMOG^2)  ;
SingVsq = sort(SingVsq);
normEf2 = sqrt(cumsum(SingVsq)) ;
normEf2 = sort(normEf2,'descend');
% hold on
% lS = log10(normEf2);
% figure(1)
% hold on
% h = plot([1:length(S)-1],lS(2:end),'r');
% xlabel('MODES')
% ylabel('LOG10 SVD ERROR')
% AAA =axis ;
%
% AAA(2)  =NMODES_SHOW ;
% axis(AAA) ;

figure(nfigure)
hold on

%dbstop('32')
h = plot([ 1:length(S)],[ normEf2(2:end); eSVD]/nTOTAL*100,[COLOR,'-*']');
xlabel('MODES')
ylabel([' SVD error (%)' ])
if ~isempty(NMODES_SHOW)
    AAA =axis ;
    AAA(2)  =NMODES_SHOW ;
    axis(AAA) ;
end


figure(nfigure+1000)
hold on
h2 = plot([ 1:length(S)],[ log10([normEf2(2:end); eSVD]/nTOTAL)],[COLOR,'-*']');
xlabel('MODES')
ylabel([' logarithm SVD error' ])
if ~isempty(NMODES_SHOW)
    AAA =axis ;
    AAA(2)  =NMODES_SHOW ;
    axis(AAA) ;
end

