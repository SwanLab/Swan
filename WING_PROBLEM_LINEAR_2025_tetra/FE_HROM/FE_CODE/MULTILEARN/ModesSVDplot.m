function [BasisUdef,SingVal_disp ]=...
    ModesSVDplot(DATAIN,dDOM,nfigure,LEGENDG,COLOR,DATA_REFMESH,INDICES_PROJECTS)

if nargin == 0
    load('tmp4.mat')
end
Ui_Si = [] ;
DATAIN = DefaultField(DATAIN,'NO_SVD_ON_SNAPSHOTS',0) ;
DATAIN = DefaultField(DATAIN,'NORMALIZE_COLUMNS_PROJECTS_SVD',0) ;
if DATAIN.NORMALIZE_COLUMNS_PROJECTS_SVD == 1
    for iproject =1 :length(dDOM)
        normPROJECTS(iproject) = norm(dDOM{iproject},'fro') ;
    end
    maxNORM = max(normPROJECTS) ;
    COEFF = maxNORM./normPROJECTS ;
    
    for iproject =1 :length(dDOM)
        dDOM{iproject} = dDOM{iproject}*COEFF(iproject) ;
    end
    
end

if  DATAIN.NO_SVD_ON_SNAPSHOTS ==0
    [U,S,V,h1,h2,Ui_Si] = SVD_dom(dDOM,nfigure,LEGENDG,DATAIN.NMODES_SHOW,COLOR...
        ,DATAIN,INDICES_PROJECTS ) ;
    
else
    
    U = cell2mat(dDOM) ;
    S =[] ;
    V = [] ;
    h1 = [] ;
    h2 =[] ;
    Ui_Si = [] ;
    
    
end


if ~isempty(h1)
    hold on
    %   title(['Domains type =',num2str(itype)])
    legend([h1 ],{LEGENDG})
    legend([h2 ],{LEGENDG})
end
DATAIN  = DefaultField(DATAIN,'NMODES_TRUNCATE',[]);
if isempty(DATAIN.NMODES_TRUNCATE) % & ~isempty( DATAIN.TOL_LOC)
    nnn  =size(U,2) ;
else
    nnn = min(DATAIN.NMODES_TRUNCATE,size(U,2)) ;
end
BasisUdef  = U(:,1:nnn) ; % = [2] ;
if isempty(S)
    SingVal_disp = 0 ;
else
    SingVal_disp  = S(1:nnn) ;
end

% Plot MODES
DATAIN = DefaultField(DATAIN,'PLOT_MODES',1) ;

if DATAIN.PLOT_MODES == 1
    idom = 1;
    %  ELEMS =  DATAOUT.DOMAINVAR.ListElements{idom} ;
    CNref = DATA_REFMESH.CN  ;
    COOR = DATA_REFMESH.COOR  ;
    NAME_MODES = [DATAIN.NAME_WS_MODES(1:end-4),LEGENDG ];
    DATA= [] ;
    GidPostProcessModes_dom(COOR,CNref,DATA_REFMESH.TypeElement,BasisUdef,DATA_REFMESH.posgp,...
        NAME_MODES,DATA,LEGENDG);
end

