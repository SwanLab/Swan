function [BasisUdef,SingVal_disp] = BasisUdef_computation(DATA_TRAINING,DATAIN)


disp('--------------------------')
disp('Collecting disp. snapshots')
disp('--------------------------')
%  MAXd = [] ;
DATAIN = DefaultField(DATAIN,'DOMAINS_TO_INCLUDE_TRAINING',[]) ;
dDOM = cell(1,length(DATA_TRAINING)) ;  % Displacement snapshots
% Loop over FE projects

for iproject = 1:length(DATA_TRAINING)
    disp('------------------------------------')
    [~, bbb] = fileparts(DATA_TRAINING{iproject})
    disp(['PROJECT = ',bbb]) ;
    disp('------------------------------------')
    % Reading displacements vectors of each domain (stored in dDOM)
    [dDOM{iproject},DATAOUT ]= ExtractDisplMatrix(DATA_TRAINING{iproject},DATAIN) ;
    % DomdDOMains to be included in the SVD
    if ~isempty(DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject})
        
        selected_columnsLOC = DATAIN.DOMAINS_TO_INCLUDE_TRAINING{iproject} ;
        selected_columnsLOC = intersect(1:size(dDOM{iproject},2),selected_columnsLOC) ;
        
        dDOM{iproject} = dDOM{iproject}(:,selected_columnsLOC) ;
    end
end
% % -----------------------------------------------------
%% Determination of deformational modes modes (SVD of dDOM)
%------------------------------------------------------
% dbstop('111')

DATAIN.LEGEND_GRAPHS = 'Displacements' ;
nTYPErve  =1
nfigure = 1;
LEGENDG = 'Disp. ' ;
COLOR = 'r' ;

DATAIN = DefaultField(DATAIN,'NMODES_PROJECT_DISPLACEMENTS',[]) ;
Ui_Si = [] ;

DATAIN.NMODES_PROJECT_LOC = DATAIN.NMODES_PROJECT_DISPLACEMENTS ;
[U,S,V,h1,h2,Ui_Si] = SVD_dom(dDOM,nfigure,LEGENDG,NMODES_SHOW,COLOR...
    ,DATAIN ) ;


if ~isempty(h1)
    hold on
    %   title(['Domains type =',num2str(itype)])
    legend([h1 ],{'DISPLACEMENTS'})
    legend([h2 ],{'DISPLACEMENTS'})
end

if isempty(DATAIN.NMODES_TRUNCATE) % & ~isempty(DATAIN.TOL_LOC)
    nnn  =size(U,2) ;
else
    nnn = min(DATA.NMODES_TRUNCATE,size(U,2)) ;
end
BasisUdef  = U(:,1:nnn) ; % = [2] ;
if isempty(S)
    SingVal_disp{itype} = 0 ;
else
    SingVal_disp  = S(1:nnn) ;
end

% Plot MODES
idom = 1;
ELEMS =  DATAOUT.DOMAINVAR.ListElements{idom} ;
CNref = DATAOUT.CN(ELEMS,:) ;
NODES = DATAOUT.DOMAINVAR.ListNodesDom{idom} ;
COOR = DATAOUT.COOR(NODES,:) ;
NAME_MODES = DATAIN.NAME_WS_MODES(1:end-4) ;

GidPostProcessModes_dom(COOR,CNref,DATAOUT.TypeElement,BasisUdef,DATAOUT.posgp,...
    NAME_MODES,DATA,[]);