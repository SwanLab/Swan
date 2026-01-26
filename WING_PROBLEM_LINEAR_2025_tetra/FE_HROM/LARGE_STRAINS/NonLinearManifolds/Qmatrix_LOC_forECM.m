function Q = Qmatrix_LOC_forECM(SNAPredFINT_nw,DATA,wSTs,DATAoffline,sqrt_wST)
if nargin == 0
    load('tmp1.mat')
end
%[Q] = QbasisMatrixIntegrand_givenSNAPr(SNAPredFINT_nw,DATA,wSTs,DATAoffline)
nsnap  = length(SNAPredFINT_nw) ;

if isempty(DATAoffline.Exponent_Function_relating_global_local_TOL_fint)
    % Global SVD
    SNAPredFINT_nw  = bsxfun(@times,cell2mat(SNAPredFINT_nw),sqrt_wST) ;
    if isempty(DATAoffline.errorFINT_reactions)
        
        if ~isempty(DATA.Matrix_Include_SVD)
            Matrix_Include_SVD_w  = bsxfun(@times,DATA.Matrix_Include_SVD,sqrt_wST) ;
            DATAsvdLOC.HIDE_OUTPUT =  0;
            [Q,S,V] = SRSVD({Matrix_Include_SVD_w,SNAPredFINT_nw},[0,DATAoffline.errorFINT],DATAsvdLOC) ;
            
            %             eFROB = Matrix_Include_SVD_w-Q*(Q'*Matrix_Include_SVD_w) ;
            %             eFROB_1 = norm(eFROB,'fro')/norm(Matrix_Include_SVD_w,'fro')
            %
            %              eFROB = SNAPredFINT_nw-Q*(Q'*SNAPredFINT_nw) ;
            %             eFROB_2 = norm(eFROB,'fro')/norm(SNAPredFINT_nw,'fro') ;
            %             DATAsvdLOC.RELATIVE_SVD = 1;
            %             [Qon,Son,V] = SVDT( SNAPredFINT_nw ,DATAoffline.errorFINT,DATAsvdLOC) ;
            
            
            
        else
            
            DATAsvdLOC.HIDE_OUTPUT =  0;
            [Q,S,V] = SRSVD(SNAPredFINT_nw,DATAoffline.errorFINT,DATAsvdLOC) ;
        end
        
    else
        
        
        
        IndLOC = DATA.IndexMatrixECM.fint ; %  IndLOC =1
        
        IndLOC = repmat(IndLOC,nsnap,1) ;
        nentriesLOC = length(DATA.IndexMatrixECM.fint) + length(DATA.IndexMatrixECM.react) ;
        FactorsM = (1:nentriesLOC:nentriesLOC*nsnap)-1 ;
        IndLOC = bsxfun(@plus,IndLOC,FactorsM');
        IndFINT = IndLOC(:) ;
        
        
        IndLOC = DATA.IndexMatrixECM.react ;  % Ind_react = 2:5
        IndLOC = bsxfun(@plus,IndLOC,FactorsM');
        IndREACT = IndLOC(:) ;
        
        
        
        
        if ~isempty(DATA.Matrix_Include_SVD)
            Matrix_Include_SVD_w  = bsxfun(@times,DATA.Matrix_Include_SVD,sqrt_wST) ;
            DATAsvdLOC.HIDE_OUTPUT =  0;
            [Q,S,V] = SRSVD({Matrix_Include_SVD_w,SNAPredFINT_nw(:,IndFINT),SNAPredFINT_nw(:,IndREACT)},...
                [0,DATAoffline.errorFINT,DATAoffline.errorFINT_reactions],DATAsvdLOC) ;
        else
            
            DATAsvdLOC.HIDE_OUTPUT =  0;
            [Q,S,V] = SRSVD({SNAPredFINT_nw(:,IndFINT),SNAPredFINT_nw(:,IndREACT)},[DATAoffline.errorFINT,DATAoffline.errorFINT_reactions],DATAsvdLOC) ;
            
            % [Q,S,V] = SRSVD([SNAPredFINT_nw,SNAPreact_nw],[DATAoffline.errorFINT],DATAsvdLOC) ;
        end
        
    end
    
else
    if ~isempty(DATAoffline.errorFINT_reactions)
        error('Option only available for unified treatment of internal forces and reactions')
    end
    [Q,S,V] =  Qmatrix_for_ECM_exponentialLAWtolerances(SNAPredFINT_nw,sqrt_wST,DATAoffline) ;
    
end




%else
%    [Q,S,V,eSVD] = SVDT(SNAPredFINT,DATAoffline.errorFINT,DATAsvd) ;