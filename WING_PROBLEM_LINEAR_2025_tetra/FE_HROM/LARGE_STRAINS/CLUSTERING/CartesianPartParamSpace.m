function [COOR,CN,NeighboringElements] = CartesianPartParamSpace(DATAoffline,xMIN,xMAX)
% Limits cartesian domain
xMIN = min(MATRIX_POINTS_SPACE_PARAMETER) ;
xMAX = max(MATRIX_POINTS_SPACE_PARAMETER) ;
 DATAoffline = DefaultField(DATAoffline,'PARTITION_PARAMETER',[]) ;
        if  isempty(DATAoffline.PARTITION_PARAMETER)
            nclusters = DATAoffline.NCLUSTERS_BASIS_DISP ;
            ndiv = ceil((nclusters)^(1/ninp))*ones(1,ninp) ;
            Xinp = cell(1,ninp) ;
            for iparaminp = 1:ninp
                Xinp{iparaminp} = linspace(xMIN(iparaminp),xMAX(iparaminp),ndiv) ;
            end
        else
            Xinp = cell(1,ninp) ;
            for iparaminp = 1:ninp
                Xinp{iparaminp} = xMIN(iparaminp) + DATAoffline.PARTITION_PARAMETER{iparaminp}*(xMAX(iparaminp)-xMIN(iparaminp)) ;
            end
            [aaaa,ndiv] = cellfun(@size,Xinp)  ;
        end
        
        
        
        
        XGRID = cell(size(Xinp)) ;
        [XGRID{:}]  = ndgrid(Xinp{:}) ;
        
        % Coordinates
        % ************
        COOR = cell(1,ninp) ;
        for iparaminp = 1:ninp
            COOR{iparaminp} = XGRID{iparaminp}(:) ;
        end
        COOR = cell2mat(COOR) ;
        
        % Connectivities
        % **************
        nnodes = size(COOR,1) ;
        if ninp == 1
            CN = [1:nnodes; 2:nnodes+1]' ;
        elseif ninp == 2
            nelem = (ndiv(1)-1)*(ndiv(2)-1) ;
            nnodeE = 4 ;
            CN = zeros(nelem, nnodeE) ;
            for yy = 1:(ndiv(2)-1)
                for xx = 1:(ndiv(1)-1)
                    NODO_1 = ndiv(1)*(yy-1)+ xx ;
                    NODO_2 = NODO_1 + 1 ;
                    NODO_3 = NODO_1 + ndiv(1) ;
                    NODO_4 = NODO_3 + 1;
                    e = (ndiv(1)-1)*(yy-1)+ xx ;
                    CN(e,:) = [NODO_1,NODO_2,NODO_3,NODO_4] ;
                end
            end
            [InvCNmatrix, ElemNode, NeighboringElements, ElemShared]= InverseConnectMatrix(CN) ;
        else
            error('Option not impleemnted')
        end
        
        