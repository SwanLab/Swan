function  [Relem]= ComputeRelemBoundALL(COOR,CNb,TypeElementB) ;
%%%%
% See ALEX_THESIS_mine.pdf
if nargin == 0
    load('tmp2.mat')
    COOR = [-1 -1; 1 -1 ; 1 1 ; -1 1] ;
    CNb = [1 2; 2 3 ; 3 4; 4 1] ; 
    TypeElementB = 'Linear' ;
    ndim = 2 ; 
end
nnode = size(COOR,1); ndim = size(COOR,2)  ; nelemB = size(CNb,1); nnodeEb = size(CNb,2) ;
% Shape function routines (for calculating shape functions and derivatives)
TypeIntegrand = 'RHS';
[weig,posgp,shapef,dershapef] = ComputeElementShapeFun(TypeElementB,nnodeEb,TypeIntegrand) ;
ngaus = length(weig) ;

%% Elemental coordinates (boundary elements)
XeALL= COORallELEM(ndim,nelemB,nnodeEb,CNb,COOR) ;

%  Unit vectors normal to each boundary element (nelem*dim x 1)
[SeALL normals]= ChangeCoordBndVect_andNORMALS(XeALL,ndim);


if ndim == 3
    nstrain = 6;
else
    nstrain = 3;
end
% Matrix with all elemental Rs, including weights and Jacobians
% ---------------------------
Relem = zeros(nelemB*nstrain,ndim*nnodeEb) ;

if ndim == 3
 n1 = normals(1:ndim:(ndim*nelemB)) ;
 n2 = normals(2:ndim:(ndim*nelemB)) ;
 n3 = normals(3:ndim:(ndim*nelemB)) ;
else
 n1 = normals(1:ndim:(ndim*nelemB)) ;
 n2 = normals(2:ndim:(ndim*nelemB)) ;
end



for  g = 1:ngaus
    % Matrix of shape functions and  derivatives at Gauss point "g"
    Nbar =  shapef(g,:) ;
    BeXi = dershapef(:,:,g) ;
    % ----------------------------------------
    % Coordinate change (vectorized form)
    %%%%%
    % Jacobian Matrix for the g-th G. point of all elements %
    JeALL = SeALL*BeXi' ;
    %%%%%%%%%
    % JAcobian
    detJeALL= determinantVECTORIZE(JeALL,ndim-1) ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for inode = 1: nnodeEb
        iCOL = (inode-1)*ndim+1   ;
        if ndim == 3
            % ----------
            for istrain = 1:nstrain ;
                % ----------
                ROWS = istrain:nstrain:nelemB*nstrain ;
                if istrain == 1
                    Relem(ROWS,iCOL) = Relem(ROWS,iCOL) + n1*Nbar(inode).*detJeALL*weig(g) ;
                elseif istrain == 2
                    Relem(ROWS,iCOL+1) = Relem(ROWS,iCOL+1) + n2*Nbar(inode).*detJeALL*weig(g) ;
                elseif istrain == 3
                    Relem(ROWS,iCOL+2) = Relem(ROWS,iCOL+2) + n3*Nbar(inode).*detJeALL*weig(g) ;
                elseif istrain == 4
                    Relem(ROWS,iCOL+1) = Relem(ROWS,iCOL+1) + n3*Nbar(inode).*detJeALL*weig(g) ;
                    Relem(ROWS,iCOL+2) = Relem(ROWS,iCOL+2) + n2*Nbar(inode).*detJeALL*weig(g) ;
                elseif istrain == 5
                    Relem(ROWS,iCOL) = Relem(ROWS,iCOL) + n3*Nbar(inode).*detJeALL*weig(g) ;
                    Relem(ROWS,iCOL+2) = Relem(ROWS,iCOL+2) + n1*Nbar(inode).*detJeALL*weig(g) ;
                elseif istrain ==6
                    Relem(ROWS,iCOL) = Relem(ROWS,iCOL) + n2*Nbar(inode).*detJeALL*weig(g) ;
                    Relem(ROWS,iCOL+1) = Relem(ROWS,iCOL+1) + n1*Nbar(inode).*detJeALL*weig(g) ;
                end
                
            end
            
            
        else
            for istrain = 1:nstrain ;
                % ----------
                ROWS = istrain:nstrain:nelemB*nstrain ;
                if istrain == 1
                    Relem(ROWS,iCOL) = Relem(ROWS,iCOL) + n1*Nbar(inode).*detJeALL*weig(g) ;
                elseif istrain == 2
                    Relem(ROWS,iCOL+1) = Relem(ROWS,iCOL+1) + n2*Nbar(inode).*detJeALL*weig(g) ;
                elseif istrain == 3
                    Relem(ROWS,iCOL) = Relem(ROWS,iCOL) + n2*Nbar(inode).*detJeALL*weig(g) ;
                    Relem(ROWS,iCOL+1) = Relem(ROWS,iCOL+1) + n1*Nbar(inode).*detJeALL*weig(g) ;
                end
                
            end
        end
        
        
    end
    
end


end



%
%
% % Therefore, the sought-after Nelem matrix can be computed by simply making  nelem tiling
% % copies of Ne
% NelemB = repmat(shapef,nelemB,1);
% % -------------------------------
% %%%% Computation of the boundary weights wSTb
% % -------------------------------
% wSTb = zeros(nelemB*ngaus,1) ;  % Vector containinig the product of weights and Jacobians at all gauss points
% % COORDINATE MATRIX (for boundary elements) arranged in a nelemB*ndim x nnodeEB matrix
% % Change of coordinates
% % Let us define a matrix ROWSgauss such that   Belem(ROWSgauss(g,:),:) returns
% % the B-matrices of the g-th points of all elements
% weigREP = repmat(weig',nelemB,1)  ; % nelem x 1 tiling copies of weig
%
% for  g = 1:ngaus
%     % Matrix of derivatives for Gauss point "g"
%     BeXi = dershapef(:,:,g) ;
%     % ----------------------------------------
%     % Coordinate change (vectorized form)
%
%
%     %%%%%
%     % Jacobian Matrix for the g-th G. point of all elements %
%     JeALL = SeALL*BeXi' ;
%     %%%%%%%%%
%     % JAcobian
%     detJeALL= determinantVECTORIZE(JeALL,ndimB) ;
%     % Weight vectors
%     % --------------
%     wLOCa = detJeALL.*weigREP(g:ngaus:nelemB*ngaus) ;
%     wSTb(g:ngaus:nelemB*ngaus) =  wLOCa ;
% end


