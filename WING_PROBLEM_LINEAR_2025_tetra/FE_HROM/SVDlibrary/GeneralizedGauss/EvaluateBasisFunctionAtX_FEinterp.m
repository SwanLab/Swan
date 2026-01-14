function  [PHIk_y,dPHIk_y,POLYINFO]=     EvaluateBasisFunctionAtX_FEinterp(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO)
% This function returns the value of the basis functions for BAsisF at
% xNEW, as well as their derivatives. VAR_SMOOTH_FE contains the
% information required for the spatial search, as well as  the  value of
% the integrand. The integrand is assumed of the form BdomRED^T*BasisS
% The relationship between the integrand and the basis functions is
% provided by the singular values and right vectors of BasisF
% ---------------------------------------------------------------------
% Variables at our disposal
% -----------------------------------------------
% Reduced B-matrix at nodes
% VAR_SMOOTH_FE.BdomRED_nodes = BdomRED_nodes;
% % Stress basis matrix at nodes
% VAR_SMOOTH_FE.BasisS_nodes = BasisS_nodes;
% % Inverse of the singular values of BasisF times its right-singular vectors
% % --------------------------------------------------------------------------
% VAR_SMOOTH_FE.invSVsingular_F = bsxfun(@times,VrightVal_F',1./SingVal_F)' ;
% % Matrix of connectivities
% % ---------------------------
% VAR_SMOOTH_FE.CN = DATA_REFMESH.CN ;
% % Matrix of coordinates
% % ---------------------------
% VAR_SMOOTH_FE.CN = DATA_REFMESH.COOR ;
% % Type of element
% % --------------------------
% VAR_SMOOTH_FE.TypeElement = DATA_REFMESH.TypeElement ;
% % Dela. triangulation
% % --------------------------
% VAR_SMOOTH_FE.DELTRIANG = delaunayTriangulation(DATA_REFMESH.COOR);
% % ORDER POLYNOMIALS
% % ----------------------------Ç
%  nnodeE = size(DATA_REFMESH.CN,2) ;
% switch VAR_SMOOTH_FE.TypeElement
%    case 'Quadrilateral'
%           IND_POLYG_ELEMENT = [1 2 3 4 1] ;
%             if nnodeE == 9
%               ORDER_POLYNOMIALS  =[2 2] ;
%
%             elseif nnodeE == 4
%                   ORDER_POLYNOMIALS  =[1 1] ;
%             else
%                 error('element not implemented')
%             end
%     otherwise
%          error('element not implemented')
% end
% VAR_SMOOTH_FE.ORDER_POLYNOMIALS = ORDER_POLYNOMIALS;
% VAR_SMOOTH_FE.IND_POLYG_ELEMENT = IND_POLYG_ELEMENT;   % Local numbering of corner nodes (polygon)
%
%
%
%
% JAHO, 18th-April-2020, Saturday (36th day of confinment due to the  COVID19)
if nargin == 0
    load('tmp1.mat')
end


DATA = DefaultField(DATA,'OnlyCheckIfIsInside',0) ;

npoints = size(xNEW,1) ; % Number of points at which to evaluate the function
% Number of functions
nfun = length(DATALOC.ExactIntegral) ;
nsingVAL = size(VAR_SMOOTH_FE.invSVsingular_F ,2) ; % Number of singular values
%   nsingVAL    not being equal to the number of functions means that the constrained
% of exact integration of volume has been imposed. If
%  VAR_SMOOTH_FE.IMPOSE_VOLUME_CONSTRAINT =3, then, the first function is
%  constant and equal to 1
if VAR_SMOOTH_FE.IMPOSE_VOLUME_CONSTRAINT == 3   && nfun ~=nsingVAL
    CONSTANT_FUNCTION = VAR_SMOOTH_FE.ExactIntegral(1)/sum(VAR_SMOOTH_FE.wSTs);
    INDICES_FUNCTIONS = 2:nfun ;
elseif VAR_SMOOTH_FE.IMPOSE_VOLUME_CONSTRAINT ~= 3   && nfun ~=nsingVAL
    error('Incompatible options')
else
    CONSTANT_FUNCTION = 1;
    INDICES_FUNCTIONS = 1:nfun ;
end



PHIk_y = CONSTANT_FUNCTION*ones(npoints,nfun) ;  % Value of the function at xNEW. Initialization with ones
ndim =  size(VAR_SMOOTH_FE.COOR,2) ;
dPHIk_y = cell(1,ndim) ;
for idim = 1:ndim
    dPHIk_y{idim} = zeros(npoints,nfun) ;   % Initialization with zeros
end

%
% dPHIk_y{1} = zeros(npoints,nfun) ;  % Derivative x at xNEW
% dPHIk_y{2} = zeros(npoints,nfun) ;  % Devirative y at XNEW
nelem = size(VAR_SMOOTH_FE.CN,1)  ;
nnodeE = size(VAR_SMOOTH_FE.CN,2)  ;
ndim = size(VAR_SMOOTH_FE.COOR,2)  ;


POLYINFO = DefaultField(POLYINFO,'COEFFSpolynomial',cell(nelem,1)) ;  % Coefficients of the interpolation polynomial for each element



% find nodes which are closer to the nodes of the FE discretization
INDEXES_NEAR = nearestNeighbor(VAR_SMOOTH_FE.DELTRIANG, xNEW);
ELEMENTS_CONTAINING_xNEW = zeros((npoints),1) ; % Indices of the elements containing the points
IND_POLYG = VAR_SMOOTH_FE.IND_POLYG_ELEMENT  ;  % Local numbering of corner nodes

% Loop over new points
%IND_gausspointsCLOSE = zeros(npoints,1) ;  % This weill be used for later purposes

INDfin = 0 ;
inew = 1;

POLYINFO = DefaultField(POLYINFO,'SCALING_VARIABLES',[]) ;
POLYINFO.SCALING_VARIABLES = DefaultField(POLYINFO.SCALING_VARIABLES,'LENGTH',zeros(nelem,1)) ;
POLYINFO.SCALING_VARIABLES = DefaultField(POLYINFO.SCALING_VARIABLES,'coorREF',zeros(nelem,ndim)) ;

SCALING_ACTIVE = 1;
BdomRED_interp = cell(npoints,1) ;

while inew <=npoints
    xLOC = xNEW(inew,:); % Coordinate of the point at which the function is to be evaluated
    INDnear = INDEXES_NEAR(inew) ;  % Nearest FE mesh  node to xLOC
    [ELEMnear, aaaa ]= find(VAR_SMOOTH_FE.CN == INDnear) ;  % Elements sharing INDnear.
    ielem = 1;
    elemCONTAINER = [] ;
    % Searching for the element containing xLOC --> elemCONTAINER
    while ielem <= length(ELEMnear)
        elemLOC = ELEMnear(ielem) ;
        
        [inELEM,onELEM] = IsInsideGeneral(xLOC,VAR_SMOOTH_FE.COOR,VAR_SMOOTH_FE.CN,elemLOC,IND_POLYG) ;
        
        if inELEM == 1 || onELEM == 1
            elemCONTAINER  = elemLOC ;
            break
        end
        ielem = ielem + 1;
    end
    if isempty(elemCONTAINER)
        disp(['Point outside the domain...'])
        PHIk_y = [];
        break
    else
        
        ELEMENTS_CONTAINING_xNEW(inew ) = elemCONTAINER ;  % Element containing the point
        
    end
    
    
    if DATA.OnlyCheckIfIsInside == 0
        
        % Interpolation
        % --------------
        if ~isempty(VAR_SMOOTH_FE.BasisS_nodes)
            % interpolation using  values at nodes (smoothed)
            NODESloc = VAR_SMOOTH_FE.CN(elemCONTAINER,:) ;
        else
            % Interpolation using values at GAUSS POINTS
            NODESloc = small2large(elemCONTAINER,VAR_SMOOTH_FE.ngausE) ;
        end
        
        
        if ~isempty(POLYINFO.COEFFSpolynomial{elemCONTAINER})
            %  N = POLYINFO.SHAPEFUNCTIONelem{elemCONTAINER}
            COEFFSpol = POLYINFO.COEFFSpolynomial{elemCONTAINER} ;
            LelemSCALING = POLYINFO.SCALING_VARIABLES.LENGTH(elemCONTAINER) ;
            coorREF = POLYINFO.SCALING_VARIABLES.coorREF(elemCONTAINER,:) ;
        else
            % -------------------------------------------------
            % Now we have the element defined by the nodes NODESloc = CN(elemCONTAINER,:)
            % with coordinates COOR(NODESloc,:), which contains xLOC
            %     indGAUSS = small2large(elemCONTAINER,ngausELEM) ;
            %     IND_gausspointsCLOSE(inew) = indGAUSS(1) ;
            if ~isempty(VAR_SMOOTH_FE.BasisS_nodes)
                COORelement = VAR_SMOOTH_FE.COOR(NODESloc,:) ; % Nodal coordinates of the element
            else
                COORelement = VAR_SMOOTH_FE.COORg(NODESloc,:) ; %   coordinates of the Gauss points
            end
            
            % Scaling (required to avoid that P is badly scaled)
            LelemSCALING = norm(COORelement(1,:)-COORelement(2,:)) ;
            coorREF  = COORelement(1,:) ;
            POLYINFO.SCALING_VARIABLES.LENGTH(elemCONTAINER) = LelemSCALING ;
            POLYINFO.SCALING_VARIABLES.coorREF(elemCONTAINER,:) = coorREF ;
            
            if SCALING_ACTIVE == 1
                COORelement = bsxfun(@minus,COORelement',coorREF')'/LelemSCALING   ;
            end
            P = CoordinateMatrixPolynomial(COORelement,VAR_SMOOTH_FE.ORDER_POLYNOMIALS)  ;
            
            
            COEFFSpol =  inv(P) ; % P\eye(nnodeE);  % Coefficients of the polynomials
            POLYINFO.COEFFSpolynomial{elemCONTAINER} = COEFFSpol ;
        end
        % Evaluating the shape functions a point xLOC
        factorDERIVATIVE = 1;
        if SCALING_ACTIVE == 1
            xLOC = (xLOC-coorREF)/LelemSCALING   ;
            factorDERIVATIVE = 1/LelemSCALING ;
        end
        
        
        [Pevaluate PevaluateDER]= CoordinateMatrixPolynomial(xLOC,VAR_SMOOTH_FE.ORDER_POLYNOMIALS)  ;
        N = Pevaluate*COEFFSpol ;  % Shape functions at the given points COORevaluate
        derN = zeros(length(PevaluateDER),nnodeE) ;
        for idim = 1:length(PevaluateDER)
            derN(idim,:) =   PevaluateDER{idim}*COEFFSpol*factorDERIVATIVE ;
        end
        ndim = length(PevaluateDER) ;
        % Computing matrices employed in the interpolation
        % -----------------------------------------------
        M = N'*N ;
        G = cell(1,ndim) ;
        for idim = 1:ndim
            G{idim} = derN(idim,:)'*N + N'*derN(idim,:) ;
        end
        
        % Computing B^T BasisS  (in analogy to BasisFfromStress.m)
        % ---------------------
        nstrain = VAR_SMOOTH_FE.nstrain ;
        DOFs = small2large(NODESloc,nstrain) ;
        if ~isempty(VAR_SMOOTH_FE.BasisS_nodes)
            BasisS = VAR_SMOOTH_FE.BasisS_nodes(DOFs,:) ;
            BdomRED = VAR_SMOOTH_FE.BdomRED_nodes(DOFs,:) ;
        else
            BasisS = VAR_SMOOTH_FE.BasisS_gauss(DOFs,:) ;
            BdomRED = VAR_SMOOTH_FE.BdomRED_gauss(DOFs,:) ;
            
        end
        
        
        nmodesS  = size(BasisS,2) ;         nmodesU  = size(BdomRED,2) ;
        nmodesF = nmodesS*nmodesU;
        A = zeros(1,nmodesF) ;
        derA = zeros(ndim,nmodesF) ;
        BdomRED_interp{inew} = zeros(nstrain,nmodesU);
        
        I =1 ;
        for j = 1:nmodesU
            for  k =1:nmodesS
                Aloc = 0;
                derAloc = zeros(ndim,1) ;
                for  l = 1:nstrain
                    BdomRED_interp{inew}(l:nstrain:end,j) = N*BdomRED(l:nstrain:end,j) ;
                    Aloc = Aloc + BdomRED(l:nstrain:end,j)'*M*BasisS(l:nstrain:end,k) ;
                    for idim = 1:ndim
                        derAloc(idim,:) =  derAloc(idim,:)  + BdomRED(l:nstrain:end,j)'*G{idim}*BasisS(l:nstrain:end,k) ;
                    end
                end
                A(I) = Aloc;
                derA(:,I) = derAloc;
                I = I+1;
            end
        end
        
        %
        PHIk_y(inew,INDICES_FUNCTIONS) = A*VAR_SMOOTH_FE.invSVsingular_F ;
        for idim = 1:ndim
            dPHIk_y{idim}(inew,INDICES_FUNCTIONS) = derA(idim,:)*VAR_SMOOTH_FE.invSVsingular_F ;
        end
        
    end
    
    inew = inew +1 ;
end


POLYINFO.ELEMENTS_CONTAINING_xNEW = ELEMENTS_CONTAINING_xNEW ; %(inew ) = elemCONTAINER ;  % Element containing the point
POLYINFO.BdomRED_interp = BdomRED_interp ;


%
% % COMPUTATION OF THE LOCAL  B-MATRIX OF THIS POINT
% % -------------------------------------------------
% % -------------------------------------------------
% % Now we have the element defined by the nodes NODESloc = CN(elemCONTAINER,:)
% % with coordinates COOR(NODESloc,:), which contains xLOC
% indGAUSS = small2large(elemCONTAINER,ngausELEM) ;
% IND_gausspointsCLOSE(inew) = indGAUSS(1) ;
% NODESloc = CN(elemCONTAINER,:) ;
% [N,  BeTILDE] = ComputeBmatrixDirectly(COOR(NODESloc,:),xLOC,TypeElement) ;
% ndim = size(COOR,2) ;
% % number of stain components
%
% Bpoint = QtransfB(BeTILDE,ndim) ;
% if size(Bpoint,1) ~= nstrain
%     % Plane strain
%     Bpoint = [Bpoint; zeros(1,size(Bpoint,2))] ;
% end
%
% % Contribution to the strain at this point
%
% DOFSloc = small2large(NODESloc,ndim) ;
% Bpoint_hyper = Bpoint*BasisU(DOFSloc,:) ;
% INDINI = INDfin +1;
% INDfin = INDINI +nstrain-1;
% BstHYPER(INDINI:INDfin,:) = Bpoint_hyper ;
% inew = inew +1 ;
% end
