function [PHIk_y,dPHIk_y,POLYINFO]=  EvaluateBasisFunctionDIRECTFITopt(xNEW,DATA,VAR_SMOOTH_FE,POLYINFO)
% See
%  Similar to EvaluateBasisFunctionDIRECTFIT.m, but without evaluating the
%  function at those points in which the points remain stationary
%  (specified by
% POLYINFO.PREVIOUS_STEP.IND_POINTS_CHANGE_POSITION


if nargin == 0
    load('tmp.mat')
end

IndPointsChangePosition = POLYINFO.PREVIOUS_STEP.IND_POINTS_CHANGE_POSITION ;


DATA = DefaultField(DATA,'OnlyCheckIfIsInside',0) ;

npoints = size(xNEW,1) ; % Number of points at which to evaluate the function
% Number of functions to be integrated
nfun = length(VAR_SMOOTH_FE.ExactIntegral) ;
% Initialization outputs
PHIk_y =zeros(npoints,nfun) ;  % Value of the integrand function at xNEW.
ndim =  size(VAR_SMOOTH_FE.COOR,2) ; % Number of spatial dimensions
dPHIk_y = cell(1,ndim) ;  % Value of the gradient of the integrand functions at xNEW
for idim = 1:ndim
    dPHIk_y{idim} = zeros(npoints,nfun) ;   % Initialization with zeros
end
nelem = size(VAR_SMOOTH_FE.CN,1)  ;  % Total number of FE elements
nnodeE = size(VAR_SMOOTH_FE.CN,2)  ;  % Number of nodes per element

% Coefficients of the interpolation polynomial for each element
% At the beginning, they are empty. As the search algorithm proceeds,
% coefficients are calculated and stored in this cell array.
POLYINFO = DefaultField(POLYINFO,'COEFFSpolynomial',cell(nelem,1)) ;

%
% Find the nodes of the FE discretization which are closer to xNEW
INDEXES_NEAR = nearestNeighbor(VAR_SMOOTH_FE.DELTRIANG, xNEW);
ELEMENTS_CONTAINING_xNEW = zeros((npoints),1) ; % Indices of the elements containing the points
IND_POLYG = VAR_SMOOTH_FE.IND_POLYG_ELEMENT  ;  % Local numbering of corner nodes
%
%INDfin = 0 ;
inew = 1;

% POLYINFO = DefaultField(POLYINFO,'SCALING_VARIABLES',[]) ;
% POLYINFO.SCALING_VARIABLES = DefaultField(POLYINFO.SCALING_VARIABLES,'LENGTH',zeros(nelem,1)) ;
% POLYINFO.SCALING_VARIABLES = DefaultField(POLYINFO.SCALING_VARIABLES,'coorREF',zeros(nelem,ndim)) ;

%SCALING_ACTIVE = 1;
BdomRED_interp = cell(npoints,1) ;  % Interpolation of the reduced B-matrix

while inew <=npoints
      
    III = find(inew == IndPointsChangePosition) ; 
    
    if   isempty(III) 
        ELEMENTS_CONTAINING_xNEW(inew ) = POLYINFO.PREVIOUS_STEP.ELEMENTS_CONTAINING_xNEW(inew) ; 
         PHIk_y(inew,:)  =    POLYINFO.PREVIOUS_STEP.PHIk_y(inew,:) ; 
         for idim = 1:ndim 
              dPHIk_y{idim}(inew,:) =  POLYINFO.PREVIOUS_STEP.dPHIk_y{idim}(inew,:) ; 
         end
    else
        
        
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
                %             if elemCONTAINER == 99
                %                 disp('Borrar esto')
                %             end
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
            
            
            % INDEXES ASSOCIATED GAUSS POINTS OF THE ELEMENT
            NODESloc = small2large(elemCONTAINER,VAR_SMOOTH_FE.ngausE) ;
            if ~isempty(POLYINFO.COEFFSpolynomial{elemCONTAINER})
                % THE COEFFICIENTS HAVE BEEN ALREADY CALCULATED
                COEFFSpol = POLYINFO.COEFFSpolynomial{elemCONTAINER} ;
                LelemSCALING = POLYINFO.SCALING_VARIABLES.LENGTH(elemCONTAINER) ;
                coorREF = POLYINFO.SCALING_VARIABLES.coorREF(elemCONTAINER,:) ;
            else
                % cOORDINATES OF THE GAUSS POINTS OF THE ELEMENT CONTAINING THE POINT
                % UNDER CONSIDERATION
                COORelement = VAR_SMOOTH_FE.COORg(NODESloc,:) ; %   coordinates of the Gauss points
                % Scaling (required to avoid that P is badly scaled)
                LelemSCALING = norm(COORelement(1,:)-COORelement(2,:)) ; % Typical length
                coorREF  = COORelement(1,:) ;  % Coordinates reference point
                POLYINFO.SCALING_VARIABLES.LENGTH(elemCONTAINER) = LelemSCALING ;
                POLYINFO.SCALING_VARIABLES.coorREF(elemCONTAINER,:) = coorREF ;
                % Transformed coordinates (reltive to coorREF), and scaled by   LelemSCALING
                COORelement = bsxfun(@minus,COORelement',coorREF')'/LelemSCALING   ;
                % Coefficientes of  the polynomial
                P = CoordinateMatrixPolynomial(COORelement,VAR_SMOOTH_FE.ORDER_POLYNOMIALS)  ;
                
                COEFFSpol =  inv(P) ; % P\eye(nnodeE);  % Coefficients of the polynomials
                POLYINFO.COEFFSpolynomial{elemCONTAINER} = COEFFSpol ;
            end
            % Evaluating the shape functions (and derivatives) at  point xLOC
            % -------------------------------------------
            xLOC = (xLOC-coorREF)/LelemSCALING   ;  % Transformed coordinate of the point under study
            factorDERIVATIVE = 1/LelemSCALING ;    % This for computing the derivatives
            % -------------------------------------------------------------------------------------------
            % Computatio  of the shape functions (recall that the "nodes" here of such functions are the Gauss points)
            [Pevaluate, PevaluateDER]= CoordinateMatrixPolynomial(xLOC,VAR_SMOOTH_FE.ORDER_POLYNOMIALS)  ;
            N = Pevaluate*COEFFSpol ;  % Shape functions at the given points COORevaluate
            PHIk_y(inew,:)  =   N*VAR_SMOOTH_FE.BasisIntegrand(NODESloc,:) ;
            for idim = 1:length(PevaluateDER)
                dN = PevaluateDER{idim}*COEFFSpol ;
                dPHIk_y{idim}(inew,:) = factorDERIVATIVE*dN*VAR_SMOOTH_FE.BasisIntegrand(NODESloc,:) ;
            end
            
            
            
            
            
            
            
            %  ndim = length(PevaluateDER) ;
            %         % Computing matrices employed in the interpolation
            %         % -----------------------------------------------
            %         M = N'*N ;
            %         G = cell(1,ndim) ;
            %         for idim = 1:ndim
            %             G{idim} = derN(idim,:)'*N + N'*derN(idim,:) ;
            %         end
            
            %         % Computing  projected B matrix  (this should be done outside in future versions, 16-Dec-2021)
            %         % ---------------------
            %         nstrain = VAR_SMOOTH_FE.nstrain ; % Number of strain entries
            %         DOFs = small2large(NODESloc,nstrain) ;  % INdexes associated to the Gauss points of the elements
            %         BdomRED = VAR_SMOOTH_FE.BdomRED_gauss(DOFs,:) ;  % reduced B matrix  of the element
            %         nmodesU  = size(BdomRED,2) ;
            %         BdomRED_interp{inew} = zeros(nstrain,nmodesU);
            
            
            
            
            
            
            %         PHIk_y(inew,INDICES_FUNCTIONS) = A*VAR_SMOOTH_FE.invSVsingular_F ;
            %         for idim = 1:ndim
            %             dPHIk_y{idim}(inew,INDICES_FUNCTIONS) = derA(idim,:)*VAR_SMOOTH_FE.invSVsingular_F ;
            %         end
            
            
            %
            %         I =1 ;
            %         for j = 1:nmodesU
            %             for  k =1:nmodesS
            %                 Aloc = 0;
            %                 derAloc = zeros(ndim,1) ;
            %                 for  l = 1:nstrain
            %                     BdomRED_interp{inew}(l:nstrain:end,j) = N*BdomRED(l:nstrain:end,j) ;
            %                     Aloc = Aloc + BdomRED(l:nstrain:end,j)'*M*BasisS(l:nstrain:end,k) ;
            %                     for idim = 1:ndim
            %                         derAloc(idim,:) =  derAloc(idim,:)  + BdomRED(l:nstrain:end,j)'*G{idim}*BasisS(l:nstrain:end,k) ;
            %                     end
            %                 end
            %                 A(I) = Aloc;
            %                 derA(:,I) = derAloc;
            %                 I = I+1;
            %             end
            %         end
            
            %
            %         PHIk_y(inew,INDICES_FUNCTIONS) = A*VAR_SMOOTH_FE.invSVsingular_F ;
            %         for idim = 1:ndim
            %             dPHIk_y{idim}(inew,INDICES_FUNCTIONS) = derA(idim,:)*VAR_SMOOTH_FE.invSVsingular_F ;
            %         end
            
        end
        
    end
    
    inew = inew +1 ;
end


POLYINFO.ELEMENTS_CONTAINING_xNEW = ELEMENTS_CONTAINING_xNEW ; %(inew ) = elemCONTAINER ;  % Element containing the point
%POLYINFO.BdomRED_interp = BdomRED_interp ;
