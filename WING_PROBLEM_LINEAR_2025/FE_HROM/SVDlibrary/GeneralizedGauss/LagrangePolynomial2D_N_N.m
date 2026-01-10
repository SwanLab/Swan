function  [A, dA_x,dA_y,DATAOUT ]=    LagrangePolynomial2D_N_N(xLIM,P,xMAT,DATA)
% See comments in LagrangePolynomial2D_N_N_aux.mlx
if nargin == 0
    load('tmp.mat')
end
DATAOUT = [] ; dA_x = [] ; dA_y = [] ; 
DATA = DefaultField(DATA,'EVALUATE_GRADIENT',0) ;
DATA = DefaultField(DATA,'EVALUATE_SHAPE_FUNCTIONS',0) ;

EVALUATE_GRADIENT = DATA.EVALUATE_GRADIENT ;
 
% 1) Nshape = M points x  nmonomials, i.e., number of shape functions
[Nshape, B_x,B_y,DATAOUTlagrag] =    LagrangePolynomial2D(xLIM,P,xMAT) ;

DATAOUT.Nshape = Nshape ; 

 

% 
%2) Matrix of derivatives  (with respect to xi1_1 and xi_2 )
B_xi ={B_x ;B_y} ;  %
% 
% 3) Matrix of second derivatives (with respect to xi1_1 and xi_2 )
H_xi = {DATAOUTlagrag.H_xx,DATAOUTlagrag.H_xy
    DATAOUTlagrag.H_xy,DATAOUTlagrag.H_yy} ;


% 3) Parameterization in terms of alpha
% --------------------------------------
X = CoordElementsFromShapeParameters(DATAOUTlagrag,DATA) ;

if DATA.EVALUATE_SHAPE_FUNCTIONS == 1
    
     Xdistorted  = Nshape*X(1,:)' ; 
     Ydistorted  = Nshape*X(2,:)' ; 
     DATAOUT.COORdistorted = [Xdistorted,Ydistorted] ; 
end

DATAOUT.Xcoor = X ; 


% 4) Jacobian matrix
% ------------------------
%
% Jacobian matrix (vectorized)
ndim = 2;
npoints =size(xMAT,1) ;
J =zeros(ndim*npoints,ndim) ;
npointsE =  size(X,2) ;
ipoints = 1:npoints ;
for  idim = 1:ndim
    idimLOC = (ipoints-1)*ndim+idim ;
    for Ipoints = 1:npointsE
        for jdim = 1:ndim
            J(idimLOC,jdim) = J(idimLOC,jdim) + X(idim,Ipoints)*B_xi{jdim}(ipoints,Ipoints) ;
        end
    end
end

% Determinant
% ------------------------------
detJ= determinantVECTORIZE(J,ndim) ;
DATAOUT.detJ = detJ ; 

% Inverse
% ----------
invJ= inverseVECTORIZED(J,ndim,detJ) ;

% % B matrix (with respect to x1 and x2)
% % ------------------------------------
% % See LagrangePolynomial2D_B_B_aux.mlx
% B = {zeros(size(B_xi{1})),zeros(size(B_xi{2}))} ; % Initialization
% ipoints = 1:npoints ;
% for idim = 1:ndim
%     for Inodes = 1:npointsE
%         for hdim = 1:ndim
%             hdimGLO = (ipoints-1)*ndim+hdim ;
%             B{idim}(:,Inodes) = B{idim}(:,Inodes) + B_xi{hdim}(:,Inodes).*invJ(hdimGLO,idim) ;
%         end
%     end
% end

if EVALUATE_GRADIENT == 1
    
    % Derivative Jacobian
    % -------------------
 R = cell(1,ndim) ;
  S = cell(1,ndim) ;
    for n_dim = 1:ndim
        S{n_dim} = zeros(size(J)) ;
        for idim = 1:ndim
            idimGLO = (ipoints-1)*ndim+idim ;
            for kdim = 1:ndim
                for Jnode = 1:npointsE
                    S{n_dim}(idimGLO,kdim) =  S{n_dim}(idimGLO,kdim)  + X(idim,Jnode).*H_xi{kdim,n_dim}(:,Jnode) ;
                end
            end
        end
        R{n_dim} =zeros(size(detJ)) ;
        
        for  kdim = 1:ndim
            kdimGLO = (ipoints-1)*ndim+kdim ;
            for idim = 1:ndim
                idimGLO = (ipoints-1)*ndim+idim ;
                R{n_dim} =  R{n_dim} + J(kdimGLO,idim).*S{n_dim}(idimGLO,kdim);
            end
        end
        R{n_dim} =   R{n_dim}.*detJ ;
        
    end
    
    %%%
%     E = cell(1,ndim) ;
%     for n_dim = 1:ndim
%         E{n_dim} = zeros(size(J)) ; 
%         for hdim = 1:ndim 
%             hdimGLO = (ipoints-1)*ndim+hdim ;
%             for idim = 1:ndim 
%                 for kdim = 1:ndim 
%                      kdimGLO = (ipoints-1)*ndim+kdim ;
%                     for mdim  = 1:ndim 
%                         mdimGLO = (ipoints-1)*ndim+mdim ;
%                         E{n_dim}(hdimGLO,idim) =  E{n_dim}(hdimGLO,idim) ...
%                                                   -invJ(hdimGLO,kdim).*S{n_dim}(kdimGLO,mdim).*invJ(mdimGLO,idim); 
%                     end
%                 end
%             end
%         end 
%     end


%     
%     
%     T = cell(ndim,ndim) ; 
%     
%         for  n_dim = 1:ndim
%             for idim = 1:ndim
%                 T{n_dim,idim} = zeros(size(B{n_dim}));
%                 for Inode = 1:npointsE
%                     for hdim = 1:ndim
%                         hdimGLO = (ipoints-1)*ndim+hdim ;
%                         T{n_dim,idim}(:,Inode) = T{n_dim,idim}(:,Inode)  + H_xi{n_dim,hdim}(:,Inode).*invJ(hdimGLO,idim) + B_xi{hdim}(:,Inode).*E{n_dim}(hdimGLO,idim) ;
%                     end
%                 end
%             end
%         end
    
    
    
end


% Functions to be integrated
% --------------------------
n_funct =  (npointsE)^2;
A = zeros(npoints,n_funct) ;
kacum = 1;

if EVALUATE_GRADIENT == 1
    dA{1} = A;
    dA{2} = A;
else
    dA{1} = [] ; dA{2} = [];
end

%for idim = 1:ndim
    for Inode = 1:npointsE
       % for jdim = 1:ndim
            for Jnode = 1:npointsE
                A(:,kacum) = Nshape(:,Inode).*Nshape(:,Jnode).*detJ ;
                
                if  EVALUATE_GRADIENT == 1
                    
                    for n_dim = 1:ndim
                        dA{n_dim}(:,kacum) =   Nshape(:,Inode).*Nshape(:,Jnode).*R{n_dim} ; % G_C
                    end
                    
                end
                
                
                kacum = kacum + 1;
            end
       % end
    end
%end

dA_x  = dA{1} ; 
dA_y = dA{2} ;




%
%
%
% alpha = DATA.IRREGULAR_SHAPE_PARAMETRIZATION_ALPHA ;
%
% % COordinates irregular element
%
% X= [-1, 1, +1-alpha, -1
%     -1, -1, 1, 1  ] ;
%
% N =@(xi,eta) (0.25*[(1-xi).*(1-eta), (1+xi).*(1-eta), (1+xi).*(1+eta), (1-xi).*(1+eta) ]);
%
% DATAOUT.Nshape = N ;
% DATAOUT.Xcoor =X ;
%
%
% B_xi = @(xi,eta)  (0.25*[ -(1-eta),  (1-eta),  (1+eta), -(1+eta) ]) ;
% B_eta = @(xi,eta)  (0.25*[ -(1-xi) , -(1+xi) , (1+xi), (1-xi) ]) ;
%
% B = {B_xi(xMAT(:,1),xMAT(:,2))
%     B_eta(xMAT(:,1),xMAT(:,2))} ;
%
%
%

%
% DATAOUT.detJ = detJ ;
%
% % Therefore, the function we wish to integrate is given by A = f_L*detJ
%
%
% [f, df_x,df_y] =    LagrangePolynomial2D(xLIM,P,xMAT) ;
% A = bsxfun(@times,f,detJ) ;
%
%
%
% % The derivative, on the other hand, is given by
% % dA_x = df_x*detJ + f*der{(detJ)}{x}
% % dA_y = df_y*detJ + f*der{(detJ)}{y}
% dA_x = bsxfun(@times,df_x,detJ) ;
% dA_y = bsxfun(@times,df_y,detJ) ;
%
% % Now we have to sum the contribution (See comments in
% % LagrangePolynomial2D_irregular_aux.mlx) C_x = f*der{(detJ)}{x}, C_y= f*der{(detJ)}{y}
% %-----------------------------
% % der{(detJ)}{x} = detJ*TRACE(inv(J)*der{J}{x})
% % der{(detJ)}{y} = detJ*TRACE(inv(J)*der{J}{y})
%
% % Let us calculate inv(J)
%   invJ= inverseVECTORIZED(J,ndim,detJ) ;
%
% % -----------------------------------------------------------------------------------
% % Derivative J with respect to x (\xi_1, or \xi, in LagrangePolynomial2D_irregular_aux.mlx)
% %---------------------------------------------------------------------------------
%
% % B_xi = @(xi,eta)  (0.25*[ -(1-eta),  (1-eta),  (1+eta), -(1+eta) ]) ;
% der_B_xi_xi =  [ 0,  0,  0, 0 ]  ;
% % B_eta = @(xi,eta)  (0.25*[ -(1-xi) , -(1+xi) , (1+xi), (1-xi) ]) ;
% der_B_eta_xi =  0.25*[ +1 ,-1 , +1, -1  ] ;
%
% derB_x = { repmat(der_B_xi_xi,npoints,1),repmat(der_B_eta_xi,npoints,1) } ;
%
% derJ_x =zeros(ndim*npoints,ndim) ;
% ipoints = 1:npoints ;
% for  idim = 1:ndim
%     idimLOC = (ipoints-1)*ndim+idim ;
%     for jpoints = 1:npointsE
%         for kdim = 1:ndim
%             derJ_x(idimLOC,kdim) = derJ_x(idimLOC,kdim) + X(idim,jpoints)*derB_x{kdim}(ipoints,jpoints) ;
%         end
%     end
% end
%
%   % -----------------------------------------------------------------------------------
% % Derivative J with respect to y (\xi_2, or \eta in LagrangePolynomial2D_irregular_aux.mlx)
% %---------------------------------------------------------------------------------
%
% % B_xi = @(xi,eta)  (0.25*[ -(1-eta),  (1-eta),  (1+eta), -(1+eta) ]) ;
% der_B_xi_eta =  0.25*[ +1,  -1,  +1, -1]  ;
% % B_eta = @(xi,eta)  (0.25*[ -(1-xi) , -(1+xi) , (1+xi), (1-xi) ]) ;
% der_B_eta_eta =  [ 0 ,0 , 0,0 ]  ;
%
% derB_y = { repmat(der_B_xi_eta,npoints,1),repmat(der_B_eta_eta,npoints,1) } ;
%
% derJ_y =zeros(ndim*npoints,ndim) ;
% ipoints = 1:npoints ;
% for  idim = 1:ndim
%     idimLOC = (ipoints-1)*ndim+idim ;
%     for jpoints = 1:npointsE
%         for kdim = 1:ndim
%             derJ_y(idimLOC,kdim) = derJ_y(idimLOC,kdim) + X(idim,jpoints)*derB_y{kdim}(ipoints,jpoints) ;
%         end
%     end
% end
%
% derJ_x = invJ.*derJ_x ;
% derJ_y = invJ.*derJ_y ;
%
% if ndim == 2
%     TRACE_x = derJ_x(1:ndim:end,1) + derJ_x(2:ndim:end,2) ;
%     TRACE_y = derJ_y(1:ndim:end,1) + derJ_y(2:ndim:end,2) ;
%     % C_x = detJ*TRACE(inv(J)*der{J}{x})
% % C_y = detJ*TRACE(inv(J)*der{J}{y})
% TRACE_x = detJ.*TRACE_x ;
% TRACE_y = detJ.*TRACE_y ;
%
%
% else
%     error('option not implemented')
% end
%
% dA_x = dA_x + f.*TRACE_x ;
% dA_y = dA_y + f.*TRACE_y ;
%
% % INTEGRATE_VOLUME = 0 ;
% %  % In order to integrate exactly the surface of the reference element
% %  if INTEGRATE_VOLUME == 1
% % A = [ones(npoints,1),A] ;
% % dA_x = [zeros(npoints,1),dA_x] ;
% % dA_y= [zeros(npoints,1),dA_y] ;
% %  end
%
%
% % C_x = detJ*TRACE(inv(J)*der{J}{x})
% % C_y = detJ*TRACE(inv(J)*der{J}{y})
%
%
% %end
%
% %
% % xi = xMAT(:,1) ;
% % eta = xMAT(:,2) ;
% %
% % B_xi = 0.25*[ -(1-eta),  (1-eta),  (1+eta), -(1+eta)] ;
% % B_eta =    0.25*[   -(1-xi) , -(1+xi) , (1+xi), (1-xi)   ] ;
% %
% % X_x = X(1,:) ;
% % X_y = X(2,:) ;
% %
% % IND_local = 1:ndim:npoints;
% % npointsELEM = size(X,2) ;
% %
% % for ipoints = 1:npointsELEM
% % end
%
%
%
%
% %     iini = (ipoints-1)*ndim+1;
% %     ifin = ipoints*ndim ;
% %
% %     J(iini:ifin,1) =
%
%
%
%
% % Derivatives  shape functions at points xMAT  (matrix B)
%
% %  [Ne BeXi] = Quadrilateral4N(xiV) ;
% % % Shape functions and derivatives for 4-node quadrilateral
% % xi = xiV(1) ; eta = xiV(2) ;
% % % Matrix of shape functions
% % Ne =0.25*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta) ];
% % % Matrix of the gradient of shape functions
% % BeXi = 0.25*[ -(1-eta),  (1-eta),  (1+eta), -(1+eta) ;
% %              -(1-xi) , -(1+xi) , (1+xi), (1-xi)   ] ;
%
%
