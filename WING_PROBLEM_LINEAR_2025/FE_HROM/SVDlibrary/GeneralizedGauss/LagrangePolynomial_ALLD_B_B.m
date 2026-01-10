function  [A, dA,DATAOUT ]=    LagrangePolynomial_ALLD_B_B(xLIM,P,xMAT,DATA)
% See comments in LagrangePolynomial2D_B_B_aux.mlx
% This is a generalization to any spatial dimension of LagrangePolynomial2D_B_B.m
if nargin == 0
    load('tmp.mat')
end
ndim = size(xMAT,2) ;
DATAOUT = [] ; % dA = cell(1,ndim); % = [] ; dA_y = [] ;
DATA = DefaultField(DATA,'EVALUATE_GRADIENT',0) ;
DATA = DefaultField(DATA,'EVALUATE_SHAPE_FUNCTIONS',0) ;

EVALUATE_GRADIENT = DATA.EVALUATE_GRADIENT ;

% 1) Nshape = M points x  nmonomials, i.e., number of shape functions
DATALOC.EVALUATE_GRADIENT = 1;
if ndim == 2
    
    [Nshape, B_x,B_y,DATAOUTlagrag] =    LagrangePolynomial2D(xLIM,P,xMAT,DATALOC) ;
    % 2) Matrix of derivatives  (with respect to xi1_1 and xi_2 )
    B_xi ={B_x ;B_y} ;  %
    % 3) Matrix of second derivatives (with respect to xi1_1 and xi_2 )
    H_xi = {DATAOUTlagrag.H_xx,DATAOUTlagrag.H_xy
        DATAOUTlagrag.H_xy,DATAOUTlagrag.H_yy} ;
elseif ndim == 3
    [Nshape, B_x,B_y,B_z,DATAOUTlagrag] =    LagrangePolynomial3D(xLIM,P,xMAT,DATALOC) ;
    B_xi ={B_x ;B_y; B_z} ;  %
    H_xi =     {DATAOUTlagrag.H_xx,DATAOUTlagrag.H_xy,DATAOUTlagrag.H_xz
        DATAOUTlagrag.H_xy,DATAOUTlagrag.H_yy,DATAOUTlagrag.H_yz
        DATAOUTlagrag.H_xz,DATAOUTlagrag.H_yz,DATAOUTlagrag.H_zz} ;
else
    error('Option not implemented')
end

DATAOUT.Nshape = Nshape ;

% 3) Parameterization in terms of alpha
% --------------------------------------
if ndim == 2
    X = CoordElementsFromShapeParameters(DATAOUTlagrag,DATA) ;
elseif ndim==3
    X = CoordElementsFromShapeParameters3D(DATAOUTlagrag,DATA) ;
else
    error('Option not implemented')
end

if DATA.EVALUATE_SHAPE_FUNCTIONS == 1
    
    %      Xdistorted  = Nshape*X(1,:)' ;
    %      Ydistorted  = Nshape*X(2,:)' ;
    %      DATAOUT.COORdistorted = zeros(size(Nshape,1),ndim) ;
    
    for idim = 1:ndim
        DATAOUT.COORdistorted(:,idim) = Nshape*X(idim,:)' ;
    end
end

DATAOUT.Xcoor = X ;


% 4) Jacobian matrix
% ------------------------
%
% Jacobian matrix (vectorized)
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

% B matrix (with respect to x1 and x2)
% ------------------------------------
% See LagrangePolynomial2D_B_B_aux.mlx
B = cell(1,ndim) ;
B(:) = {zeros(size(B_xi{1}))} ;
ipoints = 1:npoints ;
for idim = 1:ndim
    for Inodes = 1:npointsE
        for hdim = 1:ndim
            hdimGLO = (ipoints-1)*ndim+hdim ;
            B{idim}(:,Inodes) = B{idim}(:,Inodes) + B_xi{hdim}(:,Inodes).*invJ(hdimGLO,idim) ;
        end
    end
end

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
                %  R{n_dim} =  R{n_dim} +
                %  J(kdimGLO,idim).*S{n_dim}(idimGLO,kdim);  % Before
                %  29-Oct-2021
                R{n_dim} =  R{n_dim} + invJ(kdimGLO,idim).*S{n_dim}(idimGLO,kdim);
                
            end
        end
        R{n_dim} =   R{n_dim}.*detJ ;
        
    end
    
    %%%
    E = cell(1,ndim) ;
    for n_dim = 1:ndim
        E{n_dim} = zeros(size(J)) ;
        for hdim = 1:ndim
            hdimGLO = (ipoints-1)*ndim+hdim ;
            for idim = 1:ndim
                for kdim = 1:ndim
                    kdimGLO = (ipoints-1)*ndim+kdim ;
                    for mdim  = 1:ndim
                        mdimGLO = (ipoints-1)*ndim+mdim ;
                        E{n_dim}(hdimGLO,idim) =  E{n_dim}(hdimGLO,idim) ...
                            -invJ(hdimGLO,kdim).*S{n_dim}(kdimGLO,mdim).*invJ(mdimGLO,idim);
                    end
                end
            end
        end
    end
    
    
    
    
    T = cell(ndim,ndim) ;
    
    for  n_dim = 1:ndim
        for idim = 1:ndim
            T{n_dim,idim} = zeros(size(B{n_dim}));
            for Inode = 1:npointsE
                for hdim = 1:ndim
                    hdimGLO = (ipoints-1)*ndim+hdim ;
                    T{n_dim,idim}(:,Inode) = T{n_dim,idim}(:,Inode)  + H_xi{n_dim,hdim}(:,Inode).*invJ(hdimGLO,idim) + B_xi{hdim}(:,Inode).*E{n_dim}(hdimGLO,idim) ;
                end
            end
        end
    end
    
    
    
end


% Functions to be integrated
% --------------------------
n_funct = ndim^2*(npointsE)^2;
dA = cell(1,ndim) ;


DATA = DefaultField(DATA,'PROCESS_MATRIX_SNAPSHOTS_BY_BLOCKS',0) ;
DATA = DefaultField(DATA,'NAMEWS_BLOCK_MATRICES',['DATAWS/tmp_mat']) ;
DATA = DefaultField(DATA,'PROCESS_MATRIX_SNAPSHOTS_BY_BLOCKS_STORE_IN_MEMORY',0) ;

if DATA.PROCESS_MATRIX_SNAPSHOTS_BY_BLOCKS_STORE_IN_MEMORY == 1
    error('Option not implemented correctly yet...')
end




if     DATA.PROCESS_MATRIX_SNAPSHOTS_BY_BLOCKS == 0
    kacum = 1;
    
    A = zeros(npoints,n_funct) ;
    if  EVALUATE_GRADIENT == 1
        dA(:) = {A} ;
    end
    
    for idim = 1:ndim
        for Inode = 1:npointsE
            for jdim = 1:ndim
                for Jnode = 1:npointsE
                    A(:,kacum) = B{idim}(:,Inode).*B{jdim}(:,Jnode).*detJ ;
                    if  EVALUATE_GRADIENT == 1
                        for n_dim = 1:ndim
                            dA{n_dim}(:,kacum) = T{n_dim,idim}(:,Inode).*B{jdim}(:,Jnode).*detJ ...   % G_A
                                + B{idim}(:,Inode).*T{n_dim,jdim}(:,Jnode).*detJ ... % G_B
                                + B{idim}(:,Inode).*B{jdim}(:,Jnode).*R{n_dim} ; % G_C
                        end
                    end
                    kacum = kacum + 1;
                end
            end
        end
    end
    
else
    % PARTITIONED BY BLOCKS
    % ---------------------
    A = cell(1, ndim^2) ;
    dA  = cell(ndim,ndim^2) ; 
    kacumBLOCK = 1;
    n_functBLOCK = (npointsE)^2;
    
    
    for idim = 1:ndim
        for jdim = 1:ndim
            % LOOP OVER SHAPE FUNCTIONS
            % ----------------------
            kacum = 1;
            Aloc = zeros(npoints,n_functBLOCK) ;
            dAloc = cell(ndim,1) ;
            if  EVALUATE_GRADIENT == 1
                dAloc(:) = {Aloc} ;
            end
            
            for Inode = 1:npointsE
                for Jnode = 1:npointsE
                    Aloc(:,kacum) = B{idim}(:,Inode).*B{jdim}(:,Jnode).*detJ ;
                    if  EVALUATE_GRADIENT == 1
                        for n_dim = 1:ndim
                            dAloc{n_dim}(:,kacum) = T{n_dim,idim}(:,Inode).*B{jdim}(:,Jnode).*detJ ...   % G_A
                                + B{idim}(:,Inode).*T{n_dim,jdim}(:,Jnode).*detJ ... % G_B
                                + B{idim}(:,Inode).*B{jdim}(:,Jnode).*R{n_dim} ; % G_C
                        end
                    end
                    
                    
                    kacum = kacum + 1;
                end
            end
            % END LOOP OVER SHAPE FUNCTIONS
            if DATA.PROCESS_MATRIX_SNAPSHOTS_BY_BLOCKS_STORE_IN_MEMORY == 1
                error('Option not implemented yet...')
                nameLOC = [DATA.NAMEWS_BLOCK_MATRICES,num2str(kacumBLOCK),'.mat'] ;
                disp(['Saving in memoring block idim = ',num2str(idim),' jdim =',num2str(jdim)])
                save(nameLOC,'Aloc') ;
                disp('...Done')
                A{kacumBLOCK} = nameLOC ;
            else
                A{kacumBLOCK} = Aloc ;
                for idimLOC = 1:ndim
                    dA{idimLOC,kacumBLOCK} = dAloc{idimLOC}  ; 
                end
            end
            kacumBLOCK = kacumBLOCK + 1 ;
        end
    end
    
    
    
    
end


 