function Lbool = Lbool_vectorized(CN,nnode,ngaus,ndim)
% %% GOAL:  Given d (global displacements), and U  (displacements at the Gauss points of  elements CN),
% we wish to contruct two matrices Nall and Lbool such that 
% U = diag(Nall)*Lbool*d  
% Here we focus on Lbool
% All elements are assumed to have the same number of nodes and Gauss
% points
% Non-vectorized version in Lbool_nonvectorized.m 
% 
% JAHO, 29-Jun-2021, Cartagena
% --------------------------------------

if nargin == 0
    % See FLOAT_IMPLE.pdf, KW:Lbool_nonvectorized_scalar
    CN =[1 2 5 6
        2 3 5 4] ;
    
    nnode = 6 ;
    ngaus = 2 ;
    ndim = 3;
    
end


nelem = size(CN,1) ;
nnodeE = size(CN,2) ;




%%%%%%
m = nelem*ndim*ngaus*nnodeE ;  % Number of rows 
n = nnode*ndim ;  % Number of columns 

nzmaxLOC = m*nnodeE*ndim ;   % Maximum number of zeros 

Lbool =  sparse([],[],[],m,n,nzmaxLOC); % Allocating memory for Lbool

e = 1:nelem  ;
 s = ones(size(e)) ; 
% for  igaus = 1:ngaus
%     igausGLO =  (e-1)*ngaus + igaus ; %   GAUSSglo_e(igaus)  ;
%     
%     for innode = 1:nnodeE
%         NODELOC = CN(e,innode) ;
%         
%         for idim = 1:ndim
%             
%             irow = (igausGLO-1)*ndim + idim ;
%             for jdim =1:ndim
%                 
%                 jrow = (NODELOC-1)*ndim + jdim ;
%                 if idim == jdim
%                    
%                     Lbool   =  Lbool  + sparse(irow,jrow,s,m,n,m);
%                 end
%             end
%         end
%     end
% end



    for igausLOC = 1:ngaus
        igausGLO =  (e-1)*ngaus + igausLOC ;
        for inodeLOC = 1:nnodeE
            inodeGLO = CN(e,inodeLOC) ;
            irowGLO = (igausGLO-1)*nnodeE + inodeLOC ;
            for idofLOC = 1:ndim
                idofGLO = (irowGLO-1)*ndim +  idofLOC ;
                for jdofLOC =1:ndim
                    jdofGLO =  (inodeGLO-1)*ndim +jdofLOC ;
                    if idofLOC == jdofLOC
                       % Lbool_vector_gauss(idofGLO,jdofGLO) = 1;
                        
                         Lbool   =  Lbool  + sparse(idofGLO,jdofGLO,s,m,n,m);
                    end
                end
            end
        end
    end


 
  

 

