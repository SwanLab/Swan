function Lbool_vector_gauss = Lbool_nonvectorized(CN,nnode,ngaus,ndim)

if nargin == 0
    % See FLOAT_IMPLE.pdf, KW:Lbool_nonvectorized_scalar
    CN =[1 2 4 3
         3 4 5 6] ; 
    
    nnode = 6 ;
    ngaus = 2 ; 
    ndim = 1; 
    
end


nelem = size(CN,1) ; 
nnodeE = size(CN,2) ; 
nrows = nelem ; 
ncols = nnode ; 

Lbool_scalar= zeros(nelem*nnodeE,ncols) ;  

for e = 1:nelem    
    for inodeLOC = 1:nnodeE
        inodeGLO = CN(e,inodeLOC) ; 
        irowGLO = (e-1)*nnodeE + inodeLOC ; 
        Lbool_scalar(irowGLO,inodeGLO) = 1;
    end
end

%%%%%% 
nrows = nelem*ndim*nnodeE ; 
ncols = nnode*ndim ; 
Lbool_vector= zeros(nrows,ncols) ;  

for e = 1:nelem    
     for inodeLOC = 1:nnodeE
        inodeGLO = CN(e,inodeLOC) ;   
        irowGLO = (e-1)*nnodeE + inodeLOC ; 
        for idofLOC = 1:ndim
            idofGLO = (irowGLO-1)*ndim +  idofLOC ; 
            for jdofLOC =1:ndim
                jdofGLO =  (inodeGLO-1)*ndim +jdofLOC ; 
                if idofLOC == jdofLOC
                Lbool_vector(idofGLO,jdofGLO) = 1;
                end
            end
        end
    end
end

%%% ngaus 


%%%%%% 
nrows = nelem*ndim*ngaus*nnodeE ; 
ncols = nnode*ndim ; 
Lbool_vector_gauss= zeros(nrows,ncols) ;


for e = 1:nelem
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
                        Lbool_vector_gauss(idofGLO,jdofGLO) = 1;
                    end
                end
            end
        end
    end
end









disp('')









% for e = 1:nelem
%     %GAUSSglo_e = small2large(e,ngaus) ;
%     for  igaus = 1:ngaus
%         igausGLO =  (e-1)*ngaus + igaus ; %   GAUSSglo_e(igaus)  ;
%         %   ROWLOC = small2large(igausGLO,ndim) ;
%         for innode = 1:nnodeE
%             NODELOC = CN(e,innode) ;
%             %   COLLOC = small2large(NODELOC,ndim) ;
%             for idim = 1:ndim
%                 %  irow = ROWLOC(idim) ;
%                 irow = (igausGLO-1)*ndim + idim ;
%                 for jdim =1:ndim
%                     %  jrow =  COLLOC(jdim) ;
%                     jrow = (NODELOC-1)*ndim + jdim ;
%                     if idim == jdim
%                         Lbool_vector_gauss(irow,jrow) = 1;
%                     end
%                 end
%             end
%         end
%     end
% end


  

 

