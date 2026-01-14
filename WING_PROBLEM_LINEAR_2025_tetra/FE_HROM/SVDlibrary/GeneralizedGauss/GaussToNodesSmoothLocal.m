function VARNODES= GaussToNodesSmoothLocal(VARGAUSS,COOR,CN,...
        TypeElement,posgp,NameFileMesh,DATAINloc,W,SMOOTH_OPERATOR) 
  %  TRANSPORTING  VARIABLES FROM GAUSS POINTS TO NODES 
  % JAHO, 17-Apr-2020  (35th of confinement because of COVID19)
 % ---------------------------------------------------------------- 
 if nargin == 0
     load('tmp1.mat')
     SMOOTH_OPERATOR = [] ; 
 end
 
 %
 nstrain = DATAINloc.nstrain ; 
 nmodes = size(VARGAUSS,2) ;
 nelem = size(CN,1) ; nnodeE = size(CN,2) ; ndim = size(COOR,2) ; nnode = size(COOR,1) ; 
 TypeIntegrand = 'K'; 
 
 [weig,posgpOUT,N,dershapef] = ComputeElementShapeFun(TypeElement,nnodeE,TypeIntegrand) ; 
 ngausE = size(posgpOUT,2) ; 
 
 if  norm(posgp-posgpOUT) ~= 0
     error('Input and output GAuss points position should coincide')
 end
 
 % Notice that u = N*d --> Given a variable defined at nodes (d), N*d
 % returns the values the m gauss points.  
 if rank(N)~= min(size(N))
     error('Option valid only when N is full rank') 
 end
 
 % CONSTRUCTION OF SMOOTHING OPERATOR 
 % ----------------------------------
 if ~isempty(SMOOTH_OPERATOR)
     % BOOLEAN MATRIX NstBOOL --> Similar to Nst, but assuming that 
 end
 
 
 
 
 VARNODES  =zeros(nnode*nstrain,nmodes) ; 
 
 for  imode  = 1:nmodes 
     % Loop over number of nodes      
     for istrain = 1:nstrain
              % Loop over number of strain/stress components 
         VARGLOC = VARGAUSS(istrain:nstrain:end,imode) ; 
         % VARGLOC is a scalar variable defined at Gauss points
         VARGLOC = reshape(VARGLOC,ngausE,[]) ;   % ngausE x nelem 
         % Equivalent nodal values 
         VARNLOC = N\VARGLOC;   % nnode x NELEM  
         VARNODES(istrain:nstrain:end,imode) = SMOOTH_OPERATOR*VARNLOC  ; 
         
         
     end
     
 end
 
 
 
    

 