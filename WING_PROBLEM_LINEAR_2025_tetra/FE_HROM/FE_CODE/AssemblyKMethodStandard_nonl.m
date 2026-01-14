function   Kstiff = AssemblyKMethodStandard_nonl(ASSEMBLY_INFO,Cglo) ; 
%%  Assembly  method in which the stiffness matrix is computed by first determining
% the product B'(e,xi)*C(e,xi)*B(e,xi), and then performing the assembly in
% the standard, vectorized fashion
% Joaquín A. Hernández (jhortega@cimne.upc.edu), 28-Oct-2015
% Adaption: 24-june-2019
% ---------------------------------------------------------------------------
%dbstop('9')
%error('Method not IMPLEMENTED')
if nargin == 0
  load('tmp2.mat')
end


















% 
% Belem = ASSEMBLY_INFO.Belem ;      % nstrain*ngaus*nelem* x ndof 
% nelem = size(ASSEMBLY_INFO.CN,1) ; 
% ngaus = ASSEMBLY_INFO.ngaus ;
% nstrain = ASSEMBLY_INFO.nstrain ;
% ndofT = size(Belem,2) ; 
% % -----------------------------------------
% % Reshaping Belem 
% BelemR = reshape(Belem,nstrain*ngaus,ndofT,[]) ; 
% BelemR = reshape(BelemR,nstrain,ndofT,ngaus,nelem) ;   % nstrain2 x ndofT x ngaus x nelem 
% 
% % Reshaping Cglo
% CgloR = reshape(Cglo,nstrain,nstrain,ngaus, nelem) ;  % nstrain1 x nstrain2 x ngaus x nelem
% 
% % Matrix CgloR*BelemR   % ---- > nstrain1 x ndofT x ngaus x nelem  
% % -------------------
% % BelemRperm = permute(BelemR,[5,2,3,4,1])    ; %  1 x ndofT x ngaus x nelem x nstrain2 
% % CgloRperm = permute(CgloR,[1,5,3,4,2])    ;      %  nstrain1 x 1 x ngaus x nelem x nstrain2
	
% % BC  = sum(BCext,5) ;                              % nstrain1 x ndofT x ngaus x nelem  
% TIME.BC= tic ; 
% BC = sum(bsxfun(@times, permute(CgloR,[1,5,3,4,2]), permute(BelemR,[5,2,3,4,1])) ,5) ; 
% TIME.BC=  toc(TIME.BC) ; 
% disp(['Time  BC product = ',  num2str(TIME.BC)]) ; 
% % -------------------------
% % Matrix B^T*C*B 
% % ---------------
% BelemTperm = permute(BelemR,[2,5,3,4,1])    ; %  ndofT x 1 x ngaus x nelem x nstrain  
% BCperm = permute(BC,[5,2,3,4,1])    ;         %  1     x ndofT x ngaus x nelem x nstrain 
% BtCBperm = bsxfun(@times,BelemTperm,BCperm) ; %  ndofT x ndofT x ngaus x nelem x nstrain 
% BtCB = sum(BtCBperm,5) ;                      %  ndofT x ndofT x ngaus x nelem  



% 
% 
% %ngaus = length(wSTs)/nelem ;
% % disp('Global elasticity matrix  ...')
% % % Matrix with all elasticity tensors (at all gauss points)
% % Cglo = DefineElastMatGLO(celasglo,ngaus)  ; % Cglo =[C(e=1,g=1);C(e=1,g=2) ; ...
% %%% 
%  
% ndim = ASSEMBLY_INFO.ndim ; 
% nnode =length(unique(ASSEMBLY_INFO.CN)) ; 
% 
% nelem = size(ASSEMBLY_INFO.CN,1) ; 
% nnodeE = size(ASSEMBLY_INFO.CN,2) ; 
% 
% B = ASSEMBLY_INFO.Belem ; 
% 
% 
% disp('Point-wise product  C*B')
% tic
% CB = ProducMatrBlock(Cglo,ASSEMBLY_INFO.Belem) ;
% toc
% % disp('Matrix W*(B)^T') 
% % % Product of the vector of weights times Belem
% % BelemT = bsxfun(@times,ASSEMBLY_INFO.Belem,ASSEMBLY_INFO.wST) ; 
% % Transponse
% nelem = ASSEMBLY_INFO.ngaus; 
% ngaus = 1; 
% BelemT =   TransposeVectorizeGen(ASSEMBLY_INFO.Belem,ngaus*nelem)  ; 
% disp('Point-wise product ...  Bt*C*B')
% BCB = ProducMatrBlockGen(BelemT,CB,ngaus*nelem) ;
% % 
% disp('Matrix of element stiffness matrices')
% Kelem = ElementKMatrices(nelem,nnodeE,ndim,BCB,ngaus) ; 
% disp('Assembly of K ...')
%  
% K = AssemblyKGlobal(Kelem,nelem,nnodeE,ndim,CN,nnode) ;
%  
%  