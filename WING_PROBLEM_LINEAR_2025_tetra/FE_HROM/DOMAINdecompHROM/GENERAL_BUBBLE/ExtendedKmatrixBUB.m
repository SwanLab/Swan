function K = ExtendedKmatrixBUB(Kbb,Kb0,K0b,K00,ndim,nBUB)
% See /home/joaquin/Desktop/CURRENT_TASKS/MATLAB_CODES/TESTING_PROBLEMS_FEHROM/104_EIFEM_plast2D/19_ExactLinearStiff.mlx
% JAHO, M6-March-2024, Balmes 185, Barcelona
if nargin == 0
    load('tmp1.mat')
    nBUB = 1;
    Kb0 = Kb0(:,1);
    K0b = K0b(1,:);
    K00 = K00(1,1) ;
end

if ndim ==nBUB
    % Nothing is to be done
    K = [Kbb,Kb0;K0b,K00 ] ;
elseif nBUB >  ndim
    % Number of bubble modes greater than the number spatial dimensions
    % Boundary interscale B-matrix is to be expanded with zeros
    nnodeE_b = size(Kbb,2)/ndim ;  % Number of boundary nodes
    nDOFStotB = nnodeE_b*nBUB;   % Total number of DOFs, including the "ghost" dofs
    Kbb_exp = zeros(nDOFStotB,nDOFStotB) ;
    Kb0_exp = zeros(nDOFStotB,size(Kb0,2)) ;
    K0b_exp = zeros(size(K0b,1),nDOFStotB ) ;
    for idim = 1:nnodeE_b
        iiniNEW = (idim-1)*nBUB + 1;
        ifinNEW = (idim-1)*nBUB + ndim ;
        iiniOLD = (idim-1)*ndim + 1;
        ifinOLD = (idim-1)*ndim + ndim ;
        Kb0_exp(iiniNEW:ifinNEW,:) = Kb0(iiniOLD:ifinOLD,:) ;
        for  jdim = 1:nnodeE_b
            jiniNEW = (jdim-1)*nBUB + 1;
            jfinNEW = (jdim-1)*nBUB + ndim ;
            jiniOLD = (jdim-1)*ndim + 1;
            jfinOLD = (jdim-1)*ndim + ndim ;
            Kbb_exp(iiniNEW:ifinNEW,jiniNEW:jfinNEW) = Kbb(iiniOLD:ifinOLD,jiniOLD:jfinOLD) ;
            if idim == 1
                K0b_exp(:,jiniNEW:jfinNEW) = K0b(:,jiniOLD:jfinOLD) ;
            end
        end
    end
    K = [Kbb_exp,Kb0_exp;K0b_exp,K00 ] ;
else
    % Number of bubble modes less than number of spatial dimensions
    %  nnodeE_b = size(BmatB,2)/ndim ;  % Number of boundary nodes
    if  nBUB == 0
        K = Kbb; 
    else
    diffZEROS = ndim-nBUB ;
    nbound = size(Kbb,1) ;
    cZ = zeros(nbound,diffZEROS) ;
    Kb0 = [Kb0,cZ]  ;
    K0b =  [K0b;cZ'] ;
    
    cZ = zeros(nBUB,diffZEROS) ;
    K00 =[K00,cZ];
    cZ = zeros(nBUB,size(K00,2)) ;
    K00 =[K00;cZ];
    
    K = [Kbb ,Kb0 ;K0b ,K00 ] ;
    end
    
end