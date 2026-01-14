function Me = ComputeMe_EIFEnv(EIFEoper,NeALL,weig) ;

if nargin == 0
    load('tmp.mat')
end

Me = zeros(size(NeALL,2)) ;
ngaus  = length(weig) ; 
ndim = size(EIFEoper.MESH.COOR,2) ; 
for  g = 1:ngaus
    % Jacobian matrix at point g  (transformation of coordinates )
    % Matrix of derivatives for Gauss point "g" (polynomial shape functions)
    %  BeXi = dershapef(:,:,g) ;
%     BeXi_transf = zeros(ndim,nnodeE) ;
%     for idim = 1:ndim
%         BeXi_transf(idim,:) =  EIFEoper.BodyForces.Bgrad_transf{idim}(g,:) ;
%     end
%     % Jacobian Matrix
%     Je = Xe*BeXi_transf' ;
%     % JAcobian
%     detJe = det(Je) ;
%     if g == 1
%         detJe_ref = detJe;
%         % Rotation V matrix (modes interfaces)
%         %    Vall = EIFEoper.MODES.Vall ;
% 
% %         Vall_rot = zeros(size(EIFEoper.MODES.Vall)) ;
% %         for imodes = 1:size(EIFEoper.MODES.Vall,2)
% %             LOCV = (Je'\reshape(EIFEoper.MODES.Vall(:,imodes),ndim,[])) ;
% %             Vrot(:,imodes) = LOCV(:) ;
% %         end
%     else
%         if abs(detJe-detJe_ref)/detJe_ref >= 1e-10
%             error('Only constant Jacobian transformations allowed')
%         end
%     end
    
%     Ginv_PhiRBfT_m_HrbHdefINV_PsiDEFfT = EIFEoper.OPER.Ginv_PhiRBfT_m_HrbHdefINV_PsiDEFfT  ;
%         HdefINV_PsiDEFfT = EIFEoper.OPER.HdefINV_PsiDEFfT  ;
%         
%         

    % Rotation/deformation PhiDEFe
    %------------------------------
    
    
    idof = small2large(g,ndim ) ;
    
%     NmatDEFred = EIFEoper.BodyForces.NmatDEFred(idof,:);
%     NmatRBred = EIFEoper.BodyForces.NmatRBred(idof,:);
%     
%     Vrot = EIFEoper.MODES.Vall; 
    
   % Ne = NmatDEFred*HdefINV_PsiDEFfT*Vrot  +  NmatRBred*Ginv_PhiRBfT_m_HrbHdefINV_PsiDEFfT*Vrot ;
   Ne = NeALL(idof,:); 
    dens = EIFEoper.BodyForces.MATPRO.dens(g,:) ;
    Me = Me + weig(g)*(Ne'*dens*Ne) ;
end