function [BBnwLOC detJ] = BmatQuad_locBBAR(X,dN_dXi,weightsE)

 if nargin ==0
    %     X = [-1 -1;
    %         1 -1;
    %         1  1
    %         -1  1];
    %     [N dN_dXi] = QuadElem_N_derN  ;
    load('tmp5.mat') ;
end

[dN_dx detJ] = QuadElem_derNX(X,dN_dXi) ;
ndime = size(X,2);
ngaus = size(dN_dx,1)/ndime ;
nstrain = 4 ;
nnode = size(dN_dXi,2) ;

BBnwLOC  = zeros(nstrain*ngaus,nnode*ndime);


%%%% Modified dilatational part
vol = sum(detJ.*weightsE) ;
c1 = 1:ndime:ngaus*ndime ; c2 = 2:ndime:ngaus*ndime ;
B1 = dN_dx(c1,:) ;
B2 = dN_dx(c2,:) ;
Bbar1 = bsxfun(@times,B1,(detJ.*weightsE)') ;
Bbar1 = sum(Bbar1,1)/vol ;

Bbar2 = bsxfun(@times,B2,(detJ.*weightsE)') ;
Bbar2 = sum(Bbar2,1)/vol ;
%%%%



for g = 1:ngaus
    iniROW = (g-1)*nstrain+1;     finROW = g*nstrain ;
    iniROWN = (g-1)*ndime+1;     finROWN = g*ndime ;
    for a = 1:nnode
        iniCOL = (a-1)*ndime+1;         finCOL = a*ndime  ;
        dN_dxL = dN_dx(iniROWN:finROWN,a)  ;
%         Bloc = [dN_dxL(1)  0 ;
%             0     dN_dxL(2)
%             dN_dxL(2) dN_dxL(1)
%             0   0];
        Bloc =   BBarMATRIXquad(dN_dxL(1),dN_dxL(2),Bbar1(a),Bbar2(a));
        BBnwLOC(iniROW:finROW,iniCOL:finCOL) = Bloc ;
    end
end



