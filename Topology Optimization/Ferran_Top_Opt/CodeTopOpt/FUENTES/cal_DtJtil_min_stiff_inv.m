function [ DtJtil ] = cal_DtJtil_min_stiff_inv(smoothDtC,matCh,DtC,dim,ngaus,weights,hnorm,Msmooth,coordinatesn,coordinatesa,element,problembsc,phifunct,phigp)
                                                
%strain0 = alpha;
%weights = alpha'*alpha;
npnod=dim.npnod; ndime=dim.ndime; nelem=dim.nelem; nnode=dim.nnode; 
nunkn = dim.nunkn; nstre = dim.nstre;

weights_inv = inv(matCh)*weights*inv(matCh);
DtJtil = [];


[DtC] = smooth_DtC(DtC,nstre,npnod,ngaus,nelem,ndime,nnode,Msmooth,element,coordinatesn,coordinatesa,problembsc);
DtC1 = zeros(npnod,1);
DtJtil = zeros(npnod,1);
for a=1:nstre
    for b=1:nstre
        DtC1(:) = DtC(a,b,:);
        DtJtil = DtJtil - weights_inv(a,b)*DtC1(:)/abs(hnorm);
    end
end

% if (smoothDtC==1)
%     [DtC] = smooth_DtC(DtC,nstre,npnod,ngaus,nelem,ndime,nnode,Msmooth,element,coordinatesn,coordinatesa,problembsc);
%     DtC1 = zeros(npnod,1);
%     DtJtil = zeros(npnod,1);
%     for a=1:3
%         for b=1:3
%         DtC1(:) = DtC(a,b,:);
%         DtJtil = DtJtil - weights_inv(a,b)*DtC1(:)/abs(hnorm);
%         end
%     end
%     
%     
% 
%     
% elseif (smoothDtC==0)
%     DtC1 = zeros(ngaus,nelem);
%     DtJtil = zeros(ngaus,nelem);
%     for igaus=1:ngaus
%         for a=1:3
%             for b=1:3
%                 DtC1(igaus,:) = squeeze(DtC(a,b,igaus,:));
%                 DtJtil(igaus,:) = DtJtil(igaus,:) - weights_inv(a,b)*DtC1(igaus,:)/abs(hnorm);
%             end
%         end
%     end
%         
%     DtJtil2 = zeros(size(DtJtil));
%     for igaus=1:ngaus
%         DtJtil2(igaus,:)= -DtJtil(igaus,:);
%         outside= (phigp > 0);
%         if any(outside)
%             DtJtil2(igaus,outside)= DtJtil(igaus,outside);
%         end
%     end
%     DtJtil = DtJtil2;
%     [DtJtil] = smooth(nelem,npnod,ndime,nnode,DtJtil,Msmooth,element,coordinatesn,coordinatesa,problembsc);
% end



end

