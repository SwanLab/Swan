function [ StifMat] = stiff_unitary_triang( dim,...
    element,problembsc,coordinatesn,coordinatesa,nstre,nunkn,nndof)

nelem=dim.nelem;  nnode=dim.nnode;
npnod=dim.npnod; ndime = dim.ndime; 
dirichlet_data = zeros(nnode,nelem);
idx   = zeros(nnode*nunkn,nelem);
ptype = problembsc.problemtype;
ftype = 'THERMAL';
StifMat = sparse(nndof,nndof);

for i=1:nnode
    dirichlet_data(i,:)= element.conectivities(:,i);
end

etype = element.type;
for a=1:nnode
    for i=1:nunkn
        idx(nunkn*a-nunkn+i,:) = nunkn.*dirichlet_data(a,:)-nunkn+i;
    end
end
[posgp,weigp,ngaus] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);
estiff = zeros(nnode*nunkn,nnode*nunkn,nelem);


for igaus=1:ngaus
    [cartd,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,ptype);
    [Bmat] = cal_B(cartd,nstre,nnode,nunkn,nelem,etype,ptype,ftype);
    dvolu = weigp(igaus)*djacb;
 
            for iv=1:nnode*nunkn
                for jv=1:nnode*nunkn
                    for istre=1:nstre
                        %for jstre=1:nstre
                            v = squeeze(Bmat(istre,iv,:).*Bmat(istre,jv,:));
                            estiff(iv,jv,:)=squeeze(estiff(iv,jv,:)) + v(:).*dvolu(:);
                       %M end
                    end
                end
            end

end

% assembling of elemental matrix
for k=1:nnode*nunkn
    for l=1:nnode*nunkn
        vestiff = squeeze(estiff(k,l,:));
        StifMat = StifMat + sparse(idx(k,:),idx(l,:),vestiff,nndof,nndof);
    end
end


end

