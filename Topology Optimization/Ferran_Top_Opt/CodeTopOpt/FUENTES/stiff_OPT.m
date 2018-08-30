function [ StifMat,post] = stiff_OPT( dim,element,problembsc,coordinatesn,coordinatesa,gamma)

nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;

dirichlet_data = zeros(nnode,nelem);
idx   = zeros(nnode*nunkn,nelem);
ptype = problembsc.problemtype;
ftype = problembsc.phisical_type;
StifMat = sparse(nndof,nndof);

for i=1:nnode
    dirichlet_data(i,:)= element.conectivities(:,i);
end

etype = element.type;
for inode=1:nnode
    for idime=1:nunkn
        idx(nunkn*inode-nunkn+idime,:) = nunkn.*dirichlet_data(inode,:)-nunkn+idime;
    end
end
[posgp,weigp,ngaus] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);
estiff = zeros(nnode*nunkn,nnode*nunkn,nelem);
post.chi = zeros(ngaus,nelem);

for igaus=1:ngaus
    [cartd,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,ptype);
    [Bmat] = cal_B(cartd,nstre,nnode,nunkn,nelem,etype,ptype,ftype);
    dvolu = weigp(igaus)*djacb;
    
    %[chi] = cal_caracteristic_function(igaus,phifunct,element,dim,problembsc,vol_void);
    %chi = 1-vol_void;
    %chi(abs(chi)<1e-12) =  element.material.opt_epsi;

    [Ce] = compute_consitutive_law(element,problembsc,igaus,dim,gamma);

    switch element.material.subtype
        case {'PLANESTRAIN','PLANESTRES'}
            for iv=1:nnode*nunkn
                for jv=1:nnode*nunkn
                    for istre=1:nstre
                        for jstre=1:nstre
                            v = squeeze(Bmat(istre,iv,:).*Ce(istre,jstre,:).*Bmat(jstre,jv,:));
                            estiff(iv,jv,:)=squeeze(estiff(iv,jv,:)) + v(:).*dvolu(:);
                        end
                    end
                end
            end
        case 'AXI'
    end
    post.chi(igaus,:) = gamma;
end
% assembling of elemental matrix
for k=1:nnode*nunkn
    for l=1:nnode*nunkn
        vestiff = squeeze(estiff(k,l,:));
        StifMat = StifMat + sparse(idx(k,:),idx(l,:),vestiff,nndof,nndof);
    end
end

end

