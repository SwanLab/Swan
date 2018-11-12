function [ fextstra] = fext_strain(dim,element,problembsc,...
    coordinatesn,coordinatesa,gamma,stra0)

% basic dimensions
nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;

dirichlet_data = zeros(nnode,nelem);
idx   = zeros(nnode*nunkn,nelem);
ptype = problembsc.problemtype;
ftype = problembsc.phisical_type;
fextstra = zeros(nndof,1);

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
eforce = zeros(nnode*nunkn,nelem);

for igaus=1:ngaus
    [cartd,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,ptype);
    [Bmat] = cal_B(cartd,nstre,nnode,nunkn,nelem,etype,ptype,ftype);
    dvolu = weigp(igaus)*djacb;

    %[chi] = cal_caracteristic_function(igaus,gamma,element,dim,problembsc,vol_void);
    [Ce] = compute_consitutive_law(element,problembsc,igaus,dim,gamma);
        
    switch element.material.subtype
        case {'PLANESTRAIN','PLANESTRES'}
            
            stre0=zeros(nstre,nelem);
            for istre=1:nstre
                for jstre=1:nstre
                    stre0(istre,:) = stre0(istre,:) + squeeze(Ce(istre,jstre,:)*stra0(jstre))';
                end
            end
            for iv=1:nnode*nunkn
                ivBmat(1:nstre,1:nelem) = Bmat(1:nstre,iv,1:nelem);
                for istre=1:nstre
                    eforce(iv,:)=eforce(iv,:)+ivBmat(istre,:).*stre0(istre,:).*dvolu';
                end
            end
            
        case 'AXI'
            
    end
end

% assembling of elemental matrix
for k=1:nunkn*nnode
    fextstra = fextstra + sparse(idx(k,:),1,eforce(k,:),nndof,1);
end
fextstra = -fextstra;


end

