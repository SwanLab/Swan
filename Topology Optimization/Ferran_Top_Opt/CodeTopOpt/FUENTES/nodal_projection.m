function fn = nodal_projection(phi,element,problembsc,dim,coordinates)

ndime=dim.ndime; nelem=dim.nelem; nnode=dim.nnode; nunkn = 1; npnod= dim.npnod;
vol_mat = cal_vol_mat(phi,dim,element,problembsc,coordinates);

ptype=problembsc.problemtype;
[posgp,weigp,ngaus] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);

neres=0; nnode=dim.nnode; 

lnods=zeros(nnode,nelem);
for inode=1:nnode
    lnods(inode,:)=element.conectivities(:,inode);
end

idx   = zeros(nnode*nunkn,nelem);
for inode=1:nnode
    for idime=1:nunkn
        idx(nunkn*inode-nunkn+idime,:) = nunkn.*lnods(inode,:)-nunkn+idime;
    end
end

f = zeros(nnode*nunkn,nelem);
for igaus=1:ngaus
    [shape,~,~] = shape_deriv_functions(igaus,posgp,ptype,element.type,nnode,neres);
    [cartd,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinates,coordinates,ptype);
    dvolu = weigp(igaus)*djacb;
    
    for iv=1:nnode*nunkn
           f(iv,:)=f(iv,:)+ shape(inode)*vol_mat.*dvolu';
    end

end
fn = zeros(dim.npnod,1);
for k=1:nunkn*nnode
    fn = fn + sparse(idx(k,:),1,f(k,:),npnod,1);
end


end