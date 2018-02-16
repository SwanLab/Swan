function [force_ext_vol] = fext_vol(dim,...
    coordinatesn,coordinatesa,element,problembsc,coordinates)


nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;

lnods = zeros(nnode,nelem);
idx   = zeros(nnode*nunkn,nelem);
dvolu = zeros(1,nelem);
force_ext_vol  = zeros(nndof,1);
ptype = problembsc.problemtype;
ftype = problembsc.phisical_type;
etype = element.type;
neres = 0; 

for idime=1:nnode
    lnods(idime,:)= element.conectivities(:,idime);
end


[posgp,weigp,ngaus] = cal_posgp_weigp(element.type,ndime,nnode);
coord_pg = zeros(ndime,ngaus,nelem);

for inode=1:nnode
    for idime=1:nunkn
        idx(nunkn*inode-nunkn+idime,:) = nunkn.*lnods(inode,:)-nunkn+idime;
    end
end


for idime = 1:ndime
    coord_pg(idime,:,:) = interpol(coordinatesn(:,idime),element,dim,problembsc,coordinates);
end

fvol = force_volumetric(coord_pg,ndime,ftype);

for igaus=1:ngaus
    fext_vol = zeros(nnode*nunkn,nelem);
    [cartd,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,ptype);
    [shape,deriv,heslo] = shape_deriv_functions(igaus,posgp,ptype,etype,nnode,neres);
    dvolu = weigp(igaus)*djacb;
    switch element.material.subtype
        case {'PLANESTRAIN','PLANESTRES'}
                for inode=1:nnode
                    for iknown=1:nunkn
                        fext_vol(nunkn*inode-nunkn+iknown,:) =  fext_vol(nunkn*inode-nunkn+iknown,:) + shape(inode)*squeeze(fvol(iknown,igaus,:))'.*dvolu'; %fext_vol(nunkn*inode-nunkn+iknown,:) +
                        %fext_vol(iv,:)=fext_vol(iv,:)+fvol(igaus,:).*dvolu';
                    end
                end
        case 'AXI'
    end
   
end    

 % assembling of the external forces
    for k=1:nunkn*nnode
        %size(sparse(idx(k,:),1,efint(k,:),nndof,1));
        force_ext_vol = force_ext_vol + sparse(idx(k,:),1,fext_vol(k,:),nndof,1);
    end



end