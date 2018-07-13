function [volume,dvolume] = cal_omega(gamma_gp,dim,element,problembsc,coordinates)


nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;

[coordinatesn,coordinatesa] = init_coord(coordinates);
etype = element.type;
ptype = problembsc.problemtype;
[posgp,weigp,ngaus] = cal_posgp_weigp(etype,ndime,nnode,element.ngaus);
evol=zeros(ngaus,nelem);
geometric_volum=zeros(ngaus,nelem);
%vol=zeros(npnod,1);

for igaus=1:ngaus
    [~,djacb] = cal_cartd(igaus,posgp,element,ndime,nnode,nelem,coordinatesn,coordinatesa,ptype);
    evol(igaus,:) = weigp(igaus)*(djacb.*gamma_gp)';
    geometric_volum(igaus,:) = weigp(igaus)*djacb';
end

volume = sum(sum(evol))/sum(sum(geometric_volum));
dvolume = geometric_volum(igaus,:)/(sum(sum(geometric_volum)));

end






