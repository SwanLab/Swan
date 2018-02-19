function pos_nat = global_coordinate_2_natural_coordinate(global_coord,index_element,element,problembsc,dim,coordinates)

nelem=dim.nelem; nndof=dim.nndof; nnode=dim.nnode;
ndime=dim.ndime; npnod=dim.npnod; nunkn = dim.nunkn; nstre = dim.nstre;
ptype = problembsc.problemtype;
[posgp,weigp,ngaus] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);
[~,~,injacb] = cal_cartd(ngaus,posgp,element,ndime,nnode,nelem,coordinates,coordinates,ptype);
pos_nat = zeros(sum(index_element),ndime);
for idime = 1:ndime
        for jdime = 1:ndime
            dirichlet_data_1 = element.conectivities(index_element,1);
            pos_nat(:,idime) = pos_nat(:,idime)+squeeze(injacb(jdime,idime,index_element)).*(global_coord(:,jdime) - coordinates(dirichlet_data_1,jdime));
        end
end




end