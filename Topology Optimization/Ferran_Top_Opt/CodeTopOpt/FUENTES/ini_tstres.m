function tstres = ini_tstres(dim,element)

ndime=dim.ndime; nelem=dim.nelem;
nnode=dim.nnode; nstre = dim.nstre;
[~,~,ngaus] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);

tstres = zeros(nstre,ngaus,nstre,nelem);

end

