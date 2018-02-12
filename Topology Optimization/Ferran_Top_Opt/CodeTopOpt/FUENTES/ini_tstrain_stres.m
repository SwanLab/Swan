function [tstrain,tstres] = ini_tstrain_stres(dim,element,problembsc)

ndime=dim.ndime; nelem=dim.nelem;
nnode=dim.nnode; nstre = dim.nstre;
[~,~,ngaus] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);
switch problembsc.phisical_type;
    case {'ELASTIC'}
        switch element.material.subtype
            case 'PLANESTRAIN'
                tstrain = zeros(nstre,ngaus,nstre,nelem);
                tstres = zeros(nstre,ngaus,nstre+1,nelem);
                
            case 'PLANESTRES'
                tstrain = zeros(nstre,ngaus,nstre+1,nelem);
                tstres = zeros(nstre,ngaus,nstre,nelem);
        end
        
    case {'THERMAL'}
        tstrain = zeros(nstre,ngaus,nstre,nelem);
        tstres = zeros(nstre,ngaus,nstre,nelem);
        
end