function [g_gp] = compute_topological_derivative(structural_values,structural_values_adjoint,gamma,problembsc,hnorm,element,Msmooth,coordinates,dim)
                   
npnod=dim.npnod; nndof=dim.nndof; ndime=dim.ndime; nelem=dim.nelem;
nnode=dim.nnode; neleq=dim.neleq; nunkn = dim.nunkn; nstre = dim.nstre;

%[coordinatesn,coordinatesa] = init_coord(coordinates);

[P_interp] = derivative_constitutive_tensor(element,dim,gamma);

global post_info

etype = element.type;
ptype = problembsc.problemtype;
[~,~,ngaus] = cal_posgp_weigp(etype,ndime,nnode,element.ngaus);


switch problembsc.TYPE
    case 'MACRO'
        strain = structural_values_adjoint.strain;
        stres = structural_values.stres;
        g_gp = zeros(dim.nelem,ngaus);

        for igaus=1:ngaus
            for istre=1:dim.nstre
                for jstre = 1:dim.nstre
                    g_gp(:,igaus) = g_gp(:,igaus) + (squeeze(stres(igaus,istre,:)).*squeeze(P_interp(istre,jstre,:)).*squeeze(strain(igaus,jstre,:)));
                end
            end
        end
        
    case 'MICRO'
        strain = structural_values_adjoint.strain;
        stres = structural_values.stres;
        matCh = structural_values.matCh;
        g_gp = zeros(dim.nstre,dim.nstre,ngaus,dim.nelem);
        
        for i = 1:dim.nstre
            for j = 1:dim.nstre
                for igaus=1:ngaus
                    for istre=1:dim.nstre
                        for jstre = 1:dim.nstre
                            g_gp(i,j,igaus,:) = squeeze(g_gp(i,j,igaus,:)) + (squeeze(stres(i,igaus,istre,:)).*squeeze(P_interp(istre,jstre,:)).*squeeze(strain(j,igaus,jstre,:)));
                        end
                    end
                end
            end
        end
        post_info.DtC = g_gp;
        post_info.dim = dim;
        post_info.ngaus = ngaus;
        g_gp = compliance_derivative_micro(matCh,g_gp,problembsc,dim,ngaus);
        

        
end

g_gp = g_gp/abs(hnorm);

