function [g_nodal,gfunc_til] = cal_nodgfunc(phigp,phifunct,structural_values,problembsc,weights,hnorm,Perim0,element,constr,lambda,penalty,Msmooth,coordinates,dim,StifMat)
                   

npnod=dim.npnod; nndof=dim.nndof; ndime=dim.ndime; nelem=dim.nelem;
nnode=dim.nnode; neleq=dim.neleq; nunkn = dim.nunkn; nstre = dim.nstre;
[~,~,ngaus] = cal_posgp_weigp(element.type,ndime,nnode,element.ngaus);


smoothDtC = problembsc.smoothDtC;
switch problembsc.TYPE
    case 'MICRO'
        tstres = structural_values.tstres;
        [DtC]=cal_DtC(dim,element.material.opt_epsi,tstres,ngaus,phigp,element);
        structural_values.DtC = DtC;
        matCh = structural_values.matCh;
       
        switch problembsc.costfunct
            case 'HORIZONTAL'
                [DtJtil] = cal_DtJtil_hori(smoothDtC,matCh,DtC,npnod,ngaus,nelem,ndime,nnode,Msmooth,coordinates);
            case 'BULK_MAX'
                [DtJtil] = cal_DtJtil_bulk(smoothDtC,matCh,DtC,npnod,ngaus,nelem,ndime,nnode,Msmooth,coordinates);
            case 'SHEAR_MAX'
                [DtJtil] = cal_DtJtil_shear(smoothDtC,matCh,DtC,npnod,ngaus,nelem,ndime,nnode,Msmooth,coordinates);
            case 'MIN_MINUS_STIFF'
                [DtJtil] = cal_DtJtil_min_minus_stiff(smoothDtC,matCh,DtC,npnod,ngaus,nelem,ndime,nnode,weights,hnorm,Msmooth,coordinates);
            case 'MIN_INV_STIFF'
                [DtJtil] = cal_DtJtil_min_inv_stiff(smoothDtC,matCh,DtC,dim,ngaus,weights,hnorm,Msmooth,coordinates,element,problembsc);
            case 'MIN_STIFF_INV'
                [DtJtil] = cal_DtJtil_min_stiff_inv(smoothDtC,matCh,DtC,dim,ngaus,weights,hnorm,Msmooth,coordinates,element,problembsc);
                
        end
        
    case 'MACRO'
        DtJtil = structural_values.DtJtil/abs(hnorm); 
end



dcontr_domega = 1/element.material.Vfrac;

tipo = problembsc.lagrange_update;

switch tipo
    case 'AUGMENTED'
        %DtC = max(L*dcontr_domega + element.material.penalty*constr*dcontr_domega,0);
        DtCons = lambda*dcontr_domega + 1/penalty*constr*dcontr_domega;
    case 'POTENCIAL'
        n = problembsc.potential_constraint;
        DtCons = lambda*n*constr^(n-1)*dcontr_domega;

end

g_minus = DtJtil - DtCons;
g_plus  = DtJtil + DtCons; 


if (smoothDtC==1)
    g_nodal = -g_minus;
    outside = (phifunct > 0);
    if any(outside)
        g_nodal(outside)= g_plus(outside);
    end
    
elseif (smoothDtC==0)
    gfunc = zeros(ngaus,nelem);
    for igaus=1:ngaus
        gfunc(igaus,:) = -g_minus(igaus,:);
        outside= (phigp(igaus,:) > 0);
        if any(outside)
            gfunc(igaus,outside)= g_plus(igaus,outside);
        end
    end
    
    [coordinatesn,coordinatesa] = init_coord(coordinates);
    [g_nodal] = smooth(nelem,npnod,ndime,nnode,gfunc,Msmooth,element,coordinatesn,coordinatesa,problembsc);
    
end

switch problembsc.regularization_perimeter
    case 'YES'
        alpha_per = problembsc.alpha_perimeter;
        g_perimeter = compute_Perimeter_derivative(dim,element,element.epsilon,StifMat,Msmooth,phifunct);
        g_nodal = g_nodal + alpha_per*g_perimeter/Perim0;   
end

[gfunc_til] = interpol(g_nodal,element,dim,problembsc);


end