function [gfunc] = cal_gfunc_micro(phigp,lambda,DtC,matCh,problembsc,npnod,ngaus,nelem)


DtJtil = zeros(ngaus,nelem);
smoothDtC = 0;
switch problembsc.costfunct
    case 'HORIZONTAL'
        [DtJtil] = cal_DtJtil_hori(smoothDtC,matCh,DtC,DtJtil,npnod,ngaus,nelem);
    case 'BULK_MAX'
        [DtJtil] = cal_DtJtil_bulk(smoothDtC,matCh,DtC,DtJtil,npnod,ngaus,nelem);
    case 'SHEAR_MAX'
        [DtJtil] = cal_DtJtil_shear(smoothDtC,matCh,DtC,DtJtil,npnod,ngaus,nelem);
    case 'MIN_MINUS_STIFF'
        [DtJtil] = cal_DtJtil_min_minus_stiff(smoothDtC,matCh,DtC,DtJtil,npnod,ngaus,nelem,weights,hnorm);
    case 'MIN_INV_STIFF'
        [DtJtil] = cal_DtJtil_min_inv_stiff(smoothDtC,matCh,DtC,DtJtil,npnod,ngaus,nelem,weights,hnorm);
    case 'MIN_STIFF_INV'
        [DtJtil] = cal_DtJtil_min_stiff_inv(smoothDtC,matCh,DtC,DtJtil,npnod,ngaus,nelem,weights,hnorm);
        
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

gfunc = zeros(ngaus,nelem);
for igaus=1:ngaus
    gfunc(igaus,:) = -g_minus(igaus,:);
    outside= (phigp(igaus,:) > 0);
    if any(outside)
        gfunc(igaus,outside)= g_plus(igaus,outside);
    end
end



end