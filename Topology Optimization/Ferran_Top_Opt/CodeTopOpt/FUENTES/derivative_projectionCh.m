function derivative_projection_Ch = derivative_projectionCh (inv_matCh,alpha,beta,g_gp,dim,ngaus)

weights = alpha*beta';
        
weights_inv = inv_matCh*weights*inv_matCh;

DtC1 = zeros(ngaus,dim.nelem);
DtC = zeros(ngaus,dim.nelem);
for igaus=1:ngaus
    for a=1:dim.nstre
        for b=1:dim.nstre
            DtC1(igaus,:) = squeeze(g_gp(a,b,igaus,:));
            DtC(igaus,:) = DtC(igaus,:) - weights_inv(a,b)*DtC1(igaus,:);
        end
    end
end
derivative_projection_Ch = DtC;

end