function [ DtC ] = compliance_derivative_micro( matCh,g_gp,problembsc,dim,ngaus)

switch problembsc.costfunct
    case 'alpha'
        alpha = problembsc.alpha;
        beta = problembsc.beta;
        inv_matCh = inv(matCh);
        DtC = derivative_projectionCh (inv_matCh,alpha,beta,g_gp,dim,ngaus);
        
    case 'fraction'
        alpha = problembsc.alpha;
        beta = problembsc.beta;
        inv_matCh = inv(matCh);
        
        projection_alpha_beta = projectionCh(inv_matCh,alpha,beta);
        projection_beta_alpha = projectionCh(inv_matCh,beta,alpha);
        projection_alpha_alpha = projectionCh(inv_matCh,alpha,alpha);
        projection_beta_beta = projectionCh(inv_matCh,beta,beta);
        
        beta1 = projection_alpha_alpha*beta - projection_alpha_beta*alpha;
        beta2 = projection_beta_beta*alpha - projection_beta_alpha*beta;
        
        DtC_num1 = derivative_projectionCh (inv_matCh,alpha,beta1,g_gp,dim,ngaus);
        DtC_num2 = derivative_projectionCh (inv_matCh,beta,beta2,g_gp,dim,ngaus);
        
        DtC = DtC_num1/(projection_alpha_alpha)^2 + DtC_num2/(projection_beta_beta)^2;
        
    case 'C-C*'
        C_C = matCh - problembsc.Ch_star;
        C_C = problembsc.selectiveC_Cstar.*C_C;
        DtC1 = zeros(ngaus,dim.nelem);
        DtC = zeros(ngaus,dim.nelem);
        for igaus=1:ngaus
            for a=1:dim.nstre
                for b=1:dim.nstre
                    DtC1(igaus,:) = squeeze(g_gp(a,b,igaus,:));
                    DtC(igaus,:) = DtC(igaus,:) + 2*C_C(a,b)*DtC1(igaus,:);
                end
            end
        end
        
    case 'linear_components'
        C_C = problembsc.selectiveC_Cstar;
        DtC1 = zeros(ngaus,dim.nelem);
        DtC = zeros(ngaus,dim.nelem);
        for igaus=1:ngaus
            for a=1:dim.nstre
                for b=1:dim.nstre
                    DtC1(igaus,:) = squeeze(g_gp(a,b,igaus,:));
                    DtC(igaus,:) = DtC(igaus,:) + C_C(a,b)*DtC1(igaus,:);
                end
            end
        end

end


end

