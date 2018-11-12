function [costfunc] = cal_cost_funct(d_u,fext_adjoint,hnorm,structural_values,problembsc,Msmooth,StifMat,dim,element,gamma_gp,emass)

switch problembsc.TYPE
    case 'MACRO'
        costfunc  = d_u'*fext_adjoint;
        
    case 'MICRO'      
        switch problembsc.costfunct
            case 'alpha'
                alpha = problembsc.alpha;
                beta = problembsc.beta;
                inv_matCh = inv(structural_values.matCh);
                costfunc = projectionCh(inv_matCh,alpha,beta);
            case 'fraction'
                alpha = problembsc.alpha;
                beta = problembsc.beta;
                inv_matCh = inv(structural_values.matCh);
                costfunc_alpha_beta = projectionCh(inv_matCh,alpha,beta);
                costfunc_beta_alpha = projectionCh(inv_matCh,beta,alpha);
                costfunc_alpha_alpha = projectionCh(inv_matCh,alpha,alpha);
                costfunc_beta_beta = projectionCh(inv_matCh,beta,beta);
                costfunc = costfunc_alpha_beta/costfunc_alpha_alpha + costfunc_beta_alpha/costfunc_beta_beta;
            case 'C-C*'
                C_C = structural_values.matCh - problembsc.Ch_star;
                C_C = problembsc.selectiveC_Cstar.*C_C;
                costfunc = sum(bsxfun(@times,C_C(:),C_C(:)));
            case 'linear_components'
                C_C = structural_values.matCh - problembsc.Ch_star;
                C_C = problembsc.selectiveC_Cstar.*C_C;
                costfunc = sum(C_C(:));
        end


end
costfunc = costfunc/abs(hnorm);

end

