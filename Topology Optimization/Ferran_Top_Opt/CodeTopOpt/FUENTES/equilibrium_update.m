function [cost,g_gp,structural_values] = equilibrium_update(gamma_gp,element,problembsc,fixnodes,coordinates,fext,fext_adjoint,dim,Msmooth,M0,Stiff_smooth,emass,h_C_0)
global iter

[structural_values] = module_M(gamma_gp,element,fixnodes,problembsc,coordinates,fext,dim);
iter = iter + 1;

switch problembsc.TYPE
    
    case 'MACRO'
        if isempty(fext_adjoint.pointload)
            structural_values_adjoint.strain = -structural_values.strain;
            structural_values_adjoint.d_u = -structural_values.d_u;
            structural_values_adjoint.fext = structural_values.fext;
        else
            fext_adjoint.pointload(:,3) = -fext_adjoint.pointload(:,3);    
            [structural_values_adjoint] = module_M(gamma_gp,element,fixnodes,problembsc,coordinates,fext_adjoint,dim);
            structural_values_adjoint.fext = -structural_values_adjoint.fext;
        end
        structural_values.matCh = zeros(3);
        
    case 'MICRO'
        structural_values_adjoint.strain = structural_values.tstrain;
        structural_values.stres = structural_values.tstres;
        structural_values_adjoint.d_u = -structural_values.d_u;
        structural_values_adjoint.fext = structural_values.fext;

        
    otherwise
        error('Type not valid (MACRO/MICRO)');
        
end

[cost] = cal_cost_funct(structural_values.d_u,structural_values_adjoint.fext,h_C_0,structural_values,problembsc,Msmooth,Stiff_smooth,dim,element,gamma_gp,emass);

[g_gp] = compute_topological_derivative(structural_values,structural_values_adjoint,gamma_gp,problembsc,h_C_0,element,Msmooth,coordinates,dim);

structural_values.fobj = cost;


end