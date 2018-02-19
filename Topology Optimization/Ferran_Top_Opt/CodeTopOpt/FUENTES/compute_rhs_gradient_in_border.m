function rhs = compute_rhs_gradient_in_border(gradient,phi,element,dim,problembsc,coordinates)

        [~,x_gp_border,perimeter,index_border] = cal_vol_mat(phi,dim,element,problembsc,coordinates);
       pos_gp_border_nat = global_coordinate_2_natural_coordinate(x_gp_border',index_border,element,problembsc,dim,coordinates);
        shape = shape_deriv_functions_triangles_special_point(pos_gp_border_nat');
        
        norm_grad_phi = compute_norm_gradient_phi(phi,element,dim,problembsc,coordinates,index_border);
        rhs = zeros(dim.npnod,1);
        for inode = 1:dim.nnode
            dirichlet_data = element.conectivities(index_border,inode);
            rhs(dirichlet_data,1) = rhs(dirichlet_data,1) + 1./norm_grad_phi.*gradient(index_border).*(perimeter.*shape(inode,:))';
        end

end