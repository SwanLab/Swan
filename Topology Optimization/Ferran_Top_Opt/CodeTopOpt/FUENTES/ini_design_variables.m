function [gamma0,gamma0_nodal,design_variable0] = ini_design_variables(dim,element,coordinates,Msmooth,problembsc,Stiff_smooth,P_operator,algorithm,case_domain,scalar_product,init_design_variable)

%% Special cases
switch problembsc.TYPE 
    case 'MACRO'
        switch case_domain
            case {'Square','SquareFine'}
                element.initial_holes = coordinates(:,1).^2 + coordinates(:,2).^2 - 0.025^2 < 0;
                fracc = 1;
            otherwise
                element.initial_holes = false(size(coordinates,1),1);
                fracc = 1;
        end

end

%% Select initial holes
switch init_design_variable
    case 'circle'
        width = max(coordinates(:,1)) - min(coordinates(:,1));
        height = max(coordinates(:,2)) - min(coordinates(:,2));
        center_x = 0.5*(max(coordinates(:,1)) + min(coordinates(:,1)));
        center_y = 0.5*(max(coordinates(:,2)) + min(coordinates(:,2)));
        radius = 0.2*min([width,height]);

        element.initial_holes = (coordinates(:,1)-center_x).^2 + (coordinates(:,2)-center_y).^2 - radius^2 < 0;
        fracc = 1;

    case 'horiz'
        element.initial_holes = coordinates(:,2) > 0.6 | coordinates(:,2) < 0.4;
        fracc = 1;

    case 'rand'
        element.initial_holes = rand(size(coordinates,1),1) > 0.1;
        fracc = 1;

    case 'square'
        width = max(coordinates(:,1)) - min(coordinates(:,1));
        height = max(coordinates(:,2)) - min(coordinates(:,2));
        center_x = 0.5*(max(coordinates(:,1)) + min(coordinates(:,1)));
        center_y = 0.5*(max(coordinates(:,2)) + min(coordinates(:,2)));

        offset_x = 0.2*width;
        offset_y = 0.2*height;

        xrange = coordinates(:,1) < (center_x+offset_x) & coordinates(:,1) > (center_x-offset_x);
        yrange = coordinates(:,2) < (center_y+offset_y) & coordinates(:,2) > (center_y-offset_y);

        element.initial_holes = and(xrange,yrange);
        fracc = 1;
    case 'full'
        element.initial_holes = false(size(coordinates,1),1);
        fracc = 1;

     case 'feasible'
        element.initial_holes = false(size(coordinates,1),1);
        fracc = min(1,element.Vfrac);
        
    otherwise
        error('Initialize design variable case not detected.');
end

%% Initialize design variable
switch algorithm
    case 'level_set'
        phi = init_phi(dim,element,coordinates,Msmooth,problembsc);
        phi0 = phi/sqrt(scalar_product(phi,phi));
%         [gamma0] = regularization_diff_reaction_equation(phi0,estimate_mesh_size(element,coordinates),'P1_kernel',Msmooth,Stiff_smooth,element.regularization,coordinates,element,problembsc,dim,P_operator,0,'gamma',0,1,algorithm);
        gamma0_nodal = phi0 < 0;
        gamma0 = interpol(gamma0_nodal,element,dim,problembsc,coordinates);
        design_variable0 = phi0;

    otherwise %case 'Projected_gradient'
        gamma0 = fracc*ones(dim.npnod,1);
        gamma0(element.initial_holes) = 0;
        design_variable0 = gamma0;
        gamma0_nodal = gamma0;
        gamma0 = interpol(gamma0_nodal,element,dim,problembsc,coordinates);
end
end