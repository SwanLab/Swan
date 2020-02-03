classdef DOF_Stokes < DOF
    %DOF_Stokes Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = DOF_Stokes(filename,mesh,geometry,interp,nFields) % Replace mesh for pdim
            nunkn_u = 2;
            nunkn_p = 1;
            nFields = numel(interp);
            obj.nunkn = [nunkn_u nunkn_p];
            [dirichlet_data,neumann_data,full_dirichlet_data] = Preprocess.getBC_fluids(filename,mesh,geometry,interp);
            obj.getDOFconditions(nFields,dirichlet_data,neumann_data,full_dirichlet_data);
            obj.computeDOF(interp);
        end
    
       function obj = computeDOF(obj,interp)
            nfields = numel(interp);
            for ifield = 1:nfields
                nunkn = obj.nunkn(ifield);
                nnode = interp{ifield}.nnode;
                npnod = interp{ifield}.npnod;
                obj.in_elem{ifield} = obj.compute_idx(interp{ifield}.T,nunkn,nnode);
                obj.ndof(ifield) = nunkn*npnod;
                obj.constrained{ifield} = obj.compute_constrained_dof(ifield);
                obj.free{ifield} = obj.compute_free_dof(ifield);
            end
        end
    end
end

