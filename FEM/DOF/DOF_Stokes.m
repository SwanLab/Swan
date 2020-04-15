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
            obj.computeDOF(mesh,interp);
        end
        
        function obj = computeDOF(obj,mesh,interp)
            nfields = numel(interp);
            for ifield = 1:nfields
                nunkn = obj.nunkn(ifield);
                int = interp{ifield};
              %  npnod = int.npnod;

                %T = interp{ifield}.T;
                [T,npnod] = obj.computeConnec(mesh,int);
                nnode = size(T,2);
                
                obj.in_elem{ifield} = obj.compute_idx(T,nunkn,nnode);
                obj.ndof(ifield) = nunkn*npnod;
                obj.constrained{ifield} = obj.compute_constrained_dof(ifield);
                obj.free{ifield} = obj.compute_free_dof(ifield);
            end
        end
    end
end

