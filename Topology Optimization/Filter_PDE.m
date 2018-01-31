classdef Filter_PDE < Filter
    properties
        dof_per
        solver
        rhs
        epsilon
        A_nodal_2_gauss
    end
    methods
        function preProcess(obj,physicalProblem)
            preProcess@Filter(obj,physicalProblem);
            switch physicalProblem.mesh.scale
                case 'MACRO'
                    obj.dof_per=DOF(physicalProblem.geometry.nnode,physicalProblem.mesh.connec,1,physicalProblem.mesh.npnod,physicalProblem.bc.fixnodes_perimeter);
                    obj.solver = Solver_Dirichlet_Conditions;
                    data.fixnodes=physicalProblem.bc.fixnodes_perimeter;
                    obj.solver.setSolverVariables(data);
                case 'MICRO'
                    obj.dof_per=DOF(physicalProblem.geometry.nnode,physicalProblem.mesh.connec,1,physicalProblem.mesh.npnod,[]);
                    obj.solver = Solver_Periodic;
                    data.pnodes = physicalProblem.bc.pnodes;
                    data.nunkn=1;
                    obj.solver.setSolverVariables(data);
            end


            obj.epsilon=0.03;
            obj.A_nodal_2_gauss=obj.computeA(physicalProblem);
        end
        function A_nodal_2_gauss=computeA(obj,physProblem)
            nelem=physProblem.mesh.nelem; nnode=physProblem.geometry.nnode;
            A_nodal_2_gauss = sparse(nelem,physProblem.mesh.npnod);
            fn=ones(1,physProblem.mesh.npnod);
            lnods=obj.connectivities';
            fe=zeros(nnode,nelem);
            for inode=1:nnode
                fe(inode,:)=fn(lnods(inode,:));
            end
            fg=zeros(physProblem.geometry.ngaus,nelem);
            shape=physProblem.geometry.shape;
            for igaus=1:physProblem.geometry.ngaus
                for inode=1:nnode
                    fg(igaus,:) = fg(igaus,:) + shape(inode)*fe(inode,:);
                    
                    A_nodal_2_gauss = A_nodal_2_gauss + sparse([1:nelem],[lnods(inode,:)],ones(nelem,1)*shape(inode),nelem,physProblem.mesh.npnod);
                    % B_nodal_2_gauss = B_nodal_2_gauss + sparse([1:nelem],[lnods(inode,:)],dvolu*shape(inode),nelem,dim.npnod);
                end
            end
        end
    end
end