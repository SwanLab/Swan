classdef Filter_P1 < Filter
    properties
    end
    methods
        function preProcess(obj,physicalProblem)
            preProcess@Filter(obj,physicalProblem)
            obj.P_operator=obj.computePoperator(obj.Msmooth,physicalProblem);
        end
    end
    methods (Static)
        function P_operator=computePoperator(Msmooth,physicalProblem)
            nelem=physicalProblem.mesh.nelem;
            nnode=physicalProblem.geometry.nnode;
            npnod=physicalProblem.mesh.npnod;
            
            lnods=zeros(nnode,nelem);
            for inode=1:nnode
                lnods(inode,:)=physicalProblem.mesh.connec(:,inode);
            end
            
            T_nodal_2_gauss = sparse(nelem,npnod);
            
            for inode=1:nnode
                T_nodal_2_gauss = T_nodal_2_gauss + sparse([1:nelem],[lnods(inode,:)],ones(nelem,1),nelem,npnod);
            end
            
            m = T_nodal_2_gauss*sum(Msmooth,2);
            P_operator = diag(m)\T_nodal_2_gauss;
        end
    end
end