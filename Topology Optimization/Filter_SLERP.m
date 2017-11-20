classdef Filter_SLERP < handle
    properties
        
    end
    methods (Static)
        function x_gp = getP0fromP1(x,coordinates,conectivities,P_operator)
            M2 = faireF2(coordinates',conectivities',x);
            x_gp = P_operator*M2;
        end
        function x_reg = getP1fromP0(x,M0,P_operator)
            x_reg = P_operator'*M0*x;
        end
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