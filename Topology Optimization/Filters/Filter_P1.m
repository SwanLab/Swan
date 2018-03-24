classdef Filter_P1 < Filter
    properties
        P_operator
    end
    
    methods
        function obj = Filter_P1(problemID,scale)
            obj@Filter(problemID,scale);
        end
        
        function preProcess(obj)
            preProcess@Filter(obj)
            obj.P_operator = obj.computePoperator(obj.diffReacProb.element.M);
        end
        
        function P_operator=computePoperator(obj,Msmooth)
            dirichlet_data=zeros(obj.nnode,obj.nelem);
            for inode=1:obj.nnode
                dirichlet_data(inode,:)=obj.connectivities(:,inode);
            end
            
            T_nodal_2_gauss = sparse(obj.nelem,obj.npnod);
            
            for inode=1:obj.nnode
                T_nodal_2_gauss = T_nodal_2_gauss + sparse(1:obj.nelem,dirichlet_data(inode,:),ones(obj.nelem,1),obj.nelem,obj.npnod);
            end
            
            m = T_nodal_2_gauss*sum(Msmooth,2);
            P_operator = diag(m)\T_nodal_2_gauss;
        end
    end
end