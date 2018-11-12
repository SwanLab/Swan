classdef Filter < handle
    properties
        diffReacProb
        M0
        coordinates
        connectivities
        nnode
        nelem
        npnod
        ngaus
        shape
        x
        x_reg
        P_operator
    end
    
    methods
        function obj = Filter(problemID,scale)
            switch scale
                case 'MACRO'
                    obj.diffReacProb = DiffReact_Problem(problemID);
                case 'MICRO'
                    obj.diffReacProb = DiffReact_Problem_Micro(problemID);
            end
        end        
        function preProcess(obj)
            obj.diffReacProb.preProcess;
            quadrature = Quadrature.set(obj.diffReacProb.geometry.type);
            quadrature.computeQuadrature('LINEAR');
            obj.diffReacProb.element.interpolation_u.computeShapeDeriv(quadrature.posgp)
            obj.diffReacProb.geometry.computeGeometry(quadrature,obj.diffReacProb.element.interpolation_u);
            
            for igauss = 1:quadrature.ngaus
                obj.M0{igauss} = sparse(1:obj.diffReacProb.geometry.interpolation.nelem,1:obj.diffReacProb.geometry.interpolation.nelem,...
                    obj.diffReacProb.geometry.dvolu(:,igauss));
            end
            
            obj.coordinates = obj.diffReacProb.mesh.coord;
            obj.connectivities = obj.diffReacProb.mesh.connec;
            obj.nelem = obj.diffReacProb.geometry.interpolation.nelem;
            obj.nnode = obj.diffReacProb.geometry.interpolation.nnode;
            obj.npnod = obj.diffReacProb.geometry.interpolation.npnod;
            obj.ngaus = quadrature.ngaus;
            obj.shape = obj.diffReacProb.element.interpolation_u.shape;
        end
        
        function A_nodal_2_gauss = computeA(obj)
            A_nodal_2_gauss = sparse(obj.nelem,obj.npnod);
            fn = ones(1,obj.npnod);
            
            dirichlet_data = obj.connectivities';
            fe = zeros(obj.nnode,obj.nelem);
            
            fg = zeros(obj.ngaus,obj.nelem);
            
            for igaus = 1:obj.ngaus
                for inode = 1:obj.nnode
                    fe(inode,:) = fn(dirichlet_data(inode,:));
                    fg(igaus,:) = fg(igaus,:) + obj.shape(inode,igaus)*fe(inode,:);
                    A_nodal_2_gauss = A_nodal_2_gauss + sparse(1:obj.nelem,dirichlet_data(inode,:),ones(obj.nelem,1)*obj.shape(inode,igaus),obj.nelem,obj.npnod);
                end
            end
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
    methods (Static)
        function obj = create(settings)
            switch settings.filter
                case 'P1'
                    switch settings.optimizer
                        case {'MMA','PROJECTED GRADIENT','IPOPT'}
                            obj = Filter_P1_Density(settings.filename,settings.ptype);
                        case {'SLERP','HAMILTON-JACOBI','PROJECTED SLERP'}
                            switch settings.pdim
                                case '2D'
                                    obj = Filter_P1_LevelSet_2D(settings.filename,settings.ptype,settings.unfitted_mesh_algorithm);
                                case '3D'
                                    obj = Filter_P1_LevelSet_3D(settings.filename,settings.ptype,settings.unfitted_mesh_algorithm);
                            end
                    end
                case 'PDE'
                    switch settings.optimizer
                        case {'MMA','PROJECTED GRADIENT','IPOPT'}
                            obj = Filter_PDE_Density(settings.filename,settings.ptype);
                        case {'SLERP','HAMILTON-JACOBI','PROJECTED SLERP'}
                            switch settings.pdim
                                case '2D'
                                    obj = Filter_PDE_LevelSet_2D(settings.filename,settings.ptype,settings.unfitted_mesh_algorithm);
                                case '3D'
                                    obj = Filter_PDE_LevelSet_3D(settings.filename,settings.ptype,settings.unfitted_mesh_algorithm);
                            end
                    end
            end
        end
    end
end