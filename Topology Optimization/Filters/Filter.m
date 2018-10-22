classdef Filter < handle
    properties
        diffReacProb
        M0
        mesh
        nnode
        nelem
        npnod
        ngaus
        shape
        x
        x_reg
        P_operator
        
        geometry
        quadrature
        interpolation
    end
    
    methods (Static, Access = public)
        function obj = create(settings)
            switch settings.filter
                case 'P1'
                    switch settings.optimizer
                        case {'MMA','PROJECTED GRADIENT','IPOPT'}
                            obj = Filter_P1_Density;
                        case {'SLERP','HAMILTON-JACOBI','PROJECTED SLERP'}
                            switch settings.pdim
                                case '2D'
                                    obj = Filter_P1_LevelSet_2D_Interior;
                                case '3D'
                                    obj = Filter_P1_LevelSet_3D_Interior;
                            end
                    end
                case 'PDE'
                    switch settings.optimizer
                        case {'MMA','PROJECTED GRADIENT','IPOPT'}
                            obj = Filter_PDE_Density;
                        case {'SLERP','HAMILTON-JACOBI','PROJECTED SLERP'}
                            switch settings.pdim
                                case '2D'
                                    obj = Filter_PDE_LevelSet_2D_Interior;
                                case '3D'
                                    obj = Filter_PDE_LevelSet_3D_Interior;
                            end
                    end
            end
        end
    end
    
    methods (Access = public)
        function preProcess(obj)
            obj.diffReacProb.preProcess;
            obj.mesh = obj.diffReacProb.mesh;
            
            obj.setQuadrature;
            obj.setInterpolation;
            
            obj.interpolation.computeShapeDeriv(obj.quadrature.posgp)
            obj.computeGeometry;
            obj.loadParams;
            
            for igauss = 1:obj.quadrature.ngaus
                obj.M0{igauss} = sparse(1:obj.geometry.interpolation.nelem,1:obj.geometry.interpolation.nelem,...
                    obj.geometry.dvolu(:,igauss));
            end
        end
        
        function obj = setupFromMesh(obj,mesh,scale)
            obj.setDiffusionReactionProblem(scale);
            obj.diffReacProb.setupFromMesh(mesh);
        end
        
        function obj = setupFromGiDFile(obj,problemID,scale)
            obj.setDiffusionReactionProblem(scale);
            obj.diffReacProb.setupFromGiDFile(problemID);
        end
        
        function A_nodal_2_gauss = computeA(obj)
            A_nodal_2_gauss = sparse(obj.nelem,obj.npnod);
            fn = ones(1,obj.npnod);
            
            dirichlet_data = obj.mesh.connec';
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
        
        function P_operator = computePoperator(obj,Msmooth)
            dirichlet_data=zeros(obj.nnode,obj.nelem);
            for inode=1:obj.nnode
                dirichlet_data(inode,:)=obj.mesh.connec(:,inode);
            end
            
            T_nodal_2_gauss = sparse(obj.nelem,obj.npnod);
            
            for inode=1:obj.nnode
                T_nodal_2_gauss = T_nodal_2_gauss + sparse(1:obj.nelem,dirichlet_data(inode,:),ones(obj.nelem,1),obj.nelem,obj.npnod);
            end
            
            m = T_nodal_2_gauss*sum(Msmooth,2);
            P_operator = diag(m)\T_nodal_2_gauss;
        end
    end
    
    methods (Access = private)
        function setDiffusionReactionProblem(obj,scale)
            switch scale
                case 'MACRO'
                    obj.diffReacProb = DiffReact_Problem;
                case 'MICRO'
                    obj.diffReacProb = DiffReact_Problem_Micro;
            end
        end
        
        function computeGeometry(obj)
            obj.geometry = Geometry(obj.mesh,'LINEAR');
            obj.geometry.interpolation.computeShapeDeriv(obj.quadrature.posgp);
            obj.geometry.computeGeometry(obj.quadrature,obj.geometry.interpolation);
        end
        
        function setQuadrature(obj)
            obj.quadrature = Quadrature.set(obj.mesh.geometryType);
            obj.quadrature.computeQuadrature('LINEAR');
        end
        
        function setInterpolation(obj)
            obj.interpolation = Interpolation.create(obj.mesh,obj.quadrature.order);
        end
        
        function loadParams(obj)
            obj.nelem = obj.geometry.interpolation.nelem;
            obj.nnode = obj.geometry.interpolation.nnode;
            obj.npnod = obj.geometry.interpolation.npnod;
            obj.ngaus = obj.quadrature.ngaus;
            obj.shape = obj.interpolation.shape;
        end
    end
end