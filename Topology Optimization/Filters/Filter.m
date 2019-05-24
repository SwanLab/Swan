classdef Filter < handle
    
    properties (GetAccess = public, SetAccess = private)
        diffReacProb
        ngaus
        nelem
    end
    
    properties (Access = protected)
        x
        x_reg
    end
    
    properties (GetAccess = protected, SetAccess = private)
        P_operator
        
        geometry
        quadrature
        interpolation
        
        mesh
        nnode
        npnod
        shape
    end
    
    properties (Access = private)
        quadratureOrder
    end
    
    methods (Access = public)
        
        function obj = Filter(cParams)
            obj.createDiffReacProblem(cParams);
            obj.mesh = cParams.designVar.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
        end
        
        function preProcess(obj)
            obj.diffReacProb.preProcess();
            
            obj.setQuadrature();
            obj.setInterpolation();
            
            obj.interpolation.computeShapeDeriv(obj.quadrature.posgp)
            obj.computeGeometry();
            obj.storeParams();
            
            obj.P_operator = obj.computePoperator(obj.diffReacProb.element.M);
        end
        
        function obj = createDiffReacProblem(obj,cParams)
            s = cParams.femSettings;
            switch cParams.designVar.mesh.scale
                case 'MACRO'
                    obj.diffReacProb = DiffReact_Problem(s);
                case 'MICRO'
                    obj.diffReacProb = DiffReact_Problem_Micro(s);
            end
        end
        
    end
    
    methods (Access = protected)
        
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
            dirichlet_data = zeros(obj.nnode,obj.nelem);
            for inode = 1:obj.nnode
                dirichlet_data(inode,:)=obj.mesh.connec(:,inode);
            end
            
            T_nodal_2_gauss = sparse(obj.nelem,obj.npnod);
            
            for inode = 1:obj.nnode
                T_nodal_2_gauss = T_nodal_2_gauss + sparse(1:obj.nelem,dirichlet_data(inode,:),ones(obj.nelem,1),obj.nelem,obj.npnod);
            end
            
            m = T_nodal_2_gauss*sum(Msmooth,2);
            P_operator = diag(m)\T_nodal_2_gauss;
        end
        
        function itHas = xHasChanged(obj,x)
            itHas = ~isequal(x,obj.x);
        end
        
        function updateStoredValues(obj,x,x0)
            obj.x = x;
            obj.x_reg = x0;
        end
        
    end
    
    methods (Access = private)
        
        function computeGeometry(obj)
            obj.geometry = Geometry(obj.mesh,'LINEAR');
            obj.geometry.interpolation.computeShapeDeriv(obj.quadrature.posgp);
            obj.geometry.computeGeometry(obj.quadrature,obj.geometry.interpolation);
        end
        
        function setQuadrature(obj)
            obj.quadrature = Quadrature.set(obj.mesh.geometryType);
            obj.quadrature.computeQuadrature(obj.quadratureOrder);
        end
        
        function setInterpolation(obj)
            obj.interpolation = Interpolation.create(obj.mesh,'LINEAR');
        end
        
        function storeParams(obj)
            obj.nelem = obj.geometry.interpolation.nelem;
            obj.nnode = obj.geometry.interpolation.nnode;
            obj.npnod = obj.geometry.interpolation.npnod;
            obj.ngaus = obj.quadrature.ngaus;
            obj.shape = obj.interpolation.shape;
        end
        
        function computeElementalMassMatrix(obj)
            nel = obj.geometry.interpolation.nelem;
            for igauss = 1:obj.quadrature.ngaus
                dvolu = obj.geometry.dvolu(:,igauss);
                obj.M0{igauss} = sparse(1:nel,1:nel,dvolu);
            end
        end
        
    end
    
end