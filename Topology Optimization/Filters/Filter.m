classdef Filter < handle
    
    properties (GetAccess = public, SetAccess = protected)
        diffReacProb
        ngaus
        nelem
    end
    
    properties (Access = protected)
        x
        x_reg
    end
    
    properties (GetAccess = protected, SetAccess = protected)
        P_operator
        
        geometry
        quadrature
        %interpolation
        
        mesh
        nnode
        npnod
        shape
        
        quadratureOrder
        
    end
    
    properties (Access = private)
        interp
    end
    
    methods (Access = public, Static)
       
        function obj = create(cParams)
           f = FilterFactory();
           obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
        
        function preProcess(obj)
            obj.diffReacProb.preProcess();
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
            obj.storeParams();            
        end
        
        function obj = createDiffReacProblem(obj,cParams)
            s = cParams.femSettings;
            if isprop(cParams,'mesh')
                s.mesh = cParams.mesh;
            end            
            switch s.scale
                case 'MACRO'
                    obj.diffReacProb = DiffReact_Problem(s);
                case 'MICRO'
                    obj.diffReacProb = DiffReact_Problem_Micro(s);
            end
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.createDiffReacProblem(cParams);
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;                        
        end
        
        function A_nodal_2_gauss = computeA(obj)
            A0 = sparse(obj.nelem,obj.npnod);
            A_nodal_2_gauss = cell(obj.ngaus,1);
            fn = ones(1,obj.npnod);
            
            nodes = obj.mesh.connec;
            fN = zeros(obj.nnode,obj.nelem);
            fg = zeros(obj.ngaus,obj.nelem);
            
            for igaus = 1:obj.ngaus
                A_nodal_2_gauss{igaus} = A0;
                for inode = 1:obj.nnode
                    node   = nodes(:,inode);
                    fN     = fn(node);
                    shapeN = obj.shape(inode,igaus);
                    fg(igaus,:) = fg(igaus,:) + shapeN*fN;
                    Ni = ones(obj.nelem,1)*shapeN;
                    A  = sparse(1:obj.nelem,node,Ni,obj.nelem,obj.npnod);
                    A_nodal_2_gauss{igaus} = A_nodal_2_gauss{igaus} + A;
                end
            end
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
        
        function createInterpolation(obj)
            obj.interp = Interpolation.create(obj.mesh,'LINEAR');    
        end                
        
        function createGeometry(obj)
            s.mesh = obj.mesh;
            obj.geometry = Geometry.create(s);
            obj.geometry.computeGeometry(obj.quadrature,obj.interp);   
        end
        
        function createQuadrature(obj)
            obj.quadrature = Quadrature.set(obj.mesh.type);
            obj.quadrature.computeQuadrature(obj.quadratureOrder);
        end
        
        function storeParams(obj)
            obj.nelem = obj.mesh.nelem;
            obj.nnode = obj.mesh.nnode;
            obj.npnod = obj.mesh.npnod;
            obj.ngaus = obj.quadrature.ngaus;
            obj.shape = obj.interp.shape;
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