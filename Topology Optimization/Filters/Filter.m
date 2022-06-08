classdef Filter < handle
    
    properties (GetAccess = public, SetAccess = protected)
        diffReacProb
        ngaus
        nelem
    end
    
    properties (Access = protected)
        x
        x_reg
        M
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
        field
    end
    
    methods (Access = public, Static)
       
        function obj = create(cParams)
           f = FilterFactory();
           obj = f.create(cParams);
        end
        
    end
    
    methods (Access = public)
        
        function preProcess(obj)
            obj.createQuadrature();
            obj.createInterpolation();
            obj.createGeometry();
            obj.storeParams();
        end
        
        function obj = createDiffReacProblem(obj,cParams)
            s = cParams.femSettings;
            if isfield(cParams.femSettings,'LHStype')
                s.LHStype = cParams.femSettings.LHStype;
            else
                s.LHStype = 'DiffReactNeumann';
            end
            if isprop(cParams,'mesh')
                s.mesh = cParams.mesh;
            end
            s.type = 'DIFF-REACT';
            obj.diffReacProb = FEM.create(s);
        end
        
    end
    
    methods (Access = protected)
        
        function init(obj,cParams)
            obj.createDiffReacProblem(cParams);
            obj.mesh = cParams.mesh;
            obj.quadratureOrder = cParams.quadratureOrder;
            obj.createField();
            obj.computeMassMatrix();
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
        
        function createField(obj)
            s.mesh               = obj.mesh;
            s.ndimf              = 1;
            s.interpolationOrder = 'LINEAR';
            s.quadratureOrder    = 'QUADRATICMASS';
            obj.field = Field(s);
        end

        function computeMassMatrix(obj)
            s.type  = 'MassMatrix';
            s.mesh  = obj.mesh;
            s.field = obj.field;
            LHS = LHSintegrator.create(s);
            obj.M = LHS.compute();
        end
        
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
            obj.nnode = obj.mesh.nnodeElem;
            obj.npnod = obj.mesh.nnodes;
            obj.ngaus = obj.quadrature.ngaus;
            obj.shape = obj.interp.shape;
        end
        
%         function computeElementalMassMatrix(obj)
%             nel = obj.geometry.interpolation.nelem;
%             for igauss = 1:obj.quadrature.ngaus
%                 dvolu = obj.geometry.dvolu(:,igauss);
%                 obj.M0{igauss} = sparse(1:nel,1:nel,dvolu);
%             end
%         end
        
    end
    
end