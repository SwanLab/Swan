classdef IntegratorCutMesh < Integrator
    
    properties (Access = private)
        backgroundMeshType
        xCoordsIso
        cellContainingSubcell
    end
    
    methods (Access = public)
        
        function obj = IntegratorCutMesh(cParams)
            obj.init(cParams);
            obj.globalConnec          = cParams.globalConnec;
            obj.xCoordsIso            = cParams.xCoordsIso;
            obj.cellContainingSubcell = cParams.cellContainingSubcell;
            obj.backgroundMeshType    = cParams.backgroundMeshType;
        end
        
        function rhs = integrate(obj,fNodal)
            s.type      = 'CutMesh';
            s.mesh      = obj.mesh;
            s.backgroundMeshType  = obj.backgroundMeshType; % meshType
            s.fType     = 'Nodal';
            s.fNodal    = fNodal;
            s.quadOrder = 'LINEAR';
            s.npnod     = obj.npnod;
            s.xCoordsIso = obj.xCoordsIso;
            s.globalConnec = obj.globalConnec;
            s.cellContainingSubcell = obj.cellContainingSubcell;
            RHS = RHSintegrator.create(s);
            rhs = RHS.compute();
        end
        
    end
    
end

