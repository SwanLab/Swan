classdef UnfittedMesh_AbstractBuilder < handle
    properties (GetAccess = public, SetAccess = private, Abstract)
        type
        max_subcells
        nnodes_subcell
        
        subcells_Mesher
        cutPoints_Calculator
    end
    
    properties (GetAccess = public, SetAccess = private)
        mesh_unfitted
    end
    
    methods (Access = public)   
        function obj = UnfittedMesh_AbstractBuilder(mesh_background,interpolation_background)
            obj.createBasicUnfittedMesh(mesh_background,interpolation_background);
        end
        
        function mesh_unfitted = buildMesh(obj)
            obj.mesh_unfitted.type = obj.type;
            obj.mesh_unfitted.max_subcells = obj.max_subcells;
            obj.mesh_unfitted.nnodes_subcell = obj.nnodes_subcell;
            obj.mesh_unfitted.subcells_Mesher =	obj.subcells_Mesher;
            obj.mesh_unfitted.cutPoints_Calculator = obj.cutPoints_Calculator;
            
            mesh_unfitted = obj.mesh_unfitted;
        end
    end
    
        methods (Access = private)
        function createBasicUnfittedMesh(obj,mesh_background,interpolation_background)
            obj.mesh_unfitted = Mesh_Unfitted(mesh_background,interpolation_background);
        end
    end
end

