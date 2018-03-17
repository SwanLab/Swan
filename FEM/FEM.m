classdef FEM < handle
    %FEM Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = public)
        problemID
        geometry
        mesh
        dof
        element
        Fext
        LHS
        variables
    end
    
    %% Restricted properties definition ===================================
    properties (GetAccess = ?Postprocess, SetAccess = private)
    end
    
    %% Private properties definition ======================================
    properties (Access = protected)
        solver
    end
    
    %% Public methods definition ==========================================
    methods (Static, Access = public)
        function obj = create(problemID)
            mesh = Mesh(problemID); % Mesh defined twice, but almost free
            switch mesh.ptype
                case 'ELASTIC'
                    obj = Elastic_Problem(problemID);
                case 'THERMAL'
                    obj = Thermal_Problem(problemID);
                case 'DIFF-REACT'
                    obj = DiffReact_Problem(problemID);
                case 'HYPERELASTIC'
                    obj = Hyperelastic_Problem(problemID);
                case 'Stokes'
                    obj = Stokes_Problem(problemID);
            end
        end
    end
    
    methods (Abstract)
        %% !! OBLIGATED METHODS !!
    end 
end
