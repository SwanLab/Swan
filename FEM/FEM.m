classdef FEM < handle
    %FEM Summary of this class goes here
    %   Detailed explanation goes here
    
    %% Public GetAccess properties definition =============================
    properties (GetAccess = public, SetAccess = public)
        geometry
        Fext
        LHS
    end
    
    %% Restricted properties definition ===================================
    properties (GetAccess = ?Postprocess, SetAccess = private)
    end
    
    %% Private properties definition ======================================
    properties (Access = protected)
        solver
    end
    
    %% Public methods definition ==========================================
    methods (Access = public)
    end
end
