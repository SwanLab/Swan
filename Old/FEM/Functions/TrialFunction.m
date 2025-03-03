classdef TrialFunction < LagrangianFunction
    
    % properties (Access = public)
    %     fValues
    % end
    % 
    % properties (Access = private)
    % 
    % end
    % 
    % properties (Access = private)
    %     mesh
    %     ndimf
    %     order
    % end
    
    methods (Access = public)
        
        % function obj = TrialFunction(mesh, ndimf, order)
        %     obj.mesh  = mesh;
        %     obj.ndimf = ndimf;
        %     obj.order = order;
        %     fun = LagrangianFunction.create(obj.mesh, obj.ndimf, obj.order);
        %     obj.fValues = ones(size(fun.fValues));
        % end
        % 
        % function fxV = evaluate(obj, xV)
        %     fun = LagrangianFunction.create(obj.mesh, obj.ndimf, obj.order);
        %     fun.fValues = ones(size(fun.fValues));
        %     fxV = fun.evaluate(xV);
        % end
        
    end
    
    methods (Access = private)
        
    end
    
end