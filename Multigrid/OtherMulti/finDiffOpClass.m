% Class 'finDiffOpClass' implements a finite difference operator that
% approximates an differential operator in left part of elliptic equation:
% - d/dx (a(x, y) du/dx) - d/dy (b(x, y) du/dy) + q(x, y) u = f(x, y)
classdef finDiffOpClass < handle
    properties (Access = private)
        % coefficients of equation 
        a_coeff = [];
        b_coeff = [];
        q_coeff = [];       
    end
        
    methods (Access = public)
        function obj = finDiffOpClass(a_coeff, b_coeff, q_coeff)
            obj.a_coeff = a_coeff;
            obj.b_coeff = b_coeff;
            obj.q_coeff = q_coeff;            
        end
        
        % func returns matrix of finite difference operator,
        % approximation by finite volumes
        function A_mat = getMatrix(obj, grid_m) 
            h = grid_m.h;
            N = grid_m.N;
            x = grid_m.x;
            y = grid_m.y;
            
            a = (obj.a_coeff(x(2:N+1)+h/2, y(2:N+1))+...
                 obj.a_coeff(x(2:N+1)-h/2, y(2:N+1)))/h^2+...
                (obj.b_coeff(x(2:N+1),       y(2:N+1)+h/2)+...
                 obj.b_coeff(x(2:N+1),       y(2:N+1)-h/2))/h^2+...
                 obj.q_coeff(x(2:N+1),       y(2:N+1));

            b = -obj.a_coeff(x(2:N+1)+h/2, y(2:N+1))/h^2;

            c = -obj.b_coeff(x(2:N+1),       y(2:N+1)-h/2)/h^2;

            d = -obj.a_coeff(x(2:N+1)-h/2, y(2:N+1))/h^2;

            e = -obj.b_coeff(x(2:N+1),       y(2:N+1)+h/2)/h^2;
            
            A_mat = matOfFinDiffOpClass(a, b, c, d, e);
        end
    end    
end

