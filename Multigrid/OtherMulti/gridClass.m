% Class 'gridClass' implements a grid in square area of plane
classdef gridClass < handle
    properties (Access = public)        
        % central point
        center = [];
        % length of edge
        edge_len = [];
        
        % a square can be set by two point of plane:
        % a - a left and down vertex,
        % b - a right and up  vertex
        a = [];
        b = [];
        
        % number of inner nodes by ONE direction
        N = 0;
        % step of grid
        h = 0;
        
        % vectors of grid nodes
        x = [];
        y = [];
    end
    
    methods (Access = public)
        % PARAMS IN:
        %   * c   - center of square,
        %   * e_l - edge length,
        %   * N   - number of inner nodes by ONE direction
        function obj = gridClass(c, e_l, N)
            obj.center = c;
            obj.edge_len = e_l;            
            
            obj.N = 2^(ceil(log2(N+1))) - 1;
            
            obj.a = [obj.center(1), obj.center(2)] - obj.edge_len/2;
            obj.b = [obj.center(1), obj.center(2)] + obj.edge_len/2;
                        
            obj.h = e_l/(N+1);            
            
            obj.x = obj.a(1):obj.h:obj.b(1);            
            obj.y = obj.a(2):obj.h:obj.b(2);
        end
        
        % func returns points on left edge
        function [x, y] = leftEdge(obj)
            x = obj.a(1);
            y = obj.y;
        end
        
        % func returns points on right edge
        function [x, y] = rightEdge(obj)
            x = obj.b(1);
            y = obj.y;
        end
        
        % func returns points on bottom edge
        function [x, y] = bottomEdge(obj)
            x = obj.x;
            y = obj.a(2);
        end
        
        % func returns points on top edge
        function [x, y] = topEdge(obj)
            x = obj.x;
            y = obj.b(2);
        end
        
        % func returns a grid which has a double step
        function grid_2h = getGrid2h(obj)
            N_half = round((obj.N+1)/2) - 1;
            if N_half >= 0
                grid_2h = gridClass(obj.center, obj.edge_len, N_half);   
            else
                disp('gridClass::getGrid2h::ERROR: wrong size of grid')
                grid_2h = [];
            end
        end
        
        % func returns a grid which has step = 2^(level-1)*h
        function grid_m = getSubgrid(obj, level)
            % first level is init grid
            grid_m = obj;
            for i = 2:level
                grid_m = grid_m.getGrid2h();
            end            
        end
    end
end
