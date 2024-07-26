classdef Circle < handle
    properties (Access = public)
        bx
        by
        fx
        fy
    end

    properties (Access = private)
        cx
        cy
        color
        r
    end

    methods (Access = public)
        function self = Circle(prop)
            self.cx = prop.cx;
            self.cy = prop.cy;
            self.color = prop.color;
            self.r = prop.r;
            self.compute_b_f();            
        end

        function plot(self)
            x1 = self.cx + self.r * 2^0.5 * cos(225*pi/180);
            y1 = self.cy + self.r * 2^0.5 * sin(225*pi/180);
            pos = [x1 y1 2*self.r 2*self.r];
            rectangle('Position',pos,'Curvature',[1 1],'FaceColor',self.color)
        end    
    end

    methods (Access = private)
        function compute_b_f(self)
            self.bx = self.cx - self.r;
            self.by = self.cy;
            self.fx = self.cx + self.r;
            self.fy = self.cy;
        end
    end
end