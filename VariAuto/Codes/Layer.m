classdef Layer < handle 
    properties (Access = public)
        theta
    end

    properties (Dependent)
        W
        b
    end

    properties (Access = private)
        prevL
        nextL
    end

    methods (Access = public)
        function self = Layer(thv,prevL,nextL)
            self.theta = thv;
            self.prevL = prevL;
            self.nextL = nextL;
        end
    end

    methods
        function value = get.W(self)
            aux = self.theta(1:self.prevL*self.nextL);
            value = reshape(aux,[self.prevL,self.nextL]);
        end
        function value = get.b(self)
            value = self.theta(self.prevL*self.nextL+1:end);
        end
    end
end