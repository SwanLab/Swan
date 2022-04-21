classdef LHSintegrator_Bending < LHSintegrator
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        
    end
    
    properties (Access = private)
        
    end
    
    methods (Access = public)
        
        function obj = LHSintegrator_Bending(cParams)
            obj.init(cParams)
            
        end

    end

    methods (Access = protected)

        function lhs = computeElementalLHS(obj)
            d = obj.dim;
            E = obj.youngModulus;
            I = obj.inertiaMoment;
            nElem = d.nelem;
            Edof = d.nnode*d.nstre;
            Be = zeros(Edof ,Edof ,nElem);
            l = obj.computeLength();
            for iElem = 1:nElem
                length = l(iElem);
                [c1,c2,c3,c4] = obj.coeffsBending(length,E,I);
                Be(1,1:4,iElem) = c1*[c2 c3 -c2 c3];
                Be(2,1:4,iElem) = c1*[c3 c4 -c3 c4/2];
                Be(3,1:4,iElem) = c1*[-c2 -c3 c2 -c3];
                Be(4,1:4,iElem) = c1*[c3 c4/2 -c3 c4];
            end
            lhs  = Be;
        end


    end

    methods (Access = private)

        
    end
    
end