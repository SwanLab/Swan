classdef SingularitiesComputer < handle

   properties (GetAccess = public, SetAccess = private)  
        nSing
        isElemSingular
    end      

    properties (Access = private)
        orientation
        mesh
    end

    methods (Access = public)

        function obj = SingularitiesComputer(cParams)
            obj.init(cParams)
        end

        function compute(obj)
            obj.computeSingularElements();
            obj.computeNumberOfSingularities();
        end

        function plot(obj)
            obj.isElemSingular.plot();
        end

    end

    methods (Access = private)

       function init(obj,cParams)
            obj.mesh        = cParams.mesh;
            obj.orientation = cParams.orientation;
        end

        function computeSingularElements(obj)
            aC = obj.orientation;
            aD = aC.project('P1D');
            aD = permute(aD.fValues, [1 3 2]);
            a1 = zeros(3,obj.mesh.nelem);
            a2 = zeros(3,obj.mesh.nelem);
            a3 = zeros(3,obj.mesh.nelem);
            a1(1:2,:) = aD(:,:,1);
            a2(1:2,:) = aD(:,:,2);
            a3(1:2,:) = aD(:,:,3);
                      
            a1a2 = dot(a1,a2);
            a1a3 = dot(a1,a3);
            a2a3 = dot(a2,a3);
            
            isS = sign(a1a2.*a1a3.*a2a3)';
            s.fValues = isS<0;
            s.mesh    = obj.mesh;
            f = P0Function(s);
            obj.isElemSingular = f;
        end

        function computeNumberOfSingularities(obj)
            obj.nSing = sum(obj.isElemSingular.fValues);
        end

    end

 

end