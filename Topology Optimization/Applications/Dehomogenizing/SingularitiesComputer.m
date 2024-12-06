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

        function f = compute(obj)
            aC = obj.orientation;
            aD = project(aC,'P1D');
            aD = permute(aD.getFvaluesByElem(), [1 3 2]);
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
            s.order   = 'P0';
            s.mesh    = obj.mesh;
            f = LagrangianFunction(s);
        end



    end

    methods (Access = private)

       function init(obj,cParams)
            obj.mesh        = cParams.mesh;
            obj.orientation = cParams.orientation;
        end

    end

 

end