classdef SingularitiesComputer < handle

    properties (Access = private)
        meshDisc
        isElemSingular
        singularitiesCoord
    end    

    properties (Access = private)
        orientation
        mesh
    end

    methods (Access = public)

        function obj = SingularitiesComputer(cParams)
            obj.init(cParams)
        end

        function sC = compute(obj)
            obj.computeSingularElements();
            obj.computeSingluaritiesCoord();
            sC = obj.singularitiesCoord;
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
            aD = permute(aD.fValues, [3 1 2]);
            a1 = aD(:,:,1);
            a2 = aD(:,:,2);
            a3 = aD(:,:,3);
                      
            a1a2 = obj.scalarProduct(a1,a2);
            a1a3 = obj.scalarProduct(a1,a3);
            a2a3 = obj.scalarProduct(a2,a3);

            aV(1,:,1) = a2a3;
            aV(1,:,2) = a1a3;
            aV(1,:,3) = a1a2;
            
            s.fValues = aV;
            s.mesh    = obj.mesh;
            f = P1DiscontinuousFunction(s);
            f.plot()

            isS = sign(a1a2.*a1a3.*a2a3);

            s.fValues = isS<0;
            s.mesh    = obj.mesh;
            f = P0Function(s);
            obj.isElemSingular = f;
        end

        function computeSingluaritiesCoord(obj)
            isS = obj.isElemSingular.fValues;
            coord = obj.mesh.computeBaricenter();
            coord = transpose(coord);
            sC    = coord(isS,:); 
            obj.singularitiesCoord = sC;
        end

    end

    methods (Access = private, Static)

        function ab = scalarProduct(a,b)
            ab = a(:,1).*b(:,1) + a(:,2).*b(:,2);
        end

    end

end