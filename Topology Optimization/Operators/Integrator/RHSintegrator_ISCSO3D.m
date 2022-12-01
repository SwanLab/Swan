classdef RHSintegrator_ISCSO3D < handle

    properties (Access = private)
        ndofs
        neumann
    end
    
    methods (Access = public)

        function obj = RHSintegrator_ISCSO(cParams)
            obj.init(cParams);
        end

        function Fext = compute(obj)
            Fext = obj.computePunctualFext();
        end

    end

    methods (Access = private)

        function init(obj, cParams)
            obj.ndofs   = cParams.ndofs;
            obj.neumann = cParams.neumann;
        end

        function Fp = computePunctualFext(obj)
            %Compute Global Puntual Forces (Not well-posed in FEM)
            Fdata = obj.convertDataToDof();
            Fdof = Fdata(:,1);
            Fval = Fdata(:,2);
            Fp = zeros(obj.ndofs,1);
            if ~isempty(obj.neumann)
                Fp(Fdof) = Fval;
            end
        end

        function Fdata = convertDataToDof(obj)
            Fdata = [];
            for i = 1:height(obj.neumann)
                data = obj.neumann(i,:);
                inod = data(1);
                iunk = data(2);
                idof = obj.nod2dof(inod,iunk);
                Fdata = [Fdata; idof, data(3)];
            end

        end

        function idof = nod2dof(obj, inode, iunkn)
            ndimf = 6;
            idof(:,1)= ndimf*(inode - 1) + iunkn;
        end

    end
    
end