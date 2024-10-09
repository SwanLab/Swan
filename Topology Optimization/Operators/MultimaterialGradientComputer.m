classdef MultimaterialGradientComputer < handle

    properties (Access = private)
        mesh
        designVariable
    end

    methods (Access = public)
        function obj = MultimaterialGradientComputer(cParams)
            obj.init(cParams);
        end

        function dt = compute(obj,topDers)
            x = obj.designVariable;
            TD = topDers;

            charfun = x.obtainDomainFunction();
            [~,tfiFun] = charfun.computeAtNodesAndElements();
            tfi = tfiFun.fValues';
            psi2 = x.levelSets{1,2}.fun.fValues;
            psi3 = x.levelSets{1,3}.fun.fValues;

            t = obj.mesh.connec';
            p = obj.mesh.coord';
            [tXi2,~] = integ_exact(t,p,psi2); chi2 = (1 - tXi2); %- Mixed formulation method
            [tXi3,~] = integ_exact(t,p,psi3); chi3 = (1 - tXi3); %- Mixed formulation method
            %     fi = (pdeintrp(p,t,fi)).'; % interpolation at gauss point - P1 projection method
            %     chi2 = (pdeintrp(p,t,(psi(:,2)<0))).'; %- P1 projection method

            dt = [];
            dt(1,:) = - tfi(1,:).*TD{1,end} - tfi(2,:).*TD{2,end} - tfi(3,:).*TD{3,end} ...
                + tfi(4,:).*( (1-chi2).*TD{end,1} + (1-chi3).*chi2.*TD{end,2} + chi2.*chi3.*TD{end,3} );

            dt(2,:) = - tfi(2,:).*TD{2,1} - tfi(3,:).*TD{3,1} + tfi(1,:).*( (1-chi3).*TD{1,2} + chi3.*TD{1,3} );

            dt(3,:) = tfi(2,:).*TD{2,3} - tfi(3,:).*TD{3,2};
        end
    end

    methods (Access = private)
        function init(obj,cParams)
            obj.mesh           = cParams.mesh;
            obj.designVariable = cParams.designVariable;
        end
    end
end