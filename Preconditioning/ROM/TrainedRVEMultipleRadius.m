classdef TrainedRVEMultipleRadius < handle

    properties (Access = public)
        ndimf
        Kcoarse
        Udef
        Urb
        T
        Mcoarse
        PhiDef
        PhiRb
        Grb
    end

    properties (Access = private)


    end

    properties (Access = private)

    end

    methods (Access = public)

        function obj = TrainedRVEMultipleRadius(radius,mesh)
            obj.init(radius,mesh)

        end

    end

    methods (Access = private)

        function init(obj,radius,mesh)
            K = cell(size(radius)); M = K; T = M;
            for i = 1:size(radius,1)
                for j=1:size(radius,2)
                    [K{i,j}, M{i,j}, T{i,j}] = RebuildKMTData.compute('A', radius(i,j), mesh);
                end

            end

            obj.Kcoarse = K;
            obj.T = T;
            obj.Mcoarse = M;

            obj.ndimf   = 2;

        end

    end

end