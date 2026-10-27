classdef TrainedRVE < handle

    properties (GetAccess = public , SetAccess = private)
        ndimf
        Kcoarse
        Tdef
        Trb
        T
        PhiDef
        PhiRb
        Grb
        dKcoarse
        dUrb
        dUdef
        dU
    end

    properties (Access = private)


    end

    properties (Access = private)

    end

    methods (Access = public)

        function obj = TrainedRVE(filename)
            obj.init(filename)

        end

    end

    methods (Access = private)

        function init(obj,data)
            if ischar(data)
                load(data);
                obj.Kcoarse = EIFEoper.Kcoarse;
                obj.Tdef    = EIFEoper.Tdef;
                obj.Trb     = EIFEoper.Trb;
                if isfield(EIFEoper,'U')
                    obj.T = EIFEoper.T;
                end
                if isfield(EIFEoper,'dKcoarse')
                    obj.dKcoarse = EIFEoper.dKcoarse;
                end
                if isfield(EIFEoper,'dUdef')
                    obj.dUdef = EIFEoper.dUdef;
                end
                if isfield(EIFEoper,'dUrb')
                    obj.dUrb = EIFEoper.dUrb;
                end
                if isfield(EIFEoper,'dU')
                    obj.dU = EIFEoper.dU;
                end
            else
                obj.Kcoarse = data.Kcoarse;
                obj.Tdef    = data.Tdef;
                obj.Trb     = data.Trb;
            end

            obj.ndimf   = 2;
            %             obj.PhiDef  = EIFEoper.RECONSTRUCTION.DEF_DISP.BASIS;
            %             obj.PhiRb   = EIFEoper.RECONSTRUCTION.RB_DISP.BASIS;
            %             obj.Grb     = EIFEoper.Grb;
        end

    end

end