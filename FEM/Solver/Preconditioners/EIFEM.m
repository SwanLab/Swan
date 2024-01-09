classdef EIFEM < handle
    
    properties (Access = public)
        
    end
    
    properties (Access = private)
        RVE
        mesh
        
    end
    
    properties (Access = private)
        LHS
        Kel
    end
    
    methods (Access = public)
        
        function obj = EIFEM(cParams)
            obj.init(cParams)
            LHS = obj.computeLHS();
            obj.LHS = LHS;

        end

        function x = apply(obj,r)
            x=obj.D\r;
        end
        
    end
    
    methods (Access = private)
        
        function init(obj,cParams)
            obj.mesh = cParams.mesh;
            obj.RVE  = cParams.RVE;
            obj.Kel  = repmat(obj.RVE.Kcoarse,[1,1,obj.mesh.nelem]);
        end

        function dim = getDims(obj)
            d.ndimf     = obj.RVE.ndimf;
            d.nnodes    = size(obj.mesh.coord, 1);
            d.ndofs     = d.nnodes*d.ndimf;
            d.nnodeElem = obj.mesh.nnodeElem; % should come from interp..
            d.ndofsElem = d.nnodeElem*d.ndimf;
            dim = d;
        end

        function LHS = computeLHS(obj)
            s.dim          = obj.getDims();
            s.nnodeEl      = obj.mesh.nnodeElem;
            s.globalConnec = obj.mesh.connec;
            assembler = Assembler(s);
            LHS = assembler.assemble(obj.Kel);
        end
        
    end
    
end