classdef Element_Elastic_Micro < Element_Elastic
    %Element_Elastic Summary of this class goes here
    %   Detailed explanation goes here
    % !! CONSIDER TO IMPLEMENT A CONSTRUCTOR THAT DEFINES B & C DIMENS AT
    % THE PRE-PROCESS !!
    
    properties
    end
    
    methods (Access = ?Physical_Problem)
        function obj = computeRHS(obj,nunkn,nelem,nnode,bc,idx)  
            computeRHS@Element(obj,nunkn,nelem,nnode,bc,idx);
            RHSStrain = obj.computeStrainRHS(nunkn,nelem,nnode,bc,idx,vstrain);
            obj.RHS = obj.RHS + RHSStrain;
            %             RHSPuntual = obj.computePuntualRHS(nunkn,nelem,nnode,bc,idx);
            %             RHSSuperficial  = obj.computeSuperficialRHS(nunkn,nelem,nnode,bc,idx);
            %             RHSVolumetric  = obj.computeVolumetricRHS(nunkn,nelem,nnode,bc,idx);
            %             RHSStrain = obj.computeStrainRHS(nunkn,nelem,nnode,bc,idx);
            %             obj.RHS = RHSSuperficial + RHSVolumetric + RHSPuntual + RHSStrain;
        end
    end
    
    methods (Access = private)
        function F = computeStrainRHS(nunkn,nelem,nnode,bc,idx,vstrain)
            nstre=3;
            F = zeros(nnode*nunkn,1,nelem);
            tstrain = zeros(nstre,ngaus,nstre,nelem);
            tstres = zeros(nstre,ngaus,nstre+1,nelem);
            stre0=zeros(nstre,nelem);
            for istre=1:nstre
                for jstre=1:nstre
                    stre0(istre,:) = stre0(istre,:) + squeeze(Ce(istre,jstre,:)*stra0(jstre))';
                end
            end
            for iv=1:nnode*nunkn
                ivBmat(1:nstre,1:nelem) = Bmat(1:nstre,iv,1:nelem);
                for istre=1:nstre
                    eforce(iv,:)=eforce(iv,:)+ivBmat(istre,:).*stre0(istre,:).*dvolu';
                end
            end
            F= eforce;
        end
    end
    
end
