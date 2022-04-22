classdef PST4VoigtFourthOrderTensorSymbolically < PST4VoigtFourthOrderTensor
    
    properties (Access = protected)
        TensorInPlaneStress
    end
    
    properties (Access = private)
        tensorVoigt
        
        Dim
        VoigtDim
        IndexTransformer
        
        OutOfPlaneIndex
        InPlaneIndex
        
        SymbolicStrains
        InPlaneSymbolicStrains
        OutOfPlaneSymbolicStrains
        
        Strains
        VoigtStrains
        
        Stresses
        InPlaneStresses
        OutOfPlaneStresses
    end
    
    methods (Access = public)
        
        function obj = PST4VoigtFourthOrderTensorSymbolically(tensorVoigt)
            obj.createPlaneStressTensor();
            obj.createVariables(tensorVoigt)
            obj.computeOutOfPlaneStrains()
            obj.updateVariables()
            obj.computeTensorVoigtInPlaneStress()
            obj.storePsTensor()
        end
        
        function t = getValue(obj)
            t = obj.TensorInPlaneStress;
        end
    end
    
    methods (Access = private)
        
        function createVariables(obj,tensorVoigt)
            obj.tensorVoigt = tensorVoigt;
            obj.IndexTransformer = TensorVoigtIndexTransformer3D();
            obj.VoigtDim = 6;
            obj.Dim = 3;
            obj.createTensorVoigtInPlaneStress()
            obj.createOutOfPlaneIndex()
            obj.createInPlaneIndex()
            obj.createStrains()
            obj.createSymbolicStrains()
            obj.computeVoigtStrains()
            obj.computeStresses()
        end
        
        function createTensorVoigtInPlaneStress(obj)
            obj.TensorInPlaneStress = sym(zeros(obj.Dim,obj.Dim));
        end
        
        function createOutOfPlaneIndex(obj)
            obj.OutOfPlaneIndex(1) = obj.IndexTransformer.tensor2Voigt(3,3);
            obj.OutOfPlaneIndex(2) = obj.IndexTransformer.tensor2Voigt(2,3);
            obj.OutOfPlaneIndex(3) = obj.IndexTransformer.tensor2Voigt(1,3);
        end
        
        function createInPlaneIndex(obj)
            obj.InPlaneIndex(1) = obj.IndexTransformer.tensor2Voigt(1,1);
            obj.InPlaneIndex(2) = obj.IndexTransformer.tensor2Voigt(2,2);
            obj.InPlaneIndex(3) = obj.IndexTransformer.tensor2Voigt(1,2);
        end
        
        function  createStrains(obj)
            obj.Strains =  sym('eps',[obj.VoigtDim 1],'real');
        end
        
        function createSymbolicStrains(obj)
            obj.SymbolicStrains           = obj.Strains;
            obj.InPlaneSymbolicStrains    = obj.SymbolicStrains(obj.InPlaneIndex);
            obj.OutOfPlaneSymbolicStrains = obj.SymbolicStrains(obj.OutOfPlaneIndex);
        end
        
        function computeStrains(obj)
            obj.Strains(obj.InPlaneIndex,:)    = obj.InPlaneSymbolicStrains;
            obj.Strains(obj.OutOfPlaneIndex,:) = obj.OutOfPlaneSymbolicStrains;
        end
        
        function updateVariables(obj)
            obj.computeStrains()
            obj.computeVoigtStrains()
            obj.computeStresses()
        end
        
        function computeVoigtStrains(obj)
            HydroStaticIndex = 1:3;
            ShearIndex = 4:6;
            obj.VoigtStrains = sym(zeros(obj.VoigtDim,1));
            obj.VoigtStrains(HydroStaticIndex) =   obj.Strains(HydroStaticIndex);
            obj.VoigtStrains(ShearIndex)       = 2*obj.Strains(ShearIndex);
        end
        
        function computeStresses(obj)
            Ctensor = obj.tensorVoigt;
            strain  = obj.VoigtStrains;
            stress  = Ctensor*strain;
            obj.Stresses = simplify(stress);
            obj.InPlaneStresses    = obj.Stresses(obj.InPlaneIndex);
            obj.OutOfPlaneStresses = obj.Stresses(obj.OutOfPlaneIndex,1);
        end
        
        function computeOutOfPlaneStrains(obj)
            OutStresses = obj.OutOfPlaneStresses;
            OutStrains  = obj.OutOfPlaneSymbolicStrains;
            OutStrains  = solve(OutStresses == 0,OutStrains);
            oS = struct2cell(OutStrains);
            obj.OutOfPlaneSymbolicStrains = [oS{:}];
        end
        
        function computeTensorVoigtInPlaneStress(obj)
            for iStress = 1:length(obj.InPlaneIndex)
                for iStrain = 1:length(obj.InPlaneIndex)
                    Stress = obj.InPlaneStresses(iStress);
                    Strain = obj.InPlaneSymbolicStrains(iStrain);
                    [TensorValue,StrainValue] = coeffs(Stress,Strain);
                    obj.fillTensorComponent(iStress,iStrain,TensorValue,StrainValue,Strain)
                end
            end
        end
        
        function fillTensorComponent(obj,i,j,TensorValue,StrainValue,Strain)
            
            for icoef = 1:length(StrainValue)
                if StrainValue(icoef) == Strain
                    if ~isempty(TensorValue(1))
                        if obj.isThirdComponent(j)
                            obj.TensorInPlaneStress(i,j) = 0.5*TensorValue(icoef);
                        else
                            obj.TensorInPlaneStress(i,j) = TensorValue(icoef);
                        end
                    end
                end
            end
        end
        
        function storePsTensor(obj)
            t = obj.TensorInPlaneStress;
            if ~obj.isInputNotSymbollic()
               t = double(t);
            end
            obj.psTensor.setValue(t);
        end
        
        function itIsNotSym = isInputNotSymbollic(obj)
            itIsNotSym = isa(obj.tensorVoigt,'sym');
        end
    end
    
    methods (Access = private, Static)
        
        function isThirdComponent = isThirdComponent(Component)
            isThirdComponent = Component == 3;
        end
        
    end
    

    
end

