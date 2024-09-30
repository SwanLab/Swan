classdef EIFEMtesting < handle

    properties (Access = public)

    end

    properties (Access = private)

    end

    properties (Access = private)
        EIFEMfilename
        meshDomain
        boundaryConditions
        bcApplier
        material
        meshReference
        interfaceMeshReference
        meshSubDomain
        ninterfaces
        interfaceMeshSubDomain
        globalMeshConnecSubDomain
        interfaceMeshConnecSubDomain
        subDomainContact
        cornerNodes
        quad
        interfaceConnec
        locGlobConnec
        localGlobalDofConnec
        interfaceDof
        interfaceDom
        weight


        eigenModes
        Kmodal
        MmodalPrecond
        displacementFun
        LHS
        RHS

        EIFEMsolver
        refLHS
        KeifemContinuous
        EIFEMprojection
    end

    methods (Access = public)

       


        function obj = EIFEMtesting()
     
%             obj.createModelPreconditioning();
%             u = obj.solver2(LHS,RHS,refLHS);
 
   
            gP = GeneralPreconditioner();


            LHSf = @(x) LHS*x;
            RHSf = RHS;            

            Mid = @(r) r;
            MidOrth = @(r,A,z) z+0.3*(r-A(z));

            Meifem = @(r) gP.solveEIFEM(r);
            Milu = @(r) gP.applyILU(r);
            MiluCG = @(r) gP.ILUCG(r,LHSf);
            MgaussSeidel = @(r) gP.applyGaussSeidel(r);

         %  LHSf = @(x) P*LHS*x;            
            %  RHSf = P*RHS;
            tol = 1e-8;
            P = @(r) Mid(r); %obj.multiplePrec(r,Mid,Mid,LHSf);
            tic
            x0 = zeros(size(RHSf));
            [uCG,residualCG,errCG,errAnormCG] = PCG.solve(LHSf,RHSf,x0,P,tol,Usol);
            toc
            %[uCG,residualCG,errCG,errAnormCG] = RichardsonSolver.solve(LHSf,RHSf,x0,P,tol,0.1,Usol);
            
            
            M = MiluCG;%Milu_m;%Meifem; %Milu %Pm
            M2 = Meifem;
            M3 = MiluCG;
%             [uPCG,residualPCG,errPCG,errAnormPCG] = obj.solverTestEifem(LHSf,RHSf,Usol,M);
            tol = 1e-8;
            P = @(r) gP.multiplePrec(r,M,M2,M3,LHSf);  
%            P = Milu;
%              P = @(r) obj.additivePrec(r,Mid,Mmodal,LHSf);  
            tic
            x0 = zeros(size(RHSf));            
            [uPCG,residualPCG,errPCG,errAnormPCG] = PCG.solve(LHSf,RHSf,x0,P,tol,Usol);
            toc
            %[uCG,residualCG,errCG,errAnormCG] = RichardsonSolver.solve(LHSf,RHSf,x0,P,tol,0.1,Usol);

            figure
            plot(residualPCG,'linewidth',2)
            hold on
            plot(residualCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend({'CG + ILU-EIFEM-ILU','CG'},'FontSize',12)
            xlabel('Iteration')
            ylabel('Residual')

             figure
            plot(errPCG,'linewidth',2)
            hold on
            plot(errCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend('CG + EIFEM+ ILU(CG-90%-L2)','CG')
            xlabel('Iteration')
            ylabel('||error||_{L2}')

            figure
            plot(errAnormPCG,'linewidth',2)
            hold on
            plot(errAnormCG,'linewidth',2)
            set(gca, 'YScale', 'log')
            legend('CG + EIFEM+ ILU(CG-90%-L2)','CG')
            xlabel('Iteration')
            ylabel('Energy norm')





        end

 
    end

end
