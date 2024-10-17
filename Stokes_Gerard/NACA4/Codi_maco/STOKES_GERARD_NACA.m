clear
close all

O = 1;

% for MM = 0.5:1:9.5
    H = 1;
%     M = MM/100;
%     for pp = 2.5:1:7.5
%         p = pp/10;

        % % INPUT DATA
        m = QuadMesh(10,4,100,40); % MESH

        % % NACA 4
        M=9/100;
        p=8/10;
        t=12/100;

        % Biga (posada en el centre de màx t):
        alt = 0.11;
        ampl = 0.095;
        x_pos = 0.3; %Centre de la biga respecte el LE

        AOAd = 10; %deg
        x_centr = 3.5;
        y_centr = 2;

        %% Airfoil creation

        % [t,Q] = findthickness(x_pos,ampl,alt,M,p);


        fH = Find_fH_circles(M,p,t,x_centr,y_centr,AOAd);

        % fH = Find_fH_points(M,p,t,x_centr,y_centr,AOAd);

        %% Create mesh

        s.fHandle = fH;
        s.type='Given';
        g = GeometricalFunction(s);
        lsFun = g.computeLevelSetFunction(m); %D'aquí surt la malla de quadrats sense el forat
        sUm.backgroundMesh = m;
        sUm.boundaryMesh = m.createBoundaryMesh(); %sUm.boundaryMesh conté les mesh de les quatre fronteres del voltant. No té res del forat
        uMesh = UnfittedMesh(sUm);
        uMesh.compute(lsFun.fValues); % uMesh.boundaryCutMesh.mesh  és el forat
        mesh = uMesh.createInnerMesh();
        % figure
        % plot(uMesh)
        % hold on
        % scatter(punts_rot(1,:,1),punts_rot(2,:,1))
        %figure
        %plot(lsFun)
        e.type  = 'STOKES';
        e.nelem = mesh.nelem;
        material = Material.create(e);
        dtime = Inf; %Estacionari

        % VELOCITY AND PRESSURE FUNCTIONS
        velocityFun = LagrangianFunction.create(mesh, 2, 'P2');
        pressureFun = LagrangianFunction.create(mesh, 1, 'P1');
        n_dofs = velocityFun.nDofs + pressureFun.nDofs;
        
        %% Boudary conditions

        [forcesFormula,dirichlet,dir_dofs,nodespresscyl] = boundary_conditions(mesh,uMesh,velocityFun,pressureFun);
        %[forcesFormula,dirichlet,dir_dofs,nodespresscyl] = boundary_conditions_parab(mesh,uMesh,velocityFun,pressureFun);


        %% Solver

        [velocityFun,pressureFun] = solver_stokesG(forcesFormula,dirichlet,dir_dofs,velocityFun,pressureFun,dtime,mesh,material,n_dofs);


        %% PLOT RESULTS
        velocityFun.plot()
        pressureFun.plot()
        caxis([-50 50]);


        %% Lift and drag

        [Li,Di] = aero_forces(nodespresscyl,pressureFun,mesh);

        L(H,O) = Li;
        D(H,O) = Di;
        tt(H,O) = t;

%         clearvars('-except', 'time','H','D','L','O','M','p','MM','pp','tt');
%         H=H+1;
% 
%     end
%     O=O+1;
%     disp(O);
% end