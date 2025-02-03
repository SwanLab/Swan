clear
close all

O = 1;

% for MM = 0.5:1:9.5
    H = 1;
%     M = MM/100;
%     for pp = 2.5:1:7.5
%         p = pp/10;

        % % INPUT DATA
        m = TriangleMesh(2,1,150,75);
        %QuadMesh(10,4,100,40); % MESH de referència sense cap objecte, quadrilateral

        % % NACA 4
        M = 0.02;
        %9/100;
        p = 0.5;
        %8/10;
        t=12/100;

        % Biga (posada en el centre de màx t):
        % alt = 0.11;
        % ampl = 0.095;
        % x_pos = 0.3; %Centre de la biga respecte el LE

        AOAd = 0;
        %10; %deg
        x_centr = 1;
        %3.5;
        y_centr = 0.5;
        %2;

        %% Airfoil creation

        % [t,Q] = findthickness(x_pos,ampl,alt,M,p);


        fH = Find_fH_circles(M,p,t,x_centr,y_centr,AOAd); % Troba la funció 

        % fH = Find_fH_points(M,p,t,x_centr,y_centr,AOAd);

        %% Create mesh

        s.fHandle = fH;
        s.type='Given';
        g = GeometricalFunction(s); %fabrica de diferents casos, cicles, en aquest cas donada l'equació
        lsFun = g.computeLevelSetFunction(m); %D'aquí surt la funció de superficie
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
        material = Material.create(e); % definir les propietats físiques del fluid, com la viscositat dinàmica
        dtime = Inf; %Estacionari 

        % VELOCITY AND PRESSURE FUNCTIONS
        velocityFun = LagrangianFunction.create(mesh, 2, 'P2'); % interpreto que aqui el que fa es crear trial functions, és a dir, crear el shape function amb l'ordre que toca per solucions
        pressureFun = LagrangianFunction.create(mesh, 1, 'P1');
        n_dofs = velocityFun.nDofs + pressureFun.nDofs; % HE SALTAT AIXO

        % Començo per aquí:
        
        %% Boudary conditions

        [forcesFormula,dirichlet,dir_dofs,nodespresscyl] = boundary_conditions(mesh,uMesh,velocityFun,pressureFun);
        %[forcesFormula,dirichlet,dir_dofs,nodespresscyl] = boundary_conditions_parab(mesh,uMesh,velocityFun,pressureFun);


        %% Solver

        [velocityFun,pressureFun] = solver_stokesG(forcesFormula,dirichlet,dir_dofs,velocityFun,pressureFun,dtime,mesh,material,n_dofs);


        %% PLOT RESULTS
        velocityFun.plot() %Surt error al realitzar l'últim pas d'edge amv step in
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

close all;

save("datas.mat");