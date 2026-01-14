function [COORtransf,CentrF2,dP_1_GLO,dANG,R_G1_GLO,R_12_GLO]= ...
    DivideSpline3D(COORref,f,df,L,h,COORrefDUMMY,sKNOTS,DATAIN,ANG_TORS,DATAINloc)
% INPUTS
% ----------------------
% L - length slice, also chordal distance
% f  % f=f(x)... Curve function
% df  % df = df/dx
% nelem % Plotting curve
% xmin  % domain x
% xmax
% h  : Semi-height slice
% JAHO 16-Oct-2020. Copy of DivideSpline3D
% ----------------------------------
if nargin == 0
    load('tmp2.mat')
elseif nargin == 9
    DATAINloc = [] ;
end


figure(1)
hold on
xlabel('x')
ylabel('y')
ylabel('z')


% S1) HOW TO DECIDE HOW TO MAKE THE DISCRETIZATION ? Function f(s) is
% defined  in terms of a parameter s, such that f(0) gives the coordinates
% x and y of the initial point.  --- sKNOTS(i)

% INTERFACEcoorREL = [0 h
%     0  0
%     0 -h]' ;
%POINTS DEFINING THE FACE
dANG = cell(1,length(sKNOTS)-1) ; % RElative angle [dtheta_x,dtheta_y,dtheta_z]
dP_1_GLO =zeros(3,length(sKNOTS)-1) ;    % Relative translations [dX,dY,dZ]
R_G1_GLO = cell(length(sKNOTS)-1,1) ; % Rotation of FACE 1 of each slice
R_12_GLO = cell(length(sKNOTS)-1,1) ;
COORtransf = cell(1,length(sKNOTS)-1);

DATAIN = DefaultField(DATAIN,'ROTATION_INITIAL_PLANE','SVD_GIVEN') ;   % AXIS_Z
DATAINloc = DefaultField(DATAINloc,'TYPEslices',[]) ; 
DATAIN.TYPEslices = DATAINloc.TYPEslices ; 

ROTATION_INITIAL_PLANE = DATAIN.ROTATION_INITIAL_PLANE ;
R_G1_prev = [] ;
% LOOP OVER POINTS OF THE SPLINE CURVVE
for i = 1:length(sKNOTS)-1
    disp('*******************************+')
    disp(['Domain ',num2str(i)])
    disp('*******************************+')
    s = sKNOTS(i) ;
    P1_G = f(s) ;  % Initial point (point 1, or centroid face 1, in GLOBAL coordinates )
    t1_G = df(s) ;   % Vector tangent to the curve  (Global Coordinates )
    nt1_G = t1_G/norm(t1_G) ;  % Unit vector
    
    ISPLANE =0;
    
    
    
    if i==1  ||  ISPLANE == 1  % TEMPORARILY
        %%%% method for determining the orientation of the other two local
        %%%% axes ---the first one is always tangent to the curve
        switch ROTATION_INITIAL_PLANE
            case 'SVD_GIVEN'
                [R_G1,SingV,Vright]=svd(nt1_G) ;
            case 'AXIS_Z'
                R_G1_3 = [0,0,1]' ;
                R_G2_3 = cross(R_G1_3,nt1_G) ;
                R_G1 = [nt1_G,R_G2_3,R_G1_3] ;
                
            otherwise
                error('Option not implemented yet')
        end
    end
    %  R_G1
    
    
    R_G1_GLO{i} = R_G1;  %It maps vectors expressed in the
    
    % RS intrinsic to point 1 to its counterpart in the global system
    % -----------------------------------------------------------------------------------
    %%%% Second point ---Centroid of face 2
    s = sKNOTS(i+1) ;
    P2_G = f(s) ; % Coordingates
    t2_G = df(s) ; % Vector tangent to the curve at the end point
    nt2_G = t2_G/norm(t2_G) ;  % Unit vector
    % Expression of the nt2_G in the reference system (RS) intrinsic to
    % point 1
    nt2_1 = R_G1'*nt2_G ;
    dP_G = P2_G-P1_G ;
    dP_1 = R_G1'*dP_G ;   % Coordinates centroid 2 in the RS intrinsic to 1
    dP_1_GLO(:,i) = dP_1 ;
    
    
    % Matrix R12 that maps a vector expressed in the RS intrinsic to point
    % 2 into RS of point 1.
    METHOD_ROTATION_PARTIAL = 'Rotation_z_y';
    switch METHOD_ROTATION_PARTIAL
        case 'RODRIGUES_FORMULA'
            error('Option not available')
            %The rotation takes place around the axis
            eAXIS = cross([1,0,0]',ntFIN_1) ;
            % and the angle is equal to
            thetaROT = acos(ntFIN_1'*[1,0,0]');
            % Rodrigues formula ()
            R12 = RodriguesFormulaRotation(thetaROT,eAXIS);
        case 'Rotation_z_y'
            % ntFIN_1 -- > See SymbolicRotations.m
            %    R_12 = Ry*Rz ;
            %
            % [cos(ay)*cos(az), -cos(ay)*sin(az), -sin(ay)]
            % [        sin(az),          cos(az),        0]
            % [cos(az)*sin(ay), -sin(ay)*sin(az),  cos(ay)]
            % nt2_1 = R_12*nt2_2 = R_12*[1;0;0] =
            % cos(ay)*cos(az)
            %        sin(az)
            % cos(az)*sin(ay)
            
            %             method_1 = 0 ;
            %             if  method_1 == 1
            %                 az = asin(nt2_1(2)) ;  % Incremental rotation, local z axis
            %                 cosay = nt2_1(1)/cos(az) ;
            %                 ay = real(acos(cosay)) ;  % Incremental rotation, local y axis
            %                 ax = 0 ;   % Incremental rotation, local x axis, zero by default (but torsion rotation may be seamlessly included here)
            %                 dANG{i} = [ax;ay;az] ;
            %                 R_12 = [cos(ay)*cos(az), -cos(ay)*sin(az), -sin(ay)
            %                     sin(az),          cos(az),        0
            %                     cos(az)*sin(ay), -sin(ay)*sin(az),  cos(ay)] ;
            %             else
            
             
            sin_az = nt2_1(2)  ;
            %cos_az = sqrt(nt2_1(1)^2 + nt2_1(3)^2)  ;
            
            cos_az = sqrt(1-sin_az^2) ;
            
            sin_ay =  nt2_1(3)/cos_az ;
            %   cos_ay =  nt2_1(1)/cos_az ;
            cos_ay = sqrt(1-sin_ay^2) ;
            
            az = real(acosd(cos_az)) ;
            ay = real(acosd(cos_ay)) ;
            ax = 0 ;
            dANG{i} = [ax;ay;az] ;
            R_12 = [cos_ay*cos_az, -cos_ay*sin_az, -sin_ay
                sin_az,          cos_az,        0
                cos_az*sin_ay, -sin_ay*sin_az,  cos_ay] ;
            
            
            disp('Checking orthogonality')
            nnnn = norm(R_12'*R_12 - eye(3)) ;
            disp(['norm(R_12^T R_12 - ident )= ',num2str(nnnn)]) ;
            
            %      end
            
            
    end
    
    
    if ~isempty(ANG_TORS)
        % TORSION
        ax =  ANG_TORS(i+1)-ANG_TORS(i) ;
        dANG{i}(1) = ax ;
        Rx  = [1 0  0
            0 cos(ax) -sin(ax)
            0  sin(ax) cos(ax)] ;
        R_12 = R_12*Rx ;
    else
        ax = [] ;
    end
    
    
    R_12_GLO{i} = R_12 ;
    %     plot3([xini,xfin,zfin],[yini,yfin,zfin],'r*')
    %
    %     if i == 1
    %         R = [cos(ANGini) -sin(ANGini); sin(ANGini), cos(ANGini)] ;
    %
    %         INTERFACEcoor = R*INTERFACEcoorREL ;
    %         INTERFACEcoor(1,:) =  INTERFACEcoor(1,:) + xini ;
    %         INTERFACEcoor(2,:) =  INTERFACEcoor(2,:) + yini ;
    %         plot(INTERFACEcoor(1,:),INTERFACEcoor(2,:),'r','LineWidth',2)
    %
    %     end
    
    %     R = [cos(ANGfin) -sin(ANGfin); sin(ANGfin), cos(ANGfin)] ;
    %     INTERFACEcoor = R*INTERFACEcoorREL ;
    %     INTERFACEcoor(1,:) =  INTERFACEcoor(1,:) + xfin ;
    %     INTERFACEcoor(2,:) =  INTERFACEcoor(2,:) + yfin ;
    %     plot(INTERFACEcoor(1,:),INTERFACEcoor(2,:),'r','LineWidth',2)
    %RADIO_CURVATURA(i+1) = RADIUS(xfin) ;
    % Transformation that has to undergo the reference domain
    
    if isempty(DATAIN.TYPEslices)
          COORrefLOC = COORref ; 
    else
        COORrefLOC = COORref{DATAIN.TYPEslices(i)} ; 
    end
    
    %%% VARIABLE CROSS-SECTION AREA 
    % -------------------------------------
    VARCROSS = DATAINloc.VARIABLE_CROSS_SECTION ; 
    VARCROSS = DefaultField(VARCROSS,'YP',@(s,smax) (1)) ; 
    VARCROSS = DefaultField(VARCROSS,'YN',@(s,smax) (1)) ; 
    VARCROSS = DefaultField(VARCROSS,'ZP',@(s,smax) (1)) ; 
    VARCROSS = DefaultField(VARCROSS,'ZN',@(s,smax) (1)) ;
    
    smax = sKNOTS(end) ; 
     
    efact(1).yp =  VARCROSS.YP(sKNOTS(i),smax) ; 
    efact(2).yp =  VARCROSS.YP(sKNOTS(i+1),smax) ; 
    efact(1).yn =  VARCROSS.YN(sKNOTS(i),smax) ; 
    efact(2).yn =  VARCROSS.YN(sKNOTS(i+1),smax) ; 
    efact(1).zp =  VARCROSS.ZP(sKNOTS(i),smax) ; 
    efact(2).zp =  VARCROSS.ZP(sKNOTS(i+1),smax) ; 
    efact(1).zn =  VARCROSS.ZN(sKNOTS(i),smax) ; 
    efact(2).zn =  VARCROSS.ZN(sKNOTS(i+1),smax) ; 
    
    DATAcubic.efact = efact ; 

    
    [COORtransf{i},~,CentrF2{i}] = ...
        CubicTransfSpline3D(P1_G,P2_G,R_G1,R_12,COORrefLOC,COORrefDUMMY,ax,DATAcubic) ;
    % Therefore
    R_G2 = R_G1*R_12 ;   % See Rotation12.xoj
    R_G1 = R_G2
    
    %   [R_G1,SS,VV] =  svd(R_G1) ;
    
    %%%%%%%%%%
    
    
    %  Pini = Pfin ;
    
end


axis equal


figure(2)
hold on
subplot(2,1,1)
hold on
xlabel('slice')
ylabel('Bending angle Y (degrees)')
ANGULOS = cell2mat(dANG)
bar(ANGULOS(2,:))
subplot(2,1,2)
hold on
xlabel('slice')
ylabel('Bending angle Z (degrees)')
bar(ANGULOS(3,:))

figure(3)
hold on
subplot(3,1,1)
hold on
xlabel('slice')
ylabel('\Delta X (mm)')
bar(dP_1_GLO(1,:))
subplot(3,1,2)
hold on
xlabel('slice')
ylabel('\Delta Y (mm)')
bar(dP_1_GLO(2,:))

subplot(3,1,3)
hold on
xlabel('slice')
ylabel('\Delta Z (mm)')
bar(dP_1_GLO(3,:))


%
% figure(3)
% hold on
% xlabel('slice')
% ylabel('dX')
% bar(dX)
%
% figure(4)
% hold on
% xlabel('slice')
% ylabel('dY')
% bar(dY)

% figure(3)
% hold on
% xlabel('x')
% ylabel('Radius of curvature')
% plot(xVECT,RADIO_CURVATURA,'k*')
% xMED = (xVECT(2:end)+xVECT(1:end-1))/2 ;
% bar(xMED,abs(restimated))
% xGAUSS = 0.5*(x(2:end)+x(1:end-1)) ;
% W = (x(2:end)-x(1:end-1))'  ; % Full-order weights
% y = a*cos(b*xGAUSS) ;
% figure(1)
% plot(xGAUSS,y)
% xlabel('x')
% ylabel('y')
% dy =  -a*b*sin(b*xGAUSS) ;
% f = sqrt(1 + dy.^2 ) ;
% figure(2)
% hold on
% plot(xGAUSS,f)
% % Integral of f ...
% s = W'*f'