function SpatialSVD_3D(xx,PHI,INDICES_ECM_POINTS,TOL_SVD,LABELloc,NFIG)

if nargin == 0
    load('tmp.mat')
end

% % PHI IS A VECTOR THAT ARISES FROM reshaping the points of a cartesian mesh
% % NX x NY x NZ. By default, the cartesian points are stored such that the
% % tridimensional array is   Y-X-Z. Thefore, in reshaping, to convert it
% % into a vector, the order is as explained in the following.
% % Suppose that you have  a matrix
% ny = 3; nx = 2; nz = 2;
% A = zeros(ny,nx,nz) ;
% iacum = 0 ;
% for ix = 1:nx
%     for iy = 1:ny
%         for iz = 1:nz
%             STRING = [num2str(iy),num2str(ix),num2str(iz)] ;
%             A(iy,ix,iz)  = str2num(STRING) ;
%         end
%     end
% end
% % If I now make
% Av = reshape(A,nx*ny*nz,[]) ;
% % I get
% %     111   211   311   121   221   321   112   212   312   122   222   322
% % If I now make
% Av_y_xz = reshape(Av,ny,[]) ; % we get
% % %
% % Av_y_xz =
% %
% %    111   121   112   122
% %    211   221   212   222
% %    311   321   312   322
% %  Therefore, an SVD will gives Av_y_xz = U(y)*S*V^T(xz)
% % Notice that V(xz) is a linear combination of rows. Hence, they have the
% % structure
% % B_xz = [ y11
% %   y21
% %   y12
% %   y22]
% %
% % A reshape(B_xz,nx,[]) will produce a matrix
% %  [y11  y12
% %   y21  y22]  = U2(x)*S2*V2(z)^T


% We have to rephape PHI into a matrix PHI_x_y,
% so that, upon application of the SVD, we get
%  U(x)*S*V^T(y) = PHI_x_y
PHI_y_xz = reshape(PHI,length(xx{2}),[]) ;
DATAsvd.RELATIVE_SVD = 1;
[Uy,S_y_xz,Vxz] = RSVDT(PHI_y_xz,TOL_SVD,[],0,DATAsvd) ;
 figure(NFIG)
subplot(3,1,1)
title(LABELloc)
hold on
xlabel('y')
ylabel('U(y)')
grid on
%
subplot(3,1,2)
hold on
xlabel('x')
ylabel('V(x)')
grid on
%
subplot(3,1,3)
hold on
xlabel('z')
ylabel('V(z)')
grid on

hy = [] ; LEGENDy = [] ;
hx = [] ; LEGENDx = [] ;


COLORSLOCy = rand(length(S_y_xz),3) ;
 iacum = 0; 
for ixz = 1:length(S_y_xz)
    
    subplot(3,1,1)
    
    idim = 2;
    setINDEX = INDICES_ECM_POINTS(:,idim) ;
    xxLOC = xx{idim}(setINDEX) ;
    hy(ixz) = plot(xx{idim},Uy(:,ixz),'Marker','.','Color',COLORSLOCy(ixz,:)) ;
    LEGENDy{ixz} = ['S_',num2str(ixz),'=',num2str(S_y_xz(ixz))] ;
    yyLOC = Uy(setINDEX,ixz) ;
    plot(xxLOC,yyLOC,'rx','MarkerSize',6) ;
    
    M_xz = reshape(Vxz(:,ixz),length(xx{1}),[]) ;
    [Ux,S_x_z,Vz]= RSVDT(M_xz,TOL_SVD,[],0,DATAsvd) ;
    COLORSLOCxz  = rand( length(S_x_z),3) ;
   
    
    for i_xz = 1:length(S_x_z)
        iacum = iacum +1 ; 
        subplot(3,1,2)
        
        idim = 1;
        setINDEX = INDICES_ECM_POINTS(:,idim) ;
        xxLOC = xx{idim}(setINDEX) ;
        hx(iacum) = plot(xx{idim},Ux(:,i_xz),'Marker','.','Color',COLORSLOCy(ixz,:)) ;
        LEGENDx{iacum} = ['S(',num2str(ixz),',',num2str(i_xz),')','=',num2str(S_y_xz(ixz)*S_x_z(i_xz))] ;
        yyLOC = Ux(setINDEX,i_xz) ;
      %  plot(xxLOC,yyLOC,'rx','MarkerSize',6) ;
        
         subplot(3,1,3)
        
        idim = 3;
        setINDEX = INDICES_ECM_POINTS(:,idim) ;
        xxLOC = xx{idim}(setINDEX) ;
        hz(iacum) = plot(xx{idim},Vz(:,i_xz),'Marker','.','Color',COLORSLOCy(ixz,:)) ;
        LEGENDz{iacum} = ['S(',num2str(ixz),',',num2str(i_xz),')','=',num2str(S_y_xz(ixz)*S_x_z(i_xz))] ;
        yyLOC = Vz(setINDEX,i_xz) ;
     %   plot(xxLOC,yyLOC,'rx','MarkerSize',6) ;
        
        
    end
    
    
    
end
legend(hy,LEGENDy)
legend(hx,LEGENDx)
legend(hz,LEGENDz)


%
% PHI_x_y = PHI_x_y' ;
% DATAsvd.RELATIVE_SVD = 1;
% figure(NFIG)
%
%
% subplot(2,1,1)
% hold on
% title(LABELloc)
% xlabel('x')
% ylabel('U(x)')
% grid on
% subplot(2,1,2)
% hold on
% xlabel('y')
% ylabel('V(y)')
% grid on
% COLORSLOC = rand(length(S),3) ;
% hx = [] ;    LEGENDx={} ;
% subplot(2,1,1)
% hold on
%
%
%
%
%
%
%
%
% legend(hx,LEGENDx) ;
% hy = [] ;    LEGENDy={} ;
% subplot(2,1,2)
% hold on
% idim = 2;
% setINDEX = INDICES_ECM_POINTS(:,idim) ;
% xxLOC = xx{idim}(setINDEX) ;
% for i = 1:length(S)
%     hy(i) = plot(xx{2},Vy(:,i),'Marker','.','Color',COLORSLOC(i,:)) ;
%     LEGENDy{i} = ['S_',num2str(i) ] ;
%     yyLOC = Vy(setINDEX,i) ;
%     plot(xxLOC,yyLOC,'rx','MarkerSize',6) ;
% end
% legend(hy,LEGENDy) ;
%
%
