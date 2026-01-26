function [GeneralizedForces]=...
    ForcesMoments_diagrams(DATAROM,MESH1D,rDEF,rRB,DATA_REFMESH,ndimINTF,DATAIN) ;


if nargin ==0
    load('tmp1.mat')
end

DATAIN  =DefaultField(DATAIN,'TIME_STEP_LOCAL',[]) ; 

%V= V(:,1:6) ;
%
% if   length(DATAROM) == 1
%     GeneralizedForces = generalized_forces_onlyslices(DATA_REFMESH{1},DATAROM{1},rDEF,rRB,ndimINTF)  ;
% else
[GeneralizedForces,RotationMatrixFace2 ]= generalized_forces_jslices(DATA_REFMESH,DATAROM,rDEF,rRB,MESH1D,DATAIN)  ;
%end

if ~isempty(RotationMatrixFace2)
    % Polar coordinates (angle )
    X = atan2(MESH1D.COOR(:,2),MESH1D.COOR(:,1))*180/pi ; 
    XLL = '\theta  (Deg. )' ; 
    %%%%% ERROR: There is some mistake here. We disable this option for
    %%%%% curved elements. 
    PLOT_DIAGRAMS = 0 ;
else
    X = MESH1D.COOR(:,1) ; 
    XLL = 'x' ; 
    PLOT_DIAGRAMS = 1 ;
end

if PLOT_DIAGRAMS == 1

nmodes = size(rRB{1},1) ; 
if nmodes == 6
NAMES = {'Axial force','Shear force -y','Shear force -z','Torsion','Bending Moment -y','Bending Moment -z'};
else
  NAMES = {'Axial force','Shear force -y', 'Bending Moment -z'};

end

figureREF = 1999 ;

for iforce = 1:length(NAMES)
    figure(figureREF+iforce) ;
    hold on
    xlabel(XLL)
    ylabel(NAMES{iforce}) ;
    
    hhh=plot(X,GeneralizedForces(iforce,:));
    
    VALMIN = min(GeneralizedForces(iforce,:)) ;
    VALMAX = max(GeneralizedForces(iforce,:)) ;
    
    title(['MAX VAL = ',num2str(VALMAX),' , ', 'MIN VAL = ',num2str(VALMIN)])
    
    if ~isempty(DATAIN.TIME_STEP_LOCAL)
        legend(hhh,['t = ',num2str(DATAIN.TIME_STEP_LOCAL)])
    end
    grid on
end


end



