function IdealIndexCluster =  PlotParametricSpaceClusters(DATAoffline,MATRIX_POINTS_SPACE_PARAMETER,IND,DISTANCE_CENT,BasisU_cluster,...
    LOWERBOUND_matrix,TransMatrix,UPPERBOUND_matrix,CENT,BasisU_ALL,TESTING_SNAPSHOTS,TESTING_PARAM,MESH,NAME_INPUTS)

if nargin == 0
    load('tmp1.mat')
  %  close all
    
end

DATAoffline= DefaultField(DATAoffline,'PLOT_CLUSTERS_parameter_space',1) ; 
DATAoffline= DefaultField(DATAoffline,'OffSetShowNumberModes',0) ; 
DATAoffline= DefaultField(DATAoffline,'AxisYShow',[]) ; 

 
if DATAoffline.PLOT_CLUSTERS_parameter_space == 1
    
    
    ifig = 590; 
    figure(ifig)
    hold on
    xlabel('Parameter 1')
    if size(MATRIX_POINTS_SPACE_PARAMETER,2) == 2
    ylabel('Parameter 2')
    end
    title('Distribution clusters in parameter space')
    h = []  ;
    LEGGG = {} ;
    %  colores = ColoresMatrix(DATAoffline.NCLUSTERS_BASIS_DISP);
    for iii = 1:length(BasisU_cluster)
        if ~iscell(IND)
        INDLOC = find(IND ==iii) ;
        else
            INDLOC = IND{iii} ; 
        end
       
        
        if size(MATRIX_POINTS_SPACE_PARAMETER,2) == 2  
        h(iii) =   plot(MATRIX_POINTS_SPACE_PARAMETER(INDLOC,1),MATRIX_POINTS_SPACE_PARAMETER(INDLOC,2),'.','Markersize',12) ;
        else
            h(iii) =   plot(MATRIX_POINTS_SPACE_PARAMETER(INDLOC,1),zeros(size(MATRIX_POINTS_SPACE_PARAMETER(INDLOC,1))),'.','Markersize',12) ;
     
        end
        LEGGG{iii} = ['Cluster = ',num2str(iii),'; n_u = ',num2str(size(BasisU_cluster{iii},2))] ;
        
        % CENTROID
        if isempty(CENT)
                PLOT_CENTROID = 0 ;
        else 
        PLOT_CENTROID = 1 ;
        end
        if PLOT_CENTROID ==1 
        DATAoffline = DefaultField(DATAoffline,'CentroidIndex',[]) ;
        if isempty(DATAoffline.CentroidIndex)
            [MINDDDDD,CCC] =   min(DISTANCE_CENT(:,iii)) ;
            
        else
            CCC =DATAoffline.CentroidIndex(iii)  ; 
        end
        
        
        XXX = MATRIX_POINTS_SPACE_PARAMETER(CCC,:) ;
            if length(XXX) == 1
                XXX(2) = 0 ;
            end
       
       plot(XXX(1),XXX(2),'x','MarkerSize',14)
       text(XXX(1),XXX(2)+DATAoffline.OffSetShowNumberModes,num2str(iii),'FontSize',16) ; 
        end
        
    end
    
    legend(h,LEGGG)
    
    if ~isempty(DATAoffline.AxisYShow)
        %DATAoffline= DefaultField(DATAoffline,'MaxAxisYShow',[]) ; 
        aaa = axis; 
        aaa([3,4]) = DATAoffline.AxisYShow ; 
        axis(aaa); 
    end
    
    
    % PLOTTING TESTING TRAJECTORY 
    DATAoffline = DefaultField(DATAoffline,'TESTING_TRAJ',[]) ;
    DATAoffline.TESTING_TRAJ = DefaultField(DATAoffline.TESTING_TRAJ,'TESTING_TRAJ',[]) ;

    NAMEsnap_base = DATAoffline.TESTING_TRAJ.NAMEsnap_base ; CASES = [] ;
    if ~isempty(NAMEsnap_base)
        [BasisU,SNAPdisp,MATRIX_POINTS_SPACE_PARAMETER,MESH,TIME_STEPS,IdealIndexCluster ]= ...
            FindOutClusterTransitionFE(BasisU_cluster,NAMEsnap_base,DATAoffline,ifig,LOWERBOUND_matrix,TransMatrix,...
            UPPERBOUND_matrix,CENT,BasisU_ALL,TESTING_SNAPSHOTS,TESTING_PARAM,MESH,NAME_INPUTS) ;
        
        
           
        
    end
    
    

    
    
    
    
    
end