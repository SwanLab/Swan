function DATAoffline =  PlotDeformedShapeClusters(DATAoffline,IND,MESH,CENT,SNAPdisp,DISTANCE_CENT)



if DATAoffline.PLOT_CLUSTERS_DEFORMED_SHAPE ==1
    
    DATAoffline = DefaultField(DATAoffline,'FREQ_STEPS',1) ;
    DATAoffline = DefaultField(DATAoffline,'SPECIFIC_SET_POINTS_PLOT',[]) ;
    DATAoffline = DefaultField(DATAoffline,'COLOR_plot_def','r.') ;
    DATAoffline = DefaultField(DATAoffline,'SizeFont',10) ;
    DATAoffline = DefaultField(DATAoffline,'LINEWIDTH_plot_def',1) ;
     DATAoffline = DefaultField(DATAoffline,'nfig_DEFORMED',1090) ;
         
 
    
    
    %     if ~isempty(DATAoffline.SPECIFIC_SET_POINTS_PLOT)
    %         figure(345)
    %         hold on
    %         xlabel('Snapshot')
    %         ylabel('Cluster')
    %         bar(IND)
    %     end
    ifig = DATAoffline.nfig_DEFORMED  ;
    figure(ifig)
    hold on
    
    
    
    if  ~isempty(DATAoffline.SPECIFIC_SET_POINTS_PLOT)
        load(DATAoffline.SPECIFIC_SET_POINTS_PLOT) ;
        
        NodesBoundary = load(DATAoffline.SPECIFIC_SET_POINTS_PLOT) ;
        NodesBoundary = NodesBoundary(:,1) ;
        
    else
        NodesBoundary = unique(MESH.CNb(:)) ;
    end
    
    if size(CENT,1) > size(CENT,2)
        CENT = CENT' ;
    end
    
    INDexam = 1:DATAoffline.FREQ_STEPS:size(CENT,1) ;
    
    CENT = CENT(INDexam,:) ;
    
    
    for i=1:size(CENT,1)
        if size(MESH.COOR,2) == 2
            % 2D plot
            INDx = 2*NodesBoundary-1; % x-indexes
            INDy = 2*NodesBoundary ; % y-indexes
            
            % Deformed shape corresponding to each centroid
            CENTx = CENT(i,INDx) + MESH.COOR(NodesBoundary,1)' ;
            CENTy = CENT(i,INDy) + MESH.COOR(NodesBoundary,2)' ;
            %           for idim = 1:size(MESH.COOR,2)
            %    IND = idim:size(MESH.COOR,2):size(SNAPdisp,1);
            %
            %    SNAPdisp(IND,:) = bsxfun(@plus,SNAPdisp(IND,:),MESH.COOR(:,idim));
            %
            % end
            nC = CENTx.^2 + CENTy.^2 ;
            [~,III] = max(nC) ;
            %
            %
            x = CENTx(III) ;
            y = CENTy(III);
            
            
            plot(CENTx,CENTy,DATAoffline.COLOR_plot_def,'LineWidth',DATAoffline.LINEWIDTH_plot_def)
            %
            if DATAoffline.PLOT_CLUSTERS_DEFORMED_SHAPE_ONLY_CENTROIDS ==0
                nsteps = find(IND == i) ;
                maxN = max(nsteps) ;
                minN = min(nsteps) ;
                text(x,y,['i = ',num2str(i),'(n  = ',num2str(minN),':',num2str(maxN),')']);
                SNAPdispx = SNAPdisp(INDx,maxN) + MESH.COOR(NodesBoundary,1) ;
                SNAPdispy = SNAPdisp(INDy,maxN) + MESH.COOR(NodesBoundary,2) ;
                plot(SNAPdispx,SNAPdispy,'y.')
            else
                text(x,y,[num2str(INDexam(i))],'FontSize',DATAoffline.SizeFont);
                
            end
            
        else
            error('Option not available')
        end
        
    end
    
end