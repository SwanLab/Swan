function PlotInterfaceVariablesMOD(DATAIN,a,ndimINTF,MESH1D)

DATAIN = DefaultField(DATAIN,'PlotInterfaceVariables',0) ; 



if DATAIN.PlotInterfaceVariables == 1
    
    DATAIN = DefaultField(DATAIN,'ScalingFactorsInterfaceVariables',ones(1,ndimINTF)) ; 
    
    figure(5000)
    hold on 
    xlabel('x') 
    ylabel('a_i')
    title(['Intensity of each mode along the x-coordinate'])
    aPLOT = reshape(a,ndimINTF,[]) ; 
    COLOR = ColoresMatrix(ndimINTF) ; 
    LGG = {} ; 
    for iii = 1:size(aPLOT,1)
        ScalingFactor = DATAIN.ScalingFactorsInterfaceVariables(iii) ; 
        hpp(iii) = plot(MESH1D.COOR(:,1),aPLOT(iii,:)*ScalingFactor,'Color',COLOR(iii,:));
        if ScalingFactor == 1
        LGG{iii} = ['Mode ',num2str(iii)] ; 
        else
                LGG{iii} = ['Mode ',num2str(iii),' x ',num2str(ScalingFactor)] ;   
        end
    end
    legend(hpp,LGG) ; 
end
