function SpatialSVD_2D(xx,PHIloc,INDICES_ECM_POINTS,TOL_SVD,LABELloc,NFIG)


 % We have to rephape PHIloc into a matrix PHIloc_x_y,
            % so that, upon application of the SVD, we get
            %  U(x)*S*V^T(y) = PHIloc_x_y
            PHIloc_x_y = reshape(PHIloc,length(xx{2}),[]) ;
            PHIloc_x_y = PHIloc_x_y' ;
            DATAsvd.RELATIVE_SVD = 1;
            [Ux,S,Vy] = RSVDT(PHIloc_x_y,TOL_SVD,[],0,DATAsvd) ;
            figure(NFIG)
            
            
            subplot(2,1,1)
            hold on
            title(LABELloc)
            xlabel('x')
            ylabel('U(x)')
            grid on
            subplot(2,1,2)
            hold on
            xlabel('y')
            ylabel('V(y)')
            grid on
            COLORSLOC = rand(length(S),3) ;
            hx = [] ;    LEGENDx={} ;
            subplot(2,1,1)
            hold on
            
            idim = 1;
            setINDEX = INDICES_ECM_POINTS(:,idim) ;
            xxLOC = xx{idim}(setINDEX) ;
            for i = 1:length(S)
                hx(i) = plot(xx{1},Ux(:,i),'Marker','.','Color',COLORSLOC(i,:)) ;
                LEGENDx{i} = ['S_',num2str(i),'=',num2str(S(i))] ;
                yyLOC = Ux(setINDEX,i) ;
                plot(xxLOC,yyLOC,'rx','MarkerSize',6) ;
            end
            
            
            
            
            
            
            legend(hx,LEGENDx) ;
            hy = [] ;    LEGENDy={} ;
            subplot(2,1,2)
            hold on
            idim = 2;
            setINDEX = INDICES_ECM_POINTS(:,idim) ;
            xxLOC = xx{idim}(setINDEX) ;
            for i = 1:length(S)
                hy(i) = plot(xx{2},Vy(:,i),'Marker','.','Color',COLORSLOC(i,:)) ;
                LEGENDy{i} = ['S_',num2str(i) ] ;
                yyLOC = Vy(setINDEX,i) ;
                plot(xxLOC,yyLOC,'rx','MarkerSize',6) ;
            end
            legend(hy,LEGENDy) ;
            
            
            