function [coordinatesa,pressure] = actual_df(actividad,free_df,coordinatesa,pressure,problembsc,element,fixnodes,d_u,contactb )

% Actualiza grados de libertad restringidos

    switch element.type
        case {'TRIANGLE','QUAD'}
            switch actividad % Actualizamos las coordenas actuales con el desplazamiento impuesto
                case 'BC'
                    for i = 1 : size(fixnodes,1)
                        coordinatesa(fixnodes(i,1),fixnodes(i,2))=coordinatesa(fixnodes(i,1),fixnodes(i,2))+(1/problembsc.loadsteps)*fixnodes(i,3);
                    end

                case 'SOL'
                    for i=1: size(free_df,2)
                        if mod(free_df(i),2)==0
                            coordinatesa(ceil(free_df(i)/2),2)=coordinatesa(ceil(free_df(i)/2),2)+d_u(free_df(i),1);
                        else
                            coordinatesa(ceil(free_df(i)/2),1)=coordinatesa(ceil(free_df(i)/2),1)+d_u(free_df(i),1);
                        end
                    end
                otherwise
            end
  
        case {'HEXAHEDRA'}
            switch actividad
                case 'BC'
                    for i = 1 : size(fixnodes,1)   
                        coordinatesa(fixnodes(i,1),fixnodes(i,2)) = coordinatesa(fixnodes(i,1),fixnodes(i,2)) + (1/problembsc.loadsteps)*fixnodes(i,3);
                    end

                case 'SOL'
                    for i=1: size(free_df,2)
                        if mod(free_df(i)+3,3)==1
                            coordinatesa((free_df(i)+2)/3,1)=coordinatesa((free_df(i)+2)/3,1)+d_u(free_df(i),1);
                        end
                        if mod(free_df(i)+3,3)==2
                            coordinatesa((free_df(i)+1)/3,2)=coordinatesa((free_df(i)+1)/3,2)+d_u(free_df(i),1);
                        end
                        
                    end

                otherwise
            end
        otherwise
            error('No existe el elemento')
    end

end