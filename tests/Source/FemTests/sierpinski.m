function carpet = sierpinski(levels,classname)
            if nargin == 1
                classname = 'single';
            end
            
            msize = 3^levels;
            carpet = ones(msize,classname);
            
            cutCarpet(1,1,msize,levels) % Begin recursion
            
            function cutCarpet(x,y,s,cl)
                if cl
                    ss = s/3; % Define subsize
                    for lx = 0:2
                        for ly = 0:2
                            if lx == 1 && ly == 1
                                % Remove center square
                                carpet(x+ss:x+2*ss-1,y+ss:y+2*ss-1) = 0;
                            else
                                % Recurse
                                cutCarpet(x+lx*ss,y+ly*ss,ss,cl-1)
                            end
                        end
                    end
                end
            end
        end