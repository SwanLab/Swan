classdef FileComparator < handle

    properties(Access = private)
        tol = 1e-4;
    end

    methods (Access = public)

        function theyAre = areFilesDifferent(obj,fileComputed,fileStored)

            A = obj.getContent(fileComputed);
            B = obj.getContent(fileStored);

            is_equal = obj.isContentEqual(A,B);

            theyAre = ~is_equal;
        end

    end

    methods (Access = private)

        function itIs = isContentEqual(obj,A,B)
            nLines      = size(A,1);
            checkVector = false(size(A));
            for i = 1:nLines
                s.Ai  = A(i);
                s.Bi  = B(i);
                s.tol = obj.tol;
                checkVector(i) = obj.isLineEqual(s);
            end
            itIs = all(checkVector);
        end

    end

    methods (Access = private, Static)

        function a = getContent(file)
            fid = fopen(file);

            a = textscan(fid,'%s','delimiter','\n');
            a = a{1};

            fclose(fid);
        end

        function check = isLineEqual(cParams)
            Ai  = cParams.Ai;
            Bi  = cParams.Bi;
            tol = cParams.tol;
            if isempty(cell2mat(Ai)) && isempty(cell2mat(Bi))
                check = true;
            else
                a = cell2mat(textscan(char(Ai),'%f'));
                if isempty(a)
                    check = strcmp(Ai,Bi);
                else
                    b = cell2mat(textscan(char(Bi),'%f'));
                    if norm(b) == 0
                        err = norm(a - b);
                    else
                        err = norm(a - b)/norm(b);
                    end
                    if err<tol
                        check = true;
                    else
                        check = false;
                    end
                end
            end
        end

    end

end