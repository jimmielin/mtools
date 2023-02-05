function TropLev = TropLevFromTROPP(fourth_idx, IM, JM, LM, pedge, TropP)
	if fourth_idx > 0
        TropLev = zeros(IM, JM);
        for I = 1:IM
        for J = 1:JM
            for L = 1:LM
                if pedge(I,J,L+1) <= TropP(I,J,fourth_idx)
                    TropLev(I,J) = L;
                    break
                end
            end
        end
        end
    else
        TropLev = zeros(IM, JM, 12);
        for M = 1:12
            for I = 1:IM
            for J = 1:JM
                for L = 1:LM
                    if pedge(I,J,L+1,M) <= TropP(I,J,M)
                        %fprintf("found tropp l %d %f %f\n", L, pedge_left(I,J,L+1), TropP_left(I,J,M))
                        TropLev(I,J,M) = L;
                        break
                    end
                end
            end
            end
        end
    end
end