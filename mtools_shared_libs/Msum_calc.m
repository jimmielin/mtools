function Msum = Msum_calc(spec_current, AD, pedge, VMRtoMMR, ...
                          IM, JM, LM, fourth_idx, ...
                          lm_range_option, LM_range, TropLev)    
    Msum = 0.0;
    if lm_range_option == "manual"
        for L = LM_range
            % sum by levels ..                 vv note this is element-wise mult.
            if fourth_idx > 0
                by_layer = sum(sum(AD(:,:,L) .* spec_current(:,:,L,fourth_idx) * VMRtoMMR));
            elseif fourth_idx == -1 % assuming fourth_idx spans from 1:12 ... yearly sum
                by_layer = 0; % accum.
                for mo = 1:12
                    by_layer = by_layer + sum(sum(AD(:,:,L,mo) .* spec_current(:,:,L,mo) / 12 * VMRtoMMR));
                end
            end
            Msum = Msum + by_layer;
        end
    else
        for I = 1:IM
        for J = 1:JM
            if fourth_idx > 0
                TropLev = TropLev(I,J);
                if lm_range_option == "strat"
                    LM_range = (TropLev+1):LM;
                elseif lm_range_option == "trop"
                    LM_range = 1:TropLev;
                end
                %fprintf("[L] i=%d,j=%d,troplev=%d", I, J, TropLev);
    
                for L = LM_range
                    by_box = AD(I,J,L) .* spec_current(I,J,L,fourth_idx) * VMRtoMMR;
                    Msum = Msum + by_box;
                end
            elseif fourth_idx == -1
                for M = 1:12
                    TropLev = TropLev(I,J,M);
                    if lm_range_option == "strat"
                        LM_range = (TropLev+1):LM;
                    elseif lm_range_option == "trop"
                        LM_range = 1:TropLev;
                    end
                    %fprintf("[L] i=%d,j=%d,troplev=%d", I, J, TropLev);
        
                    for L = LM_range
                        % sum by levels ..                 vv note this is element-wise mult.
                        % assuming fourth_idx spans from 1:12 ... yearly sum
                        by_box = AD(I,J,L,M) .* spec_current(I,J,L,M) / 12 * VMRtoMMR;
                        Msum = Msum + by_box;
                    end
                end
            end
        end
        end
    end
end