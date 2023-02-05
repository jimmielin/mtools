
% AD: calculates air mass for each level
% dry air mass [kg]
function AD = ad_calc(DELP_DRY, AREA_M2)
    % from physconstants.F90 (G-C)
    G0 = 9.80665;
    % G0_100 = 100.0 / G0;

    IM = size(DELP_DRY, 1);
    JM = size(DELP_DRY, 2);
    LM = size(DELP_DRY, 3);

    AD = zeros(IM, JM, LM);
    for L = 1:LM
        %            note element-wise mult vv
        AD(:,:,L) = DELP_DRY(:,:,L) .* AREA_M2(:,:) / G0; % [Pa] calc.
        % if using [hPa], need to use G0_100 like in G-C
    end
end

 
