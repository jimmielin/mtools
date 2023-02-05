

% delp_dry: give PEDGE_DRY, returns deltas for each grid box.
function delp_dry = delp_dry_calc(PEDGE)
    IM = size(PEDGE, 1);
    JM = size(PEDGE, 2);
    LM = size(PEDGE, 3) - 1;

    delp_dry = zeros(IM, JM, LM);
    for L = 1:LM
        delp_dry(:,:,L) = PEDGE(:,:,L) - PEDGE(:,:,L+1);
    end
end
