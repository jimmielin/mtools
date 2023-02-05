
% bxheight_approx: gives approximation of box height [m] for each
% grid box. this is an approximation because we do not use
% moist air in this calculation, instead assuming q = 0.
% virtual temperature should be used instead, but Q/H2O is not always
% available.
function bxheight_approx = bxheight_approx(T, PEDGE)
    IM = size(PEDGE, 1);
    JM = size(PEDGE, 2);
    LM = size(PEDGE, 3) - 1;

    bxheight_approx = zeros(IM, JM, LM);
    for L = 1:LM
        bxheight_approx(:,:,L) = (287.0 / 9.80665) * T(:,:,L) .* log(PEDGE(:,:,L) ./ PEDGE(:,:,L+1));
    end
end
 
