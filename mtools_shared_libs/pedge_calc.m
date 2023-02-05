
% pedge: calculates EDGE pressure levels given hyai, hybi, PSFC
% for the entire grid. [Pa]
function pedge = pedge_calc(hyai, hybi, PSFC)
    % sizes
    IM = size(PSFC, 1);
    JM = size(PSFC, 2);
    LM = size(hyai, 1) - 1;
    %DT = size(PSFC, 3);

    % assumes PSFC is in [Pa] (but outputs just follow PSFC units anyway)
    % if you want to use this for delp_dry calc and AD (kg) calc though
    % must give Pa
    % pedge = zeros(IM, JM, LM+1, DT);
    pedge = zeros(IM, JM, LM+1);

    fprintf("\n=> Computing PEDGE ... IM=%d,JM=%d,LM=%d\n",IM,JM,LM)

    %for T = 1:DT
    for L = 1:LM+1
        % pedge(:,:,L,T) = hyai(L) + hybi(L) .* PSFC(:,:,T);

        % GEOS-Chem format
        %pedge(:,:,L) = hyai(L) + (hybi(L) * PSFC(:,:));

        % CAM-chem format https://www.ncl.ucar.edu/Document/Functions/Built-in/pres_hybrid_ccm.shtml
        pedge(:,:,L) = hyai(L) * 100000.0 + hybi(L) * PSFC(:,:);
    end
    %end
end
 
