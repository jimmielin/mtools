% area_m2: calculates M2 area of grid box given CENTER lons, lats
% which have been expanded to 2-D sizes lons=lats=(IM,JM)
%
% assumes CENTER because that is what we work with out of COARDS/CF
% conventions. cf. https://cfconventions.org/cf-conventions/cf-conventions.html
% "Bounds for 2-D coordinate variables with 4-sided cells"
%
% currently only works on cartesian grid.
%
% assumes lats are from south to north
function area_m2 = area_m2_calc(lons, lats)
    % Re: Radius of Earth [m]
    Re = 6.3710072e6;

    % assuming these are center quantities.
    IM = size(lons, 1);
    JM = size(lats, 2);

    % create area_m2:
    % fixme: may not be exactly accurate for curvilinear
    % due to the calculation of dx is
    % based on adjacent center cells.
    area_m2 = zeros(IM, JM);

    % loop over the grid
    for J = 1:JM
        for I = 1:IM
            % compute the edge radian of the y-coordinate to north and
            % south. given the centers, we can do 1/2
            %if lons(I,J) == 90.0
            %    YEdgeRad_N = pi/2;
            %elseif lons(I,J) == -90.0
            %    YEdgeRad_N = -pi/2;
            %else
           if J == 1
               YEdge_N = (lats(I,J+1) + lats(I,J))/2;
               YEdge_S = -90.0;
           elseif J == JM
               YEdge_N =  90.0;
               YEdge_S = (lats(I,J-1) + lats(I,J))/2;
           else
               YEdge_N = (lats(I,J+1) + lats(I,J))/2;
               YEdge_S = (lats(I,J-1) + lats(I,J))/2;
           end

           YEdgeRad_N = YEdge_N * pi / 180;
           YEdgeRad_S = YEdge_S * pi / 180;

           YEdgeSin_N = sin(YEdgeRad_N);
           YEdgeSin_S = sin(YEdgeRad_S);

           % get adjacent grid box dx (assuming rectilinearity here...)
           if I ~= IM
               DX = lons(I+1,J) - lons(I,J);
           else
               DX = lons(I,J)   - lons(I-1,J);
           end

           area_m2(I,J) = (DX * pi / 180) * (Re^2) * (YEdgeSin_N - YEdgeSin_S);
        end
    end


end
