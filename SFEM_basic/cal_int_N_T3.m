function [int_N] = cal_int_N_T3(x_coord, y_coord, area, x_subSD, y_subSD)

%---------------------------------------------------------------------------------------------------------
% Calculate the integral of N on one sub-smoothing domain
% 'area' is the area of this T3 element
% x_subSD is a vector that contains the x-coordinates of the 4 vertices of
% the sub-smoothing domain, similar for y_subSD, they are in the
% counter-clockwise (important)
%---------------------------------------------------------------------------------------------------------
int_N = zeros(1,3);

x_gauss = [-sqrt(3)/3; sqrt(3)/3; sqrt(3)/3; -sqrt(3)/3];
y_gauss = [-sqrt(3)/3; -sqrt(3)/3; sqrt(3)/3; sqrt(3)/3];
weight = ones(4,1);

nnel = 4;  % Number of nodes per element

for i = 1:4   % here "4" is not nnel, but the number of Gauss integration points
    % One gauss integration point in the real domain
    [shapeq4, dNdrq4, dNdsq4]=get_shape_Q4(x_gauss(i),y_gauss(i)); % shapeq4 is a column vector
    detJ = det(cal_jacob_2D(nnel, dNdrq4, dNdsq4, x_subSD, y_subSD));   % determinant of the Jacobian mapping matrix
    
    % One gauss integration point in the real domain
    x = shapeq4' *  x_subSD;
    y = shapeq4' *  y_subSD;

    N1 = 1/(2*area)*(x_coord(2)*y_coord(3) - x_coord(3)*y_coord(2) + ...
        (y_coord(2) - y_coord(3))*x + (x_coord(3) - x_coord(2))*y);
    N2 = 1/(2*area)*(x_coord(3)*y_coord(1) - x_coord(1)*y_coord(3) + ...
        (y_coord(3) - y_coord(1))*x + (x_coord(1) - x_coord(3))*y);
    N3 = 1/(2*area)*(x_coord(1)*y_coord(2) - x_coord(2)*y_coord(1) + ...
        (y_coord(1) - y_coord(2))*x + (x_coord(2) - x_coord(1))*y);

% Evaluate N^T*N
int_N = int_N + weight(i) * [N1, N2, N3]*detJ;

end

% The difficulty is too find the vertices of the sub-smoothing domain