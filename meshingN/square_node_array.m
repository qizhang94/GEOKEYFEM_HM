function node=square_node_array(pt1,pt2,pt3,pt4,nnx,nny)

% Generates a quadratleral array of nodes between the counterclockwise 
% ordering of nodes pt1 - pt4,given number of elements in x and y direction 


if ( nargin < 6 )
   disp('Not enough parameters specified for quare_node_array function')

end

% get node spacing along u direction

  xi_pts=linspace(-1,1,nnx);


% get node spacing along v direction

  eta_pts=linspace(-1,1,nny);


x_pts=[pt1(1),pt2(1),pt3(1),pt4(1)];
y_pts=[pt1(2),pt2(2),pt3(2),pt4(2)];

for r=1:nny
  eta=eta_pts(r);
  for c=1:nnx
    xi=xi_pts(c);
    % get interpolation basis at xi, eta
    N=shape_func('Q4',[xi,eta]);
    N=N(:,1);
    node((r-1)*nnx+c,:)=[x_pts*N,y_pts*N];
  end
end
end

