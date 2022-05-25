function [area_nod, area_T3] = cal_area_nod_T3(nod_adjele)

%------------------------------------------------------------------------
%  Purpose:
%     Compute the area of SD associated with node and element areas
%
%  Synopsis:
%     [area_nod,area_T3] = cal_area_nod_T3(nod_adjele)  
%
%  Variable Description:
%   area_nod : vector containing the area of SD associated with node
%   area_T3  : vector containing the area of elements
%  nod_adjele: matrix containing adjacent elements of each node
%--------------------------------------------------------------------------
% Coded by Dr. Nguyen Thoi Trung (Nguyen-Thoi T or Nguyen T.T)            %
% University of Science - Vietnam National University HCMC, Vietnam     %
% National University of Singapore (NUS)                                  %
% email: thoitrung76@gmail.com                                            %
% Last modified: December 2009                                            %
%--------------------------------------------------------------------------

% Important note: The authors decided to release the source codes free of charge with the hope that the S-FEM technique 
% can be applied to more problems and can be further developed into even more powerful methods. 
% The authors are not be able to provide any services or user-guides, but appreciate very much your feedback on errors and suggestions, 
% so that we can improve these codes/methods in the future. If the idea, method, and any part of these codes are used in anyway,
% the users are required to cite the book and the following related original papers of the authors:

% Liu GR, The Finite element method a practical course, Elsevier (BH), UK.  2003.
% Liu, G. R. and Nguyen Thoi Trung, Smoothed Finite Element Method, CRC press, Boca Raton, USA, 2010.
% Liu GR, Dai KY, Nguyen-Thoi T (2007) A smoothed finite element method for mechanics problems. Computational Mechanics; 39: 859-877.
% Liu GR, Nguyen-Thoi T, Nguyen-Xuan H, Lam KY (2009) A node-based smoothed finite element method (NS-FEM) for upper bound solutions to solid mechanics problems. Computers and Structures; 87: 14-26.
% Nguyen-Thoi T, Liu GR, Nguyen-Xuan H (2009) Additional properties of the node-based smoothed finite element method (NS-FEM) for solid mechanics problems. International Journal of Computational Methods, 6(4): 633-666.
% Nguyen-Thoi T, Liu GR, Nguyen-Xuan H, Nguyen-Tran C (2010) Adaptive analysis using the node-based smoothed finite element method (NS-FEM). International Journal for Numerical Methods in Biomedical Engineering; in press, doi: 10.1002/cnm.1291.
%--------------------------------------------------------------------------

global nel nnode gcoord ele_nods

%Compute the area of each triangle
for i = 1:nel
    x1 = gcoord(ele_nods(i,1),1);
    y1 = gcoord(ele_nods(i,1),2);
    x2 = gcoord(ele_nods(i,2),1);
    y2 = gcoord(ele_nods(i,2),2);
    x3 = gcoord(ele_nods(i,3),1);
    y3 = gcoord(ele_nods(i,3),2);
    area_T3(i) = 0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2); % area of triangle;
end

%----------------------------------------------
% Compute the area of SDs associated with nodes
%----------------------------------------------

area_nod = zeros(nnode,1);
for ino=1:nnode
    for j=1:length(nod_adjele{ino})
        area_nod(ino) = area_nod(ino) + area_T3(nod_adjele{ino}(j))/3;        
    end
end

% or
%area_nod=sparse(ele_nods,ones(nel,3),1/3*area_T3'*[1,1,1],nnode,1);

return
