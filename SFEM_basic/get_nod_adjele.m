function [nod_adjele]=get_nod_adjele

%------------------------------------------------------------------------
%  Purpose:
%     find adjacent elements of each node
%
%  Synopsis:
%     [nod_adjele]=get_nod_adjele  
%
%  Variable Description:
%    nod_adjele - matrix containing adjacent elements of each node
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
% Liu GR, Nguyen-Thoi T, Nguyen-Xuan H, Lam KY (2009) A node-based smoothed finite element method (NS-FEM) for upper bound solutions to solid mechanics problems. Computers and Structures; 87: 14-26.
% Nguyen-Thoi T, Liu GR, Nguyen-Xuan H (2009) Additional properties of the node-based smoothed finite element method (NS-FEM) for solid mechanics problems. International Journal of Computational Methods, 6(4): 633-666.
% Nguyen-Thoi T, Liu GR, Nguyen-Xuan H, Nguyen-Tran C (2010) Adaptive analysis using the node-based smoothed finite element method (NS-FEM). International Journal for Numerical Methods in Biomedical Engineering; in press, doi: 10.1002/cnm.1291.
%--------------------------------------------------------------------------

global nnode ele_nods nel

nod_adjele = cell(nnode,1);
for ino=1:nnode   % loop for all nodes   
    ind=0;
    for iel=1:nel  % loop for all elements
       if (find(ino==ele_nods(iel,:)))>=1
           ind = ind+1;
           nod_adjele{ino}(ind) = iel;
       end
    end
end