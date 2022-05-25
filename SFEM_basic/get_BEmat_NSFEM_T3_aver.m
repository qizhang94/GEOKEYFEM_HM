function [nodB, B_ino, E_ino] = get_BEmat_NSFEM_T3_aver(ino_adjele, area_ino, area_T3, mat_Be, mat_Ee)

%------------------------------------------------------------------------
%  Purpose:
%     the strain-displacement matrix of (ino)-th SD and nodes contributing
%     (ino)-th smoothing domain
%
%  Synopsis:
%     [nodB,B_ino] = get_Bmat_NSFEM_T3_aver(ino_adjele,area_ino,area_T3,mat_Be) 
%
%  Variable Description:
%   nodB       = nodes contributing (ino)-th smoothing domain
%   B_ino      = strain-displacement matrix of (ino)-th smoothing domain (SD)
%   ino_adjele = adjacent elements of (ino)-th node
%   area_ino   = area of (ino)-th SD 
%   area_T3    = vector containing the area of elements
%   mat_Be     = matrix containing strain-displacement matrices of elements
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

global ele_nods

ne = length(ino_adjele); % number of adjacent elements of node 

for ie = 1:ne  % loop for adjacent elements of node 'ino'
    ie_nods = ele_nods(ino_adjele(ie),:);  % nodes of (ie)-th adjacent element, a row vector of length 3
    nn_ie = length(ie_nods);  % number of nodes of (ie)-th adjacent element
        if ie==1   %  initialization of B_ino for 1st adjacent element
             nodB = ie_nods;  % extract nodes of 1st adjacent element
             nn=nn_ie;           % current total number of nodes of 'nodB'
             B_ino =1/3*area_T3(ino_adjele(ie))*mat_Be{ino_adjele(ie)};
             E_ino = 1/3*area_T3(ino_adjele(ie))*mat_Ee{ino_adjele(ie)};
        else    % assemble into matrix 'B_ino' from other adjacent elements
            i0=0;   % number of added new nodes from each adjacent element
            for ino=1:nn_ie  % loop for nodes of (ie)-th adjacent element
                nod= ie_nods(ino);  % extract (ino)-th node of (ie)-th adjacent element => give variable 'nod'
                flag=0;      % flag=1: 'nod' belongs to the set 'nodB'  
                             % flag=0: 'nod' doesn't belongs to the set 'nodB'
                for j=1:nn   % loop for current total number of nodes of 'nodB'
                   if nodB(j)==nod    % if 'nod' is the j-th member of 'nodB'    
                                      % => assemble into 'B_ino' at the corresponding current positions
                      B_ino(:, 2*j-1:2*j)=B_ino(:, 2*j-1:2*j) + 1/3*area_T3(ino_adjele(ie))*mat_Be{ino_adjele(ie)}(:, 2*ino-1:2*ino);
                      E_ino(:, j)=E_ino(:, j) + 1/3*area_T3(ino_adjele(ie))*mat_Ee{ino_adjele(ie)}(:, ino);
                      flag=1;  break  % after determining 'nod' to be the j-th member of 'nodB'
                                      % => exit the current inmost loop 'for' to consider the next node of
                                      % (ie)-th adjacent element
                    end
                end
                if flag==0   % 'nod' doesn't belongs to the set 'nodB' => new node
                    i0=i0+1;          
                    nodB(nn+i0)=nod;   % add 'nod' into 'nodB'
                    % add strain-disp matrix of 'nod' into 'B_ino' at next new positions 
                    B_ino(:, 2*(nn+i0)-1:2*(nn+i0)) = 1/3*area_T3(ino_adjele(ie))*mat_Be{ino_adjele(ie)}(:, 2*ino-1:2*ino);
                    E_ino(:, nn+i0)=1/3*area_T3(ino_adjele(ie))*mat_Ee{ino_adjele(ie)}(:, ino);
                end
            end
            nn=nn+i0;   % current total number of nodes of 'nodB', after loop over one adjacent element, so nn is always updated to be the size of nodB
        end
end
   
B_ino = B_ino/area_ino;   % smooth (or average) the matrix B_ino
E_ino = E_ino/area_ino;

return
    