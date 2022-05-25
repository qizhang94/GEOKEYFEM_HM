function [index]=get_eledof(nod,nnel,ndof)
%----------------------------------------------------------
%  Purpose:
%     Extract system dofs associated with each element or nod
%
%  Synopsis:
%     [index]=get_eledof(nod,nnel,ndof)
%
%  Variable Description:
%     index : system dof vector associated with element 
%     nnel  : number of nodes per element
%     ndof  : number of dofs per node 
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
%--------------------------------------------------------------------------

%  edof = nnel*ndof;   % number of dofs of the element
k=0;
for i=1:nnel      % loop for nnel nodes of the element
    start = (nod(i)-1)*ndof;
    for j=1:ndof
        k=k+1;
        index(k)=start+j;
    end
end

 
