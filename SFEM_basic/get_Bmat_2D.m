function [Be]=get_Bmat_2D(nnel,dNdx,dNdy)

%------------------------------------------------------------------------
%  Purpose:
%     determine the strain-displacement matrix for 2D solids
%
%  Synopsis:
%     [Be]=get_Bmat_2D(nnel,dNdx,dNdy) 
%
%  Variable Description:
%       Be  : strain-displacement matrix of element
%     nnel  : number of nodes per element
%     dNdx  : derivatives of shape functions with respect to x   
%     dNdy  : derivatives of shape functions with respect to y
%------------------------------------------------------------------------
%      [dNdx(1) 0       dNdx(2) 0       ... dNdx(nnel) 0         ]
% Be = [0       dNdy(1) 0       dNdy(2) ... 0          dNdy(nnel)]
%      [dNdy(1) dNdx(1) dNdy(2) dNdx(2) ... dNdy(nnel) dNdx(nnel)]
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

 for i=1:nnel
    i1=(i-1)*2+1;  
    i2=i1+1;
    Be(1,i1)=dNdx(i);
    Be(2,i2)=dNdy(i);
    Be(3,i1)=dNdy(i);
    Be(3,i2)=dNdx(i);
 end
