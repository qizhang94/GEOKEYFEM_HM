function [shapeq4,dNdrq4,dNdsq4]=get_shape_Q4(rvalue,svalue)

%------------------------------------------------------------------------
%  Purpose:
%     compute isoparametric four-node quadilateral shape functions
%     and their derivatves at the selected (integration) point
%     in terms of the natural coordinate 
%
%  Synopsis:
%     [shapeq4,dNdrq4,dNdsq4]=get_shape_Q4(rvalue,svalue)  
%
%  Variable Description:
%     shapeq4 - shape functions for four-node element
%     dNdrq4 - derivatives of the shape functions w.r.t. r
%     dNdsq4 - derivatives of the shape functions w.r.t. s
%     rvalue - r coordinate value of the selected point   
%     svalue - s coordinate value of the selected point
%
%  Notes:
%     1st node at (-1,-1), 2nd node at (1,-1)
%     3rd node at (1,1), 4th node at (-1,1)
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

% shape functions
 shapeq4  = zeros(4,1);
 shapeq4(1)=0.25*(1-rvalue)*(1-svalue);
 shapeq4(2)=0.25*(1+rvalue)*(1-svalue);
 shapeq4(3)=0.25*(1+rvalue)*(1+svalue);
 shapeq4(4)=0.25*(1-rvalue)*(1+svalue);

% derivatives

 dNdrq4(1)=-0.25*(1-svalue);
 dNdrq4(2)=0.25*(1-svalue);
 dNdrq4(3)=0.25*(1+svalue);
 dNdrq4(4)=-0.25*(1+svalue);

 dNdsq4(1)=-0.25*(1-rvalue);
 dNdsq4(2)=-0.25*(1+rvalue);
 dNdsq4(3)=0.25*(1+rvalue);
 dNdsq4(4)=0.25*(1-rvalue);
