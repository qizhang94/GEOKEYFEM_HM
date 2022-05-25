function [jacob2D]=cal_jacob_2D(nnel,dNdr,dNds,xco,yco)

%------------------------------------------------------------------------
%  Purpose:
%     determine the Jacobian matrix for two-dimensional mapping
%
%  Synopsis:
%     [jacob2D]=cal_jacob_2D(nnel,dNdr,dNds,xco,yco) 
%
%  Variable Description:
%     jacob2D - Jacobian for two-dimensional mapping
%     nnel - number of nodes per element   
%     dNdr - derivative of shape functions w.r.t. natural coordinate r
%     dNds - derivative of shape functions w.r.t. natural coordinate s
%     xco - x axis coordinate values of nodes
%     yco - y axis coordinate values of nodes
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

 jacob2D=zeros(2,2);

 for i=1:nnel
    jacob2D(1,1)=jacob2D(1,1)+dNdr(i)*xco(i);
    jacob2D(1,2)=jacob2D(1,2)+dNdr(i)*yco(i);
    jacob2D(2,1)=jacob2D(2,1)+dNds(i)*xco(i);
    jacob2D(2,2)=jacob2D(2,2)+dNds(i)*yco(i);
 end
