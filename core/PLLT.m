function [A,theta,Gamma,CL,CDi,alfind] = PLLT(ynd,c,b,AR,a,alf,alfZL,V)
%PLLT  Analyze wings using Prandtl's lifting-line theory
% Inputs:
%    ynd: dimensionless grid (-1,1) on which distributions are defined
%      c: chord distribution (scalar or column vector)
%      b: span
%      a: lift-curve slope [1/rad] (scalar or column vector)
%     AR: aspect ratio
%    alf: angle of attack distribution [rad] (column vector)
%  alfZL: zero-lift angle of attack distribution [rad] (column vector)
%      V: velocity
% Outputs:
%      A: Fourier coefficients
%  theta: grid spacing in polar domain [rad]
%  Gamma: circulation distribution
%     CL: wing lift coefficient
%    CDi: wing induced drag coefficient
% alfind: induced angle of attack [rad]

n = 1:2:length(ynd);
% generate theta stencil
theta = acos(ynd);
% compute cofactors in the Fourier sum
cofactors = sin(n.*theta) + c/(4*b).*a.*n.*sin(n.*theta)./sin(theta);
RHS = c/(4*b).*a.*(alf - alfZL);
A = cofactors \ RHS; % solve Fourier coefficients
Gamma = 2*b*V * sin(n.*theta)*A;
CL = pi*AR*A(1);
CDi = pi*AR*n*(A.^2);
alfind = alf - alfZL - 2*Gamma./(a.*c*V);
