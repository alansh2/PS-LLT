function [Gsol,E] = solveGamma(G0,alpha,vij,un,ua,zeta,a0,alf_ZL)
%SOLVEGAMMA  Solves for the gammas that converge the two lift equations
% Inputs:
%     G0: initial guess for vortex strengths
%  alpha: angle of attack (deg)
%    vij: dimensionless induced velocity influence matrix
%     un: normal unit vector
%     ua: chordwise unit vector
%   zeta: dimensionless directed differential vortex length vector
%     a0: local section lift-curve slope (1/rad)
% alf_ZL: local section zero-lift angle of attack (deg)
% Outputs:
%   Gsol: gamma solution vector
%      E: residuals at each spanwise station
N = length(G0);

uinf = [cosd(alpha) 0 sind(alpha)];

[Gsol,E] = fsolve(@objfun,G0(1:N/2));
Gsol = [Gsol;flipud(Gsol)];

    function fval = objfun(Ghalf)
        G = [Ghalf;flipud(Ghalf)];
        v = uinf + reshape(G.'*reshape(vij,[N 3*N]),[N 3]);
        alf = atan2d(dot(v,un,2),dot(v,ua,2));
        Cl = a0.*(alf - alf_ZL)*pi/180;
        fval = 2*vecnorm(cross(v(1:N/2,:),zeta(1:N/2,:),2),2,2).*G(1:N/2) - Cl(1:N/2);
    end
end