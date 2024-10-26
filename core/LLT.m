function [C_L,Cl] = LLT(geom,alpha)
%LLT  Analyze wing using lifting-line method
% Inputs:
%  geom: wing object with fields
%         -vert: N+1-by-3 vector of horseshoe vortex nodes
%         -ctrl: N-by-3 vector of control points
%         -chrd: local section chord length
%         -ua: chordwise unit vector
%         -un: normal unit vector
%         -us: spanwise unit vector
%         -a0: local section lift-curve slope (1/rad)
%         -alf_ZL: local section zero-lift angle of attack (deg)
% alpha: angle of attack (deg)
% Outputs:
%   C_L: wing lift coefficient
%    Cl: sectional lift coefficient distribution
%
% Reference:
% Phillips, W. F., and Snyder, D. O., "Modern Adaptation of Prandtl's Classic
% Lifting-Line Theory," Journal of Aircraft, Vol. 37, No. 4, 2000, pp. 662-670.
N = size(geom.ctrl,1);

uinf = [cosd(alpha) 0 sind(alpha)];

dl = diff(geom.vert);
dA = diff(geom.vert(:,2)).*geom.chrd; % panel areas
zeta = dl./dA;

%--------------------------------------------------- build influence matrix
vij = zeros(N,N,3);
for i = 1:N
    for j = 1:N
        r1 = geom.ctrl(j,:) - geom.vert(i,:);
        r2 = geom.ctrl(j,:) - geom.vert(i+1,:);
        r_1 = sqrt(dot(r1,r1));
        r_2 = sqrt(dot(r2,r2));
        if i == j
            % straight vortex segment induces no downwash along its own length
            vij(i,j,:) = (cross(uinf,r2)/(r_2*(r_2-dot(uinf,r2))) ...
                -cross(uinf,r1)/(r_1*(r_1-dot(uinf,r1))))/(4*pi);
        else
            vij(i,j,:) = (cross(uinf,r2)/(r_2*(r_2-dot(uinf,r2))) ...
                +(r_1+r_2)*cross(r1,r2)/(r_1*r_2*(r_1*r_2+dot(r1,r2))) ...
                -cross(uinf,r1)/(r_1*(r_1-dot(uinf,r1))))/(4*pi);
        end
    end
end

%----------------------------- solve linear problem to use as initial guess
A = diag(2*vecnorm(cross(repmat(uinf,N,1),zeta,2),2,2),0) ...
    -geom.a0.*(sum(vij.*reshape(geom.un,[1 N 3]),3).');
RHS = geom.a0.*(sum(uinf.*geom.un,2) - geom.alf_ZL*pi/180);
G0 = A \ RHS;

%--------------------------------------------------- solve the full problem
[G,E] = solveGamma(G0,alpha,vij,geom.un,geom.ua,zeta,geom.a0,geom.alf_ZL);

Cl = 2*vecnorm(cross(uinf + reshape(G.'*reshape(vij,[N 3*N]),[N 3]), zeta, 2),2,2).*G;

%Vj = uinf + reshape(G.'*reshape(vij,[N 3*N]),[N 3]);
%alf = atan2(dot(Vj,geom.un,2),dot(Vj,geom.ua,2));
%Cl = 2*pi*(alf - alf_ZL); % should be the same as other calc

C_L = Cl.'*dA/sum(dA);