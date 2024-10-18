function geom = designWing(ys,Sref,bref,taper,LEswp,dih,twist,a0,alfZL)
%DESIGNWING  Create wing object usable in LLT analysis
% Inputs:
%    ys: dimensionless spanwise coordinate
%  Sref: reference area
%  bref: reference wingspan
% taper: taper ratio
% LEswp: leading-edge sweep angle (deg)
%   dih: dihedral angle (deg)
% twist: twist distribution (deg) evaluated at midpoints between ys
%    a0: local section lift-curve slope (1/rad) specified as a scalar or vector
% alfZL: local section zero-lift angle of attack (deg)
% Output:
%  geom: wing object (see LLT)
%
% See also LLT, IMPORTWING.
N = length(ys) - 1;

% Guarantee appropriate variable sizes
ys    = reshape(ys,N+1,1);
twist = reshape(twist,N,1);
a0    = reshape(a0,[],1) + zeros(N,1);
alfZL = reshape(alfZL,[],1) + zeros(N,1);

cr = 2*Sref/bref/(1+taper);
ct = cr*taper;

y = 0.5*bref*ys;
xLE = tand(LEswp)*y;
xTE = xLE + cr + (ct-cr)*ys;
z = tand(dih)*y;

c1 = xTE(1:N) - xLE(1:N);
c2 = xTE(2:N+1) - xLE(2:N+1);

geom.vert = [0.75*xLE+0.25*xTE y z]; % nodes of hshoe vortices on c/4 line
geom.ctrl = 0.5*(geom.vert(1:N,:) + geom.vert(2:N+1,:));
geom.chrd = 0.5*(c1 + c2);

dy = diff(geom.vert(:,2));
dz = diff(geom.vert(:,3));

% Aligned unit vectors
ua = [cosd(twist) zeros(N,1) -sind(twist)];
ua(:,2:3) = ua(:,3).*[-sin(dih) cos(dih)];
us = [zeros(N,1) dy dz]./sqrt(dy.^2 + dz.^2);
un = cross(ua,us,2);

% Mirror for full span
geom.vert = [flipud(geom.vert);geom.vert(2:N+1,:)];
geom.vert(1:N,2) = -geom.vert(1:N,2);
geom.ctrl = [flipud(geom.ctrl);geom.ctrl];
geom.ctrl(1:N,2) = -geom.ctrl(1:N,2);
geom.chrd = [flipud(geom.chrd);geom.chrd];
geom.ua = [flipud(ua);ua];
geom.ua(1:N,2) = -geom.ua(1:N,2);
geom.us = [flipud(us);us];
geom.us(1:N,2) = -geom.us(1:N,2);
geom.un = [flipud(un);un];
geom.un(1:N,2) = -geom.un(1:N,2);
geom.a0 = [flipud(a0);a0];
geom.alf_ZL = [flipud(alfZL);alfZL];