function [gamma,alf_ind,cn,CDi,CDp,CDp2] = DragOpt(vert,CL_des,Re)
%DRAGOPT
% Inputs:
%   vert: (N+1)-by-4 array containing [xle xte y z] coordinates at nodal half-span positions
% CL_des: design lift coefficient
%     Re: aircraft Reynolds number

%---------------------------------------------- basic checks
if size(vert,2) ~= 4
    error('DragOpt:invalidInput','Invalid vertex definition.')
end
if (~isscalar(CL_des) || ...
    ~isscalar(Re))
    error('DragOpt:invalidInput','CL_des and Re must be scalar.')
end

%-------------------------------------------- parse geometry
N = size(vert,1) - 1; %Number of panels
bref = 2*vert(N+1,3); %Reference span

dy = diff(vert(:,3));
dz = diff(vert(:,4));
xle = 0.5*(vert(1:N,1) + vert(2:N+1,1));
xte = 0.5*(vert(1:N,2) + vert(2:N+1,2));
y = vert(1:N,3) + 0.5*dy;
z = vert(1:N,4) + 0.5*dz;
sp = 0.5*sqrt(dy.^2 + dz.^2);
s = 2*sp/bref;
c = xte - xle;
x = xle + 0.25*c;
theta = atan2(dz,dy);

%--------------------------------------- calculate constants
sref = c.' * dy; %Reference area
cavg = sref/bref; %Mean geometric chord
cmac = 2/sref*trapz((vert(:,2)-vert(:,1)).^2,vert(:,3));
Re = Re*c/cmac;

%----------------------------------------- do inverse design
a = zeros(N,N);
abar = zeros(N+1,N+1);
b = zeros(N,1);

for i=1:N
    for j=1:N
        yp = (y(i)-y(j))*cos(theta(j))+(z(i)-z(j))*sin(theta(j));
        zp = -(y(i)-y(j))*sin(theta(j))+(z(i)-z(j))*cos(theta(j));
        r1 = zp^2+(yp-sp(j))^2;
        r2 = zp^2+(yp+sp(j))^2;

        a1 = ((yp-sp(j))/r1-(yp+sp(j))/r2)*cos(theta(i)-theta(j))+(zp/r1-zp/r2)*sin(theta(i)-theta(j));

        %Assuming symmetry[
        yp = (y(i)+y(j))*cos(-theta(j))+(z(i)-z(j))*sin(-theta(j));
        zp = -(y(i)+y(j))*sin(-theta(j))+(z(i)-z(j))*cos(-theta(j));
        r1 = zp^2+(yp-sp(j))^2;
        r2 = zp^2+(yp+sp(j))^2;
        a2 = ((yp-sp(j))/r1-(yp+sp(j))/r2)*cos(theta(i)+theta(j))+(zp/r1-zp/r2)*sin(theta(i)+theta(j));
        %]

        a(i,j) = -cavg/(4.*pi)*(a1 + a2);
    end
end

for i=1:N
    for j=1:N
        abar(i,j) = a(i,j)*s(i) + a(j,i)*s(j);
    end
    abar(i,j) = abar(i,j) + (3.15183643E-20*Re(i)^2 - 1.07237288E-11*Re(i) + 6.00306257E-03)+(0.00172387*(cavg/c(i))^2);

    b(i) = 0;
end

for i=1:N
    abar(i,N+1) = s(i)*cos(theta(i));
end
for j=1:N
    abar(N+1,j) = s(j)*cos(theta(j));
end
abar(N+1,N+1) = 0;
%Assume symmetry [
b(N+1) = CL_des/2;
%]

gamma=abar\b;

%------------------------------------------- evaluate result
CL = 0;
CDi = 0;
CDp = 0;
CDp2 = 0;
cn = zeros(N,1);
cdp = zeros(N,1);
cdp2 = zeros(N,1);
for i=1:N
    %Assume symmetry[
    CL = CL+2.*gamma(i)*s(i)*cos(theta(i));
    %]

    cn(i) = gamma(i)*cavg/c(i);
    cdp(i) = ((3.15183643E-20*Re(i)^2 - 1.07237288E-11*Re(i) + 6.00306257E-03));%+(0.00172387*((gamma(i)*cavg/c(i)))^2));
    cdp2(i) = ((3.15183643E-20*Re(i)^2 - 1.07237288E-11*Re(i) + 6.00306257E-03)+(0.00172387*((gamma(i)*cavg/c(i)))^2));
    for j=1:N
        %Assume symmetry[
        CDi = CDi+gamma(i)*gamma(j)*s(i)*a(i,j);
        %]
    end
    CDp = CDp+((3.15183643E-20*Re(i)^2 - 1.07237288E-11*Re(i) + 6.00306257E-03))*s(i)*2;%+(0.00172387*((gamma(i)*cavg/c(i)))^2))*c(i)/(sref/2);
    CDp2 = CDp2+((3.15183643E-20*Re(i)^2 - 1.07237288E-11*Re(i) + 6.00306257E-03)+(0.00172387*((gamma(i)*cavg/c(i)))^2))*s(i)*2;
end

alf_ind = gamma(1)*cavg/4/bref;

% Add graph step
plot(y,cn)
xlabel('Half-span, ft')
ylabel('cl')
title('Spanwise Lift Distribution')
