%DRAGOPT
clear
close all

N = 200;    %Number of panels
CL_des = 0.46;   %Design wing CL
sref = 2160.53; %Reference area
bref = 135.70;
cavg = sref/bref;   %Mean geometric chord

%Enter end points of wing surface
xc1 = [35.46429217 68.526752 64.36102 64.36102];
xc2 = [35.46429217 68.526752 74.87219 58.16147070];
yc = [0 68.21432 68.21432 0];
zc = [0 4 4 0]; %Why 4?
xcg = 60;
cp = 0.25;
rho = 0.000675954;
mu = 2.99135e-7;
nu = mu/rho;
Minf = 0.78;
sos = 968.076;
V = Minf*sos;

a = zeros(N,N);
abar = zeros(N+1,N+1);
b = zeros(N,1);


%CL_des = 0.6;   %Design wing CL
%sref = 86.1112888448; %Reference area
%bref = 26.24672;
%cavg = sref/bref;   %Mean geometric chord
%xc1 = [0 0 cavg cavg];
%xc2 = [0 0 cavg cavg];
%yc = [0 bref/2 bref/2 0];
%zc = [0 0 0 0];
%V = 262.4672;


%Calculate reference span
bref = sref/cavg;

%-------------------------------------------------------------- parse input
TE = [xc1(4) 0 xc2(3);yc(4) 0 yc(3)];
DNOM = xc1(3) - xc2(3) + xc2(4) - xc1(4);
if DNOM == 0;
    TE(:,2) = [];
else
    TE(4) = ((xc1(3)-xc2(3))*yc(4)+(xc2(4)-xc1(4))*yc(3))/DNOM;
    TE(3) = (xc2(4)*(xc1(3)-xc1(4))-xc1(4)*(xc2(3)-xc2(4)))/DNOM;
end

spacing = sin(pi/2*(0:N).'/N);
vert = zeros(N+1,4);
vert(:,[1 3 4]) = [xc1(1) yc(1) zc(1)] + spacing.*[xc1(2)-xc1(1) yc(2)-yc(1) zc(2)-zc(1)];
vert(:,2) = interp1(TE(2,:),TE(1,:),vert(:,3),'linear','extrap');

dy = diff(vert(:,3));
dz = diff(vert(:,4));
xle = 0.5*(vert(1:N,1) + vert(2:N+1,1));
xte = 0.5*(vert(1:N,2) + vert(2:N+1,2));
y = vert(1:N,3) + 0.5*dy;
z = vert(1:N,4) + 0.5*dz;
sp = 0.5*sqrt(dy.^2 + dz.^2);
s = 2*sp/bref;

c = xte - xle;
Re = c*V/nu;
x = xle + 0.25*c;

theta = atan2(dz,dy);

%-------------------------------------------------------- do inverse design
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

M = N+1;
lda = N+1;
ldb = lda;
nb = 1;

gamma=abar\b;

%---------------------------------------------------------- evaluate result
CL = 0;
CM = 0;
CDi = 0;
CDp = 0;
CDp2 = 0;
cn = zeros(N,1);
cdp = zeros(N,1);
cdp2 = zeros(N,1);
for i=1:N
    %Assume symmetry[
    CL = CL+2.*gamma(i)*s(i)*cos(theta(i));
    CM = CM+2.*gamma(i)*s(i)*cos(theta(i))*(xcg-(xle(i)+cp*c(i)))/cavg;
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

ar = bref^2/sref;
e = CL^2/(pi*ar*CDi);

% Add graph step
plot(y,cn)
xlabel('Half-span, ft')
ylabel('cl')
title('Spanwise Lift Distribution')

%------------------------------------------------------------ export config
M = spdiags(zeros(2*N+1,2)+[-1 1],-1:0,2*N+1,2*N);
ys = sin(pi/2*(-N:N).'/N);
ysc = 0.5*(ys(1:2*N) + ys(2:2*N+1));
G = [gamma(N:-1:1);gamma(1:N)]*V;
w = 1/pi/bref*sum(M*G./(ysc.'-ys),1).'; %downwash velocity

alf_ind = w(N+1:2*N)/V;
%alf_ind = zeros(N,1) + gamma(1)/(2*bref);
a0 = zeros(N,1) + 2*pi;
alf_ZL = zeros(N,1);
twist = (cn./a0 + alf_ind + alf_ZL*pi/180)*180/pi;

% Write to file
%fid = fopen('lltwing.dat','w');
%fprintf(fid,'%13s %13s %13s %13s %13s %13s %13s\r\n','xle','xte','y','z','twist(deg)','a0','alfZL(deg)');
%fprintf(fid,'%13.8f %13.8f %13.8f %13.8f %13.8f %13.8f %13.8f\r\n',[vert [twist a0 alf_ZL;0 0 0]].');
%fid = fclose(fid);