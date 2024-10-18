%DRAGOPT
clear
close all

N = 50;    %Number of panels
CL_des = 0.46;   %Design wing CL
sref = 2160.53; %Reference area
bref = 135.70;
cavg = sref/bref;   %Mean geometric chord

%Enter end points of wing surface
%xc1 = [0.7294 190.2256 90 90];
%xc2 = [64.442 90.5 95.6 86.422];
%yc = [0 59 59 0];
%zc = [0 4 4 0];
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

theta = zeros(N,1);
s = zeros(N,1);
xle = zeros(N,1);
xte = zeros(N,1);
xvle = zeros(N+1,1);
xvte = zeros(N+1,1);
c = zeros(N,1);
x = zeros(N,1);
a = zeros(N,N);
abar = zeros(N+1,N+1);
b = zeros(N,1);
Re = zeros(N,1);

%Calculate reference span
bref = sref/cavg;

spacing = sin(pi/2*(0:N).'/N);
midsp = 0.5*(spacing(1:N) + spacing(2:N+1));
yv = spacing*(yc(2)-yc(1));
zv = spacing*(zc(2)-zc(1));
y = 0.5*(yv(1:N) + yv(2:N+1));
z = 0.5*(zv(1:N) + zv(2:N+1));
sp = 0.5*sqrt(diff(yv).^2 + diff(zv).^2);
xvle(1) = xc1(1);
xvte(1) = xc1(4);
for i=1:N
    theta(i) = atan2(zc(2)-zc(1),yc(2)-yc(1));
    s(i) = 2.*sp(i)/bref;
    if y(i)<=25.306992  %Yehudi break <- Why defined here?
        xle(i) = xc1(1)+midsp(i)*(xc1(2)-xc1(1));
        xte(i) = xc1(4)+midsp(i)*(xc1(3)-xc1(4));
        xvle(i+1) = xc1(1)+spacing(i+1)*(xc1(2)-xc1(1));
        xvte(i+1) = xc1(4)+spacing(i+1)*(xc1(3)-xc1(4));
    else
        xle(i) = xc2(1)+midsp(i)*(xc2(2)-xc2(1));
        xte(i) = xc2(4)+midsp(i)*(xc2(3)-xc2(4));
        xvle(i+1) = xc2(1)+spacing(i+1)*(xc2(2)-xc2(1));
        xvte(i+1) = xc2(4)+spacing(i+1)*(xc2(3)-xc2(4));
    end
    c(i) = (xte(i)-xle(i));
    Re(i) = c(i)*V/nu;
    x(i) = xle(i)+0.25*c(i);
end

for i=1:N
    for j=1:N
        yp = (y(i)-y(j))*cos(theta(j))+(z(i)-z(j))*sin(theta(j));
        zp = -(y(i)-y(j))*sin(theta(j))+(z(i)-z(j))*cos(theta(j));
        r1 = zp^2+(yp-sp(j))^2;
        r2 = zp^2+(yp+sp(j))^2;

        a1 = ((yp-sp(j))/r1-(yp+sp(j))/r2)*cos(theta(i)-theta(j))+(zp/r1-zp/r2)*sin(theta(i)-theta(j));

        %Assuming symmetry{
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

figure
plot(yv,xvle,'-o')
hold on
plot(yv,xvte,'-o')
plot([y y].',[xle xte].','k-x')

%------------------------------------------------------------ export config
% M = spdiags(zeros(N+1,2)+0.5,-1:0,N+1,N+1);
% M(1,1) = 1;
% 
% out = zeros(N+1,7);
% yv = M \ [0;y];
% out(:,[1 2 4]) = interp1(y,[xle xte z],yv,'linear','extrap');
% out(:,3) = yv;
% 
% M = spdiags(zeros(2*N+1,2)+[-1 1],-1:0,2*N+1,2*N);
% ys = [-flipud(yv);yv(2:N+1)]/yv(N+1);
% ysc = [-flipud(y);y]/yv(N+1);
% G = [gamma(N:-1:1);gamma(1:N)]*V;
% w = 1/pi/bref*sum(M*G./(ysc.'-ys),1).';
% 
% alf_ind = w(N+1:2*N)/V;
% %alf_ind = gamma(1)/(2*bref);
% a0 = zeros(N,1) + 2*pi;
% alf_ZL = zeros(N,1);
% twist = (cn./a0 + alf_ind + alf_ZL*pi/180)*180/pi;
% 
% out(1:N,5:7) = [twist a0 alf_ZL];

% Write to file
%fid = fopen('lltwing.dat','w');
%fprintf(fid,'%13s %13s %13s %13s %13s %13s %13s\r\n','xle','xte','y','z','twist(deg)','a0','alfZL(deg)');
%fprintf(fid,'%13.8f %13.8f %13.8f %13.8f %13.8f %13.8f %13.8f\r\n',out.');
%fid = fclose(fid);