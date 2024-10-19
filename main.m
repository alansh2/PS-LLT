%MAIN
clear
close all

geom = importWing('lltwing.dat');

alpha= 0;

%--------------------------------------------- Phillips-Snyder lifting-line
[C_L,Cl] = LLT(geom,alpha);

figure
plot(geom.ctrl(:,2),geom.chrd.*Cl)
xlabel('Half-span')
ylabel('cC_l')
title('Spanwise Lift Distribution')

%--------------------------------------------- classic Prandtl lifting-line
N = size(geom.ctrl,1);
bref = 2*geom.vert(N+1,2);
ynd = geom.ctrl(:,2)*2/bref;
Sref = geom.chrd.' * diff(geom.vert(:,2));
AR = bref^2/Sref;

[~,~,Gamma,C_L,~,~] = PLLT(ynd, ...
                           geom.chrd, ...
                           bref, ...
                           AR, ...
                           geom.a0, ...
                           (geom.twist+alpha)*pi/180, ...
                           geom.alf_ZL*pi/180, ...
                           1);

figure
plot(geom.ctrl(:,2),2*Gamma)
xlabel('Half-span')
ylabel('cC_l')
title('Spanwise Lift Distribution')