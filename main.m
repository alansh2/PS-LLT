%MAIN
clear
close all

geom = importWing('lltwing.dat');

alpha= 0;

[C_L,Cl] = LLT(geom,alpha);

figure
plot(geom.ctrl(:,2),geom.chrd.*Cl)
xlabel('Half-span')
ylabel('cC_l')
title('Spanwise Lift Distribution')