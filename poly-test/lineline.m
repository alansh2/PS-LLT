function [P,A] = lineline(XY1,XY2)
%LINELINE
X13 = XY1(:,1) - XY2(:,1).';
Y13 = XY1(:,2) - XY2(:,2).';
X21 = XY1(:,3) - XY1(:,1);
Y21 = XY1(:,4) - XY1(:,2);
X43 = XY2(:,3).' - XY2(:,1).';
Y43 = XY2(:,4).' - XY2(:,2).';
DNOM = Y43.*X21 - X43.*Y21;
F = (X43.*Y13 - Y43.*X13) ./ DNOM;
G = (X21.*Y13 - Y21.*X13) ./ DNOM;
X = XY1(:,1) + X21.*F;
Y = XY1(:,2) + Y21.*F;
A = (F>=0) & (F<=1) & (G>=0) & (G<=1);
P = [X(A) Y(A)];