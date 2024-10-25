function [y,xle,xte] = wingify(P)
%WINGIFY  Wing description from polygon.
N = size(P,1);

[~,k] = min(P(:,2));
kp = mod(k,N)+1;
km = N-mod(1-k,N);
if P(k,2)==P(km,2), kp=k; k=km; end
P = P([kp:N 1:k],:);
if P(1)>P(N), P=flipud(P); end %ensure CW

l1 = [diff(P(:,2))>=0;false];
l2 = ~l1;

y = unique(P(:,2),'sorted');

xle = interp1(P(l1,2),P(l1,1),y);
xte = interp1(P(l2,2),P(l2,1),y);
