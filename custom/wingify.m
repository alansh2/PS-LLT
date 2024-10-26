function [y,xle,xte] = wingify(P,StableFlag)
%WINGIFY  Wing description from polygon.
%  WINGIFY(P) executes using default parameters.
%  WINGIFY(P,STABLEFLAG) returns coordinates sampled from the
%  polygon P with a choice between two span coordinate
%  distributions. STABLEFLAG can be 1 (default) to preserve
%  the input points or 0 to sample according to a sine
%  distribution.

if nargin==1, StableFlag=1; end % default

%------------------------------------ basic input correction
%-- Coordinates from decsg2 typically result in the root
%-- section being divided into multiple segments
is_zero = find(P(:,2) == 0);
[~,ord] = sort(P(is_zero,1));
if numel(ord)~=2, P(is_zero(ord(2:end-1)),:)=[]; end

%----------------------------------------- cycle point order
N = size(P,1);

[~,k] = min(P(:,2));
kp = mod(k,N)+1;
km = N-mod(1-k,N);
if P(k,2)==P(km,2), kp=k; k=km; end
P = P([kp:N 1:k],:);
if P(1)>P(N), P=flipud(P); end %ensure CW

%------------------------------------------ define LE and TE
l1 = [diff(P(:,2))>=0;false];
l2 = ~l1;

y = unique(P(:,2),'sorted');
if StableFlag==0, y=y(end)*sin(linspace(0,pi/2,numel(y)).'); end

xle = interp1(P(l1,2),P(l1,1),y);
xte = interp1(P(l2,2),P(l2,1),y);
