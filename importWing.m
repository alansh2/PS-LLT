function geom = importWing(filename)
%IMPORTWING  Assemble LLT wing object from file
fid = fopen(filename,'r');
C = textscan(fid,'%[^\n]');
fid = fclose(fid);

nrows = size(C{1},1);
A = zeros(nrows,7);
N = 0;
for i = 1:nrows
    a = sscanf(C{1}{i},'%f');
    if ~isempty(a)
        N = N + 1;
        A(N,:) = a;
    end
end
A = num2cell(A(1:N,:),1);
[xLE,xTE,y,z,twist,a0,alf_ZL] = deal(A{:});
twist(N) = []; % remove padding from values defined at control points
a0(N) = [];
alf_ZL(N) = [];
N = N - 1; % number of panels

%--------------------------------------- construct compatible input for LLT
c1 = xTE(1:N) - xLE(1:N);
c2 = xTE(2:N+1) - xLE(2:N+1);
%cbar = 2/3*(c1.^2+c1.*c2+c2.^2)./(c1+c2);
%cbar = [flipud(cbar);cbar];

geom.vert = [0.75*xLE+0.25*xTE y z]; % nodes of hshoe vortices on c/4 line
geom.ctrl = 0.5*(geom.vert(1:N,:) + geom.vert(2:N+1,:));
geom.chrd = 0.5*(c1 + c2);

dy = diff(geom.vert(:,2));
dz = diff(geom.vert(:,3));
dih = atan2(dz,dy); % dihedral

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
geom.alf_ZL = [flipud(alf_ZL);alf_ZL];
geom.twist = [flipud(twist);twist];