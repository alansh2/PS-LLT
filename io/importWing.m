function [geom,vdo] = importWing(filename)
%IMPORTWING  Assemble LLT wing object from file

%--------------------------------- general-purpose extractor
fid = fopen(filename,'r');
C = textscan(fid,'%[^\n]');
fclose(fid);

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
vdo = A(1:N,1:4); % vertex definition accepted by DragOpt
A = num2cell(A(1:N,:),1);
[xLE,xTE,y,z,geom.twist,geom.a0,geom.alf_ZL] = deal(A{:});
%-- remove padding from values defined at control points
geom.twist(N) = [];
geom.a0(N) = [];
geom.alf_ZL(N) = [];
N = N - 1; % number of panels

%------------------------ construct compatible input for LLT
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
geom.ua = [cosd(geom.twist) zeros(N,1) -sind(geom.twist)];
geom.ua(:,2:3) = geom.ua(:,3).*[-sin(dih) cos(dih)];
geom.us = [zeros(N,1) dy dz]./sqrt(dy.^2 + dz.^2);
geom.un = cross(geom.ua,geom.us,2);

% Mirror for full span
flds = fieldnames(geom);
for i = 1:numel(flds)
    [dim1,dim2] = size(geom.(flds{i}));
    geom.(flds{i}) = [geom.(flds{i})(dim1-(0:N-1),:);geom.(flds{i})];
    if dim2 == 3
        geom.(flds{i})(1:dim1,2) = -geom.(flds{i})(1:dim1,2);
    end
end