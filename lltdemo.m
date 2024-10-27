%LLTDEMO
clear
close all

initutils();

filename = mfilename('fullpath');
filepath = fileparts( filename );

%-- use the geometry creation pipeline to design a wing at the command line
N     = 200;            % number of panels
Sref  = 1840;           % reference area
AR    = 10;             % aspect ratio
bref  = sqrt(Sref*AR);  % reference span
taper = 0.23;           % taper ratio
LEswp = 25;             % LE sweep (deg)
dih   = 0;              % dihedral (deg)

cr = 2*Sref/bref/(1+taper);
ct = cr*taper;

ynd = sin(pi/2*(0:N).'/N); % dimensionless half-span coordinate
y   = 0.5*bref*ynd;
xLE = tand(LEswp)*y;
xTE = xLE + cr + (ct-cr)*ynd;

gm1 = [2*N+2;xLE;flipud(xTE);y;flipud(y)];
gm2 = [3;20;28;28;-1;35;-1]; % -1 is a hack around limitations in decsg2

%-- the set formula produces a union between the simple wing
%-- and the yehudi addition
pgon = decsg2([gm1;gm2],'wing+yehudi',{'wing','yehudi'});
pgon(pgon(:,2)<0,2) = 0; % set the -1 due to gm2 to 0 (snap to centerline)

%-- define the LE and TE from the polygon
[y,xle,xte] = wingify(pgon,0);
z = tand(dih)*y;
N = numel(y) - 1; % updated number of panels

%--------------------------------------------------- perform inverse design
CL_des = 0.46; % design CL
Re = 20e6; % aircraft Reynolds number

%-- solve the drag minimization problem for the design CL
[~,alf_ind,cn,~,~,~] = DragOpt([xle xte y z],CL_des,Re);
a0     = zeros(N,1) + 2*pi; % section lift-curve slope
alf_ZL = zeros(N,1); % section zero-lift AoA
twist  = (cn./a0 + alf_ind + alf_ZL*pi/180)*180/pi;

%-- export this config so it can be imported for analysis
exportWing([filepath '/wing-data/wing.dat'],[xle xte y z],[twist a0 alf_ZL]);

%-------------------------------------------------- wing analysis using LLT
%-- VDO is the vertex argument for DragOpt
%-- this allows imported wings to be optimized
[geom,vdo] = importWing([filepath '/wing-data/wing.dat']);

figure
plot(vdo(:,[3 3]),vdo(:,[1 2]))
set(gca,'YDir','reverse','DataAspectRatio',[1 1 1])

alpha = 0;

%-- Phillips-Snyder lifting-line
[C_L,Cl] = LLT(geom,alpha);
disp(C_L)

figure
plot(geom.ctrl(:,2),geom.chrd.*Cl)
xlabel('Half-span')
ylabel('cC_l')
title('Spanwise Lift Distribution')

%-- classic Prandtl lifting-line
N    = size(geom.ctrl,1);
bref = 2*geom.vert(N+1,2);
ynd  = geom.ctrl(:,2)*2/bref;
Sref = geom.chrd.' * diff(geom.vert(:,2));
AR   = bref^2/Sref;

[~,~,Gamma,C_L,~,~] = PLLT(ynd, ...
                           geom.chrd, ...
                           bref, ...
                           AR, ...
                           geom.a0, ...
                           (geom.twist+alpha)*pi/180, ...
                           geom.alf_ZL*pi/180, ...
                           1);
disp(C_L)

figure
plot(geom.ctrl(:,2),2*Gamma)
xlabel('Half-span')
ylabel('cC_l')
title('Spanwise Lift Distribution')