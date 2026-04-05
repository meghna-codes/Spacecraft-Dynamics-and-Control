close all;

global mu FORWARD a tol Lp
FORWARD = 1; 
a=1;
tol = 1e-12;

%% System parameters
m_earth = 5.97219e24;
m_moon  = 7.34767309e22;
mu = m_moon / (m_earth + m_moon);

%% Choose Lagrange Point (L1 or L2)
Lp = 1;        % 1 for L1, 2 for L2
Az = 0.2;

% Choosing northern (n=-1) or southern halo (n=1)
n=-1;
% Get third-order Richardson guess for an initial condition
x0 = richardson(mu,Lp,Az,n);

[xL, yL] = Lagrange(mu, Lp);

tinitial = 0;
tfinal = 2*pi; 
[x,t]=integrate(x0,tinitial,tfinal);

% we need it to intersect the y = 0 plane twice 
% Use differential correction to get better approximation

CASE = 1; % Keep z(0) constant, let x(0) and ydot(0) vary
%CASE = 2; %Keep x(0) constant, let z(0) and ydot(0) vary

[XH,TH]=diffcor(x0,CASE); 
% This gives us the initial condition XH and half-period TH
    
% Now integrate again and plot using determined initial condition
x0=XH;
T=2*TH; % TH = half period
 
[x,t]=integrate(x0,0,T); 

[xx,tt,phi_T,PHI]=getPHI(x0,T);

%controlled orbit data
xc = squeeze(out.x);
yc = squeeze(out.y);
zc = squeeze(out.z);
tc = squeeze(out.t);

figure(1);
am=384000;
plot3(x(:,1)*am,x(:,2)*am,x(:,3)*am,'b--','linewidth',2.5);
hold on
plot3(xc*am, yc*am, zc*am, 'r-', 'LineWidth', 1.5)
plot3(xL*am, yL*am, 0,'ro','MarkerSize',6,'MarkerFaceColor','r')


Rmoon =  1737.4; %km
EMdist= 385000;   %km
% make a sphere with a nondimensional radius of the Moon
[Xm,Ym,Zm] = sphere;
Xm = Xm*Rmoon + (1-mu)*am;
Ym = Ym*Rmoon ;
Zm = Zm*Rmoon ;
surf(Xm,Ym,Zm)

title('Reference and Controlled halo orbits')
xlabel('x (km)'); 
ylabel('y (km)'); 
zlabel('z (km)')
legend('Reference', 'Controlled')
axis equal; 
grid on; view(3)
set(gca,'fontsize',12)

%% calculating the state transition matrix PHI(t,0)
function [x,t,phi_T,PHI]=getPHI(x0,tf)
 
% Gets state transition matrix, phi_T, and the trajectory (x,t) for a length 
% of time, tf (2*pi is 1 period). In particular, for periodic solutions of 
% period tf=T, one can obtain the monodromy matrix, PHI(0,T).
 
global tol
OPTIONS = odeset('RelTol',1e-11,'AbsTol',1e-12);
PHI_0(1:36) = reshape(eye(6),36,1);
PHI_0(1+36:6+36) = x0;
[t,PHI] = ode78(@(t,PHI) variational(t,PHI), [0 tf], PHI_0, OPTIONS); 
x = PHI(:,1+36:6+36); % trajectory
phi_T = reshape(PHI(length(t),1:36),6,6); % monodromy matrix, PHI(O,T)
end

function [PHIdot] = variational(t,PHI)

global mu FORWARD

mu2 = 1-mu;

x(1:6) = PHI(37:42);
phi = reshape(PHI(1:36), 6, 6);

r2 = (x(1)+mu)^2 + x(2)^2 + x(3)^2;    
R2 = (x(1)-mu2)^2 + x(2)^2 + x(3)^2;    
r3 = r2^1.5; r5= r2^2.5;
R3 = R2^1.5; R5= R2^2.5;

Uxx = 1+(mu2/r5)*(3*(x(1)+mu)^2)+(mu/R5)*(3*(x(1)-mu2)^2)-(mu2/r3+mu/R3);
Uyy = 1+(mu2/r5)*(3*x(2)^2)+(mu/R5)*(3* x(2)^2)-(mu2/r3+mu/R3);
Uzz = (mu2/r5)*(3*x(3)^2)+(mu/R5)*(3* x(3)^2)-(mu2/r3+mu/R3);
Uxy = 3*x(2)*(mu2*(x(1)+mu)/r5+mu*(x(1)-mu2)/R5); 
Uxz = 3*x(3)*(mu2*(x(1)+mu)/r5+mu*(x(1)-mu2)/R5); 
Uyz = 3*x(2)*x(3)*(mu2/r5+mu/R5);

% Df is the jacobian matrix of the vector field f
Df  = [  0     0     0     1     0     0 ; 
         0     0     0     0     1     0 ; 
         0     0     0     0     0     1 ;
         Uxx   Uxy   Uxz   0     2     0 ; 
         Uxy   Uyy   Uyz   -2    0     0 ;
         Uxz   Uyz   Uzz   0     0     0  ];

phidot = Df * phi;
PHIdot        = zeros(42,1);
PHIdot(1:36)  = reshape(phidot, 36, 1);
PHIdot(37)    = x(4);
PHIdot(38)    = x(5);
PHIdot(39)    = x(6);
PHIdot(40)    = x(1)-(mu2*(x(1)+mu)/r3) -(mu*(x(1)-mu2)/R3) + 2*x(5);
PHIdot(41)    = x(2)-(mu2*x(2)/r3) - (mu*x(2)/R3) - 2*x(4);
PHIdot(42)    =     -(mu2*x(3)/r3) - (mu*x(3)/R3);
PHIdot(37:42) = PHIdot(37:42)*FORWARD;
end

%% Find an x-axis crossing

function [t1_z,x1_z]=find0(x0)

% From the given state x0, this finds the next x-crossing of the halo orbit
global t0_z x0_z x1_zgl
tolzero = 1.e-12;
options = optimset('DiffMaxChange',tolzero,'TolCon',tolzero,'TolFun',tolzero);

t0_z = pi/2 - 0.15;
t0 = 0;
[xx,~] = integrate(x0,t0,t0_z);
x0_z = xx(end,:);
t1_z = fzero(@haloy,t0_z,options);
x1_z = x1_zgl;

clear global t0_z x0_z x1_zgl

end

%% Integrating from t0 to tf with ICs
function [x,t] = integrate(x0,t0,tf)

%  Integrate x w/ initial condition x0 from time t0 to tf
 
global FORWARD tol
OPTIONS = odeset('RelTol',1e-11,'AbsTol',1e-12);

% 3D equations of motion integrated
[t,x] = ode78(@(t,x)  cr3bp(t,x), [t0 tf], x0, OPTIONS) ;
 
t=FORWARD*t;    
end

%% Circular restricted 3 body-problem equations
function xdot = cr3bp(t,x)

global mu FORWARD

mu1 = 1-mu; % mass of larger  primary 
mu2 = mu;   % mass of smaller primary 
r3 = ((x(1)+mu2)^2 + x(2)^2 + x(3)^2)^1.5; % r: distance to m1, LARGER MASS    
R3 = ((x(1)-mu1)^2 + x(2)^2 + x(3)^2)^1.5; % R: distance to m2, smaller mass  

xdot = zeros(6,1);
xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);
xdot(4) = x(1)-(mu1*(x(1)+mu2)/r3) -(mu2*(x(1)-mu1)/R3) + 2*x(5);
xdot(5) = x(2)-(mu1* x(2)/r3) -(mu2*x(2)/R3) - 2*x(4);
xdot(6) =     -(mu1* x(3)/r3) - (mu2*x(3)/R3);
xdot    = FORWARD * xdot;

end

%% for computing the y-position of the halo orbit
function y1=haloy(t1)

global t0_z x0_z x1_zgl
if t1==t0_z
  x1_zgl = x0_z ; 
else
  [xx,tt]=integrate(x0_z,t0_z,t1);
  x1_zgl=xx(end,:);
end
y1 =x1_zgl(2);
end
 
%% Differential correction step
function [XH,TH]=diffcor(X0,CASE)

 % Initial Values: (x0, 0, z0,  0, dy0,  0)
 % Final Values  : (x1, 0, z1, dx1, dy1, dz1)
 % Use diff corr to kill dx1, dz1 to get symmetric halos.
 % CASE 1: fix z0
 % CASE 2: fix x0

global mu tol

 mu2 = 1-mu;
 t0 = 0;
 c=['r' 'b' 'g' 'm' ]; clf; 
 hold on;
 axis('equal');

 dx1 = 1;
 maxdx1=1.e-12;
 iter= 0;
 maxiter=25;

 while abs(dx1) > maxdx1
 if iter > maxiter
        ERROR = 'OCCURRED'
 end
   N=rem(iter,3)+1;

   [t1,xx1]=find0(X0);
   
   iter=iter+1;

   x1 = xx1(1); y1 = xx1(2); z1 = xx1(3);
   dx1 = xx1(4); dy1 = xx1(5); dz1 = xx1(6);
    
   [x,t,phi,PHI]=getPHI(X0,t1);
    
   rho1 = 1/((x1+mu)^2 + y1^2 + z1^2)^1.5;
   rho2 = 1/((x1-mu2)^2 + y1^2 + z1^2)^1.5;

   Ux1 = -(mu2*(x1+mu)*rho1) -(mu*(x1-mu2)*rho2) + x1;
 
   ddotz1 = -(mu2*z1*rho1) -(mu*z1*rho2);
   ddotx1 = 2*dy1+Ux1;

 if CASE == 1
 %  Fix z0
    C1   = [phi(4,1) phi(4,5); phi(6,1) phi(6,5)];
    C2   = (C1 -(1/dy1)*[ddotx1 ddotz1]'*[phi(2,1) phi(2,5)]);
    C3   = inv(C2)*[-dx1 -dz1]';
    dx0  = C3(1);
    ddoty0 = C3(2);
    X0(1)= X0(1) + dx0;
    X0(5)= X0(5) + ddoty0;

 elseif CASE == 2
 %  Fix x0
    C1   = [phi(4,3) phi(4,5); phi(6,3) phi(6,5)];
    C2   = (C1 - (1/dy1)*[ddotx1 ddotz1]'*[phi(2,3) phi(2,5)]);
    C3   = inv(C2)*[-dx1 -dz1]';
    dz0  = C3(1);
    ddoty0 = C3(2);
    X0(3)= X0(3) + dz0;
    X0(5)= X0(5) + ddoty0;
 end

 XH=X0;
 TH=t1;
 PERIOD=2*TH;

 end

end

%% initial conditions for halo orbit using richardson's analytic approxn.
function x0 = richardson(mu,Lpt,Azlp,n)
% Gives initial state (position,velocity) for a 3D periodic halo orbit 
% centered on the specified collinear Lagrange point.
% [Uses Richardson's 3rd order model for analytically constructing a 3D 
%  periodic halo orbit about the points L1, L2, or L3]
% 
% output: x0 = (r0,v0) initial conditions ("initial guess") for the
%           desired halo orbit (in 3D CR3BP nondim. units)
%
% input: mu   = mass parameter of system [ mu = m2/(m1+m2) ]
%   Lpt  = 1,2, or 3, the number of the specified collinear Lagrange point
%   Azlp = out-of-plane (or z-amplitude) of the desired halo 
%            [in Lpt-primary distances] 
%   n    = +1 is northern halo (z>0, Class I), 
%        = -1 is southern halo (z<0, Class II)
%
%------------------------------------------------------------------------

Az = Azlp;

gamma = gammaL(mu, Lpt);
if     Lpt == 1, won = +1; primary = 1-mu;
elseif Lpt == 2, won = -1; primary = 1-mu; 
end

for N = 2:4
    c(N) = (1/gamma^3)*( (won^N)*mu + ...
	       ((-1)^N)*((primary)*gamma^(N+1))/((1+(-won)*gamma)^(N+1)));
end

polylambda = [ 1 0 (c(2)-2) 0 -(c(2)-1)*(1+2*c(2)) ];
lambda = roots(polylambda); % lambda = frequency of orbit

lambda = abs(lambda(1));

k = 2*lambda/(lambda^2 + 1 - c(2));
del = lambda^2 - c(2);

d1 = ((3*lambda^2)/k)*(k*(6*lambda^2 - 1) - 2*lambda);
d2 = ((8*lambda^2)/k)*(k*(11*lambda^2 -1) - 2*lambda);

a21 = (3*c(3)*(k^2 - 2))/(4*(1 + 2*c(2)));
a22 = 3*c(3)/(4*(1 + 2*c(2)));
a23 = -(3*c(3)*lambda/(4*k*d1))*( 3*(k^3)*lambda - 6*k*(k-lambda) + 4);
a24 = -(3*c(3)*lambda/(4*k*d1))*( 2 + 3*k*lambda );

b21 = -(3*c(3)*lambda/(2*d1))*(3*k*lambda - 4);
b22 = 3*c(3)*lambda/d1;
d21 = -c(3)/(2*lambda^2);

a31 = -(9*lambda/(4*d2))*(4*c(3)*(k*a23 - b21) + k*c(4)*(4 + k^2)) + ...
       ((9*lambda^2 + 1 - c(2))/(2*d2))*(3*c(3)*(2*a23 - k*b21) + c(4)*(2 + 3*k^2));
a32 = -(1/d2)*( (9*lambda/4)*(4*c(3)*(k*a24 - b22) + k*c(4)) + ...
       1.5*(9*lambda^2 + 1 - c(2))*( c(3)*(k*b22 + d21 - 2*a24) - c(4)));

b31 = (.375/d2)*( 8*lambda*(3*c(3)*(k*b21 - 2*a23) - c(4)*(2 + 3*k^2)) + ...
      (9*lambda^2 + 1 + 2*c(2))*(4*c(3)*(k*a23 - b21) + k*c(4)*(4 + k^2)));
b32 = (1/d2)*( 9*lambda*(c(3)*(k*b22 + d21 - 2*a24) - c(4)) + ...
      0.375*(9*lambda^2 + 1 + 2*c(2))*(4*c(3)*(k*a24 - b22) + k*c(4)));

d31 = (3/(64*lambda^2))*(4*c(3)*a24 + c(4));
d32 = (3/(64*lambda^2))*(4*c(3)*(a23 - d21) + c(4)*(4 + k^2));

s1 = (1/(2*lambda*(lambda*(1+k^2) - 2*k)))*(1.5*c(3)*(2*a21*(k^2 - 2) - ...
     a23*(k^2 + 2) - 2*k*b21) - 0.375*c(4)*(3*k^4 - 8*k^2 + 8) );
s2 = (1/(2*lambda*(lambda*(1+k^2) - 2*k)))*(1.5*c(3)*(2*a22*(k^2 - 2) + ...
     a24*(k^2 + 2) + 2*k*b22 + 5*d21) + 0.375*c(4)*(12 - k^2) );

a1 = -1.5*c(3)*(2*a21+ a23 + 5*d21) - 0.375*c(4)*(12-k^2);
a2 = 1.5*c(3)*(a24-2*a22) + 1.125*c(4);

l1 = a1 + 2*(lambda^2)*s1;
l2 = a2 + 2*(lambda^2)*s2;

tau1 = 0;
deltan = -n;

Ax = sqrt((-del - l2*Az.^2)/l1);

x = a21*Ax.^2 + a22*Az.^2 - Ax.*cos(tau1) + (a23*Ax.^2 - a24*Az.^2)*cos(2*tau1) ...
    + (a31*Ax.^3 - a32*Ax.*Az.^2)*cos(3*tau1);
y = k*Ax.*sin(tau1) + (b21*Ax.^2 - b22*Az.^2)*sin(2*tau1) + ...
    (b31*Ax.^3 - b32*Ax.*Az.^2)*sin(3*tau1);
z = deltan*Az.*cos(tau1) + deltan*d21*Ax.*Az.*(cos(2*tau1) - 3) + ...
    deltan*(d32*Az.*Ax.^2 - d31*Az.^3)*cos(3*tau1);
xdot = lambda*Ax.*sin(tau1) - 2*lambda*(a23*Ax.^2 - a24*Az.^2)*sin(2*tau1) ...
     - 3*lambda*(a31*Ax.^3 - a32*Ax.*Az.^2)*sin(3*tau1);
ydot = lambda*(k*Ax.*cos(tau1) + 2*(b21*Ax.^2 - b22*Az.^2)*cos(2*tau1) ...
     + 3*(b31*Ax.^3 - b32*Ax.*Az.^2)*cos(3*tau1));
zdot = - lambda*deltan*Az.*sin(tau1) - 2*lambda*deltan*d21*Ax.*Az.*sin(2*tau1) ...
     - 3*lambda*deltan*(d32*Az.*Ax.^2 - d31*Az.^3)*sin(3*tau1);

r0 = gamma*[ (primary + gamma*(-won +x))/gamma -y z ]';
v0 = gamma*[ xdot ydot zdot ]';
x0= [r0;v0];

end

%% to calculate non-dimensional distance between smaller primary and Lp
function [gamma]=gammaL(mu, Lpt)
% Calculate ratio of libration point distance from closest primary 
% to distance between two primaries  (gammaL1 = (E-L1)/AU)

mu2 = 1 - mu;

poly1 = [1  -1*(3-mu)  (3-2*mu)  -mu   2*mu  -mu ];
poly2 = [1   (3-mu)    (3-2*mu)  -mu  -2*mu  -mu ];
poly3 = [1   (2+mu)    (1+2*mu)  -mu2 -2*mu2 -mu2];

rt1 = roots(poly1); rt2 = roots(poly2); rt3 = roots(poly3);

for k=1:5
        if isreal(rt1(k)) GAMMAS(1)=rt1(k); end
        if isreal(rt2(k)) GAMMAS(2)=rt2(k); end
        if isreal(rt3(k)) GAMMAS(3)=rt3(k); end
end
gamma=GAMMAS(Lpt);

end

%% calculating the position of LP
function [x,y] = Lagrange(mu, Lpt)

% Calculate location of libration point Lpt
%
% see gammaL() 

mu2=1-mu;

LPTS = [(mu2-gammaL(mu,1)) (mu2+gammaL(mu,2)) (-mu-gammaL(mu,3))];
x=LPTS(Lpt);
y=0;

end

%% plotX for X projection of the orbit (data from various sources)
function plotX(X,color,v,LINEWIDTH)

% Plot planar (rotating coord frame) X (i.e., X(:,1) vs. X(:,2)) 
% in the units specified with a 
% v==1 means use the max and min values of X to set the axes

global mu Lpt a

% mu's, semimajor axes (km), and periods (seconds) of the nine Sun-planet systems
MU(1)=1.6601365196665584e-07; A(1)=5.79091350e+07; PRD(1)=7.600551833396541e+06;
MU(2)=2.4478336060627946e-06; A(2)=1.08208880e+08; PRD(2)=1.941415264633838e+07;
MU(3)=3.0034809249852780e-06; A(3)=1.49597927e+08; PRD(3)=3.155843300498465e+07;
MU(4)=3.2271504e-07         ; A(4)=2.27941040e+08; PRD(4)=5.935642579611749e+07;
MU(5)=9.5367840741663007e-04; A(5)=7.78328370e+08; PRD(5)=3.743416297595286e+08;
MU(6)=2.8568275263358739e-04; A(6)=1.42699081e+09; PRD(6)=9.300256264523928e+08;
MU(7)=4.3717116587359953e-05; A(7)=2.86958619e+09; PRD(7)=2.652599569683743e+09;
MU(8)=5.1614995409620897e-05; A(8)=4.49656230e+09; PRD(8)=5.188094805532950e+09;
MU(9)=7.4074073525377411e-09; A(9)=5.89021380e+09; PRD(9)=7.796863629391328e+09;

% Earth-Moon system
MU(10)=1.2150668300000000e-02;A(10)=3.85000000e+05;PRD(10)=2.360591000000000e+06;
MU(31)=1.2150668300000000e-02;A(31)=3.85000000e+05;PRD(31)=2.360593100000000e+06;
MU(32)=1.2150000000000000e-02;A(32)=3.85000000e+05;PRD(32)=2.360593100000000e+06;

% Sun-(Earth-Moon barycenter) system
MU(11)=3.03591e-6;A(11)=A(3);PRD(11)=PRD(3);

% Jupiter satellite system (data from _The New Solar System_ (Beatty) )
MU(51)= 4.684e-5 ; A(51)=421600  ; PRD(51)=1.52853510e+05; % Io
MU(52)= 2.523e-5 ; A(52)=670900  ; PRD(52)=3.06822050e+05; % Europa
MU(53)= 7.802e-5 ; A(53)=1070000 ; PRD(53)=6.18153390e+05; % Ganymede
MU(54)= 5.661e-5 ; A(54)=1883000 ; PRD(54)=1.44193116e+06; % Callisto

% Saturn satellite system (data from _The New Solar System_ (Beatty) )
MU(61)=6.592046007940774e-08; A(61)=186000 ; PRD(61)=8.184448297794924e+04; % Mimas
MU(62)=2.409471484174738e-04; A(62)=1221850; PRD(62)=1.377992433773748e+06; % Titan

% mu's found in Belbruno [1994] for study of ballistic capture by moon
% using Sun-Earth-Moon dynamics
MU(12)=0.01214  ;A(12)=A(10);PRD(12)=PRD(10);
MU(13)=3.0359e-6;A(13)=A(3); PRD(13)=PRD(3);

%%%%%%%%%%%%%%%%%%%%%%

if nargin<=3
    LINEWIDTH = 1 ;
    if nargin<=2
        v=0;
        if nargin==1
            color='r';
        end
    end
end

plot(a*X(:,1),a*X(:,2),color,'linewidth',LINEWIDTH); 
axis equal; axis(axis); AXIS=axis; hold on;
if     mu==MU(10) || mu==MU(32) 
    earthmun; % draw Earth and Moon if Earth-Moon system
elseif mu==MU(3)  || mu==MU(11)
    earthmun(1); % draw lunar orbit around Earth
else
    plot(a*[-mu],[0],'k*',a*[1-mu],[0],'ko');
end
% draw m1 and m2

[xL,~] = Lagrange(mu,Lpt);
plot(a*xL,zeros(size(a*xL)),'gx'); 

circle(a*(1-mu),[0 0],'k:'); circle(a*(mu),[0 0],'k:'); 
zoom on; grid on;axis equal; 


[MM,NN]=size(X);

if v==1
    AXIS=[min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))];
    axis normal; axis(AXIS); axis(axis); axis equal;
end

% Units
if     a>1,units='AU';if a>1000,units='km'; end
elseif a==1 
if     mu==MU(3) || mu==MU(11), units='AU'; 
elseif mu==MU(5), units='Sun-Jupiter distances'; 
elseif mu==MU(10)|| mu==MU(32), units='Earth-Moon distances'; frame=' Earth-Moon';
else, units='nondimensional units'; end;end

% Frames
if     mu==MU(3) || mu==MU(11), frame=' Sun-Earth';
elseif mu==MU(5), frame=' Sun-Jupiter';
elseif mu==MU(10)|| mu==MU(32), frame=' Earth-Moon';
else, frame=''; end

xlabel(sprintf('x (%s,%s rotating frame)',units,frame));
ylabel(sprintf('y (%s,%s rotating frame)',units,frame));
title(sprintf('mu = %.3e, a = %.3e %s (%s rotating frame )',mu,a,units,frame));

end


%% circle - customised 
function circle(r,c,color,FILL,LineWidth)
% circle(r,c,color,FILL,LineWidth);
%
% Make circle of with specified 'color' of radius r 
% at center c=[x y] for a 2D plot
%           c=[x y z] for a 3D plot
%
% if FILL=1, the fill circle

if nargin<=4, LineWidth=1; end
if nargin<=3, FILL=0; end
if nargin<=2, color ='w:';if nargin==1,c=[0 0];end; end
 
theta = 0:.0025:2*pi; theta = theta(:);

m=length(c); for k = 1:length(theta), C(k,1:m)=c; end 

if FILL==0
    if     m==2, circ = plot (C(:,1)+r*cos(theta),C(:,2)+r*sin(theta),color);
    elseif m==3, circ = plot3(C(:,1)+r*cos(theta),C(:,2)+r*sin(theta),C(:,3),color);
    set(circ,'LineWidth',LineWidth)
    end
elseif FILL==1
    fill(C(:,1)+r*cos(theta),C(:,2)+r*sin(theta),color);
end

end

%% earth-moon plotting 
function earthmun(SEM)
% earthmun(SEM)
%
% This will draw disks of the Earth and Moon in the rotating frame
% if SEM=0, and the lunar orbit in the Sun-Earth frame if SEM=1
%

global mu a

if nargin==0, SEM=0; end

Rearth=  6378.14;%km
Rmoon =  1737.4; %km
EMdist=385000;   %km
AU=149597927 ;   %km
theta = 0 : .01 : 2*pi ;

if     SEM==0, % Earth-Moon system
   if a==1, Rearth=Rearth/EMdist; Rmoon=Rmoon/EMdist;end %nondim units

   %This will draw disks of the Earth and Moon on a km plot
   fill( a*( -mu) + Rearth*cos(theta), Rearth*sin(theta),'c');
   fill( a*(1-mu) + Rmoon *cos(theta), Rmoon *sin(theta),'k');

elseif SEM==1, % Sun-Earth-Moon system
   if a==1, EMdist=EMdist/AU; Rearth=Rearth/AU; Rmoon=Rmoon/AU;end %nondim units

   fill( a*(1-mu) + Rearth*cos(theta), Rearth*sin(theta),'c');
   circle(EMdist+Rmoon,[a*(1-mu) 0],'k--');
   circle(EMdist-Rmoon,[a*(1-mu) 0],'k--');

elseif SEM==2, % Sun-Earth-Moon system (inertial frame centered on Earth)
   if a==1, EMdist=EMdist/AU; Rearth=Rearth/AU; Rmoon=Rmoon/AU;end %nondim units
 
   fill( Rearth*cos(theta), Rearth*sin(theta),'c');
   circle(EMdist+Rmoon,[0 0],'k--');
   circle(EMdist-Rmoon,[0 0],'k--');

elseif SEM==3, % Sun-Earth-Moon system (inertial frame centered on Earth)
	       % in terms of Earth-Moon distances
   if a==1, Rearth=Rearth/EMdist; Rmoon=Rmoon/EMdist; end %nondim units
   %This will draw disk of the Earth
   fill( Rearth*cos(theta), Rearth*sin(theta),'c');
   circle(1 + Rmoon,[0 0],'k--');
   circle(1 - Rmoon,[0 0],'k--');

end

end
