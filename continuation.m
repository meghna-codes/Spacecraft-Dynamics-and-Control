%% Halo - orbit families around L1 and L2

clear; clc; close all;

global mu FORWARD a tol Lp
FORWARD = 1; 
a=1;
tol = 1e-12;

%% System parameters
m_earth = 5.97219e24;
m_moon  = 7.34767309e22;
mu = m_moon / (m_earth + m_moon);

%% Choose Lagrange Point (L1 or L2)
% 1 for L1, 2 for L2 
%Lp = 1; n=1; %northern L1
%Lp = 1; n=-1; %southern L1
Lp = 2; n=-1; %northern L2
%Lp = 2; n=1; %southern L2

am = 384000; %amplitude scaling in km
Rmoon = 1737.4; %km
EMdist = 385000; %km

% make a sphere with a nondimensional radius of the Moon
[Xm,Ym,Zm] = sphere;
Xm = Xm*Rmoon + (1-mu)*am;
Ym = Ym*Rmoon ;
Zm = Zm*Rmoon ;

%% PARAMETERS
Az0 = 0.001;     % starting amplitude (small halo)
dAz = 1e-4;       % larger continuation step → wider halo family
N   = 500;         % compute more orbits to reach larger amplitudes
plot_step = 40;    % plot only every 4th orbit to reduce density

orbits_data = cell(N,1);
Az = zeros(N,1);
Az(1) = Az0;

nu1 = zeros(N,1);
nu2 = zeros(N,1);

%% FIRST ORBIT (Richardson seed)
x0 = richardson(mu,Lp,Az0,n);
 
CASE = 1; % Keep z(0) constant, let x(0) and ydot(0) vary
%CASE = 2; %Keep x(0) constant, let z(0) and ydot(0) vary

[XH,TH] = diffcor(x0,CASE);
[x,t]   = integrate(XH,0,2*TH);

orbits_data{1} = x;
Xprev = XH;

%% CONTINUATION LOOP
for i = 2:N
    
    % Predictor: use previous converged orbit
    Xguess = Xprev;
    
    % grow out-of-plane amplitude
    Xguess(3) = Xguess(3) + n*dAz;
    
    Az(i) = Az(i-1) + n*dAz;

    % Corrector
    [XH,TH] = diffcor(Xguess,CASE);
    %[x,t]   = integrate(XH,0,2*TH);
    
    [x,t,phi_T,PHI] = getPHI(XH, 2*TH);

    [nu, lambda] = stability_indices(phi_T);

    nu1(i) = real(nu(1));
    nu2(i) = real(nu(2));

    orbits_data{i} = x;
    Xprev = XH;

    fprintf('Orbit %d converged\n',i)
end

figure;
plot(Az, nu1, 'LineWidth',1.5); hold on;
plot(Az, nu2, 'LineWidth',1.5);
xlabel('A_z'); ylabel('\nu');
legend('\nu_1','\nu_2');
grid on;
title('Stability Indices vs Halo Amplitude');

[xL, yL] = Lagrange(mu, Lp)

figure(1);
clf;

% Moon
hMoon = surf(Xm,Ym,Zm);
hold on;

cmap = gray(100);
cmap = cmap(1:80,:);
colormap(cmap);

Az_norm = (Az - min(Az)) / (max(Az)-min(Az));

% Orbits
firstHaloHandle = [];   % store first halo for legend

C_all = cell(N,1);

for i = 1:plot_step:N
    x = orbits_data{i};

    % --- Jacobi constant ---
    r1 = sqrt((x(:,1)+mu).^2 + x(:,2).^2 + x(:,3).^2);
    r2 = sqrt((x(:,1)-(1-mu)).^2 + x(:,2).^2 + x(:,3).^2);

    U  = 0.5*(x(:,1).^2 + x(:,2).^2) + (1-mu)./r1 + mu./r2;
    v2 = x(:,4).^2 + x(:,5).^2 + x(:,6).^2;

    C_vals = 2*U - v2;
    C_all{i} = C_vals;

    cidx = max(1, round(1 + Az_norm(i)*(size(cmap,1)-1)));
    col = cmap(cidx,:);

    h = plot3(x(:,1)*am, x(:,2)*am, x(:,3)*am, ...
              'Color', col, 'LineWidth', 1.2);

    if isempty(firstHaloHandle)
        firstHaloHandle = h;
    end
end

cb = colorbar;
cb.Label.String = 'A_z';
clim([min(Az_norm)*am max(Az_norm)*am])

%Lagrange Point
hL = plot3(xL*am, yL*am, 0, 'ro', ...
           'MarkerFaceColor','r','MarkerSize',10);

title('Southern halo family at L1')
xlabel('x (km)'); 
ylabel('y (km)'); 
zlabel('z (km)')

axis equal
grid on
set(gca,'FontSize',14)
view(40,20)

legend([firstHaloHandle, hL, hMoon],{'Halo family','Lagrange point','Moon'},...
       'Location','northeast')

xlim([2.95e5 4.25e5]);
ylim([-5e4 5e4]);
zlim([-5e4 5e4]);

lighting gouraud
camlight headlight

%% calculating the state transition matrix PHI(t,0)
function [x,t,phi_T,PHI]=getPHI(x0,tf)
 
% Gets state transition matrix, phi_T, and the trajectory (x,t) for a length 
% of time, tf (2*pi is 1 period). In particular, for periodic solutions of 
% period tf=T, one can obtain the monodromy matrix, PHI(0,T).
 
global tol
OPTIONS = odeset('RelTol',3*tol,'AbsTol',tol);
PHI_0(1:36) = reshape(eye(6),36,1);
PHI_0(1+36:6+36) = x0;
[t,PHI] = ode45(@(t,PHI) variational(t,PHI), [0 tf], PHI_0, OPTIONS); 
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

function [t1_z,x1_z]=find0(x0,Lp)

% From the given state x0, this finds the next x-crossing of the halo orbit
global t0_z x0_z x1_zgl
tolzero = 1.e-12;
options = optimset('DiffMaxChange',tolzero,'TolCon',tolzero,'TolFun',tolzero);

if Lp==1
    t0_z = pi/2 - 0.15;
else
    t0_z = pi/2 + 0.15;   % L2 : opposite shift
end

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
OPTIONS = odeset('RelTol',3*tol,'AbsTol',tol);

% 3D equations of motion integrated
[t,x] = ode45(@(t,x)  cr3bp(t,x), [t0 tf], x0, OPTIONS) ;
 
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
 
global mu tol Lp

 mu2 = 1-mu;
 t0 = 0;
 c=['r' 'b' 'g' 'm' ]; clf; 
 hold on;
 axis('equal');

 dx1 = 1;
 maxdx1=1.e-12;
 iter= 0;
 maxiter=15;

while abs(dx1) > maxdx1
 if iter > maxiter
        ERROR = 'OCCURRED'
 end
   N=rem(iter,3)+1;

   [t1,xx1]=find0(X0,Lp);
   
   iter=iter+1

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

 XH=X0
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
deltan = n;

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

%% Stability Indices from Monodromy Matrix
function [nu, lambda] = stability_indices(phi_T)

% Eigenvalues of monodromy matrix
lambda = eig(phi_T);

tol = 1e-6;
mask = abs(lambda - 1) < tol;

% Remove ONLY two of them
remove_idx = find(mask, 2);  
lambda(remove_idx) = [];

% Now you should have 4 eigenvalues left
% Pick one from each reciprocal pair

lambda1 = lambda(1);
lambda2 = lambda(3);  % skip its reciprocal

lambda_out = [lambda1; lambda2];

% Stability indices
nu = 0.5*(lambda_out + 1./lambda_out);

end