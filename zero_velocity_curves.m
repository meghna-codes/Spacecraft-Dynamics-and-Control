clc; clear; close all;

%% Parameters (Earth–Moon)
m_earth = 5.97219e24;
m_moon  = 7.34767309e22;
mu = m_moon / (m_earth + m_moon);

%% Grid
Nx = 900; 
Ny = 900;
x = linspace(-2,2,Nx);
y = linspace(-2,2,Ny);
[X,Y] = meshgrid(x,y);

%% Distances
r1 = sqrt((X+mu).^2 + Y.^2);
r2 = sqrt((X-(1-mu)).^2 + Y.^2);

%% Pseudo-potential
Omega = (1-mu)./r1 + mu./r2 + 0.5*(X.^2 + Y.^2) + 0.5*mu*(1-mu);

%% Jacobi constants
C_levels = [3.0 3.01 3.05 3.1 3.15 3.19 3.25 3.3 3.35 3.4 3.45 3.5];

%% Plotting
figure('Color','w'); hold on; box on;

[CS,h] = contour(X,Y,2*Omega,C_levels,'LineWidth',1.2);
colormap(gray)

cb = colorbar;
ylabel(cb,'\textbf{Jacobi Constant  C}','Interpreter','latex', 'FontSize',10)

clim([min(C_levels) max(C_levels)])

L = lagrange(mu);

plot(L(:,1),L(:,2),'ro','MarkerFaceColor','r','MarkerSize', 5)
plot(1-mu,0,'ko','MarkerFaceColor','k', 'MarkerSize', 8) % Moon
plot(-mu,0,'bo','MarkerFaceColor','b', 'MarkerSize', 10)  % Earth

names = {'L1','L2','L3','L4','L5'};

for i=1:5
    text(L(i,1)-0.03,L(i,2)+0.11,names{i},'FontSize',10, 'FontWeight','bold')
end

legend('Zero velocity curves','Lagrange points','Moon','Earth', ...
    'Interpreter','latex')

xlabel('x (non-dimensional units)','Interpreter','latex', 'FontSize',14)
ylabel('y (non-dimensional units)', 'Interpreter','latex', 'FontSize',14)
title('\textbf{Zero Velocity Curves for Earth-Moon system}', ...
      'Interpreter','latex', 'FontSize',16)

%xlim([-2 2]); ylim([-2 2])
grid on;

%% to calculate lagrange points
function L = lagrange(mu)

syms x

mu1 = mu;
mu2 = 1-mu;

% Collinear points
f = x - mu2*(x+mu1)/abs(x+mu1)^3 ...
      - mu1*(x-mu2)/abs(x-mu2)^3;

L1 = double(vpasolve(f==0,x,0.8));
L2 = double(vpasolve(f==0,x,1.2));
L3 = double(vpasolve(f==0,x,-1));

% Triangular points
L4 = [0.5-mu1 ,  sqrt(3)/2];
L5 = [0.5-mu1 , -sqrt(3)/2];

L = [L1 0;
     L2 0;
     L3 0;
     L4;
     L5];
end
