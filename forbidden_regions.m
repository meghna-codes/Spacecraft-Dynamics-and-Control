clc; clear; close all;

%% Parameters (Earth–Moon)
m_earth = 5.97219e24;
m_moon  = 7.34767309e22;
mu = m_moon / (m_earth + m_moon);

%% Grid
Nx = 800; 
Ny = 800;
x = linspace(-2,2,Nx);
y = linspace(-2,2,Ny);
[X,Y] = meshgrid(x,y);

%% Distances
r1 = sqrt((X+mu).^2 + Y.^2);
r2 = sqrt((X-(1-mu)).^2 + Y.^2);

%% Pseudo-potential
Omega = (1-mu)./r1 + mu./r2 + 0.5*(X.^2 + Y.^2) + 0.5*mu*(1-mu);

Cvals = [3.24 3.19 3.16 3.01 2.80]; 
titles = { ...
'$ \mathrm{Case\ 1:}\ E < E_1$', ...
'$ \mathrm{Case\ 2:}\ E_1 < E < E_2$', ...
'$ \mathrm{Case\ 3:}\ E_2 < E < E_3$', ...
'$ \mathrm{Case\ 4:}\ E_3 < E < E_4$', ...
'$ \mathrm{Case\ 5:}\ E > E_4$'};

L = lagrange(mu);

figure('Color','w','Position',[100 100 800 600])

t = tiledlayout(2,3,'TileSpacing','none','Padding','none');

for k = 1:5
    
    C = Cvals(k);
    allowed = (2*Omega >= C);
   
    nexttile
    hold on;
    
    imagesc(x,y,allowed)
    set(gca,'YDir','normal')
    
    colormap([0.8 0.8 0.8;   % forbidden
              1 1 1])    % allowed
    clim([0 1])
    

    contour(X,Y,2*Omega,[C C],'k','LineWidth',1.8)
    
    plot(1-mu,0,'ko','MarkerFaceColor','k','MarkerSize',3)
    plot(-mu,0,'ko','MarkerFaceColor','k','MarkerSize',5)
    
    text(-mu+0.15,0,'m_1','FontSize',11,'FontWeight','bold')
    text(1-mu+0.3,0,'m_2','FontSize',11,'FontWeight','bold')
    
    title('')
    text(0,-2,titles{k}, 'HorizontalAlignment','center', ...
        'Interpreter','latex', 'FontSize',12)

    axis equal off
    xlim([-2 2]); ylim([-2 2])
end

sgtitle('\textbf{Allowed and Forbidden Regions for Different Energy Levels}', ...
    'Interpreter','latex','FontSize',15)

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
