clc; clear;

mu = 0.0121505856;

% Load reference orbit
data = readmatrix('halo_L1_one_period_nondim.csv');
t_ref = data(:,1);
X = data(:,2:7);

% Linearization point (start of orbit)
x0 = X(1,:)';   

x_ref = X(:,1);
y_ref = X(:,2);
z_ref = X(:,3);
vx_ref = X(:,4);
vy_ref = X(:,5);
vz_ref = X(:,6);

% Get linearized matrices
[A,B] = cr3bp_linearized(x0,mu);

% LQR weights
Q = diag([2e-2 2e-2 2e-2 1e-3 1e-3 1e-3]);
R = diag([1e-4 1e-4 1e-4]);
%Q = diag([1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2]);
%R = diag([5e-6,5e-6,5e-6]);

% Compute gain
K = lqr(A,B,Q,R);

disp('LQR gain K =')
disp(K)

save K.mat K


%% Linearized CR3BP
function [A,B] = cr3bp_linearized(x,mu)

X = x(1); Y = x(2); Z = x(3);
Xd = x(4); Yd = x(5); Zd = x(6);

r1 = ((X+mu)^2 + Y^2 + Z^2)^(1/2);
r2 = ((X-1+mu)^2 + Y^2 + Z^2)^(1/2);

mu1 = 1-mu;
mu2 = mu;

Uxx = 1 - mu1*(1/r1^3 - 3*(X+mu)^2/r1^5) - mu2*(1/r2^3 - 3*(X-1+mu)^2/r2^5);

Uyy = 1 - mu1*(1/r1^3 - 3*Y^2/r1^5) - mu2*(1/r2^3 - 3*Y^2/r2^5);

Uzz = - mu1*(1/r1^3 - 3*Z^2/r1^5) - mu2*(1/r2^3 - 3*Z^2/r2^5);

Uxy = 3*mu1*(X+mu)*Y/r1^5 + 3*mu2*(X-1+mu)*Y/r2^5;

Uxz = 3*mu1*(X+mu)*Z/r1^5 + 3*mu2*(X-1+mu)*Z/r2^5;

Uyz = 3*mu1*Y*Z/r1^5 + 3*mu2*Y*Z/r2^5;


Urr = [ Uxx  Uxy  Uxz
        Uxy  Uyy  Uyz
        Uxz  Uyz  Uzz ];

% State matrix
A = zeros(6,6);
A(1:3,4:6) = eye(3);
A(4:6,1:3) = Urr;
A(4:6,4:6) = [0  2  0
             -2  0  0
              0  0  0];

% Input matrix
B = [zeros(3,3); eye(3)];

end
