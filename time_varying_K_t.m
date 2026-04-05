clc; clear;

mu = 0.0121505856;

%% Load reference orbit
data = readmatrix('halo_L1_one_period_nondim.csv');
t_ref = data(:,1);
X_ref = data(:,2:7);

N = length(t_ref);
T = t_ref(end);

x0 = X_ref(1,:)';

x_ref = X_ref(:,1);
y_ref = X_ref(:,2);
z_ref = X_ref(:,3);
vx_ref = X_ref(:,4);
vy_ref = X_ref(:,5);
vz_ref = X_ref(:,6);

%% Weights
Q = diag([5e-1 5e-1 5e-1 1e-3 1e-3 1e-3]);
R = diag([1e-4 1e-4 1e-4]);
%Q = diag([5e-2 5e-2 5e-2 1e-3 1e-3 1e-3]);

%% System definition
A_fun = @(t) get_A_interp(t, t_ref, X_ref, mu);
B = [zeros(3,3); eye(3)];

%% ================= PERIODIC RICCATI =================
P_T = zeros(6);         % initial guess
max_iter = 50;
tol = 1e-6;

for iter = 1:max_iter

    % Force solution on SAME GRID (critical fix)
    t_span = flipud(t_ref);

    [t_sol, P_sol] = ode45(@(t,Pvec) riccati_ode(t, Pvec, A_fun, B, Q, R), ...
                           t_span, P_T(:));

    % Extract P(0)
    P0 = reshape(P_sol(end,:),6,6);
    P0 = (P0 + P0')/2;

    % Convergence check
    err = norm(P0 - P_T, 'fro') / max(1, norm(P_T,'fro'));
    fprintf('Iter %d: rel error = %e\n', iter, err);

    if err < tol
        disp('Converged to periodic Riccati solution');
        break;
    end

    % Update terminal condition
    P_T = P0;
end

%% Flip to forward time
t_sol = flipud(t_sol);
P_sol = flipud(P_sol);

%% ================= COMPUTE K(t) =================
K_all = zeros(3,6,N);

for k = 1:N
    P = reshape(P_sol(k,:),6,6);
    P = 0.5*(P + P');   % enforce symmetry
    
    K_all(:,:,k) = R \ (B' * P);
end

%% ================= ENFORCE PERIODIC SMOOTHNESS =================
K_avg = 0.5*(K_all(:,:,1) + K_all(:,:,end));
K_all(:,:,1)   = K_avg;
K_all(:,:,end) = K_avg;

%% Save for Simulink
t_vec = t_sol;
save('K_timevarying_periodic.mat', 'K_all', 't_vec');

disp('Periodic TV-LQR gains computed and saved')

%% ================= FUNCTIONS =================

function dPdt = riccati_ode(t, Pvec, A_fun, B, Q, R)

P = reshape(Pvec,6,6);
A = A_fun(t);

dP = -(A'*P + P*A - P*B*(R\B')*P + Q);

dPdt = dP(:);

end

function A = get_A_interp(t, t_ref, X_ref, mu)

% Periodic wrap
T = t_ref(end);
t = mod(t, T);

% Interpolate reference state
x = interp1(t_ref, X_ref, t, 'pchip')';

[A,~] = cr3bp_linearized(x, mu);

end

%% Linearized CR3BP
function [A,B] = cr3bp_linearized(x,mu)

X = x(1); Y = x(2); Z = x(3);

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

A = zeros(6,6);
A(1:3,4:6) = eye(3);
A(4:6,1:3) = Urr;
A(4:6,4:6) = [0  2  0
             -2  0  0
              0  0  0];

B = [zeros(3,3); eye(3)];

end