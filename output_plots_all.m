clc; close all;

T_orb = 2.7499;   % nondimensional orbital period
A_z   = 0.2;      % (not used in scaling but for reference)

% CR3BP Earth-Moon scaling
L_star = 384400;      % km
T_star = 27.32*86400/(2*pi); % s
V_star = L_star/T_star;   % km/s

deltar = squeeze(out.deltar);
deltav = squeeze(out.deltav);
modu = squeeze(out.modu);

x = squeeze(out.x);
y = squeeze(out.y);
z = squeeze(out.z);
t = squeeze(out.t);   % your time vector


t_norm = out.t/T_orb;   % time in units of orbital periods
t_days = t*T_star/(24*3600);  % time in days
T_day = T_orb*T_star/(24*3600);

%% FIGURE 1: NON-DIMENSIONAL 
figure('Name','Non-Dimensional Results')
tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

% --- Position error ---
nexttile
plot(t_norm, deltar, 'Color', [0.86 0.08 0.24])
grid on
ylabel('||\deltar|| (ND)')
xlim([0 30])
title('Position Error')

% --- Velocity error ---
nexttile
plot(t_norm, deltav, 'k-', 'LineWidth',1.5)
grid on
ylabel('||\deltav|| (ND)')
xlim([0 30])
ylim([0 0.025])
title('Velocity Error')

% --- Control magnitude ---
nexttile
plot(t_norm, modu)
grid on
ylabel('||u|| (ND)')
xlabel('Time (T_{orb})')
xlim([0 30])
ylim([0 4e-4])
title('Control Effort')

nexttile
plot(t_norm, x, 'LineWidth',1.2)
hold on
plot(t_norm, y, 'LineWidth',1.2)
plot(t_norm, z, 'LineWidth',1.2)
grid on
ylabel('X(t) (ND)')
xlabel('Time (T_{orb})')
legend('x', 'y', 'z');
title('Output State')
xlim([0 30])
ylim([-1.5 1.5])

sgtitle('TV-LQR Performance (Non-Dimensional)','FontWeight','bold')

%% DIMENSIONAL CONVERSION

% Position in km
deltar_dim = deltar * L_star;

% Velocity in km/s
deltav_dim = deltav * V_star;

% Control in km/s^2 (since nondim accel)
u_dim = modu*(L_star/T_star^2);

% State in km
x_dim = x * L_star;
y_dim = y * L_star;
z_dim = z * L_star;

%% FIGURE 2: DIMENSIONAL
figure('Name','Dimensional Results')

tiledlayout(2,2,'TileSpacing','compact','Padding','compact')

% --- Position error ---
nexttile
plot(t_days, deltar_dim, 'Color', [0.86 0.08 0.24])
grid on
ylabel('||\deltar|| (km)')
xlim([0 359])
title('Position Error')

% --- Velocity error ---
nexttile
plot(t_days, deltav_dim, 'k-', 'LineWidth',1.5)
grid on
ylabel('||\deltav|| (km/s)')
xlim([0 359])
ylim([0 0.025])
title('Velocity Error')

% --- Control ---
nexttile
plot(t_days, u_dim)
grid on
ylabel('||u|| (km/s^2)')
xlabel('Time (T_{orb}) in days')
xlim([0 359])
title('Control Effort')

nexttile
plot(t_days, x_dim, 'LineWidth',1.2)
hold on
plot(t_days, y_dim, 'LineWidth',1.2)
plot(t_days, z_dim, 'LineWidth',1.2)
grid on
ylabel('X(t) (km)')
xlabel('Time (T_{orb}) in days')
xlim([0 359])
title('Output State')
legend('x', 'y', 'z');
xlim([0 359])
ylim([-3e5 4e5])

sgtitle('TV-LQR Performance (Dimensional)','FontWeight','bold')


%rms error
rms_err = sqrt(mean(deltar_dim.^2));

%max error
max_err = max(deltar_dim);

%delta V per year
deltaV_total_nd = trapz(t, modu)
num_orbits = (t(end) - t(1)) / T_orb;

deltaV_orbit_nd = deltaV_total_nd / num_orbits;
deltaV_orbit_ms = deltaV_orbit_nd * V_star * 1000;

deltaV_year_ms = deltaV_orbit_ms * (365 / 11.94);


% saturation fraction
u = squeeze(out.modu);  
umax = 5e-4;        
u_mag = vecnorm(u,2,2);
tol = 0.999;
is_sat = u_mag >= tol * umax;
sat_fraction = sum(is_sat) / length(u_mag);


fprintf('DeltaV/year : %.3f m/s\n', deltaV_year_ms);
fprintf('RMS error in position : %f km\n', rms_err);
fprintf('Maximum error in position : %f km\n', max_err);
fprintf('Saturation fraction (magnitude-based): %f\n', sat_fraction);