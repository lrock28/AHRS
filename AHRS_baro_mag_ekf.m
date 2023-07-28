%% Luke Rockwell
%  AHRS
%  
%% Set up workspace
clear;
close all;
clc;
% seed RNG
rng(4319);

%% Constants
RAD2DEG = 180/pi;
DEG2RAD = pi/180;
dt_s = 0.01;

% Noise chracteristic
% GNSS char
sigma_ne_pos = 0.1;
sigma_d_pos = 0.2;
sigma_vel_ne = 1;
sigma_vel_d = 1;

% Accel char
accel_char.sigma_bias = 0.01;
accel_char.tau = 100;
accel_char.sigma_w = 0.05;
accel_char.sigma_mu = 2 * accel_char.sigma_bias ^2 / accel_char.tau;

% Gyro char
gyro_char.sigma_bias = 0.00025;
gyro_char.tau = 50;
gyro_char.sigma_w = 0.00175;
gyro_char.sigma_mu = 2 * gyro_char.sigma_bias ^2 / gyro_char.tau;

% Magnetometer char
sigma_mag = .001;
R_mag = sigma_mag^2;

% Press char
sigma_press = 5;
R_baro = sigma_press^2;

% AHRS char
Rh = .05;
Ra_phi = .05;
Ra_theta = .05;
R_ahrs = diag([Ra_phi, Ra_theta, Rh]);

% WGS84 ellipsoid parameters
wgs84.a = 6378137;
wgs84.flattening = 1 / 298.257223563;
wgs84.e2 = wgs84.flattening * (2 - wgs84.flattening);
wgs84.e = sqrt (wgs84.e2);
wgs84.eq_rad = 6378137;
wgs84.earth_rot_rate = 7.292115*10^-5;
wgs84.g_const = 3.986004418*10^14;

% Magnetic Field data
% Needs Pulled for Lat/Lon & Date (W33.21550 degrees N87.54370 degrees)
% https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm
magvar_rad = -3.26581*DEG2RAD;
magvarx_nt = 22957.5;
magvary_nt = -1310;
magvarz_nt = 42595;

%% Parse data
name = 'malt153';
cur_dir = pwd;

% IMU file
% [sys_time_s, sys_frame_time_s, imu_new_data, g_x, g_y, g_z, a_x, a_y, a_z,...
% m_x, m_y, m_z]
imu_file = strcat(cur_dir,'\dataset\',name,'_imu.csv');
imu_file_content = readmatrix(imu_file);

% GNSS file
% [sys_time_s, gnss_new_data, gnss_lat_rad, gnss_lon_rad, gnss_alt_wgs84_m,...
% gnss_ned_vel_mps]
gnss_file = strcat(cur_dir,'\dataset\',name,'_gnss.csv');
gnss_file_content = readmatrix(gnss_file);

% Online EKF gile
online_file = strcat(cur_dir,'\dataset\',name,'_online.csv');
online_file_content = readmatrix(online_file);

% Static pressure file
static_pres_file = strcat(cur_dir,'\dataset\',name,'_static.csv');
static_pres_file_content = readmatrix(static_pres_file);

% % Magnetic Field file - Needs work
% mag_field_file = strcat(cur_dir,'\dataset\','igrfwmmData.csv');
% mag_field_file_content = readmatrix(mag_field_file);


% Parse sensor data
time_vec_s = imu_file_content(:,1);
imu_new_dat = imu_file_content(:,2);
gyro_radps = imu_file_content(:,3:5);
accel_mps2 = imu_file_content(:,6:8);
mag_ut = imu_file_content(:,9:11); %microTesla
gnss_new_dat = gnss_file_content(:,2);
lat_rad = gnss_file_content(:,5);
lon_rad = gnss_file_content(:,6);
alt_m = gnss_file_content(:,7);
ned_vel_mps = gnss_file_content(:,8:10);
num_sat = gnss_file_content(:,3);
fix = gnss_file_content(:,4);
baro_new_dat = static_pres_file_content (:,2);
static_temp_c = static_pres_file_content (:,4);
static_pres_pa = static_pres_file_content (:,3);


%% Set up filter
% initial conditions
init_a_bias = [0; 0; 0];
init_g_bias = [2.4; -1.3; 5.6] * 10^-4;

% x_hat = [lat, lon, alt, v_x, v_y, v_z, phi, theta, psi, a_bias_x,...
% a_bias_y, a_bias_z, g_bias_x, g_bias_y, g_bias_z]'
% y_hat = [lat, lon, alt, v_x, v_y, v_z, phi, theta, psi, a_bias_x,...
% a_bias_y, a_bias_z, g_bias_x, g_bias_y, g_bias_z, gnd_spd_mps,...
% gnd_track_rad]'

% Pre allocate matrices
x_hat = zeros(15,size(time_vec_s,1));
es_state = x_hat;
y_hat = zeros(17, size(time_vec_s,1));
P = zeros (15,15,size(time_vec_s,1));
x_hat_ahrs = zeros(6,size(time_vec_s,1));
P_ahrs = zeros(6,6,size(time_vec_s,1));
es_state_ahrs = x_hat_ahrs;
eul_acc_mag = zeros(3,size(time_vec_s,1));
Sk_ahrs = zeros(6,6, size(time_vec_s,1));


% Process noise covariance
S = zeros (12,12);
S(1:3,1:3) = accel_char.sigma_w ^2 * eye (3);
S(4:6,4:6) = gyro_char.sigma_w ^2 * eye (3);
S(7:9,7:9) = accel_char.sigma_mu * eye (3);
S(10:12,10:12) = gyro_char.sigma_mu * eye (3);
S_ahrs = zeros(6,6);
S_ahrs(1:3,1:3) = gyro_char.sigma_bias ^2;
S_ahrs(4:6,4:6) = gyro_char.sigma_mu;

% Observation noise variance
% R = diag([sigma_ne_pos, sigma_ne_pos, sigma_d_pos, sigma_ne_vel, ...
%     sigma_ne_vel, sigma_d_vel]) .^2;
% R = diag([sigma_ne_pos, sigma_ne_pos, sigma_d_pos]) .^2;
R_gnss = diag([sigma_ne_pos, sigma_ne_pos, sigma_d_pos, ...
     sigma_vel_ne, sigma_vel_ne, sigma_vel_d]).^2;
% R_relpos = diag([sigma_rel_ne_pos, sigma_rel_ne_pos, sigma_rel_d_pos]).^2;

% Init error state
es_state(1:3,1) = 10^2 * ones(1,3);
es_state(4:6,1) = 1^2 * ones(1,3);
es_state(7:8,1) = 0.34906^2 * ones(1,2);
es_state(9,1) = 3.14159^2;
es_state(10:12,1) =  0.981^2 * ones(1,3);
es_state(13:15,1) =  0.01745^2 * ones(1,3);
es_state_ahrs(1:2,1) = 0.34906^2 * ones(1,2);
es_state_ahrs(3,1) = 3.14159^2;
es_state_ahrs(4:6,1) =  0.01745^2 * ones(1,3);


% Init state error covariance
P(:,:,1)  =  diag (es_state(:,1));
P_ahrs(:,:,1)  =  diag (es_state_ahrs(:,1));
k = 1;
j = 1;
nav_initialized = 0;
init_index = 0;

% Alt test 
baro_alt_m = zeros(1, size(time_vec_s,1));
gps_alt_m = zeros(1, size(time_vec_s,1));
est_alt_m = zeros(1, size(time_vec_s,1));

% Yaw test
mag_yaw_deg = zeros(1, size(time_vec_s,1));

%% Run INS GNSS loop
for i = 1:1:size(time_vec_s,1)
    cur_imu = [gyro_radps(i,:), accel_mps2(i,:)];
    cur_base_gps = [lat_rad(i); lon_rad(i); alt_m(i); ned_vel_mps(i,:)'];
    if nav_initialized == 0
        if (imu_new_dat(i) == 1) && (gnss_new_dat(i) == 1) && ...
                (num_sat(i) > 7) && (fix(i) == 3)
            x_hat(1:6,i) = cur_base_gps;
            [x_hat(7,i), x_hat(8,i)] = tilt(cur_imu);
            % x_hat(9,i) = 270 * DEG2RAD;
            x_hat(9,i) = tilt_comp(mag_ut(i,:)', x_hat(7,i), x_hat(8,i));
            x_hat(10:12,i) = init_a_bias;
            x_hat(13:15,i) = init_g_bias;
            es_state(1:3,i) = 5^2 * ones(1,3);
            es_state(4:6,i) = 1^2 * ones(1,3);
            es_state(7:8,i) = 0.34906^2 * ones(1,2);
            es_state(9,i) = 3.14159^2;
            es_state(10:12,i) =  0.981^2 * ones(1,3);
            es_state(13:15,i) =  0.01745^2 * ones(1,3);
            %AHRS Init
            [x_hat_ahrs(1,i), x_hat_ahrs(2,i)] = tilt(cur_imu);
            % x_hat_ahrs(3,i) = 270 * DEG2RAD;
            x_hat_ahrs(3,i) = tilt_comp(mag_ut(i,:)', x_hat_ahrs(1,i), x_hat_ahrs(2,i));
            x_hat_ahrs(4:6,i) = init_g_bias;
            es_state_ahrs(1:2,i) = 0.34906^2 * ones(1,2);
            es_state_ahrs(3,i) = 3.14159^2;
            es_state_ahrs(4:6,i) =  0.01745^2 * ones(1,3);
            nav_initialized = 1;
            init_index = i;
            home_lla = cur_base_gps(1:3)';
            home_lla(1:2) = home_lla(1:2) * RAD2DEG;
            home_pres_pa = static_pres_pa(i);
        end
    elseif nav_initialized == 1
        if (imu_new_dat(i) == 1)
            [x_hat(:,i), P(:,:,i)] = time_update (dt_s, x_hat(:,i-1), P (:,:,i-1), cur_imu, S,wgs84, accel_char, gyro_char);
            baro_alt_m(i) = -baro_alt_est(static_pres_pa(i), home_pres_pa,...
                static_temp_c(i));
            %AHRS Update
            [x_hat_ahrs(:,i), P_ahrs(:,:,i)] = ahrs_time_update (dt_s, x_hat_ahrs(:,i-1), P_ahrs(:,:,i-1), cur_imu, S_ahrs, gyro_char);
            % Mag/accel measurement update
            [x_hat_ahrs(:,i), eul_acc_mag(:,i), P_ahrs(:,:,i), es_state_ahrs(:,i)] = ahrs_measurement_update (x_hat_ahrs(:,i), home_lla(1), home_lla(3), mag_ut(i,:)', accel_mps2(i,:), P_ahrs(:,:,i), R_ahrs, wgs84);
            gps_alt_m(i) = gps_alt_m(i-1);
            est_alt_m(i) = est_alt_m(i-1);
            mag_yaw_deg(i) = tilt_comp (mag_ut(i,:)', x_hat(7), x_hat(8))*RAD2DEG;
        end

        %Loose magnetometer/barometer filter
        % if (imu_new_dat(i) == 1)
        %     [x_hat(:,i), P(:,:,i), es_state(:,i)] = mag_measurement_update (x_hat(:,i), mag_ut(i,:)', P(:,:,i), R_mag, wgs84);
        % end
        % 
        % if (baro_new_dat(i) == 1)
        %     [x_hat(:,i), P(:,:,i), es_state(:,i)] = baro_measurement_update (x_hat(:,i), static_pres_pa(i), static_temp_c(i), static_pres_pa(1), P(:,:,i), R_baro, wgs84);
        % end

        if (gnss_new_dat(i) == 1)
            [x_hat(:,i), P(:,:,i), es_state(:,i)] = gps_measurement_update (x_hat(:,i), cur_base_gps, P(:,:,i), R_gnss, wgs84);
            % est_alt_m(i) = complementary_filt(baro_alt_m(i), cur_base_gps(3), 0.5);
            lla_point(:,k) = cur_base_gps;
            k = k + 1;
        end
    end
end

%% Plotting
end_ind = size(gnss_file_content,1);

test_ekf_lla = [x_hat(1,init_index:end_ind)' * RAD2DEG,x_hat(2,init_index:end_ind)' * RAD2DEG,x_hat(3,init_index:end_ind)'];
ned_pos_test_ekf = lla2ned(test_ekf_lla,home_lla,'ellipsoid');
gnss_lla = horzcat(lat_rad(init_index:end_ind) * RAD2DEG, lon_rad(init_index:end_ind) * RAD2DEG, alt_m (init_index:end_ind));
ned_gnss = lla2ned (gnss_lla, home_lla, 'ellipsoid');
online_lla = online_file_content(init_index:end_ind,2:4);
online_lla(:,1:2) = online_lla(:,1:2) * RAD2DEG;
ned_online = lla2ned (online_lla,home_lla,'ellipsoid');
ekf_eul = [x_hat(7,init_index:end_ind)' * RAD2DEG,x_hat(8,init_index:end_ind)' * RAD2DEG,x_hat(9,init_index:end_ind)'* RAD2DEG];
ahrs_ekf_eul = [x_hat_ahrs(1,init_index:end_ind)' * RAD2DEG,x_hat_ahrs(2,init_index:end_ind)' * RAD2DEG,x_hat_ahrs(3,init_index:end_ind)'* RAD2DEG-90];
eul_acc_mag = eul_acc_mag*RAD2DEG;

% figure(1)
% plot(ned_pos_test_ekf(:,2),ned_pos_test_ekf(:,1),'DisplayName','test ekf')
% hold on
% plot (ned_gnss(:,2),ned_gnss(:,1),'.','DisplayName','gps')
% %plot (ned_online(10:end,2), ned_online(10:end,1),'DisplayName','online')
% grid on
% grid minor
% legend



l = 1:1:(end_ind - init_index + 1);
% figure(2)
% % plot(l, mod(x_hat(9,init_index:end_ind) * RAD2DEG + 360, 360), 'DisplayName','EKF')
% % hold on
% % plot(l, mod(truth_heading_rad(init_index:end_ind) * RAD2DEG + 360,360),'DisplayName','truth')
% plot(l, unwrap(x_hat(9,init_index:end_ind)) * RAD2DEG, 'DisplayName','EKF')
% hold on
% plot(l, unwrap(truth_heading_rad(init_index:end_ind)) * RAD2DEG,'DisplayName','truth')
% plot (l, unwrap(online_file_content(init_index:end_ind,10)) * RAD2DEG, 'DisplayName','online ekf');
% grid on
% grid minor
% legend
% 
% figure(3)
% plot (l, ned_pos_test_ekf(:,3),'DisplayName','EKF')
% hold on
% plot(l, ned_rover(:,3),'DisplayName','GPS')
% plot (l(10:end), ned_online(10:end,3),'DisplayName','online EKF')
% legend
% 
% figure (4)
% plot (l, x_hat(6, init_index:end_ind),'DisplayName','EKF')
% hold on
% plot(l, gnss_file_content(init_index:end_ind,10),'DisplayName','GPS')
% plot (l, online_file_content(init_index:end_ind, 7) ,'DisplayName','online EKF')
% yline(0)
% legend
% 
% figure(5)
% plot (l, baro_alt_m(init_index:end_ind),'DisplayName','baro')
% hold on;
% plot (l, alt_m(init_index:end_ind),'DisplayName','gps')
% plot (l, est_alt_m(init_index:end_ind),'DisplayName','comp')
% plot (l, -ned_pos_test_ekf(:,3),'DisplayName','ekf')
% yline(0)
% legend
% 
% figure(6)
% plot (l, mag_yaw_deg(init_index:end_ind),'DisplayName','mag')
% hold on;
% plot (l, RAD2DEG*x_hat(9, init_index:end_ind),'DisplayName','ekf')
% yline(0)
% legend

figure(7)
plot (l, ekf_eul(:,1),'DisplayName','phi GNSS')
yline(0)
hold on
plot (l, ahrs_ekf_eul(:,1),'DisplayName','phi AHRS')
plot (l, eul_acc_mag(1,:),'DisplayName','phi Accel')
legend

figure(8)
plot (l, ekf_eul(:,2),'DisplayName','theta GNSS')
yline(0)
hold on
plot (l, ahrs_ekf_eul(:,2),'DisplayName','theta AHRS')
plot (l, eul_acc_mag(2,:),'DisplayName','theta Accel')
legend

figure(9)
plot (l, ekf_eul(:,3),'DisplayName','psi GNSS')
yline(0)
hold on
plot (l, ahrs_ekf_eul(:,3),'DisplayName','psi AHRS')
plot (l, eul_acc_mag(3,:)-90,'DisplayName','psi Accel')
% plot (l, mag_yaw_deg(init_index:end_ind)-90,'DisplayName','psi mag')
legend

%% Filter function
function [x_pred, P_pred] = time_update (dt, x, P, imu, S, ellipsoid, a_char, g_char)
    %x_hat = [lat lon alt v_x v_y v_z phi theta psi a_bias_x a_bias_y a_bias_z
    % g_bias_x g_bias_y g_bias_z]'
    % Rename state
    lat = x(1);
    lon = x(2);
    alt = x(3);
    v_x = x(4);
    v_y = x(5);
    v_z = x(6);
    phi = x(7);
    theta = x(8);
    psi = x(9);
    a_bias = x(10:12);
    g_bias = x(13:15);
    quat = bfs_euler2quaternion([psi, theta, phi]);

    % Rename mesurements and apply bias
    ins_accel_mps2 = imu(4:6)' - a_bias;
    ins_gyro_radps = imu(1:3)' - g_bias;
    
    % Delta quaternion
    delta_quat = zeros (1,4);
    delta_quat(1) = 1;
    delta_quat(2) = 0.5 * ins_gyro_radps(1) * dt;
    delta_quat(3) = 0.5 * ins_gyro_radps(2) * dt;
    delta_quat(4) = 0.5 * ins_gyro_radps(3) * dt;

    quat = quatmultiply(quat, delta_quat);
    quat = quatnormalize(quat);

    % Avoid quartenion sign flip
    if quat(1) < 0
        quat = -1 * quat;
    end
    
    %% Update attitude
    ins_euler_angle = bfs_quaternion2euler(quat);
    Cbn = bfs_quaternion2dcm(quat); % Nav to body dcm
    Cnb = Cbn';
    euler_angle_pred = [ins_euler_angle(3), ins_euler_angle(2), ins_euler_angle(1)];

    %% Update velocity
    % Calculate gravity vector
    gn = calc_gravity (lat, alt, ellipsoid); 
    
    % v_dot for low grade IMU
    v_dot = Cnb * ins_accel_mps2 + gn;
    
    % update velocity
    v_ned_pred = [v_x; v_y; v_z] + dt * v_dot;
    
    %% Update position
    %p_dot = diag ([1 / (Rn + alt), 1 / ((Re + alt) * cos(lat)), -1]) * v_ned_pred;
    p_dot = bfs_llarate(v_ned_pred, [lat, lon, alt], ellipsoid);
    p_pred = [lat; lon; alt] + dt * p_dot;
    
    %% Predict states
    x_pred = x;
    x_pred (1:3) = p_pred;
    x_pred (4:6) = v_ned_pred;
    x_pred (7:9) = euler_angle_pred;
    
    %% Propagate state error covariance
    A = zeros (15, 15);

    A (1:3, 4:6) = eye (3,3);
    
    g_mag = norm(gn);
    A(4:6, 1:3) = (g_mag / ellipsoid.a) * diag ([0, 0, -2]);
    A(4:6, 7:9) = -2 * Cnb * skew(ins_accel_mps2);
    A(4:6, 10:12) = -Cnb;

    A(7:9, 7:9) = - skew(ins_gyro_radps);
    A(7:9, 13:15) = - 0.5 * eye(3,3);
    
    A(10:12, 10:12) = - (1 / a_char.tau) * eye(3, 3);
    
    A(13:15, 13:15) = - (1 / g_char.tau) * eye(3, 3);
    
    L = zeros (15,12);
    L(4:6,1:3) = -Cnb;
    L(7:9,4:6) = -0.5 * eye(3,3);
    L(10:12, 7:9) = eye(3,3);
    L(13:15, 10:12) = eye(3,3);
    
    %% Discrete transformation
    PHI = eye(15,15) + A * dt;
    Qk = PHI * dt * L * S * L';
    Qk = 0.5 * (Qk + Qk');
    
    % Covariance time update
    P_pred = PHI * P * PHI' + Qk;
    P_pred = 0.5 * (P_pred + P_pred');
end

%AHRS Time update
function [x_ahrs_pred, P_pred] = ahrs_time_update (dt, x, P, imu, S, g_char)
    %x_ahrs = [phi theta psi g_bias_x g_bias_y g_bias_z]'
    % Rename state
    phi = x(1);
    theta = x(2);
    psi = x(3);
    g_bias = x(4:6);
    quat = bfs_euler2quaternion([psi, theta, phi]);

    % Rename mesurements and apply bias
    ins_gyro_radps = imu(1:3)' - g_bias;
    
    % Delta quaternion
    delta_quat = zeros (1,4);
    delta_quat(1) = 1;
    delta_quat(2) = 0.5 * ins_gyro_radps(1) * dt;
    delta_quat(3) = 0.5 * ins_gyro_radps(2) * dt;
    delta_quat(4) = 0.5 * ins_gyro_radps(3) * dt;

    quat = quatmultiply(quat, delta_quat);
    quat = quatnormalize(quat);

    % Avoid quartenion sign flip
    if quat(1) < 0
        quat = -1 * quat;
    end
    
    %% Update attitude
    ins_euler_angle = bfs_quaternion2euler(quat);
    Cbn = bfs_quaternion2dcm(quat); % Nav to body dcm
    Cnb = Cbn';
    euler_angle_pred = [ins_euler_angle(3), ins_euler_angle(2), ins_euler_angle(1)];

    x_ahrs_pred = x;
    x_ahrs_pred (1:3) = euler_angle_pred;
    
    %% Propagate state error covariance
    A = zeros (6,6);

    A (1:3, 4:6) = -Cnb;
    
    A (4:6, 4:6) = - (1 / g_char.tau) * eye(3, 3);
    
    L = zeros (6,6);
    L(1:3,1:3) = -Cnb;
    L(4:6, 4:6) = eye(3,3);
    
    %% Discrete transformation
    PHI = eye(6,6) + A * dt;
    Qk = PHI * dt * L * S * L';
    Qk = 0.5 * (Qk + Qk');
    
    % Covariance time update
    P_pred = PHI * P * PHI' + Qk;
    P_pred = 0.5 * (P_pred + P_pred');
end

%AHRS Measurement update 
function [x_ahrs, eul_acc_mag, P_correct, es_state] = ahrs_measurement_update (x_pred, lat, alt, mag_field, accel_mps2, P_pred, R, ellipsoid)
    %x_hat = [lat lon alt v_x v_y v_z phi theta psi a_bias_x a_bias_y a_bias_z
    % g_bias_x g_bias_y g_bias_z]'
    % Rename state
    phi = x_pred(1);
    theta = x_pred(2);
    psi = x_pred(3);
    g_bias = x_pred(4:6);
    quat = bfs_euler2quaternion([psi, theta, phi]);
    gn = calc_gravity (lat, alt, ellipsoid);
    g_mag = norm(gn);
    theta_acc = -asin(accel_mps2(1)/g_mag);
    phi_acc = -atan(accel_mps2(2)/accel_mps2(3));
    psi_mag = tilt_comp (mag_field, phi, theta); %INS phi/theta
    eul_acc_mag = [phi_acc, theta_acc, psi_mag];
    %Cbn = bfs_quaternion2dcm(quat); % Nav to body dcm
    %Cnb = Cbn';

    %moving_base_vector_nav = Cnb * moving_base_vector;
    eul_res = [phi_acc, theta_acc, psi_mag] - [phi, theta, psi];
    %rel_pos_res = rel_pos' - moving_base_vector_nav;
    state_res = eul_res';

    H = zeros (3,6);
    H (1:3, 1:3) = eye(3,3);

    % Innovation
    Sk = H * P_pred * H' + R;

    % Kalman gain
    K = P_pred * H' * inv(Sk);

    % Correct covariance matrix
    P_correct = (eye(6) - K * H) * P_pred * (eye(6) - K * H)' + K * R * K';

    % Eroor state update
    es_state = K * state_res;

    %% Nominal state correction
    % Attitude correction
    delta_quat = zeros(1,4);
    delta_quat(1) = 1;
    delta_quat(2) = es_state(1);
    delta_quat(3) = es_state(2);
    delta_quat(4) = es_state(3);

    quat = quatmultiply(quat, delta_quat);
    quat = quatnormalize(quat);
    euler_from_quat = bfs_quaternion2euler(quat);
    phi_correct = euler_from_quat(3);
    theta_correct = euler_from_quat(2);
    psi_correct = euler_from_quat(1);

    gyro_bias_correct = g_bias + es_state(4:6);

    x_ahrs = [phi_correct; theta_correct; psi_correct; gyro_bias_correct(1); gyro_bias_correct(2); gyro_bias_correct(3)];

end

% Baro update 
function [x_correct, P_correct, es_state] = baro_measurement_update (x_pred, static_pres_pa, temp_c, home_pres_pa, P_pred, R, ellipsoid)
    %x_hat = [lat lon alt v_x v_y v_z phi theta psi a_bias_x a_bias_y a_bias_z
    % g_bias_x g_bias_y g_bias_z]'
    % Rename state
    lat = x_pred(1);
    lon = x_pred(2);
    alt = x_pred(3);
    v_x = x_pred(4);
    v_y = x_pred(5);
    v_z = x_pred(6);
    phi = x_pred(7);
    theta = x_pred(8);
    psi = x_pred(9);
    a_bias = x_pred(10:12);
    g_bias = x_pred(13:15);
    quat = bfs_euler2quaternion([psi, theta, phi]);
    baro_alt = baro_alt_est (static_pres_pa, home_pres_pa, temp_c);
    %Cbn = bfs_quaternion2dcm(quat); % Nav to body dcm
    %Cnb = Cbn';

    %moving_base_vector_nav = Cnb * moving_base_vector;
    alt_res = alt - baro_alt;
    %rel_pos_res = rel_pos' - moving_base_vector_nav;
    state_res = alt_res';

    H = zeros (1,15);
    H (1, 3) = 1;

    % Innovation
    Sk = H * P_pred * H' + R;
    
    % Kalman gain
    K = P_pred * H' * inv(Sk);
    
    % Correct covariance matrix
    P_correct = (eye(15) - K * H) * P_pred * (eye(15) - K * H)' + K * R * K';
    
    % Eroor state update
    es_state = K * state_res;
    
    %% Nominal state correction
    % Position correction
    [Rn, Re] = curvature (lat, ellipsoid);
    alt_correct = alt - es_state(3);
    lat_correct = lat + es_state(1) /(Rn + alt_correct);
    lon_correct = lon + es_state(2) / ((Re + alt_correct) * cos (lat));
    %alt_correct = gps(3);

    % % Velocity correction
    v_x_correct = v_x + es_state(4);
    v_y_correct = v_y + es_state(5);
    v_z_correct = v_z + es_state(6);
    
    % % Attitude correction
    delta_quat = zeros(1,4);
    delta_quat(1) = 1;
    delta_quat(2) = es_state(7);
    delta_quat(3) = es_state(8);
    delta_quat(4) = es_state(9);

    quat = quatmultiply(quat, delta_quat);
    quat = quatnormalize(quat);
    euler_from_quat = bfs_quaternion2euler(quat);
    phi_correct = euler_from_quat(3);
    theta_correct = euler_from_quat(2);
    psi_correct = euler_from_quat(1);
    
    % % Bias correction
    accel_bias_correct = a_bias + es_state(10:12);
    gyro_bias_correct = g_bias + es_state(13:15);

    x_correct = [lat_correct; lon_correct; alt_correct; v_x_correct; v_y_correct;...
        v_z_correct; phi_correct; theta_correct; psi_correct; accel_bias_correct; ...
        gyro_bias_correct];
    
end

%Magnetometer Measurement update 
function [x_correct, P_correct, es_state] = mag_measurement_update (x_pred, mag_field, P_pred, R, ellipsoid)
    %x_hat = [lat lon alt v_x v_y v_z phi theta psi a_bias_x a_bias_y a_bias_z
    % g_bias_x g_bias_y g_bias_z]'
    % Rename state
    lat = x_pred(1);
    lon = x_pred(2);
    alt = x_pred(3);
    v_x = x_pred(4);
    v_y = x_pred(5);
    v_z = x_pred(6);
    phi = x_pred(7);
    theta = x_pred(8);
    psi = x_pred(9);
    a_bias = x_pred(10:12);
    g_bias = x_pred(13:15);
    quat = bfs_euler2quaternion([psi, theta, phi]);
    psi_mag = tilt_comp (mag_field, phi, theta);
    %Cbn = bfs_quaternion2dcm(quat); % Nav to body dcm
    %Cnb = Cbn';

    %moving_base_vector_nav = Cnb * moving_base_vector;
    eul_res = psi_mag - psi;
    %rel_pos_res = rel_pos' - moving_base_vector_nav;
    state_res = eul_res';

    H = zeros (1,15);
    H (1, 9) = 1;

    % Innovation
    Sk = H * P_pred * H' + R;

    % Kalman gain
    K = P_pred * H' * inv(Sk);

    % Correct covariance matrix
    P_correct = (eye(15) - K * H) * P_pred * (eye(15) - K * H)' + K * R * K';

    % Eroor state update
    es_state = K * state_res;

    %% Nominal state correction
    % Position correction
    [Rn, Re] = curvature (lat, ellipsoid);
    lat_correct = lat + es_state(1) /(Rn + alt);
    lon_correct = lon + es_state(2) / ((Re + alt) * cos (lat));
    %alt_correct = gps(3);
    alt_correct = alt - es_state(3);

    % Velocity correction
    v_x_correct = v_x + es_state(4);
    v_y_correct = v_y + es_state(5);
    v_z_correct = v_z + es_state(6);

    % Attitude correction
    delta_quat = zeros(1,4);
    delta_quat(1) = 1;
    delta_quat(2) = es_state(7);
    delta_quat(3) = es_state(8);
    delta_quat(4) = es_state(9);

    quat = quatmultiply(quat, delta_quat);
    quat = quatnormalize(quat);
    euler_from_quat = bfs_quaternion2euler(quat);
    phi_correct = euler_from_quat(3);
    theta_correct = euler_from_quat(2);
    psi_correct = euler_from_quat(1);
   
    % Bias correction
    accel_bias_correct = a_bias + es_state(10:12);
    gyro_bias_correct = g_bias + es_state(13:15);

    x_correct = [lat_correct; lon_correct; alt_correct; v_x_correct; v_y_correct;...
        v_z_correct; phi_correct; theta_correct; psi_correct; accel_bias_correct; ...
        gyro_bias_correct];

end

% Measurement update need checking
function [x_correct, P_correct, es_state] = gps_measurement_update (x_pred, gps, P_pred, R, ellipsoid)
    %x_hat = [lat lon alt v_x v_y v_z phi theta psi a_bias_x a_bias_y a_bias_z
    % g_bias_x g_bias_y g_bias_z]'
    % Rename state
    lat = x_pred(1);
    lon = x_pred(2);
    alt = x_pred(3);
    v_x = x_pred(4);
    v_y = x_pred(5);
    v_z = x_pred(6);
    phi = x_pred(7);
    theta = x_pred(8);
    psi = x_pred(9);
    a_bias = x_pred(10:12);
    g_bias = x_pred(13:15);
    quat = bfs_euler2quaternion([psi, theta, phi]);
    %Cbn = bfs_quaternion2dcm(quat); % Nav to body dcm
    %Cnb = Cbn';

    %moving_base_vector_nav = Cnb * moving_base_vector;
    pos_res = lla2ned ([rad2deg(gps(1)), rad2deg(gps(2)), gps(3)],...
        [rad2deg(lat), rad2deg(lon), alt], "ellipsoid");
    vel_res = [gps(4), gps(5), gps(6)] - [v_x, v_y, v_z];
    %rel_pos_res = rel_pos' - moving_base_vector_nav;
    state_res = [pos_res'; vel_res'];

    H = zeros (6,15);
    H (1:6, 1:6) = eye(6);

    % Innovation
    Sk = H * P_pred * H' + R;
    
    % Kalman gain
    K = P_pred * H' * inv(Sk);
    
    % Correct covariance matrix
    P_correct = (eye(15) - K * H) * P_pred * (eye(15) - K * H)' + K * R * K';
    
    % Eroor state update
    es_state = K * state_res;
    
    %% Nominal state correction
    % Position correction
    [Rn, Re] = curvature (lat, ellipsoid);
    lat_correct = lat + es_state(1) /(Rn + alt);
    lon_correct = lon + es_state(2) / ((Re + alt) * cos (lat));
    %alt_correct = gps(3);
    alt_correct = alt - es_state(3);

    % Velocity correction
    v_x_correct = v_x + es_state(4);
    v_y_correct = v_y + es_state(5);
    v_z_correct = v_z + es_state(6);
    
    % Attitude correction
    delta_quat = zeros(1,4);
    delta_quat(1) = 1;
    delta_quat(2) = es_state(7);
    delta_quat(3) = es_state(8);
    delta_quat(4) = es_state(9);

    quat = quatmultiply(quat, delta_quat);
    quat = quatnormalize(quat);
    euler_from_quat = bfs_quaternion2euler(quat);
    phi_correct = euler_from_quat(3);
    theta_correct = euler_from_quat(2);
    psi_correct = euler_from_quat(1);
    
    % Bias correction
    accel_bias_correct = a_bias + es_state(10:12);
    gyro_bias_correct = g_bias + es_state(13:15);

    x_correct = [lat_correct; lon_correct; alt_correct; v_x_correct; v_y_correct;...
        v_z_correct; phi_correct; theta_correct; psi_correct; accel_bias_correct; ...
        gyro_bias_correct];
    
end

function [x_correct, P_correct, es_state] = relpos_measurement_update (x_pred, rel_pos, P_pred, R, ellipsoid, moving_base_vector)
    %x_hat = [lat lon alt v_x v_y v_z phi theta psi a_bias_x a_bias_y a_bias_z
    % g_bias_x g_bias_y g_bias_z]'
    % Rename state
    lat = x_pred(1);
    lon = x_pred(2);
    alt = x_pred(3);
    v_x = x_pred(4);
    v_y = x_pred(5);
    v_z = x_pred(6);
    phi = x_pred(7);
    theta = x_pred(8);
    psi = x_pred(9);
    a_bias = x_pred(10:12);
    g_bias = x_pred(13:15);

    quat = bfs_euler2quaternion([psi, theta, phi]);
    Cbn = bfs_quaternion2dcm(quat); % Nav to body dcm
    Cnb = Cbn';

    moving_base_vector_nav = Cnb * moving_base_vector;
    rel_pos_res = rel_pos' - moving_base_vector_nav;
    state_res = rel_pos_res;

    H = zeros (3,15);
    moving_base_rev = - moving_base_vector;
    H (1:3, 7:9) = Cnb * skew(moving_base_rev);
    
    % Innovation
    Sk = H * P_pred * H' + R;
    
    % Kalman gain
    K = P_pred * H' * inv(Sk);
    
    % Correct covariance matrix
    P_correct = (eye(15) - K * H) * P_pred * (eye(15) - K * H)' + K * R * K';
    
    % Eroor state update
    es_state = K * state_res;
    
    %% Nominal state correction
    % Position correction
    [Rn, Re] = curvature (lat, ellipsoid);
    lat_correct = lat + es_state(1) /(Rn + alt);
    lon_correct = lon + es_state(2) / ((Re + alt) * cos (lat));
    %alt_correct = gps(3);
    alt_correct = alt;

    % Velocity correction
    v_x_correct = v_x + es_state(4);
    v_y_correct = v_y + es_state(5);
    v_z_correct = v_z;
    %v_z_correct = gps(6);
    
    % Attitude correction
    delta_quat = zeros(1,4);
    delta_quat(1) = 1;
    delta_quat(2) = es_state(7);
    delta_quat(3) = es_state(8);
    delta_quat(4) = es_state(9);

    quat = quatmultiply(quat, delta_quat);
    quat = quatnormalize(quat);
    euler_from_quat = quat2eul(quat);
    phi_correct = euler_from_quat(3);
    theta_correct = euler_from_quat(2);
    psi_correct = euler_from_quat(1);
    %psi_correct = gps_compass(rel_pos, moving_base_vector);
    
    % Bias correction
    accel_bias_correct = a_bias + es_state(10:12);
    gyro_bias_correct = g_bias + es_state(13:15);

    x_correct = [lat_correct; lon_correct; alt_correct; v_x_correct; v_y_correct;...
        v_z_correct; phi_correct; theta_correct; psi_correct; accel_bias_correct; ...
        gyro_bias_correct];
    
end


%% Helper functions
% Radius of curvature using WGS84
function [Rn, Re] = curvature (lat, ellipsoid)
    Rn = ellipsoid.a * (1 - ellipsoid.e^2) / (1 - ellipsoid.e^2 ...
        * sin (lat)^2)^1.5;
    Re = ellipsoid.a / sqrt(1 - ellipsoid.e^2 * sin(lat)^2);
end

% Skew matrix
function S = skew (v)
    S = [   0, -v(3),  v(2);
         v(3),     0, -v(1);
        -v(2),  v(1),     0];
end

% Calculate local gravity
function g_local = calc_gravity (lat, alt, ellipsoid)
    frac = alt/ ellipsoid.a;
    g_0 = (9.7803253359 * (1 + 0.0019311853 * (sin (lat))^2)) / ...
        sqrt (1 - ellipsoid.flattening * (2 - ellipsoid.flattening) * (sin(lat))^2);
    ch = 1 - 2 * (1 + ellipsoid.flattening + (ellipsoid.a^3 * (1 - ellipsoid.flattening) *...
        ellipsoid.earth_rot_rate^2)/ellipsoid.g_const) * frac + 3 * frac^2;
    g_local = [0; 0; 9.80665];
end


% Convert euler = [psi, theta, phi] to quaternion
function quaternion = bfs_euler2quaternion (euler)
    psi = euler (1) / 2;
    the = euler (2) / 2;
    phi = euler (3) / 2;

    quaternion = zeros (1,4);
    quaternion(1) = cos(psi) * cos(the) * cos(phi) + sin(psi) * sin(the)...
        * sin(phi);
    quaternion(2) = cos(psi) * cos(the) * sin(phi) - sin(psi) * sin(the)...
        * cos(phi);
    quaternion(3) = cos(psi) * sin(the) * cos(phi) + sin(psi) * cos(the)...
        * sin(phi);
    quaternion(4) = sin(psi) * cos(the) * cos(phi) - cos(psi) * sin(the)...
        * sin(phi);
end

% Convert quaternion to euler = [psi, theta, phi] 
function euler = bfs_quaternion2euler(quaternion)
    m11 = 2 * quaternion(1)^2 + 2 * quaternion(2)^2 - 1;
    m12 = 2 * quaternion(2) * quaternion(3) + 2 * quaternion(1) * quaternion(4);
    m13 = 2 * quaternion(2) * quaternion(4) - 2 * quaternion(1) * quaternion(3);
    m23 = 2 * quaternion(3) * quaternion(4) + 2 * quaternion(1) * quaternion(2);
    m33 = 2 * quaternion(1)^2 + 2 * quaternion(4)^2 - 1;

    psi = atan2 (m12,m11);
    the = asin(-m13);
    phi = atan2 (m23,m33);

    euler = zeros (3,1);
    euler(1) = psi;
    euler(2) = the;
    euler(3) = phi;
end

function dcm = bfs_quaternion2dcm(quaternion)
    dcm = zeros(3,3);
    dcm(1,1) = 2 * quaternion(1)^2 - 1 + 2 * quaternion(2)^2;
    dcm(2,1) = 2 * quaternion(2) * quaternion(3) - 2 * quaternion(1) * ...
        quaternion(4);
    dcm(3,1) = 2 * quaternion(2) * quaternion(4) + 2 * quaternion(1) * ...
        quaternion(3);
    dcm(1,2) = 2 * quaternion(2) * quaternion(3) + 2 * quaternion(1) * ...
        quaternion(4);
    dcm(2,2) = 2 * quaternion(1)^2 - 1 + 2 * quaternion(3)^2;
    dcm(3,2) = 2 * quaternion(3) * quaternion(4) - 2 * quaternion(1) * ...
        quaternion(2);
    dcm(1,3) = 2 * quaternion(2) * quaternion(4) - 2 * quaternion(1) * ...
        quaternion(3);
    dcm(2,3) = 2 * quaternion(3) * quaternion(4) + 2 * quaternion(1) * ...
        quaternion(2);
    dcm(3,3) = 2 * quaternion(1)^2 - 1 + 2 * quaternion(4)^2;
end

function [roll, pitch] = tilt(imu)
    a = imu(4:6);
    a_norm = norm(a);
    a = a/a_norm;
    pitch = asin(a(1));
    roll = asin(-a(2)) / cos(pitch);
end

function heading = gps_compass(rel_pos, baseline_body)
baseline_body = baseline_body' .* ones(size(rel_pos));
dot_prod = dot(rel_pos, baseline_body,2);
cross_prod = cross(baseline_body, rel_pos,2);
heading = atan2(cross_prod(:,3),dot_prod);
end

function pdot = bfs_llarate (ned_vel, lla, ellipsoid)
[Rn, Re] = curvature (lla(1), ellipsoid);
vlat = ned_vel(1) / (Rn + lla(3));
vlon = ned_vel(2) / (Re + lla(3)) / cos(lla(1));
valt = - ned_vel(3);
pdot = [vlat;vlon;valt];
end

function h_baro = baro_alt_est (static_pres_pa, ref_pres_pa, temp_c)
    R = 287.058;
    temp_K = temp_c + 273.15;
    g = 9.81;
    h_baro = (R * temp_K / g) * log(static_pres_pa/ref_pres_pa);
end

function comp_out = complementary_filt (input1, input2, alpha)
    comp_out = (1-alpha)*input1 + alpha * input2;
end

function psi_mag = tilt_comp (hb, phi, theta)
    % Magnetic Field data
    % Needs Pulled for Lat/Lon & Date (W33.21550 degrees N87.54370 degrees)
    % https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm
    DEG2RAD = pi/180;
    magvar_rad = -3.26581*DEG2RAD;
    htc1 = ones(3,3);
    htc1(1,:) = [cos(theta), 0, sin(theta)];
    htc1(2,:) = [0, 1, 0];
    htc1(3,:) = [-sin(theta), 0, cos(theta)];
    htc2(1,:) = [1, 0, 0];
    htc2(2,:) = [0, cos(phi), -sin(phi)];
    htc2(3,:) = [0, -sin(phi), cos(phi)];
    htc = htc1 * htc2 * hb;
    psi_mag = -atan(htc(2)/htc(1)) + magvar_rad;
end