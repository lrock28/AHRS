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

% Noise chracteristics

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
magvar_deg = -3.26581;
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

% Static pressure file
static_pres_file = strcat(cur_dir,'\dataset\',name,'_static.csv');
static_pres_file_content = readmatrix(static_pres_file);


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
init_g_bias = [2.4; -1.3; 5.6] * 10^-4;

% x_hat_ahrs = [phi, theta, psi, g_bias_x, g_bias_y, g_bias_z]'
% y_hat = [phi, theta, psi, g_bias_x, g_bias_y, g_bias_z]'

% Pre allocate matrices
x_hat_ahrs = zeros(6,size(time_vec_s,1));
P_ahrs = zeros(6,6,size(time_vec_s,1));
es_state_ahrs = x_hat_ahrs;
eul_acc_mag = zeros(3,size(time_vec_s,1));
Sk_ahrs = zeros(6,6, size(time_vec_s,1));


% Process noise covariance
S_ahrs(1:3,1:3) = gyro_char.sigma_bias ^2;
S_ahrs(4:6,4:6) = gyro_char.sigma_mu;

% Init error state
es_state_ahrs(1:2,1) = 0.34906^2 * ones(1,2);
es_state_ahrs(3,1) = 3.14159^2;
es_state_ahrs(4:6,1) =  0.01745^2 * ones(1,3);


% Init state error covariance
P_ahrs(:,:,1)  =  diag (es_state_ahrs(:,1));
k = 1;
j = 1;
nav_initialized = 0;
init_index = 0;

% Yaw test
mag_yaw_deg = zeros(1, size(time_vec_s,1));

%% Run INS GNSS loop
for i = 1:1:size(time_vec_s,1)
    cur_imu = [gyro_radps(i,:), accel_mps2(i,:)];
    cur_base_gps = [lat_rad(i); lon_rad(i); alt_m(i); ned_vel_mps(i,:)'];
    if nav_initialized == 0
        if (imu_new_dat(i) == 1) && (gnss_new_dat(i) == 1) && ...
                (num_sat(i) > 7) && (fix(i) == 3)
            %AHRS Init
            [x_hat_ahrs(1,i), x_hat_ahrs(2,i)] = tilt(cur_imu);
            x_hat_ahrs(3,i) = tilt_comp(mag_ut(i,:)', x_hat_ahrs(1,i), x_hat_ahrs(2,i), magvar_rad);
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
            %AHRS Update
            [x_hat_ahrs(:,i), P_ahrs(:,:,i)] = ahrs_time_update (dt_s, x_hat_ahrs(:,i-1), P_ahrs(:,:,i-1), cur_imu, S_ahrs, gyro_char);
            % Mag/accel measurement update
            [x_hat_ahrs(:,i), eul_acc_mag(:,i), P_ahrs(:,:,i), es_state_ahrs(:,i)] = ahrs_measurement_update (x_hat_ahrs(:,i), home_lla(1), home_lla(3), mag_ut(i,:)', accel_mps2(i,:), P_ahrs(:,:,i), R_ahrs, wgs84, magvar_rad);
            mag_yaw_deg(i) = tilt_comp (mag_ut(i,:)', x_hat_ahrs(1), x_hat_ahrs(2), magvar_rad)*RAD2DEG;
        end
    end
end

%% Plotting
end_ind = size(gnss_file_content,1);

ahrs_ekf_eul = [x_hat_ahrs(1,init_index:end_ind)' * RAD2DEG,x_hat_ahrs(2,init_index:end_ind)' * RAD2DEG,x_hat_ahrs(3,init_index:end_ind)'* RAD2DEG-90];
eul_acc_mag = eul_acc_mag*RAD2DEG;

l = 1:1:(end_ind - init_index + 1);

figure(1)
plot (l, ahrs_ekf_eul(:,1),'DisplayName','phi AHRS')
yline(0)
hold on
plot (l, eul_acc_mag(1,:),'DisplayName','phi Accel')
legend

figure(2)
plot (l, ahrs_ekf_eul(:,2),'DisplayName','theta AHRS')
yline(0)
hold on
plot (l, eul_acc_mag(2,:),'DisplayName','theta Accel')
legend

figure(3)
plot (l, ahrs_ekf_eul(:,3),'DisplayName','psi AHRS')
yline(0)
hold on
plot (l, eul_acc_mag(3,:)-90,'DisplayName','psi Accel')
% plot (l, mag_yaw_deg(init_index:end_ind)-90,'DisplayName','psi mag')
legend

%% Filter function

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
function [x_ahrs, eul_acc_mag, P_correct, es_state] = ahrs_measurement_update (x_pred, lat, alt, mag_field, accel_mps2, P_pred, R, ellipsoid, magvar_rad)
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
    psi_mag = tilt_comp (mag_field, phi, theta, magvar_rad); %INS phi/theta
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

%% Helper functions
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

function psi_mag = tilt_comp (hb, phi, theta, magvar_rad)
    % Magnetic Field data
    % Needs Pulled for Lat/Lon & Date (W33.21550 degrees N87.54370 degrees)
    % https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm
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