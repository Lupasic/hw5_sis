%%
% Extended Kalman filter.
% 
%% Car params
dist_btw_wheels = 97;

%% Read from files
f_gyro = fopen('gyro.txt','r');
f_motor_r = fopen('motor_left.txt','r');
f_motor_l = fopen('motor_right.txt','r');
f_sonar = fopen('sonar.txt','r');

gyro_data = fscanf(f_gyro,'%f');
motor_l_data = fscanf(f_motor_l,'%f');
motor_r_data = fscanf(f_motor_r,'%f');
sonar_data = fscanf(f_sonar,'%f');
fclose('all');
% Dist finder assumption

sonar_data = sonar_data - sonar_data(1);
%transform in mm
sonar_data = sonar_data * (-10);

%transform to rad
gyro_data = deg2rad(gyro_data);
%% Algorithm params
%Process
motor_l_var = 20;
motor_r_var = 20;
enc_theta_var = 20;
%Measurement
gyro_var = 1;
sonar_var = 10;

%% Initial data

X_EKF = zeros(3,(length(gyro_data)));
P_previous = eye(3);
%for plots
X_model = zeros(3,(length(gyro_data)));
thetas = zeros((length(gyro_data)),1);
coord_x = zeros((length(gyro_data)),1);
sonar_x = zeros((length(gyro_data)),1);

%calculate jacobians
syms dsr dsl pre_theta;
m_m = (dsr + dsl)/2;
th = (dsr - dsl) / dist_btw_wheels;
Fu = [diff(m_m * cos(pre_theta + th/2),dsr) diff(m_m * cos(pre_theta + th/2),dsl); diff(mean_motor * sin(pre_theta + th/2),dsr) diff(mean_motor * sin(pre_theta + th/2),dsl); diff(th,dsr) diff(th,dsr)];
    

%% Implementation
for i = 2:length(gyro_data)
    %Prediction
    delta_motor_l = motor_l_data(i) - motor_l_data(i-1);
    delta_motor_r = motor_r_data(i) - motor_r_data(i-1);
    mean_motor = (delta_motor_r + delta_motor_l)/2;
    theta_rad = (delta_motor_r - delta_motor_l) / dist_btw_wheels;
    
%     For plots
    thetas(i) = X_model(3,i-1) + theta_rad;
    coord_x(i) = X_model(1,i-1) + mean_motor * cos(X_model(3,i-1) + theta_rad/2);
    sonar_x(i) = sonar_data(i)*cos(gyro_data(i));

    
%     
    A = [1 0 0; 0 1 0; 0 0 1];
    Q = [motor_r_var * abs(delta_motor_r) 0 ;0 motor_l_var * abs(delta_motor_l)];
    X_prediction = A * X_EKF(:,i-1) + [mean_motor * cos(X_EKF(3,i-1) + theta_rad/2); mean_motor * sin(X_EKF(3,i-1) + theta_rad/2); theta_rad];    
    
    X_model(:,i) = A * X_model(:,i-1) + [mean_motor * cos(X_model(3,i-1) + theta_rad/2); mean_motor * sin(X_model(3,i-1) + theta_rad/2); theta_rad];

    
    dsr = delta_motor_r;
    dsl = delta_motor_l;
    pre_theta = X_EKF(3,i-1);
    Fu_cur = subs(Fu);
    Fu_cur = double(Fu_cur);
    %
    P_prediction = A * P_previous * A' + Fu_cur * Q * Fu_cur';
    %Update
    C = [1 0 0;0 0 1];
    R = [sonar_var^2 0 ;0 gyro_var ^2];
    K = P_prediction * C' / (C * P_prediction * C' + R);
    X_EKF(:,i) =  X_prediction + K *([sonar_data(i)*cos(gyro_data(i));gyro_data(i)] - C * X_prediction);
    P_cur = (eye(3) - K*C)*P_prediction;
    
    %Change vals
    P_previous = P_cur;
end

%% Plot
time = 1:(length(gyro_data));

figure
% Trajectory
subplot(2,2,1)
plot(X_EKF(1,:),X_EKF(2,:))
title('Trajectory')
xlabel('X, mm')
ylabel('Y, mm')
legend('kalman filter')

% Angle
subplot(2,2,2)
plot(X_EKF(3,:),time,gyro_data, time,thetas', time)
title('Angle')
xlabel('angle, rad')
legend('kalman filer','gyroscope','encoders')


% X coords
subplot(2,2,3)
plot(X_EKF(1,:),time,sonar_x, time, coord_x', time)
title('X coords')
xlabel('X, mm')
legend('kalman filer','sonar','encoders')