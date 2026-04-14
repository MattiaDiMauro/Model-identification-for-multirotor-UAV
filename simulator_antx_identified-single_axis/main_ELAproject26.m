%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Estimation and Learning" project 2025-2026           

% Professor: Lovera Marco
% Authors:  Boschini Marianna, Cavalli Filippo, Mattia Di Mauro                     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clearvars;
close all;
clear;
addpath('datasets','common','common/simulator-toolbox','common/simulator-toolbox/attitude_library','common/simulator-toolbox/trajectory_library');
clc;

%% Model parameters
% Initial model (state: longitudinal velocity, pitch rate, pitch angle; input: normalised pitching moment; outputs: state and longitudinal acceleration)
Xu=-0.1068;
Xq=0.1192;
Mu=-5.9755;
Mq=-2.6478;
Xd=-10.1647;
Md=450.71;
A=[Xu, Xq, -9.81; Mu, Mq, 0; 0, 1, 0];
B=[Xd; Md; 0];
C=[1, 0, 0; 0, 1, 0; 0, 0, 1; Xu, Xq, 0]; 
D=[0; 0; 0; Xd];
g = 9.81;

% Noise
noise.Enabler = 1; % 0/1
noise.pos_stand_dev = noise.Enabler * 0.0011;                            	%[m]
noise.vel_stand_dev = noise.Enabler * 0.01;                               %[m/s]
noise.acc_stand_dev = noise.Enabler * 0.5;                               %[m/s^2]
noise.attitude_stand_dev = noise.Enabler * deg2rad(0.33);                 %[rad]
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(1);                   %[rad/s]

seed.x = 1;
seed.vx = 2;
seed.theta = 3;
seed.q = 4;
seed.ax = 4;

% Delays
delay.position_filter = 1;
delay.attitude_filter = 1;
delay.mixer = 1;

% Load controller parameters
parameters_controller                    

% M injection example (sweeep: first column time vector, secondo column time history of pitching moment) 
load ExcitationM
SetPoint=[0,0];

% Values selected
t=ExcitationM(:,1);
simulation_time=t(end)-t(1);

%% TASK 1
load_system('Simulator_Single_Axis');
set_param('Simulator_Single_Axis', "FastRestart", "off");
simulation = sim('Simulator_Single_Axis', 'SrcWorkspace', 'current');

% Delete temporary files
if exist('slprj','dir')
    rmdir('slprj', 's')                                                    
end

time = 0:sample_time:simulation_time;

time_mask = (time>=65) & (time<=99);
time_filt= time(time_mask);

input = simulation.Mtot(time_mask);   % Input: normalized pitching moment
output1 = simulation.vx(time_mask);      % Output: pitch rate [rad/s]
output2 = simulation.q(time_mask);
output3 = simulation.theta(time_mask); 
output4 = simulation.ax(time_mask); % Output: Longitudinal (body) acceleration [rad/s^2]

figure;
subplot(3,1,1); 
plot(ExcitationM(:,1), ExcitationM(:,2), 'b', 'LineWidth', 1.5);
title('Input Signal: Normalized Pitching Moment', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('$\delta_{lon}$', 'Interpreter', 'latex', 'FontSize', 14);
grid on; xlim tight; set(gca, 'FontSize', 12);
subplot(3,1,2); 
plot(time_filt, output2, 'r', 'LineWidth', 1.5);
title('System Output: Pitch Rate Response', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('$q(t)$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 14);
grid on; xlim tight; set(gca, 'FontSize', 12);
subplot(3,1,3); 
plot(time_filt, output4, 'g', 'LineWidth', 1.5);
title('System Output: Longitudinal (body) acceleration Response', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('$a_x(t)\,[\mathrm{m/s^2}]$','Interpreter','latex', 'FontSize', 14)
grid on; xlim tight; set(gca, 'FontSize', 12);

% Grey-Box System Identification
C_model=[0, 1, 0; Xu, Xq, 0]; 
D_model=[0; Xd];
output_model=[output2 output4];
true_sys = ss(A, B, C_model, D_model, sample_time);  % True system 
true_tf  = tf(true_sys);
true_params = [-0.1068, 0.1192, -5.9755, -2.6478, -10.1647, 450.71]; % True parameter values

sys_id_model = iddata(output_model, input, sample_time);
% figure;
sys_ft = fft(sys_id_model);         % Frequency domain

Fs = 1/sample_time;

u  = simulation.Mtot;
yq = simulation.q;
yax = simulation.ax;

figure;
subplot(3,1,1)
plot(time,u), grid on
xlabel('Time [s]'), ylabel('u=Mtot')
title('Input')

subplot(3,1,2)
pwelch(u,[],[],[],Fs)
title('PSD input')

subplot(3,1,3)
[coh_q,f]  = mscohere(u,yq,[],[],[],Fs);
[coh_ax,~] = mscohere(u,yax,[],[],[],Fs);
plot(f,coh_q,'LineWidth',1.2); hold on
plot(f,coh_ax,'LineWidth',1.2);
grid on
xlabel('Frequency [Hz]')
ylabel('Coherence')
legend('u \rightarrow q','u \rightarrow ax')
title('Input-output coherence')

init_parameters={
    'Xu_g',  0;
    'Xq_g',   0;
    'Mu_g',  0;
    'Mq_g',  0;
    'Xd_g', 0;
    'Md_g', 0
};

sys_init = idgrey(@sys_init_creat, init_parameters, 'c');
grey_id = greyest(sys_ft, sys_init);  % Estimate parameters

cov_matrix = getcov(grey_id);  % Asymptotic covariance
%std_devs = sqrt(diag(cov_matrix));  % Standard deviation (1-sigma)

% Extract and reassign estimated values
params = getpvec(grey_id)';
Xu_g = params(1); 
Xq_g = params(2);
Mu_g = params(3); 
Mq_g = params(4);
Xd_g = params(5); 
Md_g = params(6);

% Compute absolute and percentage errors
abs_errors = abs(params - true_params);
percent_errors = 100 * abs_errors ./ abs(true_params);

% Compute RMSE
rmse = sqrt(mean((params' - true_params).^2));
rmse_percent = sqrt(mean(percent_errors.^2));

A_grey = grey_id.A;
B_grey = grey_id.B;
C_grey = grey_id.C;
D_grey = grey_id.D;

grey_idsys = ss(A_grey, B_grey, C_grey, D_grey, sample_time);
grey_idtf  = tf(grey_idsys);

% Generate the 3211 input
amplitudes = [0,  0.1, -0.1,  0.1, -0.1,  0];
time_step = 2; % [s]
durations = time_step*[1, 3, 2, 1, 1, 1];  
dt = sample_time; 
t_end = sum(durations);
t_3211 = 0:dt:t_end;
u_3211 = zeros(size(t_3211));
t_start = 0;
for i = 1:length(amplitudes)
    t_end_segment = t_start + durations(i);
    idx = (t_3211 >= t_start) & (t_3211 < t_end_segment);
    u_3211(idx) = amplitudes(i);
    t_start = t_end_segment;
end
figure;
clf;
plot(t_3211,u_3211,'LineWidth',2);
xlabel('Time [s]');
ylabel('Input');
title('3211 maneuver input');
grid on;
Excitation_3211 = [t_3211', u_3211'];
save('Excitation_3211.mat','Excitation_3211');

% Validation
load Excitation_3211
assignin('base', 'A_validation', A);
assignin('base', 'B_validation', B);
assignin('base', 'C_validation', C);
assignin('base', 'D_validation', D);
load_system('Validation');
set_param('Validation', "FastRestart", "off");
validation_sim1 = sim('Validation', 'SrcWorkspace', 'base');
true_output_validation = [validation_sim1.q, validation_sim1.ax];

assignin('base', 'A_validation', A_grey);
assignin('base', 'B_validation', B_grey);
assignin('base', 'C_validation', [1 0 0; C_grey(1,:); 0 0 1; C_grey(2,:)]);
assignin('base', 'D_validation', [0; D_grey(1,:); 0;  D_grey(2,:)]);
load_system('Validation');
set_param('Validation', "FastRestart", "off");
validation_sim2 = sim('Validation', 'SrcWorkspace', 'base');
output_validation_grey = [validation_sim2.q, validation_sim2.ax];

figure;
clf;
plot(t_3211', true_output_validation(:,1), 'k', 'LineWidth', 1); hold on;
plot(t_3211', output_validation_grey(:,1), 'r--', 'LineWidth', 1);
xlabel('Time [s]');
ylabel('Pitch Rate $q(t)$ [rad/s]', 'Interpreter', 'latex');
title('3211 Validation Pitch rate');
legend('True Output', 'Grey-Box Model');
grid on; set(gca, 'FontSize', 12);

figure;
clf;
plot(t_3211', true_output_validation(:,2), 'k', 'LineWidth', 1); hold on;
plot(t_3211', output_validation_grey(:,2), 'g--', 'LineWidth', 1);
xlabel('Time [s]');
ylabel('$a_x(t)\,[\mathrm{m/s^2}]$','Interpreter','latex', 'FontSize', 14)
title('3211 Validation Acceleration');
legend('True Output', 'Grey-Box Model');
grid on; set(gca, 'FontSize', 12);

abs_errors_valid = abs(output_validation_grey - true_output_validation);

figure;
subplot(2,1,1);
plot(t_3211', abs_errors_valid(:,1), 'k', 'LineWidth', 1); 
xlabel('Time [s]');
ylabel('Error', 'Interpreter', 'latex');
title('3211 Validation – Pitch rate error');
grid on; set(gca, 'FontSize', 12);
subplot(2,1,2);
plot(t_3211', abs_errors_valid(:,2), 'k', 'LineWidth', 1); 
xlabel('Time [s]');
ylabel('Error', 'Interpreter', 'latex');
title('3211 Validation – Acceleration error');
grid on; set(gca, 'FontSize', 12);



%% TASK 2: Noise1
noise.acc_stand_dev = noise.Enabler * 0.35;                               %[m/s^2]
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(1.3);                   %[rad/s]

load_system('Simulator_Single_Axis');
set_param('Simulator_Single_Axis', "FastRestart", "off");
simulation = sim('Simulator_Single_Axis', 'SrcWorkspace', 'current');

% Delete temporary files
if exist('slprj','dir')
    rmdir('slprj', 's')                                                    
end

input = simulation.Mtot(time_mask);   % Input: normalized pitching moment
output1_mod1 = simulation.vx(time_mask);      % Output: pitch rate [rad/s]
output2_mod1 = simulation.q(time_mask);
output3_mod1 = simulation.theta(time_mask); 
output4_mod1 = simulation.ax(time_mask); % Output: Longitudinal (body) acceleration [rad/s^2]

figure;
subplot(2,1,1); 
plot(time_filt, output2_mod1, 'r', 'LineWidth', 1.5);
title('System Output with noise modification1: Pitch Rate Response', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('$q(t)$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 14);
grid on; xlim tight; set(gca, 'FontSize', 12);
subplot(2,1,2); 
plot(time_filt, output4_mod1, 'g', 'LineWidth', 1.5);
title('System Output with noise modification1: Longitudinal (body) acceleration Response', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('$a_x(t)\,[\mathrm{m/s^2}]$','Interpreter','latex', 'FontSize', 14)
grid on; xlim tight; set(gca, 'FontSize', 12);

% Grey-Box System Identification with noise_modification1
output_model_mod1=[output2_mod1 output4_mod1];
sys_id_model_mod1 = iddata(output_model_mod1, input, sample_time);
sys_ft_mod1 = fft(sys_id_model_mod1);         % Frequency domain
grey_id_mod1 = greyest(sys_ft_mod1, sys_init);  % Estimate parameters

cov_matrix_mod1 = getcov(grey_id_mod1);  % Asymptotic covariance
%std_devs = sqrt(diag(cov_matrix));  % Standard deviation (1-sigma)

% Extract and reassign estimated values
params_mod1 = getpvec(grey_id_mod1)';
Xu_mod1 = params_mod1(1); 
Xq_mod1 = params_mod1(2);
Mu_mod1 = params_mod1(3); 
Mq_mod1 = params_mod1(4);
Xd_mod1 = params_mod1(5); 
Md_mod1 = params_mod1(6);

% Compute absolute and percentage errors
abs_errors_mod1 = abs(params_mod1 - true_params);
percent_errors_mod1 = 100 * abs_errors_mod1 ./ abs(true_params);

% Compute RMSE
rmse_mod1 = sqrt(mean((params_mod1' - true_params).^2));
rmse_percent_mod1 = sqrt(mean(percent_errors_mod1.^2));

A_grey_mod1 = grey_id_mod1.A;
B_grey_mod1 = grey_id_mod1.B;
C_grey_mod1 = grey_id_mod1.C;
D_grey_mod1 = grey_id_mod1.D;

grey_idsys_mod1 = ss(A_grey_mod1, B_grey_mod1, C_grey_mod1, D_grey_mod1, sample_time);
grey_idtf_mod1  = tf(grey_idsys_mod1);

abs_diff_output2_mod1 = abs(output2_mod1 - output2);
abs_diff_output4_mod1 = abs(output4_mod1 - output4);
figure;
subplot(2,1,1)
plot(time_filt, abs_diff_output2_mod1, 'r', 'LineWidth', 1.5);
title('Response difference with noise 1 in Pitch rate', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('$q(t)$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 14);
grid on; xlim tight; set(gca, 'FontSize', 12);
subplot(2,1,2)
plot(time_filt, abs_diff_output4_mod1, 'g', 'LineWidth', 1.5);
title('Response difference with noise 1 in Acceleration', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('$a_x(t)\,[\mathrm{m/s^2}]$','Interpreter','latex', 'FontSize', 14)
grid on; xlim tight; set(gca, 'FontSize', 12);

assignin('base', 'A_validation', A_grey_mod1);
assignin('base', 'B_validation', B_grey_mod1);
assignin('base', 'C_validation', [1 0 0; C_grey_mod1(1,:); 0 0 1; C_grey_mod1(2,:)]);
assignin('base', 'D_validation', [0; D_grey_mod1(1,:); 0;  D_grey_mod1(2,:)]);
load_system('Validation');
set_param('Validation', "FastRestart", "off");
validation_sim2 = sim('Validation', 'SrcWorkspace', 'base');
output_validation_grey_mod1 = [validation_sim2.q, validation_sim2.ax];

figure;
clf;
plot(t_3211', true_output_validation(:,1), 'k', 'LineWidth', 1); hold on;
plot(t_3211', output_validation_grey_mod1(:,1), 'r--', 'LineWidth', 1);
xlabel('Time [s]');
ylabel('Pitch Rate $q(t)$ [rad/s]', 'Interpreter', 'latex');
title('3211 Validation Pitch rate with noise modification 1');
legend('True Output', 'Grey-Box Model');
grid on; set(gca, 'FontSize', 12);

figure;
clf;
plot(t_3211', true_output_validation(:,2), 'k', 'LineWidth', 1); hold on;
plot(t_3211', output_validation_grey_mod1(:,2), 'g--', 'LineWidth', 1);
xlabel('Time [s]');
ylabel('$a_x(t)\,[\mathrm{m/s^2}]$','Interpreter','latex', 'FontSize', 14)
title('3211 Validation Acceleration with noise modification 1');
legend('True Output', 'Grey-Box Model');
grid on; set(gca, 'FontSize', 12);

abs_errors_valid_mod1 = abs(output_validation_grey_mod1 - true_output_validation);
figure;
subplot(2,1,1);
plot(t_3211', abs_errors_valid_mod1(:,1), 'k', 'LineWidth', 1); 
xlabel('Time [s]');
ylabel('Error', 'Interpreter', 'latex');
title('3211 Validation with noise modification 1 – Pitch rate error');
grid on; set(gca, 'FontSize', 12);
subplot(2,1,2);
plot(t_3211', abs_errors_valid_mod1(:,2), 'k', 'LineWidth', 1); 
xlabel('Time [s]');
ylabel('Error', 'Interpreter', 'latex');
title('3211 Validation with noise modification 1 – Acceleration error');
grid on; set(gca, 'FontSize', 12);

%% TASK 2: Noise2
noise.acc_stand_dev = noise.Enabler * 0.65;                               %[m/s^2]
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(0.7);                   %[rad/s]

load_system('Simulator_Single_Axis');
set_param('Simulator_Single_Axis', "FastRestart", "off");
simulation = sim('Simulator_Single_Axis', 'SrcWorkspace', 'current');

% Delete temporary files
if exist('slprj','dir')
    rmdir('slprj', 's')                                                    
end

input = simulation.Mtot(time_mask);   % Input: normalized pitching moment
output1_mod2 = simulation.vx(time_mask);      % Output: pitch rate [rad/s]
output2_mod2 = simulation.q(time_mask);
output3_mod2 = simulation.theta(time_mask); 
output4_mod2 = simulation.ax(time_mask); % Output: Longitudinal (body) acceleration [rad/s^2]

figure;
subplot(2,1,1); 
plot(time_filt, output2_mod2, 'r', 'LineWidth', 1.5);
title('System Output with noise modification3: Pitch Rate Response', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('$q(t)$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 14);
grid on; xlim tight; set(gca, 'FontSize', 12);
subplot(2,1,2); 
plot(time_filt, output4_mod2, 'g', 'LineWidth', 1.5);
title('System Output with noise modification3: Longitudinal (body) acceleration Response', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('$a_x(t)\,[\mathrm{m/s^2}]$','Interpreter','latex', 'FontSize', 14)
grid on; xlim tight; set(gca, 'FontSize', 12);

% Grey-Box System Identification with noise_modification2
output_model_mod2=[output2_mod2 output4_mod2];
sys_id_model_mod2 = iddata(output_model_mod2, input, sample_time);
sys_ft_mod2 = fft(sys_id_model_mod2);         % Frequency domain
grey_id_mod2 = greyest(sys_ft_mod2, sys_init);  % Estimate parameters

cov_matrix_mod2 = getcov(grey_id_mod2);  % Asymptotic covariance
%std_devs = sqrt(diag(cov_matrix));  % Standard deviation (1-sigma)

% Extract and reassign estimated values
params_mod2 = getpvec(grey_id_mod2)';
Xu_mod1 = params_mod2(1); 
Xq_mod1 = params_mod2(2);
Mu_mod1 = params_mod2(3); 
Mq_mod1 = params_mod2(4);
Xd_mod1 = params_mod2(5); 
Md_mod1 = params_mod2(6);

% Compute absolute and percentage errors
abs_errors_mod2 = abs(params_mod2 - true_params);
percent_errors_mod2 = 100 * abs_errors_mod2 ./ abs(true_params);

% Compute RMSE
rmse_mod2 = sqrt(mean((params_mod2' - true_params).^2));
rmse_percent_mod2 = sqrt(mean(percent_errors_mod2.^2));

A_grey_mod2 = grey_id_mod2.A;
B_grey_mod2 = grey_id_mod2.B;
C_grey_mod2 = grey_id_mod2.C;
D_grey_mod2 = grey_id_mod2.D;

grey_idsys_mod2 = ss(A_grey_mod2, B_grey_mod2, C_grey_mod2, D_grey_mod2, sample_time);
grey_idtf_mod2  = tf(grey_idsys_mod2);

abs_diff_output2_mod2 = abs(output2_mod2 - output2);
abs_diff_output4_mod2 = abs(output4_mod2 - output4);
figure;
subplot(2,1,1)
plot(time_filt, abs_diff_output2_mod2, 'r', 'LineWidth', 1.5);
title('Response difference with noise 3 in Pitch rate', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('$q(t)$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 14);
grid on; xlim tight; set(gca, 'FontSize', 12);
subplot(2,1,2)
plot(time_filt, abs_diff_output4_mod2, 'g', 'LineWidth', 1.5);
title('Response difference with noise 3 in Acceleration', 'FontSize', 14);
xlabel('Time [s]', 'FontSize', 12);
ylabel('$a_x(t)\,[\mathrm{m/s^2}]$','Interpreter','latex', 'FontSize', 14)
grid on; xlim tight; set(gca, 'FontSize', 12);

assignin('base', 'A_validation', A_grey_mod2);
assignin('base', 'B_validation', B_grey_mod2);
assignin('base', 'C_validation', [1 0 0; C_grey_mod2(1,:); 0 0 1; C_grey_mod2(2,:)]);
assignin('base', 'D_validation', [0; D_grey_mod2(1,:); 0;  D_grey_mod2(2,:)]);
load_system('Validation');
set_param('Validation', "FastRestart", "off");
validation_sim2 = sim('Validation', 'SrcWorkspace', 'base');
output_validation_grey_mod2 = [validation_sim2.q, validation_sim2.ax];

figure;
clf;
plot(t_3211', true_output_validation(:,1), 'k', 'LineWidth', 1); hold on;
plot(t_3211', output_validation_grey_mod2(:,1), 'r--', 'LineWidth', 1);
xlabel('Time [s]');
ylabel('Pitch Rate $q(t)$ [rad/s]', 'Interpreter', 'latex');
title('3211 Validation Pitch rate with noise modification 3');
legend('True Output', 'Grey-Box Model');
grid on; set(gca, 'FontSize', 12);

figure;
clf;
plot(t_3211', true_output_validation(:,2), 'k', 'LineWidth', 1); hold on;
plot(t_3211', output_validation_grey_mod2(:,2), 'g--', 'LineWidth', 1);
xlabel('Time [s]');
ylabel('$a_x(t)\,[\mathrm{m/s^2}]$','Interpreter','latex', 'FontSize', 14)
title('3211 Validation Acceleration with noise modification 3');
legend('True Output', 'Grey-Box Model');
grid on; set(gca, 'FontSize', 12);

abs_errors_valid_mod2 = abs(output_validation_grey_mod2 - true_output_validation);
figure;
subplot(2,1,1);
plot(t_3211', abs_errors_valid_mod2(:,1), 'k', 'LineWidth', 1); 
xlabel('Time [s]');
ylabel('Error', 'Interpreter', 'latex');
title('3211 Validation with noise modification 3 – Pitch rate error');
grid on; set(gca, 'FontSize', 12);
subplot(2,1,2);
plot(t_3211', abs_errors_valid_mod2(:,2), 'k', 'LineWidth', 1); 
xlabel('Time [s]');
ylabel('Error', 'Interpreter', 'latex');
title('3211 Validation with noise modification 3 – Acceleration error');
grid on; set(gca, 'FontSize', 12);

%%
% Comparison between the two noise configuration
comparison=abs(abs_errors_valid_mod1 - abs_errors_valid_mod2);
figure;
subplot(2,1,1);
plot(t_3211', comparison(:,1), 'k', 'LineWidth', 1); 
xlabel('Time [s]');
ylabel('Error', 'Interpreter', 'latex');
title('3211 Validation error comparison – Pitch rate');
grid on; set(gca, 'FontSize', 12);
subplot(2,1,2);
plot(t_3211', comparison(:,2), 'k', 'LineWidth', 1); 
xlabel('Time [s]');
ylabel('Error', 'Interpreter', 'latex');
title('3211 Validation error comparison – Acceleration');
grid on; set(gca, 'FontSize', 12);

%%
% Calculate FIT using 3211
% Pitch rate
y_true_mean_q = mean(true_output_validation(:,1));
FIT_GreyT1_q = 100 * (1 - ((norm(true_output_validation(:,1) - output_validation_grey(:,1))) / (norm(true_output_validation(:,1) - y_true_mean_q))));
FIT_Greymod1_q = 100 * (1 - ((norm(true_output_validation(:,1) - output_validation_grey_mod1(:,1))) / (norm(true_output_validation(:,1) - y_true_mean_q))));
FIT_Greymod2_q = 100 * (1 - ((norm(true_output_validation(:,1) - output_validation_grey_mod2(:,1))) / (norm(true_output_validation(:,1) - y_true_mean_q))));

% Acceleration
y_true_mean_ax = mean(true_output_validation(:,2));
FIT_GreyT1_ax = 100 * (1 - ((norm(true_output_validation(:,2) - output_validation_grey(:,2))) / (norm(true_output_validation(:,2) - y_true_mean_ax))));
FIT_Greymod1_ax = 100 * (1 - ((norm(true_output_validation(:,2) - output_validation_grey_mod1(:,2))) / (norm(true_output_validation(:,2) - y_true_mean_ax))));
FIT_Greymod2_ax = 100 * (1 - ((norm(true_output_validation(:,2) - output_validation_grey_mod2(:,2))) / (norm(true_output_validation(:,2) - y_true_mean_ax))));

% Calculate VAF using 3211
% Pitch rate
y_q=true_output_validation(:,1);
y_est_T1q=output_validation_grey(:,1);
y_est_mod1q=output_validation_grey_mod1(:,1);
y_est_mod2q=output_validation_grey_mod2(:,1);
vaf_GreyT1_q = max(diag(100*(eye(size(y_q,2))-cov(y_q-y_est_T1q)./cov(y_q))),0);
vaf_Greymod1_q = max(diag(100*(eye(size(y_q,2))-cov(y_q-y_est_mod1q)./cov(y_q))),0);
vaf_Greymod2_q = max(diag(100*(eye(size(y_q,2))-cov(y_q-y_est_mod2q)./cov(y_q))),0);
% Acceleration
y_ax=true_output_validation(:,2);
y_est_T1ax=output_validation_grey(:,2);
y_est_mod1ax=output_validation_grey_mod1(:,2);
y_est_mod2ax=output_validation_grey_mod2(:,2);
vaf_GreyT1_ax = max(diag(100*(eye(size(y_ax,2))-cov(y_ax-y_est_T1ax)./cov(y_ax))),0);
vaf_Greymod1_ax = max(diag(100*(eye(size(y_ax,2))-cov(y_ax-y_est_mod1ax)./cov(y_ax))),0);
vaf_Greymod2_ax = max(diag(100*(eye(size(y_ax,2))-cov(y_ax-y_est_mod2ax)./cov(y_ax))),0);

%% TASK 3
N=300;
Matrix_output1 = zeros(length(output1),N);
Matrix_output2 = zeros(length(output1),N);
Matrix_output3 = zeros(length(output1),N);
Matrix_output4 = zeros(length(output1),N);
Matrix_outputGrey_q = zeros(length(output_validation_grey),N);
Matrix_outputGrey_ax = zeros(length(output_validation_grey),N);
Matrix_poles = zeros(3,N);
Matrix_zeros_q = zeros(2,N);
Matrix_zeros_ax = zeros(3,N);
Array_Xu_Monte = zeros(1,N);
Array_Xq_Monte = zeros(1,N);
Array_Mu_Monte = zeros(1,N);
Array_Mq_Monte = zeros(1,N);
Array_Xd_Monte = zeros(1,N);
Array_Md_Monte = zeros(1,N);
Array_grey_id_tf = cell(1,N);
Arraymean_Xu_Monte = zeros(1,N);
Arraymean_Xq_Monte = zeros(1,N);
Arraymean_Mu_Monte = zeros(1,N);
Arraymean_Mq_Monte = zeros(1,N);
Arraymean_Xd_Monte = zeros(1,N);
Arraymean_Md_Monte = zeros(1,N);

for i = 1:N 
    seed.x = i;
    seed.vx = i + 1;
    seed.theta = i + 2;
    seed.q = i + 3;
    seed.ax = i + 3;
    load_system('Simulator_Single_Axis');
    set_param('Simulator_Single_Axis', "FastRestart", "off");
    simulation = sim('Simulator_Single_Axis', 'SrcWorkspace', 'current');

    if exist('slprj', 'dir')
        rmdir('slprj', 's');
    end

    output1_Montecarlo = simulation.vx(time_mask);      % Output: pitch rate [rad/s]
    output2_Montecarlo = simulation.q(time_mask);
    output3_Montecarlo = simulation.theta(time_mask); 
    output4_Montecarlo = simulation.ax(time_mask); % Output: Longitudinal (body) acceleration [rad/s^2]
    Matrix_output1(:,i) = output1_Montecarlo;
    Matrix_output2(:,i) = output2_Montecarlo;
    Matrix_output3(:,i) = output3_Montecarlo;
    Matrix_output4(:,i)= output4_Montecarlo;
    
    output_model_Montecarlo=[output2_Montecarlo output4_Montecarlo];

    sys_id_model_Montecarlo = iddata(output_model_Montecarlo, input, sample_time);
    sys_ft_Montecarlo = fft(sys_id_model_Montecarlo);         % Frequency domain
    
    sys_init_Montecarlo = idgrey(@sys_init_creat, init_parameters, 'c');
    grey_id_Montecarlo = greyest(sys_ft_Montecarlo, sys_init_Montecarlo);  % Estimate parameters
    Array_grey_id_tf{i} = ss(grey_id_Montecarlo);
    
    A_grey_Montecarlo = grey_id_Montecarlo.A;
    B_grey_Montecarlo = grey_id_Montecarlo.B;
    C_grey_Montecarlo = grey_id_Montecarlo.C;
    D_grey_Montecarlo = grey_id_Montecarlo.D;
    
    Xu_grey_Monte = A_grey_Montecarlo(1,1);
    Xq_grey_Monte = A_grey_Montecarlo(1,2);
    Mu_grey_Monte = A_grey_Montecarlo(2,1);
    Mq_grey_Monte = A_grey_Montecarlo(2,2);
    Xd_grey_Monte = B_grey_Montecarlo(1,1);
    Md_grey_Monte = B_grey_Montecarlo(2,1);
    Array_Xu_Monte(i) = Xu_grey_Monte;
    Array_Xq_Monte(i) = Xq_grey_Monte;
    Array_Mu_Monte(i) = Mu_grey_Monte;
    Array_Mq_Monte(i) = Mq_grey_Monte;
    Array_Xd_Monte(i) = Xd_grey_Monte;
    Array_Md_Monte(i) = Md_grey_Monte;    
    grey_idsys_Montecarlo = ss(A_grey_Montecarlo, B_grey_Montecarlo, C_grey_Montecarlo, D_grey_Montecarlo, sample_time);
    grey_idtf_Montecarlo  = tf(grey_idsys_Montecarlo);
    Matrix_poles(:,i)= pole(grey_idtf_Montecarlo);
    Matrix_zeros_q(:,i) = tzero(grey_idtf_Montecarlo(1));
    Matrix_zeros_ax(:,i) = tzero(grey_idtf_Montecarlo(2));

    % Validation
    assignin('base', 'A_validation', A_grey_Montecarlo);
    assignin('base', 'B_validation', B_grey_Montecarlo);
    assignin('base', 'C_validation', [1 0 0; C_grey_Montecarlo(1,:); 0 0 1; C_grey_Montecarlo(2,:)]);
    assignin('base', 'D_validation', [0; D_grey_Montecarlo(1,:); 0;  D_grey_Montecarlo(2,:)]);
    load_system('Validation');
    set_param('Validation', "FastRestart", "off");
    validation_sim2 = sim('Validation', 'SrcWorkspace', 'base');
    output_validation_grey_Monte = [validation_sim2.q, validation_sim2.ax];
    Matrix_outputGrey_q(:,i) = output_validation_grey_Monte(:,1);
    Matrix_outputGrey_ax(:,i) = output_validation_grey_Monte(:,2);
end

%%
% Validation plot
figure

for j=1:N
subplot(2,1,1);
plot(t_3211', Matrix_outputGrey_q(:,j));
hold on
xlabel('Time [s]');
ylabel('Output Pitch Rate $q(t)$ [rad/s]', 'Interpreter', 'latex');
title('Output Pitch Rate');
grid on; set(gca, 'FontSize', 12);

subplot(2,1,2);
plot(t_3211', Matrix_outputGrey_ax(:,j));
hold on
xlabel('Time [s]');
ylabel('Accel $a_x(t)$ [m/s^2]', 'Interpreter', 'latex', 'FontSize', 14);
title('Output Acceleration');
grid on; set(gca, 'FontSize', 12);
end

%% Convergence Parameters
figure;
for p = 1:N
    Arraymean_Xu_Monte(p) = mean(Array_Xu_Monte(1:p));
    Arraymean_Xq_Monte(p) = mean(Array_Xq_Monte(1:p));
    Arraymean_Mu_Monte(p) = mean(Array_Mu_Monte(1:p));
    Arraymean_Mq_Monte(p) = mean(Array_Mq_Monte(1:p));
    Arraymean_Xd_Monte(p) = mean(Array_Xd_Monte(1:p));
    Arraymean_Md_Monte(p) = mean(Array_Md_Monte(1:p));
end

subplot(3,2,1);
plot(Arraymean_Xu_Monte, 'LineWidth', 1.5);
xlabel('Run'); ylabel('Xu');
title('Running Mean - Xu');
grid on; set(gca, 'FontSize', 12);

subplot(3,2,2);
plot(Arraymean_Xq_Monte, 'LineWidth', 1.5);
xlabel('Run'); ylabel('Xq');
title('Running Mean - Xq');
grid on; set(gca, 'FontSize', 12);

subplot(3,2,3);
plot(Arraymean_Mu_Monte, 'LineWidth', 1.5);
xlabel('Run'); ylabel('Mu');
title('Running Mean - Mu');
grid on; set(gca, 'FontSize', 12);

subplot(3,2,4);
plot(Arraymean_Mq_Monte, 'LineWidth', 1.5);
xlabel('Run'); ylabel('Mq');
title('Running Mean - Mq');
grid on; set(gca, 'FontSize', 12);

subplot(3,2,5);
plot(Arraymean_Xd_Monte, 'LineWidth', 1.5);
xlabel('Run'); ylabel('Xd');
title('Running Mean - Xd');
grid on; set(gca, 'FontSize', 12);

subplot(3,2,6);
plot(Arraymean_Md_Monte, 'LineWidth', 1.5);
xlabel('Run'); ylabel('Md');
title('Running Mean - Md');
grid on; set(gca, 'FontSize', 12);

%%
% Histogram Estimated parameters
figure;
subplot(3,2,1);
h=histogram(Array_Xu_Monte);
xlabel('Xu');
ylabel('Number of occurrences', 'Interpreter', 'latex');
title('Xu estimated');
grid on; set(gca, 'FontSize', 12);

subplot(3,2,2);
h=histogram(Array_Xq_Monte);
xlabel('Xq');
ylabel('Number of occurrences', 'Interpreter', 'latex');
title('Xq estimated');
grid on; set(gca, 'FontSize', 12);

subplot(3,2,3);
h=histogram(Array_Mu_Monte);
xlabel('Mu');
ylabel('Number of occurrences', 'Interpreter', 'latex');
title('Mu estimated');
grid on; set(gca, 'FontSize', 12);

subplot(3,2,4);
h=histogram(Array_Mq_Monte);
xlabel('Mq');
ylabel('Number of occurrences', 'Interpreter', 'latex');
title('Mq estimated');
grid on; set(gca, 'FontSize', 12);

subplot(3,2,5);
h=histogram(Array_Xd_Monte);
xlabel('Xd');
ylabel('Number of occurrences', 'Interpreter', 'latex');
title('Xd estimated');
grid on; set(gca, 'FontSize', 12);

subplot(3,2,6);
h=histogram(Array_Md_Monte);
xlabel('Md');
ylabel('Number of occurrences', 'Interpreter', 'latex');
title('Md estimated');
grid on; set(gca, 'FontSize', 12);

% Histogram Poles/Zeros
figure;
subplot(3,2,1);
h=histogram(real(Matrix_poles));
xlabel('Real part of poles');
ylabel('Number of occurrences', 'Interpreter', 'latex');
title('Real part of poles');
grid on; set(gca, 'FontSize', 12);

subplot(3,2,2)
h=histogram(imag((Matrix_poles)));
xlabel('Imaginary part of poles');
ylabel('Number of occurrences', 'Interpreter', 'latex');
title('Imaginary part of poles');
grid on; set(gca, 'FontSize', 12);

subplot(3,2,3);
h=histogram(real(Matrix_zeros_q));
xlabel('Real part of zeros q');
ylabel('Number of occurrences', 'Interpreter', 'latex');
title('Real part of zeros q');
grid on; set(gca, 'FontSize', 12);

subplot(3,2,4);
h=histogram(imag(Matrix_zeros_q));
xlabel('Imaginary part of zeros q');
ylabel('Number of occurrences', 'Interpreter', 'latex');
title('Imaginary part of zeros q');
grid on; set(gca, 'FontSize', 12);

subplot(3,2,5);
h=histogram(real(Matrix_zeros_ax));
xlabel('Real part of zeros ax');
ylabel('Number of occurrences', 'Interpreter', 'latex');
title('Real part of zeros ax');
grid on; set(gca, 'FontSize', 12);

subplot(3,2,6);
h=histogram(imag(Matrix_zeros_ax));
xlabel('Imaginary part of zeros ax');
ylabel('Number of occurrences', 'Interpreter', 'latex');
title('Imaginary part of zeros ax');
grid on; set(gca, 'FontSize', 12);

% Bode
figure;
for k=1:N
bode(Array_grey_id_tf{k},logspace(-1,2,500))
hold on
title('Bode Diagrams');
end

figure;
for w=1:N
plot(real(Matrix_poles(:,w)), imag(Matrix_poles(:,w)),'x')
xlabel('Re(λ)');
ylabel('Im(λ)');
title('Pole Locations (Complex Plane)');
hold on
grid on;
end

%%
function [A_sys, B_sys, C_sys, D_sys] = sys_init_creat(Xu_sys, Xq_sys, Mu_sys, Mq_sys, Xd_sys, Md_sys, ~)

g = 9.81;

A_sys = [Xu_sys Xq_sys -g; Mu_sys Mq_sys 0; 0 1 0];

B_sys = [Xd_sys Md_sys 0]';

C_sys = [0 1 0; Xu_sys Xq_sys 0];

D_sys = [0; Xd_sys];

end

