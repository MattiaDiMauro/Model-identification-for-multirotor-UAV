%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUADROTOR PARAMETERS - Controller                                       %
% Authors:  Mattia Giurato (mattia.giurato@polimi.it)                     %
%           Paolo Gattazzo (paolo.gattazzo@polimi.it)                     %
% Date: 23/04/2018                                                        %
% Adapted to ANT-X 2DoF by Salvatore Meraglia (salvatore.meraglia@polimi.it) %
% Date: 22/12/2022                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Normalized mixer

Kt = 0.25;
Kq = 0.25;
b = 2;

mixer = [ -1/(4*Kt), -2^(1/2)/(4*Kt*b),  2^(1/2)/(4*Kt*b),  1/(4*Kq);
 -1/(4*Kt),  2^(1/2)/(4*Kt*b), -2^(1/2)/(4*Kt*b),  1/(4*Kq);
 -1/(4*Kt),  2^(1/2)/(4*Kt*b),  2^(1/2)/(4*Kt*b), -1/(4*Kq);
 -1/(4*Kt), -2^(1/2)/(4*Kt*b), -2^(1/2)/(4*Kt*b), -1/(4*Kq)];

%% Angular rate controller

%q
KF_Q = 0.0;
KP_Q = 0.09;
KI_Q = 0.21;
KD_Q = 0.0016;
M_MAX = 1;
M_MIN = -1;

N_filter_rate = 100;
sample_time = 1/250;

%% Attitude controller

%theta
KP_PITCH = 12;
Q_MAX = 10; % 10 [rad/s] ~ 600 deg/s
Q_MIN = -10;

%% Velocity controller

KF_X_DOT = 0.0;
KP_X_DOT = 0.5;
KI_X_DOT = 0.04;
KD_X_DOT = 0;

F_X_MIN = -1;
F_X_MAX = 1; 

N_filter_vel = 10;

%% Position controller

KP_X = 2;

% velocity setpoint saturation
VEL_MIN = -10;
VEL_MAX = 10;

%% Baseline thrust (normalized)

BASELINE_T = 0.4;
