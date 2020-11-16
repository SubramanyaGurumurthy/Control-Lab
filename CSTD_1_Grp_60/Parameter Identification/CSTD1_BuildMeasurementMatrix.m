    % #########################################################################
% TUHH :: Institute for Control Systems :: Control Lab
% #########################################################################
% Experiment CSTD1: Identification and Control of a Torsional Plant
%
% Copyright Herbert Werner and Hamburg University of Technology, 2014
% #########################################################################
% This file is to be completed by the student.
% The completed version is to be published using 
%   publish('CSTD1_BuildMeasurementMatrix','pdf') 
% and submitted as a pdf-file at least one week prior to the scheduled date
% for the experiment
%
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!  The gaps in the code are denoted by XXX !!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
% HINT 1:
% if you want to find out more about a certain command, just type 
% 'help command' into the matlab window
% HINT 2:
% use section evaluations (Ctrl+Enter) to run the code within a single 
% section

%----------------------------
% v.1.2 - 30-10-2017
% by Antonio Mendez G
%----------------------------

%% 1 Filter the Signals
%
% Design a filter to filter measurement noise

Ts = 1e-3; % sampling time
N = 2; %Order of the system
[B,A] = butter(N,0.14,'low');

% Filter the mesurements
Theta1_ff = filtfilt(B, A, Theta1);  % To filter out the measurement noise
Theta2_ff = filtfilt(B, A, Theta2);
Theta3_ff = filtfilt(B, A,Theta3);
Plant_input_u_ff = filtfilt(B, A,Plant_input_u);

% Calculate angular velocity from input vector in rad/s
Theta1_dot(1) = (Theta1_ff(2))/(2*Ts);
Theta2_dot(1) = (Theta2_ff(2))/(2*Ts);
Theta3_dot(1) = (Theta3_ff(2))/(2*Ts);

for i=2:length(Theta1)-1
    Theta1_dot(i,1) = (Theta1_ff(i+1) - Theta1_ff(i-1))/(2*Ts);
    Theta2_dot(i,1) = (Theta2_ff(i+1) - Theta2_ff(i-1))/(2*Ts);
    Theta3_dot(i,1) = (Theta3_ff(i+1) - Theta3_ff(i-1))/(2*Ts);
end   

% Calculate angular acceleration from input vector in rad/s/s
Theta1_ddot = (Theta1_ff(2)-(2*Theta1_ff(1)))/(Ts^2);
Theta2_ddot = (Theta2_ff(2)-(2*Theta2_ff(1)))/(Ts^2);
Theta3_ddot = (Theta3_ff(2)-(2*Theta3_ff(1)))/(Ts^2);


for i=2:length(Theta1)-1
    Theta1_ddot(i,1) = (Theta1_ff(i+1)+ Theta1_ff(i-1)-(2*Theta1_ff(i)))/(Ts^2);
    Theta2_ddot(i,1) = (Theta2_ff(i+1)+ Theta2_ff(i-1)-(2*Theta2_ff(i)))/(Ts^2);
    Theta3_ddot(i,1) = (Theta3_ff(i+1)+ Theta3_ff(i-1)-(2*Theta3_ff(i)))/(Ts^2);
end

% Filter Calculated angular velocity:
Theta1_dot_ff = filtfilt(B,A,Theta1_dot);
Theta2_dot_ff = filtfilt(B,A,Theta2_dot);
Theta3_dot_ff = filtfilt(B,A,Theta3_dot);

% Filter Calculated angular acceleration:
Theta1_ddot_ff = filtfilt(B, A, Theta1_ddot);
Theta2_ddot_ff = filtfilt(B, A, Theta2_ddot);
Theta3_ddot_ff = filtfilt(B, A, Theta3_ddot);


% Find the number of samples (ThetaX_ff/ThetaX_dot_ff/ThetaX_ddot_ff/)
Theta1_ff        = Theta1_ff(2:end-1);
Theta2_ff        = Theta2_ff(2:end-1);
Theta3_ff        = Theta3_ff(2:end-1);
Plant_input_u_ff = Plant_input_u_ff(2:end-1);

no_samples = length(Theta1_ff) - 200; %Adjustable   % 

% Redefine the new signals
Theta1_ff = Theta1_ff(1:no_samples);
Theta2_ff = Theta2_ff(1:no_samples);
Theta3_ff = Theta3_ff(1:no_samples);

Theta1_dot_ff = Theta1_dot_ff(1:no_samples);
Theta2_dot_ff = Theta2_dot_ff(1:no_samples);
Theta3_dot_ff = Theta3_dot_ff(1:no_samples);

Theta1_ddot_ff = Theta1_ddot_ff(1:no_samples);
Theta2_ddot_ff = Theta2_ddot_ff(1:no_samples);
Theta3_ddot_ff = Theta3_ddot_ff(1:no_samples);

Plant_input_u_ff = Plant_input_u_ff(1:no_samples);

% -------------------------------------------------------------- 
%% 2 Build Measurement Matrix M
% 
% Below is implemented the M matrix

M1 = [ Theta1_ddot_ff,      zeros(no_samples,1), zeros(no_samples,1) ;
       zeros(no_samples,1), Theta2_ddot_ff,      zeros(no_samples,1) ;
       zeros(no_samples,1), zeros(no_samples,1), Theta3_ddot_ff      ];

M2 = [ Theta1_dot_ff,       zeros(no_samples,1), zeros(no_samples,1) ;
       zeros(no_samples,1), Theta2_dot_ff,       zeros(no_samples,1) ;
       zeros(no_samples,1), zeros(no_samples,1), Theta3_dot_ff       ];

M3 = [ Theta1_ff - Theta2_ff, zeros(no_samples,1)   ;
       Theta2_ff - Theta1_ff, Theta2_ff - Theta3_ff ;
       zeros(no_samples,1),   Theta3_ff - Theta2_ff ];

M = [M1 M2 M3]*pi/180; %rad/s

% -------------------------------------------------------------- 
%% 3 Build Matrix Torque:
%
% The torque signal is defined

Torque = zeros(3*no_samples,1);
Torque(1:no_samples,1) =Plant_input_u_ff;