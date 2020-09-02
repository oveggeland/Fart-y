% M-script for numerical integration of the attitude dynamics of a rigid 
% body represented by unit quaternions. The MSS m-files must be on your
% Matlab path in order to run the script.
%
% System:                      .
%                              q = T(q)w
%                              .
%                            I w - S(Iw)w = tau
% Control law:
%                            tau = constant
% 
% Definitions:             
%                            I = inertia matrix (3x3)
%                            S(w) = skew-symmetric matrix (3x3)
%                            T(q) = transformation matrix (4x3)
%                            tau = control input (3x1)
%                            w = angular velocity vector (3x1)
%                            q = unit quaternion vector (4x1)
%
% Author:                   2018-08-15 Thor I. Fossen and Håkon H. Helgesen

%% USER INPUTS
h = 0.01;                     % sample time (s)
N  = 4000;                    % number of samples. Should be adjusted

% model parameters
m = 180;
r = 2;
I = m*r^2*eye(3);            % inertia matrix
I_inv = inv(I);

% constants
deg2rad = pi/180;   
rad2deg = 180/pi;

phi = -5*deg2rad;            % initial Euler angles
theta = 10*deg2rad;
psi = -20*deg2rad;

q = euler2q(phi,theta,psi);   % transform initial Euler angles to q

w = [0 0 0]';                 % initial angular rates

table = zeros(N+1,20);        % memory allocation

%% FOR-END LOOP
for i = 1:N+1
   t = (i-1)*h;                  % time
%% q reference
   phi_t = 0;
   theta_t = 15*cos(0.1*t);
   psi_t = 10*sin(0.05*t);
   q_ref_euler = [phi_t, theta_t, psi_t];  %for plot
   q_ref = euler2q(phi_t,theta_t,psi_t);
   q_bar = [q_ref(1); -q_ref(2:4)];
   q_thilde = quatprod(q_bar,q);
   
%% w reference
   Theta_dot = [0;-1.5*sin(0.1*t);0.5*cos(0.05*t)];
   T = Tzyx(phi_t,theta_t);
   w_ref = T\Theta_dot;
   w_thilde = w-w_ref;  
 
%% control
   tau = -400*eye(3)*w_thilde-20*q_thilde(2:4); 
   
   [phi,theta,psi] = q2euler(q); % transform q to Euler angles
   [J,J1,J2] = quatern(q);       % kinematic transformation matrices
   
   q_dot = J2*w;                        % quaternion kinematics
   w_dot = I_inv*(Smtrx(I*w)*w + tau);  % rigid-body kinetics
   
   table(i,:) = [t q' phi theta psi w' tau' w_ref' q_ref_euler];  % store data in table
   
   q = q + h*q_dot;	             % Euler integration
   w = w + h*w_dot;
   
   q  = q/norm(q);               % unit quaternion normalization
end 

%% PLOT FIGURES
t       = table(:,1);  
q_2     = table(:,2:5); 
phi     = rad2deg*table(:,6);
theta   = rad2deg*table(:,7);
psi     = rad2deg*table(:,8);
w       = rad2deg*table(:,9:11);
tau     = table(:,12:14);
w_thilde = rad2deg*table(:,15:17);
q_ref_table   = table(:,18:20);


figure (1); clf;
hold on;
plot(t, phi, 'b');
plot(t, theta, 'r');
plot(t, psi, 'g');
plot(t, q_ref_table(:,1), 'k');
plot(t, q_ref_table(:,2), 'y');
plot(t, q_ref_table(:,3), 'c');
hold off;
grid on;
legend('\phi', '\theta', '\psi','\phi bar','\theta bar','\psi bar');
title('Euler angles');
xlabel('time [s]'); 
ylabel('angle [deg]');

figure (2); clf;
hold on;
plot(t, w(:,1), 'b');
plot(t, w(:,2), 'r');
plot(t, w(:,3), 'g');
plot(t, w_thilde(:,1), 'k');
plot(t, w_thilde(:,2), 'y');
plot(t, w_thilde(:,3), 'c');
hold off;
grid on;
legend('x', 'y', 'z','blå ref','rød ref','grønn ref');
title('Angular velocities');
xlabel('time [s]'); 
ylabel('angular rate [deg/s]');

figure (3); clf;
hold on;
plot(t, tau(:,1), 'b');
plot(t, tau(:,2), 'r');
plot(t, tau(:,3), 'g');
hold off;
grid on;
legend('x', 'y', 'z');
title('Control input');
xlabel('time [s]'); 
ylabel('input [Nm]');