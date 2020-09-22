close all;
a1 = 2.87;
a2 = -0.65;
delta_max = 30;
e_max = 15;
kp = 2*sign(a2);
zeta_phi = 0.707;
omega_n = sqrt(abs(a2)*delta_max/e_max);
kd = (2*zeta_phi*omega_n-a1)/a2;
disp(kp); disp(kd);
g = 9.81;
Vg = 580/3.6;
ki = 0.85;
kp_x = 7.5;

k = (-3:0.1:0);
k2 = (0:0.1:10);
sys = tf([0 0 0 a2], [1 (a1+a2*kd) a2*kp 0]);

A = [[-0.322 0.052 0.028 -1.12 0.002];
    [0 0 1 -0.001 0];
    [-10.6 0 -2.87 0.46 -0.65];
    [6.87 0 -0.04 -0.32 -0.02];
    [0 0 0 0 -7.5]];

%rlocus(sys,k);
rlocus(sys,k2);

%[r, k] = rlocus(sys, k);
%[r2, k2]?= rlocus(sys,k2);

%disp(r);
%disp(k);

