pkg load symbolic

% Visualization
%img = figure(1)
%clf
%hold on

% Object O = {(x,y,z) in R | y <= -2}

len_prox = 44.6
len_inter = 26.3
len_distal = 17.4

syms angle_prox
syms angle_inter
syms angle_distal

% Restrictions TODO : Triple-check the restrictions
% -pi / 3 <= Theta_M <= pi /3
% -2 pi / 3 <= Theta_P <= 0
% -2 pi / 3 <= Theta _D <= 0
lb_angle_prox = -pi /3
ub_angle_prox = pi/3

assume(lb_angle_prox <= angle_prox)
assumeAlso(angle_prox <= ub_angle_prox) 
assumptions

%
%      R  t
% T = 
%      0  1
%

T_0_1 = zeros(4)
T_1_2 = zeros(4)
T_2_3 = zeros(4)

T_0_1(0) = cos(angle_prox)
T_0_1(1) = -sin(angle_prox)
T_0_1(3) = len_prox*cos(angle_prox)
T_0_1(4) = sin(angle_prox)
T_0_1(5) = cos(angle_prox)
T_0_1(7) = len_prox*sin(angle_prox) 

T_1_2(0) = cos(angle_inter)
T_1_2(1) = -sin(angle_inter)
T_1_2(3) = len_prox*cos(angle_inter)
T_1_2(4) = sin(angle_inter)
T_1_2(5) = cos(angle_inter)
T_1_2(7) = len_prox*sin(angle_inter)

T_2_3(0) = cos(angle_distal)
T_2_3(1) = -sin(angle_distal)
T_2_3(3) = len_prox*cos(angle_distal)
T_2_3(4) = sin(angle_distal)
T_2_3(5) = cos(angle_distal)
T_2_3(7) = len_prox*sin(angle_distal)

T = T_0_1 * T_1_2 * T_2_3

sym x
sym y

target = { x, y }

q = { angle_prox, angle_inter }

J_q = zeros(2)

J_q(0) = - len_prox * sin(angle_prox) - len_inter * sin(angle_prox * angle_inter ) - len_distal * sin( angle_prox + (5/3) * angle_inter )
J_q(1) = - len_inter * sin(angle_prox * angle_inter ) - len_distal * sin( angle_prox + (5/3) * angle_inter ) 
J_q(2) = len_prox * cos(angle_prox) + len_inter * cos(angle_prox * angle_inter ) - len_distal * cos( angle_prox + (5/3) * angle_inter )
J_q(3) = len_inter * cos(angle_prox * angle_inter ) - len_distal * cos( angle_prox + (5/3) * angle_inter )

D = J_q(0)*J_q(3) - J_q(1)*J_q(2)

inv_J_q = 1/D * J_q

% IK Method : 
% nextAnglesGuess = currentAnglesGuess + inverse( J_q(i) )*( target - currentPosition ) 
% continue until the guess is within error margin