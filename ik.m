1;

function T = getTMatrix(angle_prox, angle_inter, angle_distal, len_prox, len_inter, len_distal, lb_angle_prox, ub_angle_prox, lb_angle_inter, ub_angle_inter, lb_angle_distal, ub_angle_distal)
  T = zeros(4)
  if (angle_prox < lb_angle_prox)
    angle_prox = lb_angle_prox
  elseif (angle_prox > ub_angle_prox)
    angle_prox = ub_angle_prox
  endif
    
  if (angle_inter < lb_angle_inter)
    angle_prox = lb_angle_inter
  elseif (angle_inter > ub_angle_inter)
    angle_inter = ub_angle_inter
  endif
    
  if (angle_distal < lb_angle_distal)
    angle_distal = lb_angle_distal
  elseif (angle_distal > ub_angle_distal)
    angle_distal = ub_angle_distal
  endif
  
  T_0_1 = zeros(4)
  T_1_2 = zeros(4)
  T_2_3 = zeros(4)

  T_0_1(1) = cos(angle_prox)
  T_0_1(2) = -sin(angle_prox)
  T_0_1(4) = len_prox*cos(angle_prox)
  T_0_1(5) = sin(angle_prox)
  T_0_1(6) = cos(angle_prox)
  T_0_1(8) = len_prox*sin(angle_prox) 

  T_1_2(1) = cos(angle_inter)
  T_1_2(2) = -sin(angle_inter)
  T_1_2(4) = len_prox*cos(angle_inter)
  T_1_2(5) = sin(angle_inter)
  T_1_2(6) = cos(angle_inter)
  T_1_2(8) = len_prox*sin(angle_inter)

  T_2_3(1) = cos(angle_distal)
  T_2_3(2) = -sin(angle_distal)
  T_2_3(4) = len_prox*cos(angle_distal)
  T_2_3(5) = sin(angle_distal)
  T_2_3(6) = cos(angle_distal)
  T_2_3(8) = len_prox*sin(angle_distal)

  T = T_0_1 * T_1_2 * T_2_3

endfunction

function inv_J_q = calcInvJacobian(target, angle_prox, angle_inter, angle_distal)

  %target = { x, y }
  q = { angle_prox, angle_inter }

  J_q = zeros(2)

  J_q(1) = - len_prox * sin(angle_prox) - len_inter * sin(angle_prox * angle_inter ) - len_distal * sin( angle_prox + (5/3) * angle_inter )
  J_q(2) = - len_inter * sin(angle_prox * angle_inter ) - len_distal * sin( angle_prox + (5/3) * angle_inter ) 
  J_q(3) = len_prox * cos(angle_prox) + len_inter * cos(angle_prox * angle_inter ) - len_distal * cos( angle_prox + (5/3) * angle_inter )
  J_q(4) = len_inter * cos(angle_prox * angle_inter ) - len_distal * cos( angle_prox + (5/3) * angle_inter )

  D = J_q(1)*J_q(4) - J_q(2)*J_q(3)

  inv_J_q = 1/D * J_q

endfunction


% Visualization
%img = figure(1)
%clf
%hold on

global len_prox = 44.6
global len_inter = 26.3
global len_distal = 17.4

% Restrictions TODO : Triple-check the restrictions
% -pi / 3 <= Theta_M <= pi /3
% -2 pi / 3 <= Theta_P <= 0
% -2 pi / 3 <= Theta _D <= 0

global lb_angle_prox   = -pi /3
global ub_angle_prox   = pi/3
global lb_angle_inter  = -2*pi/3
global ub_angle_inter  = 0
global lb_angle_distal = -2*pi/3
global ub_angle_distal = 0

target = {len_prox + len_inter + len_distal, -2}

global error_margin = 0.001

getTMatrix(1,1,1,2,2,2,0,0,0,0,0,0)

% Object O = {(x,y,z) in R | y <= -2}

% IK Method : 
% nextAnglesGuess = currentAnglesGuess + inverse( J_q(i) )*( target - currentPosition ) 
% continue until the guess is within error margin