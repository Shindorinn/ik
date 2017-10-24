1;
% Clamp functions as expected
function f = clamp(angle, min, max)
  f=angle;
  if(angle<min) f=min;
  elseif(angle>max) f=max; endif
endfunction

function T = rotateJoint(angle,length,lb,ub)
  angle = clamp(angle,lb,ub);
  T = [ cos(angle), -sin(angle), 0, length*cos(angle);
        sin(angle),  cos(angle), 0, length*sin(angle);
                 0,           0, 1,                 0;
                 0,           0, 0,                 1];
endfunction

% Works as expected
function inv_j = calcInvJacobian(target, 
  len_prox, len_inter, len_distal,
  angle_prox, angle_inter, angle_distal)
  theta1 = angle_prox;
  theta2 = angle_inter;
  j = [-44.6*sin(theta1       )-26.3*sin(theta1+      theta2)-17.4*sin(theta1+(5/3)*theta2), -26.3*sin(theta1+theta2)-43.5*sin(theta1+(5/3)*theta2);
        44.6*cos(theta1       )+26.3*cos(theta1+      theta2)+17.4*cos(theta1+(5/3)*theta2),  26.3*cos(theta1+theta2)+17.4*cos(theta1+(5/3)*theta2)]
%  a = -len_prox *sin(angle_prox               ) - len_inter  * sin(angle_prox +       angle_inter) - len_distal * sin( angle_prox + (5/3)*angle_inter );
%  b = -len_inter*sin(angle_prox + angle_inter ) - len_distal * sin(angle_prox + (5/3)*angle_inter);
%  c =  len_prox *cos(angle_prox               ) + len_inter  * cos(angle_prox +       angle_inter) + len_distal * cos( angle_prox + (5/3)*angle_inter );
%  d =  len_inter*cos(angle_prox + angle_inter ) + len_distal * cos(angle_prox + (5/3)*angle_inter);

  
  D = j(1,1)* j(2,2) - j(1,2)*j(2,1);
  inv_j = (1/D) * [j(2,2), -1*j(1,2); -1*j(2,1), j(1,1)];
  %J_q = [ a, b;
  %        c, d];
  %inv_J_q = (1/(a*d-b*c))*[d,-b;-c,a];

endfunction


% Visualization
img = figure(1);
clf;
hold on;
axis([0,100,-100,100],"equal");

global len_prox = 44.6;
global len_inter = 26.3;
global len_distal = 17.4; 

global lb_angle_prox   = -pi /3;
global lb_angle_inter  = -2*pi/3;
global lb_angle_distal = -2*pi/3;

global ub_angle_prox   = pi/3;
global ub_angle_inter  = 0;
global ub_angle_distal = 0;



global error_margin = 0.1;

% first angle: 0.4708, lb_angle_inter/2, lb_angle_distal
%T01 = rotateJoint(0.4708,          len_prox,  lb_angle_prox,  ub_angle_prox);
%T12 = rotateJoint(lb_angle_inter/3,len_inter, lb_angle_inter, ub_angle_inter);
%T23 = rotateJoint(lb_angle_distal, len_distal,lb_angle_distal,ub_angle_distal);
%T = T01*T12*T23;
%
%
origin = transpose([0,0,0,1]);
%joint0 = T01*origin;
%joint1 = T01*T12*origin;
%tip    = T01*T12*T23*origin;
%
%x = tip(1)
%y = tip(2)
%
%plot([origin(1), joint0(1), joint1(1), tip(1)],
%     [origin(2), joint0(2), joint1(2), tip(2)],
%     'marker','o')

plot([0,65],[-2,-2]);

%distance = sqrt(x*x+y*y)
%len_prox
%len_inter+len_distal

target = transpose([len_prox + len_inter + 10, -2]);
plot(target(1), target(2),
     'marker','x','markersize',40,'color','r');

%currentAngleGuess = [unifrnd(lb_angle_prox,  ub_angle_prox),
%                     unifrnd(lb_angle_inter, ub_angle_inter)];
currentAngleGuess = [0.1,-0.1];

i = 0;
do
  i = i+1;
  distalJoint =  2*currentAngleGuess(2)/3 ;
  T01 = rotateJoint(currentAngleGuess(1), len_prox,  lb_angle_prox,   ub_angle_prox);
  T12 = rotateJoint(currentAngleGuess(2), len_inter, lb_angle_inter,  ub_angle_inter);
  T23 = rotateJoint(distalJoint,          len_distal,lb_angle_distal, ub_angle_distal);
  j1 = T01;
  j2 = T01*T12;
  tip= T01*T12*T23;
  
  tipsx(i) = tip(1,4);
  tipsy(i) = tip(2,4);
  
  theta1 = clamp(currentAngleGuess(1), lb_angle_prox, ub_angle_prox);
  theta2 = clamp(currentAngleGuess(2), lb_angle_inter, ub_angle_inter);
  theta3 = (2*theta2)/3

  a = tip(1,4)
  b = tip(2,4)

  x = l1*cos(theta1) + l3*cos((2*theta2)/3)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)) - l3*sin((2*theta2)/3)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) + l2*cos(theta1)*cos(theta2) - l2*sin(theta1)*sin(theta2) 
  y = l1*sin(theta1) + l3*cos((2*theta2)/3)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) + l3*sin((2*theta2)/3)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)) + l2*cos(theta1)*sin(theta2) + l2*cos(theta2)*sin(theta1)


  plot([0,j1(1,4),j2(1,4),tip(1,4)],
       [0,j1(2,4),j2(2,4),tip(2,4)],
       'marker','o','color','k');
       
  inv_J_q = calcInvJacobian(target, len_prox, len_inter, len_distal, currentAngleGuess(1), currentAngleGuess(2), distalJoint);
  currentAngleGuess = currentAngleGuess + inv_J_q*( target - [tip(1,4);tip(2,4)] );
  
  error = sqrt( (target(1)-tip(1,4))^2 + (target(2)-tip(2,4))^2 );
until((error<error_margin && tip(2,4)>-2) || i>10);

i = i+1;

distalJoint =  2*currentAngleGuess(2)/3 ;
T01 = rotateJoint(currentAngleGuess(1), len_prox,  lb_angle_prox,   ub_angle_prox);
T12 = rotateJoint(currentAngleGuess(2), len_inter, lb_angle_inter,  ub_angle_inter);
T23 = rotateJoint(distalJoint,          len_distal,lb_angle_distal, ub_angle_distal);
j1 = T01;
j2 = T01*T12;
tip= T01*T12*T23;
tipsy(i) = tip(2,4);
tipsx(i) = tip(1,4);
  plot([0,j1(1,4),j2(1,4),tip(1,4)],
       [0,j1(2,4),j2(2,4),tip(2,4)],
       'marker','o','color','b');


plot(tipsx,tipsy,'linestyle','-.','color','r');

l1 = 44.6;
l2 = 26.3;
l3 = 17.4;

theta1 = clamp(currentAngleGuess(1), lb_angle_prox, ub_angle_prox);
theta2 = clamp(currentAngleGuess(2), lb_angle_inter, ub_angle_inter);
theta3 = (2*theta2)/3

a = tip(1,4)
b = tip(2,4)

x = l1*cos(theta1) + l3*cos((2*theta2)/3)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)) - l3*sin((2*theta2)/3)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) + l2*cos(theta1)*cos(theta2) - l2*sin(theta1)*sin(theta2) 
y = l1*sin(theta1) + l3*cos((2*theta2)/3)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) + l3*sin((2*theta2)/3)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)) + l2*cos(theta1)*sin(theta2) + l2*cos(theta2)*sin(theta1)



%origin = transpose([0,0,0,1]);
%joint0 = T01*origin;
%joint1 = T01*T12*origin;
%tip    = T01*T12*T23*origin;
%
%plot([origin(1), joint0(1), joint1(1), tip(1)],
%     [origin(2), joint0(2), joint1(2), tip(2)],
%     'marker','o')

% Object O = {(x,y,z) in R | y <= -2}

% IK Method : 
% nextAnglesGuess = currentAnglesGuess + inverse( J_q(i) )*( target - currentPosition ) 
% continue until the guess is within error margin