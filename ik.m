1;

% Start of Global Script
global len_prox = 44.6;
global len_inter = 26.3;
global len_distal = 17.4; 

global lb_angle_prox   = -pi /3;
global lb_angle_inter  = -2*pi/3;
global lb_angle_distal = -2*pi/3;

global ub_angle_prox   = pi/3;
global ub_angle_inter  = 0;
global ub_angle_distal = 0;

global error_margin = 0.001;


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
function inv_j = calcInvJacobian(len_prox, len_inter, len_distal, angle_prox, angle_inter)

  j = [-len_prox*sin(angle_prox )-len_inter*sin(angle_prox + angle_inter)-len_distal*sin(angle_prox+(5/3)*angle_inter), -len_inter*sin(angle_prox+angle_inter)-len_prox*sin(angle_prox+(5/3)*angle_inter);
        len_prox*cos(angle_prox )+len_inter*cos(angle_prox + angle_inter)+len_distal*cos(angle_prox+(5/3)*angle_inter),  len_inter*cos(angle_prox+angle_inter)+len_distal*cos(angle_prox+(5/3)*angle_inter)];
  
  D = j(1,1)* j(2,2) - j(1,2)*j(2,1);
  inv_j = (1/D) * [j(2,2), -1*j(1,2); -1*j(2,1), j(1,1)];

endfunction

function [X,Y] = forwardKinematics(angle_prox, angle_inter, angle_distal)
  global len_prox;
  global len_inter;
  global len_distal; 

  global lb_angle_prox;
  global lb_angle_inter;
  global lb_angle_distal;

  global ub_angle_prox;
  global ub_angle_inter;
  global ub_angle_distal;

  global error_margin;

  T1 = rotateJoint(angle_prox,   len_prox,  lb_angle_prox,   ub_angle_prox);
  T2 = rotateJoint(angle_inter,  len_inter, lb_angle_inter,  ub_angle_inter);
  T3 = rotateJoint(angle_distal, len_distal,lb_angle_distal, ub_angle_distal);
  j1 = T1;
  j2 = T1*T2;
  tip= T1*T2*T3;
  X = [0,j1(1,4),j2(1,4),tip(1,4)];
  Y = [0,j1(2,4),j2(2,4),tip(2,4)];
endfunction

function K = assignment2()
  global len_prox;
  global len_inter;
  global len_distal; 

  global lb_angle_prox;
  global lb_angle_inter;
  global lb_angle_distal;

  global ub_angle_prox;
  global ub_angle_inter;
  global ub_angle_distal;

  global error_margin;

  % Visualization
  img = figure(2);
  clf;
  hold on;
  axis([0,100,-3,100],"equal");
  plot([0,65],[-2,-2]);
  
  %forwardKinematics(0.4708, lb_angle_inter/2, lb_angle_distal)
  [X1, Y1] = forwardKinematics(0.98425, (-2*pi)/3, 0)
  [X2, Y2] = forwardKinematics(0.4835, -1.15157, (-2*pi)/3)
  plot( X1, Y1);
  plot( X2, Y2);
  hold off;  
  
endfunction

function K = assignment3()
  global len_prox;
  global len_inter;
  global len_distal; 

  global lb_angle_prox;
  global lb_angle_inter;
  global lb_angle_distal;

  global ub_angle_prox;
  global ub_angle_inter;
  global ub_angle_distal;

  global error_margin;

  % Visualization
  img2 = figure(3);
  clf;
  hold on;
  axis([0,100,-100,100],"equal");
  
  for i=lb_angle_prox:0.1:ub_angle_prox
    [X1, Y1] = forwardKinematics(i, ub_angle_inter, ub_angle_distal);
    [X2, Y2] = forwardKinematics(i, lb_angle_inter, lb_angle_distal);
     
    plot( X1(4), Y1(4));
    plot( X2(4), Y2(4));
  endfor
  
  for i=lb_angle_inter:0.1:ub_angle_inter
    [X1, Y1] = forwardKinematics(lb_angle_prox, i, ub_angle_distal);
    [X2, Y2] = forwardKinematics(ub_angle_prox, i, lb_angle_distal);
     
    plot( X1(4), Y1(4));
    plot( X2(4), Y2(4));
  endfor
  
  for i=lb_angle_distal:0.1:ub_angle_distal
    [X1, Y1] = forwardKinematics(lb_angle_prox, lb_angle_inter, i);
    [X2, Y2] = forwardKinematics(ub_angle_prox, ub_angle_inter, i);
     
    plot( X1(4), Y1(4));
    plot( X2(4), Y2(4));
  endfor
  
  hold off;
endfunction

function [Xi,Yi] = inverseKinematics(currentAngleGuess, target)
  global len_prox;
  global len_inter;
  global len_distal; 

  global lb_angle_prox;
  global lb_angle_inter;
  global lb_angle_distal;

  global ub_angle_prox;
  global ub_angle_inter;
  global ub_angle_distal;

  global error_margin;
  
  Xi = [];
  Yi = [];

  i = 0;
  do
    i = i+1;
    distalJoint =  2*currentAngleGuess(2)/3;
    [X, Y] = forwardKinematics(currentAngleGuess(1), currentAngleGuess(2), distalJoint);
    Xi = vertcat(Xi,X);
    Yi = vertcat(Yi,Y);
    inv_J_q = calcInvJacobian(len_prox, len_inter, len_distal, currentAngleGuess(1), currentAngleGuess(2));
    currentAngleGuess = currentAngleGuess + inv_J_q*( target - [X(4);Y(4)] );
    error = sqrt( (target(1)-X(4))^2 + (target(2)-Y(4))^2 );
  until((error<=error_margin && Y(4)>-2) || i>=100);
  
  i = i+1;
  distalJoint =  2*currentAngleGuess(2)/3;
  [X, Y] = forwardKinematics(currentAngleGuess(1), currentAngleGuess(2), distalJoint);
  
  Xi = vertcat(Xi,X);
  Yi = vertcat(Yi,Y);
endfunction

function K = assignment5and6()
  global len_prox;
  global len_inter;
  global len_distal; 

  global lb_angle_prox;
  global lb_angle_inter;
  global lb_angle_distal;

  global ub_angle_prox;
  global ub_angle_inter;
  global ub_angle_distal;

  global error_margin;
  
  % Visualization
  img2 = figure(5);
  clf;
  hold on;
  axis([0,100,-3,100],"equal");
  plot([0,100],[-2,-2],'color','r');
  
  target = transpose([len_prox + len_inter + 10, -2]);
  plot(target(1), target(2),'marker','x','markersize',30,'color','r');
  
  [Xi, Yi] = inverseKinematics([0.1,-0.1],target);
  
  for i=1:rows(Xi)
    px = [Xi(i,1),Xi(i,2),Xi(i,3),Xi(i,4)];
    py = [Yi(i,1),Yi(i,2),Yi(i,3),Yi(i,4)];
    plot(px,py,'marker','o');
  endfor
  
  for i=1:rows(Xi)
    tipx(i) = Xi(i,4);
    tipy(i) = Yi(i,4);
  endfor
  plot(tipx,tipy,'color','r','linestyle','--');
  hold off;
  
  img3 = figure(6);
  clf;
  hold on;
  for i=1:rows(Xi)
    diff(i) = sqrt( (Xi(i,4)-target(1))^2 + (Yi(i,4)-target(2))^2 );
  endfor
  plot(diff);
  hold off;
endfunction

function assignment7()
  global len_prox;
  global len_inter;
  global len_distal;
  global lb_angle_prox;
  global lb_angle_inter;
  global lb_angle_distal;
  global ub_angle_prox;
  global ub_angle_inter;
  global ub_angle_distal;

  
  img7 = figure(7);
  clf;
  hold on;
  axis([0,100,-3,40],"equal");
  plot([0,100],[-2,-2],'color','r');
  
  dist  = len_prox+len_inter+len_distal;
  angle = asin(-2/dist);
  [Xf,Yf] = forwardKinematics(angle,0,0);
  plot(Xf,Yf,'marker','o');
  
  [Xt,Yt] = forwardKinematics(0,lb_angle_inter,(2/3)*lb_angle_inter);
  distq  = sqrt(Xt(4)^2+Yt(4)^2);
  angleq = asin(-2/distq);
  anglediffq = atan(Yt(4)/Xt(4));
  
  [Xc,Yc] = forwardKinematics(angleq-anglediffq,lb_angle_inter,(2/3)*lb_angle_inter);
  plot(Xc,Yc,'marker','o');
  
  hold off;
endfunction

function assignment8()
    global len_prox;
  global len_inter;
  global len_distal;
  global lb_angle_prox;
  global lb_angle_inter;
  global lb_angle_distal;
  global ub_angle_prox;
  global ub_angle_inter;
  global ub_angle_distal;

  img8 = figure(8);
  clf;
  hold on;
  axis([0,100,-3,40],"equal");
  plot([0,100],[-2,-2],'color','r');
  
  [X1,Y1] = inverseKinematics([ub_angle_prox,ub_angle_inter],[60;-2]);

  s = rows(X1);
  plot([X1(s,1),X1(s,2),X1(s,3),X1(s,4)],
       [Y1(s,1),Y1(s,2),Y1(s,3),Y1(s,4)],
       'marker','x');
  
%  prev = [TODO];
%  for i=50:10:80
%    [X,Y] = inverseKinematics(prev,[i,-2]);
%    z = rows(X);
%    for j=1:z
%      plot([X(j,1),X(j,2),X(j,3),X(j,4)],
%           [Y(j,1),Y(j,2),Y(j,3),Y(j,4)],
%           'marker','o');
%    endfor
%    prev = [X(z,4);Y(z,4)];
%  endfor
  
  hold off;
endfunction

%assignment2()
%assignment3()
%assignment5and6()
%assignment7()
assignment8()

%len_prox, len_inter, len_distal, lb_angle_prox, lb_angle_inter, lb_angle_distal, ub_angle_prox, ub_angle_inter, ub_angle_distal
%
%origin = transpose([0,0,0,1]);
%target = transpose([len_prox + len_inter + 10, -2]);
%plot(target(1), target(2),
%     'marker','x','markersize',40,'color','r');
%
%currentAngleGuess = [0.1,-0.1];
%
%i = 0;
%do
%  i = i+1;
%  T01 = rotateJoint(currentAngleGuess(1), len_prox,  lb_angle_prox,   ub_angle_prox);
%  T12 = rotateJoint(currentAngleGuess(2), len_inter, lb_angle_inter,  ub_angle_inter);
%  T23 = rotateJoint((currentAngleGuess(2)*2)/3, len_distal,lb_angle_distal, ub_angle_distal);
%  j1 = T01;
%  j2 = T01*T12;
%  tip= T01*T12*T23;
%  
%  tipsx(i) = tip(1,4);
%  tipsy(i) = tip(2,4);
%  
%  theta1 = clamp(currentAngleGuess(1), lb_angle_prox, ub_angle_prox);
%  theta2 = clamp(currentAngleGuess(2), lb_angle_inter, ub_angle_inter);
%  theta3 = (2*theta2)/3
%
%  a = tip(1,4)
%  b = tip(2,4)
%
%  plot([0,j1(1,4),j2(1,4),tip(1,4)],
%       [0,j1(2,4),j2(2,4),tip(2,4)],
%       'marker','o','color','k');
%       
%  inv_J_q = calcInvJacobian(len_prox, len_inter, len_distal, currentAngleGuess(1), currentAngleGuess(2));
%  currentAngleGuess = currentAngleGuess + inv_J_q*( target - [tip(1,4);tip(2,4)] );
%  
%  error = sqrt( (target(1)-tip(1,4))^2 + (target(2)-tip(2,4))^2 );
%until((error<error_margin && tip(2,4)>-2) || i>10);
%
%i = i+1;
%
%distalJoint =  2*currentAngleGuess(2)/3 ;
%T01 = rotateJoint(currentAngleGuess(1), len_prox,  lb_angle_prox,   ub_angle_prox);
%T12 = rotateJoint(currentAngleGuess(2), len_inter, lb_angle_inter,  ub_angle_inter);
%T23 = rotateJoint(distalJoint,          len_distal,lb_angle_distal, ub_angle_distal);
%j1 = T01;
%j2 = T01*T12;
%tip= T01*T12*T23;
%tipsy(i) = tip(2,4);
%tipsx(i) = tip(1,4);
%  plot([0,j1(1,4),j2(1,4),tip(1,4)],
%       [0,j1(2,4),j2(2,4),tip(2,4)],
%       'marker','o','color','b');
%
%
%plot(tipsx,tipsy,'linestyle','-.','color','r');
%
%l1 = 44.6;
%l2 = 26.3;
%l3 = 17.4;
%
%theta1 = clamp(currentAngleGuess(1), lb_angle_prox, ub_angle_prox);
%theta2 = clamp(currentAngleGuess(2), lb_angle_inter, ub_angle_inter);
%theta3 = (2*theta2)/3
%
%a = tip(1,4)
%b = tip(2,4)
%
%x = l1*cos(theta1) + l3*cos((2*theta2)/3)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)) - l3*sin((2*theta2)/3)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) + l2*cos(theta1)*cos(theta2) - l2*sin(theta1)*sin(theta2) 
%y = l1*sin(theta1) + l3*cos((2*theta2)/3)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) + l3*sin((2*theta2)/3)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)) + l2*cos(theta1)*sin(theta2) + l2*cos(theta2)*sin(theta1)
