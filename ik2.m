1;
function [ inv_j ] = invJacobian( q )
%INVJACOBIAN Summary of this function goes here
%   Detailed explanation goes here

theta1 = q(1,1);
theta2 = q(2,1);

j = [-44.6*sin(theta1)-26.3*sin(theta1+theta2)-17.4*sin(theta1+(5/3)*theta2), -26.3*sin(theta1+theta2)-43.5*sin(theta1+(5/3)*theta2); 44.6*cos(theta1)+26.3*cos(theta1+theta2)+ 17.4*(cos(theta1+(5/3)*theta2)), 26.3*cos(theta1+theta2)+17.4*cos(theta1+(5/3)*theta2)];
D = j(1,1)* j(2,2) - j(1,2)*j(2,1);
inv_j = (1/D) * [j(2,2), -1*j(1,2); -1*j(2,1), j(1,1)];

%inv_j = (1/D) * [26.3*cos(theta1+theta2)+ 17.4*cos(theta1+(5/2)*theta2), 26.3*sin(theta1+theta2)+43.5*sin(theta1+(5/2)*theta2); -44.6*cos(theta1)-26.3*cos(theta1+theta2)-17.4*(cos(theta1+(5/2)*theta2)), -44.6*sin(theta1)- 26.3*sin(theta1+theta2)-17.4*sin(theta1+(5/2)*theta2)]


%D = (43.5*sin(theta1 + (5/2)*theta2) + 26.3*sin(theta1 + theta2))*(17.4*cos(theta1 + (5/2)*theta2) + 26.3*cos(theta1 + theta2) + 44.6*cos(theta1)) - (43.5*cos(theta1 + (5/2)*theta2) + 26.3*cos(theta1 + theta2))*(17.4*sin(theta1 + (5/2)*theta2) + 26.3*sin(theta1 + theta2) + 44.6*sin(theta1))

endfunction

function [ manhattan_dist ] = manhattanDist( v1, v2 )
%MANHATTANDIST Summary of this function goes here
%   Detailed explanation goes here

manhattan_dist = v1-v2;

endfunction


function [ manhattan_dist  ] = manhattanDistFingerTip( i_x, i_y, theta1, theta2 )
%MANHATTANDISTFINGERTIP Summary of this function goes here
%   Detailed explanation goes here

[x, y] = tipCoords(theta1, theta2);

manhattan_dist_x = i_x - x;
manhattan_dist_y = i_y - y;

manhattan_dist = [manhattan_dist_x ; manhattan_dist_y;]
endfunction

function [ x, y ] = tipCoords( theta1, theta2 )
%FINGERCOORDS Summary of this function goes here
%   Detailed explanation goes here

l1 = 44,6;
l2 = 26,3;
l3 = 17,4;

x = l1*cos(theta1) + l3*cos((2*theta2)/3)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)) - l3*sin((2*theta2)/3)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) + l2*cos(theta1)*cos(theta2) - l2*sin(theta1)*sin(theta2); 
y = l1*sin(theta1) + l3*cos((2*theta2)/3)*(cos(theta1)*sin(theta2) + cos(theta2)*sin(theta1)) + l3*sin((2*theta2)/3)*(cos(theta1)*cos(theta2) - sin(theta1)*sin(theta2)) + l2*cos(theta1)*sin(theta2) + l2*cos(theta2)*sin(theta1);
 

endfunction

function [ xCoords, yCoords ] = invFingerCoords( theta1, theta2 )
%FINGERCOORDS Summary of this function goes here
%   Detailed explanation goes here

l1 = 44.6;
l2 = 26.3;
l3 = 17.4;

T1 = [cos(theta1), -sin(theta1), 0 , l1 * cos(theta1); sin(theta1), cos(theta1), 0, l1* sin(theta1); 0,0,1,0;0,0,0,1;];
T2 = [cos(theta2), -sin(theta2), 0 , l2 * cos(theta2); sin(theta2), cos(theta2), 0, l2* sin(theta2); 0,0,1,0;0,0,0,1;];
T3 = [cos((2*theta2)/3), -sin((2*theta2)/3), 0 , l3 * cos((2*theta2)/3); sin((2*theta2)/3), cos((2*theta2)/3), 0, l3* sin((2*theta2)/3); 0,0,1,0;0,0,0,1;];

trans_j1 = eye(4);
trans_j2 = T1;
trans_j3 = T1*T2;
trans_tip = T1 * T2 * T3;


x_j1 = 0;
x_j2 = trans_j2(1,4);
x_j3 = trans_j3(1,4);
x_tip = trans_tip(1,4);

y_j1 = 0;
y_j2 = trans_j2(2,4);
y_j3 = trans_j3(2,4);
y_tip = trans_tip(2,4);

xCoords = [x_j1,x_j2,x_j3,x_tip];
yCoords = [y_j1,y_j2,y_j3,y_tip];

endfunction


function [ r_theta1, r_theta2 ] = InverseKinematicsSolver( s_theta1, s_theta2, ix, iy )
%INVERSEKINEMATICSSOLVER Summary of this function goes here
%   Detailed explanation goes here

% Initial target was [x,y] = fingerCoords(sin(-2/88.3),0,0)


% [x,y] = orgFingerCoords(q1(1,1), q1(2,1), (2*q1(2,1))/3);
% scatter(x,y);
% plot(x,y);
% hold off;

  destinationPoints = [   85, 80,   75,    70,    65,    60,    55,    50,    45,    40,   35,  30;
                       -1.95, -1.95, -1.95, -1.95, -1.95, -1.95, -1.95, -1.95, -1.95, -1.95, -1.95, -1.95;   ]

%  destinationPoints = [ 86,    70,    60,    50,    40,    30;
%                        -2, -2, -2, -2, -2, -2;   ]
                    
  lastQ = [s_theta1 ; s_theta2;];

  allIterations = [];
    
  for destination = destinationPoints
      allTheQ = [];
      allTheQ = horzcat(allTheQ, lastQ);

      destX = destination(1,1);
      destY = destination(2,1);
      
      %figure
      hold on
      scatter(destX, destY, 'filled')

      nextQ = lastQ + (invJacobian(lastQ)*(manhattanDistFingerTip(destX, destY, lastQ(1,1), lastQ(2,1) )));

      manhattanDist = manhattanDistFingerTip(destX, destY, nextQ(1,1), nextQ(2,1) );
      manhattan_dist_x = manhattanDist(1,1);
      manhattan_dist_y = manhattanDist(2,1);

      lastQ = nextQ;
      [lastX, lastY] = tipCoords(lastQ(1,1), lastQ(2,1));
      
      iterations = 1;
      allTheQ = horzcat(allTheQ, lastQ);
      while(manhattan_dist_y > 0.05 || manhattan_dist_y < -0.05 || manhattan_dist_x > 0.05 || manhattan_dist_x < -0.05)
          manhattanDist = manhattanDistFingerTip(destX, destY, lastQ(1,1), lastQ(2,1) );
          manhattan_dist_x = manhattanDist(1,1);
          manhattan_dist_y = manhattanDist(2,1);

          nextQ = lastQ + (invJacobian(lastQ)*(manhattanDistFingerTip(destX, destY, lastQ(1,1), lastQ(2,1) )));
          lastQ = nextQ;
          [lastX, lastY] = tipCoords(lastQ(1,1), lastQ(2,1));
          iterations = iterations + 1;
          allTheQ = horzcat(allTheQ, nextQ);
      endwhile

      allIterations = horzcat(allIterations, iterations);
      iterations = iterations
      allTheQ = allTheQ

      for index = allTheQ
          [x,y] = invFingerCoords(index(1,1), index(2,1));
          hold on;
          if(y >= -2.3)
            plot(x,y);
            scatter(x,y);
          endif
      endfor
  endfor
  plot([0,100],[-2,-2]);
  allIterations = allIterations;
endfunction

[ r_theta1, r_theta2 ] = InverseKinematicsSolver( -0.1, -0.1, 0, 0 )