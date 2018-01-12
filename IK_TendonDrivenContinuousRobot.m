function redundancyDemo
%% Human Neuromechanical Control: Redundancy demo.
clear all;
close all;
clc;

%% Forward kinematics.
l = [0.31, 0.34, 0.23];
%  Redoing the forward kinematics for the new arm.
thetai = [0;130;90];
thetai = thetai*pi/180;
xi = l(1)*cos(thetai(1)) + l(2)*cos(sum(thetai(1:2))) + l(3)*cos(sum(thetai(1:3)));
yi = l(1)*sin(thetai(1)) + l(2)*sin(sum(thetai(1:2))) + l(3)*sin(sum(thetai(1:3)));

theta = [0;75;0];
theta = theta*pi/180;
xf = l(1)*cos(theta(1)) + l(2)*cos(sum(theta(1:2))) + l(3)*cos(sum(theta(1:3)));
yf = l(1)*sin(theta(1)) + l(2)*sin(sum(theta(1:2))) + l(3)*sin(sum(theta(1:3)));

fprintf('xi = %f, xf = %f, xf-xi = %f\n', xi, xf, xf-xi);
fprintf('yi = %f, yf = %f, yf-yi = %f\n', yi, yf, yf-yi);

% Trajectory.
dt = 0.001; t = 0:dt:2; T= 2;
X = [(xf - xi)*polyval([6,-15,10,0,0,0],t/T) + xi; ...
     (yf - yi)*polyval([6,-15,10,0,0,0],t/T) + yi];
X_dot = [(1/T)*(xf - xi)*polyval([30,-60,30,0,0],t/T); ...
         (1/T)*(yf - yi)*polyval([30,-60,30,0,0],t/T)];
     
%% Inverse Kinematics Solution
q_dot = zeros(3,size(X_dot,1));
q = zeros(3,size(X_dot,1));
q(:,1) = thetai;
figure;
X_gen = [X(:,1),zeros(2,size(X,2)-1)];
% Starting arm pose.
armxS = [0, l(1)*cos(q(1,1)), l(1)*cos(q(1,1)) + l(2)*cos(sum(q(1:2,1))), l(1)*cos(q(1,1)) + l(2)*cos(sum(q(1:2,1))) + l(3)*cos(sum(q(1:3,1)))];
armyS = [0, l(1)*sin(q(1,1)), l(1)*sin(q(1,1)) + l(2)*sin(sum(q(1:2,1))), l(1)*sin(q(1,1)) + l(2)*sin(sum(q(1:2,1))) + l(3)*sin(sum(q(1:3,1)))];
armxM(1,:) = armxS;
armyM(1,:) = armyS;
for i = 2:length(t)
    q_dot(:,i) = transpose(J(q(:,i-1),l))*inv(J(q(:,i-1),l)*transpose(J(q(:,i-1),l)))*X_dot(:,i);    
    q(:,i) = q(:,i-1) + q_dot(:,i)*dt;
    X_gen(:,i) = X_gen(:,i-1) + X_dot(:,i)*dt;
    
    clf;
    armxM(i,:) = [0, l(1)*cos(q(1,i)), l(1)*cos(q(1,i)) + l(2)*cos(sum(q(1:2,i))), l(1)*cos(q(1,i)) + l(2)*cos(sum(q(1:2,i))) + l(3)*cos(sum(q(1:3,i)))];
    armyM(i,:) = [0, l(1)*sin(q(1,i)), l(1)*sin(q(1,i)) + l(2)*sin(sum(q(1:2,i))), l(1)*sin(q(1,i)) + l(2)*sin(sum(q(1:2,i))) + l(3)*sin(sum(q(1:3,i)))];
    plot(armxM(i,:),armyM(i,:),'b'); hold on;
    plot(armxS,armyS,'g');
    plot(X(1,:),X(2,:),'color',[1 0.5 0.5]);
    plot(X_gen(1,1:i),X_gen(2,1:i),'color',[0.5 0.5 1]);
    plot(X(1,i),X(2,i),'ko');
    title('Dealing with redundacny through speed minimization');
    axis([-0.25, 0.75, -0.25, 0.75]);
    axis square;
    drawnow;
end;

save('speed_min','armxM','armyM','armxS','armyS');

figure;
subplot(211);
plot(t',X');
subplot(212);
plot(t',(180/pi)*q');

return;

function out = J(in,l)
% Jacobian.

out(1,1) = - l(1)*sin(in(1)) - l(2)*sin(sum(in(1:2))) - l(3)*sin(sum(in(1:3)));
out(1,2) = - l(2)*sin(sum(in(1:2))) - l(3)*sin(sum(in(1:3)));
out(1,3) = - l(3)*sin(sum(in(1:3)));
out(2,1) = l(1)*cos(in(1)) + l(2)*cos(sum(in(1:2))) + l(3)*cos(sum(in(1:3)));
out(2,2) = l(2)*cos(sum(in(1:2))) + l(3)*cos(sum(in(1:3)));
out(2,3) = l(3)*cos(sum(in(1:3)));

return;