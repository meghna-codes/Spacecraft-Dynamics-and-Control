function xdot = cr3bp_uncontrolled(x, mu)
% x = [x y z xd yd zd]

mu1 = 1-mu; % mass of larger  primary 
mu2 = mu;   % mass of smaller primary 
r1 = ((x(1)+mu2)^2 + x(2)^2 + x(3)^2)^(3/2); % r1: distance to m1, LARGER MASS    
r2 = ((x(1)-mu1)^2 + x(2)^2 + x(3)^2)^(3/2); % r2: distance to m2, smaller mass  

xdot = zeros(6,1);

xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);
xdot(4) = x(1)-(mu1*(x(1)+mu2)/r1) -(mu2*(x(1)-mu1)/r2) + 2*x(5);
xdot(5) = x(2)-(mu1* x(2)/r1) -(mu2*x(2)/r2) - 2*x(4);
xdot(6) =     -(mu1* x(3)/r1) - (mu2*x(3)/r2);

end
