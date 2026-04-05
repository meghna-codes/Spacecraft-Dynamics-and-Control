function [A,B] = cr3bp_jacobian(x, mu)
% x = [x y z xd yd zd]'
% mu = m2 / (m1 + m2)

mu1 = 1 - mu;
mu2 = mu;

x1 = x(1); y1 = x(2); z1 = x(3);

r1 = (x1 + mu2)^2 + y1^2 + z1^2;
r2 = (x1 - mu1)^2 + y1^2 + z1^2;

r13 = r1^(3/2); r15 = r1^(5/2);
r23 = r2^(3/2); r25 = r2^(5/2);

Uxx = 1 - mu1*(1/r13 - 3*(x1+mu2)^2/r15) - mu2*(1/r23 - 3*(x1-mu1)^2/r25);

Uyy = 1 - mu1*(1/r13 - 3*y1^2/r15) - mu2*(1/r23 - 3*y1^2/r25);

Uzz =   - mu1*(1/r13 - 3*z1^2/r15) - mu2*(1/r23 - 3*z1^2/r25);

Uxy = 3*y1*( mu1*(x1+mu2)/r15 + mu2*(x1-mu1)/r25 );
Uxz = 3*z1*( mu1*(x1+mu2)/r15 + mu2*(x1-mu1)/r25 );
Uyz = 3*y1*z1*( mu1/r15 + mu2/r25 );

Urr = [ Uxx  Uxy  Uxz
        Uxy  Uyy  Uyz
        Uxz  Uyz  Uzz ];

A = zeros(6,6);
A(1:3,4:6) = eye(3);
A(4:6,1:3) = Urr;
A(4:6,4:6) = [0  2  0
             -2  0  0
              0  0  0];

B = [zeros(3,3); eye(3)];
end
