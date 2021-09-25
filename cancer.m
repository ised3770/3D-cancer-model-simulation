function f = cancer(t, X)

a12 = 0.73;
a13 = 2.5;
a21 = 1.5;
a31 = 0.2;
r2 = 0.6;
r3 = 4.5;
k3 = 1;
d3 = 0.5;

x = X(1); y = X(2); z = X(3);
Y = [X(4), X(7), X(10);
    X(5), X(8), X(11);
    X(6), X(9), X(12)];
f = zeros(9,1);

% Model equations
f(1)=x - x*x - a12*x*y - a13*x*z;
f(2)=r2*y - r2*y*y - a21*x*y;
f(3)=(r3*x*z)/(x+k3) - a31*x*z - d3*z;

      
A = 1-2*x-a12*y-a13*z;
B = -a12*x;
C = -a13*x;
D = -a21*y;
E = r2-2*r2*y-a21*x;
F = (r3*z)/(x+k3)*(1-(x/(x+k3)))-a31*z;
G = (r3*x)/(x+k3)-a31*x-d3;

% Jacobian 
 Jac=[A,    B,    C;
      D,    E,    0;
      F,    0,    G];
  
% Variational form 
f(4:12) = Jac * Y;