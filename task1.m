K = 1; T = 50; h = 0.1;
A = [0 1 0; 0 -1/T -K/T; 0 0 0];
B = [0; K/T; 0];
E = [0 0; 1 0; 0 1];
C = [1; 0; 0];

O = obsv(A,C');
n = rank(O);

[~, Bd] = c2d(A,B,h);
[Ad, Ed] = c2d(A,E,h);
Cd = C;