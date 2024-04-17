function [n, m, d, q, A, f] = Setup
% Setup
% Set the parameters of M-LWE and M-SIS.
n = 4;
m = 6;
d = 128;
q = 8380417;

% Randomly sample a matrix A=[A0||In] (Note: Abar = A)
% Pick A0, we store the 4×2-dimensional polynomial vector in a n*(m-n)d matrix.
A = randi([-(q-1)/2,(q-1)/2], n, d*m);
%I1 = zeros(1,d);
%I1(1,d) = 1;
%I0 = zeros(1,d);
%I11 = [I1;I0;I0;I0;I0;I0;I0];
%I12 = [I0;I1;I0;I0;I0;I0;I0];
%I13 = [I0;I0;I1;I0;I0;I0;I0];
%I14 = [I0;I0;I0;I1;I0;I0;I0];
%I15 = [I0;I0;I0;I0;I1;I0;I0];
%I16 = [I0;I0;I0;I0;I0;I1;I0];
%I17 = [I0;I0;I0;I0;I0;I0;I1];
%A = [A0,I11,I12,I13,I14,I15,I16,I17];  %matrix A=[A0||In]

% Constructing an irreducible polynomial f = x^d+1.
f = zeros(1,d+1);
f(1,1) = 1;
f(1,d+1) = 1;
end