function [SK, PK] = KeyGen(n, m, d, q, A, f)
% KeyGen algorithm
%   The private key SK is randomly sampled and stored in a 1×md dimensional matrix.
SK = randi([-1,1],1,m*d);

pk = zeros(n,2*d-1);

% Set a zero matrix to store the quotient Q.
Q = zeros(n,d-1);
% Set a zero matrix to store the remainder R.
R = zeros(n,2*d-1);

%   Compute PK = A*SK.
    for j = 1:n
        for i = 0:n
            pk(j,:) = pk(j,:) + conv(A(j,1+d*i:d+d*i),SK(1,1+d*i:d+d*i));
        end
        pk(j,:) = pk(j,:) + [zeros(1,d-1), SK(1,(1+d*j):(d+d*j))];
        %pk(j,:) = pk(j,:) + [zeros(1,d-1),SK(1,1+d*(j+n):d+d*(j+n))];
%   Next, let's take modulus of Zq[x]/x^d+1.
        [Q(j,:),R(j,:)] = deconv(pk(j,:),f);
        PK(1,1 + d*(j-1):d+d*(j-1)) = mod(R(j,d:2*d-1),q);
    end
%disp(pk);
%disp(R);
end