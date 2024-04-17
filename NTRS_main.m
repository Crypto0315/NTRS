clc; clear;
% Setup Algorithm outputs the public parameters PP. 

[n, m, d, q, A, f] = Setup;


% Algorithm KeyGen(PP) outputs (SK,PK).

tic;
% Timing tool: Start.
% Considering the instability of a single execution of KeyGen(PP) Algorithm by a personal computer,
% we repeat the execution of KeyGen(PP) Algorithm 1000 times.
for p = 1:1000
    [SK, PK] = KeyGen(n, m, d, q, A, f);
end
fprintf('The computational overhead of KeyGen(PP) Algorithm running 1000 times is %f sec.\n',toc);
% toc;
% Timing tool: End.


% In this module, the signer enters the ring size and the message to be signed.

% Let the public key set of the ring be Lpk, and to facilitate the implementation of
% the experiment, we let the index π of the signer be 1. For the public keys of other 
% ring members we can get them by random sampling in some specified set,
% which is  
% due to the determination M-LWE problem.

N = input("Please enter the size of the ring:");
% Lpk = randi([-q,q], N-1, d*n);
Lpk = [PK;randi([-q,q], N-1, d*n)];


% Sign(PP, T, miu, SK) outputs σ = (θ, z, C).

% Note that d0 and B are not part of the ring signature in the output of the Sign
% function below. d0 and B are public and only facilitate the invocation of the 
% Verify algorithm.
tic;
% Timing tool: Start.
for p=1:1000
    [B, theta, z, c, C, d0] = Sign(n, m, d, q, A, f, N, Lpk, PK, SK);
end
% Timing tool: End.
fprintf('The computational overhead of Sign Algorithm running 1000 times is %f sec.\n',toc);


% Verify(PP, T, μ, σ) outputs 0/1.

tic;
% Timing tool: Start.
for v =1:1000
    [result] = Verify(n, m, d, q, A, f, B, Lpk, theta, z, c, C, d0);
end
% Print verification result.
if (result)
    fprintf ('Ring signature is valid! Output 1\n' );
else
    fprintf ('Ring signature is invalid! Output 0\n' );
end
fprintf('The computational overhead of Verify Algorithm running 1000 times is %f sec.\n',toc);
% toc;
% Timing tool: End.




% Trace algorithm.
tic;
% Timing tool: Start.
theta1 = randi([-q,q], 1, d*n);
theta2 = randi([-q,q], 1, d*n);
for p=1:1000
    [result] = Trace(n, d, q, N, theta1, theta2);
end
% Timing tool: End.
if (result)
    fprintf ('The signer has been tracked! Output 1\n' );
else
    fprintf ('The signer has not been traced! Output 0\n' );
end
fprintf('The computational overhead of Trace Algorithm running 1000 times is %f sec.\n',toc);


fprintf("\n");

