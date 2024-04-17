function [result] = Verify(n, m, d, q, A, f, B, Lpk, theta, z, c, C, d0)
N = size(Lpk,1);
    % Compute di = d0 + i·θ, for all i∈[N].
    % Note that all di are stored in an N × (n × d) dimensional matrix.

    for i = 1:N
        di(i,1:n*d)=mod(d0+i*theta,q);
    end

for a=1:128
% Compute R0 = A·z 
    pkc = zeros(n,2*d-1);
    pkc0 = zeros(n,2*d-1);
    sum = zeros(n,2*d-1);
    Q = zeros(n,d-1);
    R = zeros(n,2*d-1);
    for j=1:n
        for i = 0:n
            sum(j,:) = sum(j,:) + conv(A(j,1+d*i:d+d*i),z(1,1+d*i:d+d*i));
        end
        sum(j,:) = sum(j,:) + [zeros(1,d-1),z(1, (1+d*j):(d+d*j))];

        for k = 1:N-1
                pkc0(j,:) = pkc0(j,:) + conv(Lpk(k+1, 1+d*(j-1):d+d*(j-1)),C(1,:));
        end
        sum(j,:) = sum(j,:) + pkc0(j,:);
        [Q(j,:),R(j,:)] = deconv(sum(j,:),f);
        t0(j,:) = mod(R(j,d:2*d-1),q);
    end


% Compute R1 = B·z
    dic = zeros(n,2*d-1);
    dic0 = zeros(n,2*d-1);
    sum = zeros(n,2*d-1);
    Q = zeros(n,d-1);
    R = zeros(n,2*d-1);
    for j=1:n
        for i = 0:n
            sum(j,:) = sum(j,:) + conv(B(j,1+d*i:d+d*i),z(1,1+d*i:d+d*i));
        end
        sum(j,:) = sum(j,:) + [zeros(1,d-1),z(1, (1+d*j):(d+d*j))];

        for k = 1:N-1
            dic0(j,:) = dic0(j,:) + conv(di(k+1,1+d*(j-1):d+d*(j-1)),C(1,:));
        end
        sum(j,:) = sum(j,:) + dic0(j,:);

        [Q(j,:),R(j,:)] = deconv(sum(j,:),f);
        t1(j,:) = mod(R(j,d:2*d-1),q);
    end

    t0_h = reshape(t0,1,d*n);
    t1_h = reshape(t1,1,d*n);

 Verify_c = hash ([Lpk; t0_h; t1_h; theta],'SHA-512');

    hash_c = [];
    for i=1:strlength(Verify_c)
        hash_c(1,i) = mod(Verify_c(1,i),3);
        if(hash_c(1,i) == 2)
            hash_c(1,i) = -1;
        end
    end

 for i=1:m*d
        if(z(1,i) > m*d^2-d && z(1,i) < -m*d^2-d)
            break
        end
 end


 if(all(c == hash_c))
        j = 1;
    else
        j = 0;
    end

    if  (i == m*d && j == 1)
        result = 1;
    else
        result = 0;
    end
end
end