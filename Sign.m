function [B, theta, z, c, C, d0] = Sign(n, m, d, q, A, f, N, Lpk, PK, SK)
  % Generate matrix B.
  % Here we have simulated the generation of matrix B using the random
  % function randi()，it is generated in a similar manner to A (Setup).
    B = randi([-(q-1)/2,(q-1)/2], n, d*m);
   % I1 = zeros(1,d);
   % I1(1,d) = 1;
   % I0 = zeros(1,d);
   % I11 = [I1;I0;I0;I0;I0;I0;I0];
   % I12 = [I0;I1;I0;I0;I0;I0;I0];
   % I13 = [I0;I0;I1;I0;I0;I0;I0];
   % I14 = [I0;I0;I0;I1;I0;I0;I0];
   % I15 = [I0;I0;I0;I0;I1;I0;I0];
   % I16 = [I0;I0;I0;I0;I0;I1;I0];
   % I17 = [I0;I0;I0;I0;I0;I0;I1];
   % B = [B0,I11,I12,I13,I14,I15,I16,I17];
    % matrix B=[B0||In].

 %  Compute dπ = B*SK.
    pd = zeros(n,2*d-1);
    Q = zeros(n,d-1);
    R = zeros(n,2*d-1);
    
    for j = 1 :n
        for i = 0:n
            pd(j,:) = pd(j,:) + conv(B(j,1+d*i:d+d*i),SK(1,1+d*i:d+d*i));
        end
    
        pd(j,:) = pd(j,:) + [zeros(1,d-1),SK(1, (1+d*j):(d+d*j))];
        %   Next, let's take modulus of Zq[x]/x^d+1.
        [Q(j,:),R(j,:)] = deconv(pd(j,:),f);
        dpai(1,1 + d*(j-1):d+d*(j-1)) = mod(R(j,d:2*d-1),q);
        %   The generated dpai is stored in a 1 × (d×n)-dimensional matrix.
    end

%   Compute d0 = H2(T,μ)
%   To get the image in the specified range, we use the random function
%   instead of the H2 hash function.
d0 = randi([-q,q], 1, d*n);
theta = mod ((dpai - d0)*1,q);

%   Compute the rest of di,i∈[N]\{π}.
for i = 2:N
    di(i,1:d*n) = mod(d0+i*theta,q);
end


r = randi([-m*d^2,m*d^2], 1, m*d);
C = randi([-1,1],1,d);

for a=1:128

 % Compute t0 = Ar + Σpk·c
        % The comments in this part are similar to the KeyGen algorithm.
        pkc = zeros(n,2*d-1);
        pkc0 = zeros(n,2*d-1);
        sum = zeros(n,2*d-1);
        Q = zeros(n,d-1);
        R = zeros(n,2*d-1);
    for g=1:1
        for j=1:n
            for i = 0:n
                sum(j,:) = sum(j,:) + conv(A(j,1+d*i:d+d*i),r(1,1+d*i:d+d*i));
            end
            sum(j,:) = sum(j,:) + [zeros(1,d-1),r(1, (1+d*j):(d+d*j))];
            
            pkc(j,:) = pkc(j,:) + conv(PK(1,1+d*(j-1):d+d*(j-1)),C(1,:));
          
            for k = 1:N-1
                pkc0(j,:) = pkc0(j,:) + conv(Lpk(k+1, 1+d*(j-1):d+d*(j-1)),C(1,:));
            end

            sum(j,:) = sum(j,:) + pkc(j,:) + pkc0(j,:);
            [Q(j,:),R(j,:)] = deconv(sum(j,:),f);
            % Store t0 into a n × d dimensional matrix.
            t0(j,:) = mod(R(j,d:2*d-1),q);
        end

        
 % Compute t1 = Br + Σdi·c
        dic = zeros(n,2*d-1);
        dic0 = zeros(n,2*d-1);
        sum = zeros(n,2*d-1);
        Q = zeros(n,d-1);
        R = zeros(n,2*d-1);
        for j=1:n
            for i = 0:n
                sum(j,:) = sum(j,:) +conv(B(j,1+d*i:d+d*i),r(1,1+d*i:d+d*i));
            end
            sum(j,:) = sum(j,:) + [zeros(1,d-1),r(1, (1+d*j):(d+d*j))];
            
            dic(j,:) = dic(j,:) + conv(dpai(1,1+d*(j-1):d+d*(j-1)),C(1,:));
          
            for k = 1:N-1
                dic0(j,:) = dic0(j,:) + conv(di(k+1,1+d*(j-1):d+d*(j-1)),C(1,:));
            end
            sum(j,:) = sum(j,:) + dic(j,:) + dic0(j,:);
            [Q(j,:),R(j,:)] = deconv(sum(j,:),f);
            % store t1 into a n × d dimensional matrix.
            t1(j,:) = mod(R(j,d:2*d-1),q);
        end
    end
    
t0_h = reshape(t0,1,d*n);
t1_h = reshape(t1,1,d*n);

cc = hash ([Lpk; t0_h; t1_h; theta],'SHA-512');
        % The result of the hash output (character vector) is transformed
        % into a 1*d vector c composed of its ASCII code, which needs to
        % be assigned to c after ASCII mod 3 first, and then change 2 to -1.
        c = [];
        for i=1:strlength(cc)
            c(1,i) = mod(cc(1,i),3);
            if(c(1,i) == 2)
                c(1,i) = -1;
            end
        end

        sc1 = zeros(m,2*d-1);
        Q = zeros(m,d-1);
        R = zeros(m,2*d-1);
        % Rejection Sampling z = sc + r.
        for i=1:m
            sc1(i,:) = conv(C,SK(1,1+d*(i-1):d+d*(i-1)));
            %modx^d+1
            [Q(i,:),R(i,:)] = deconv(sc1(i,:),f);
            z(1,1+d*(i-1):d+d*(i-1))= R(i,d:2*d-1)+r(1,1+d*(i-1):d+d*(i-1));
        end

end
end