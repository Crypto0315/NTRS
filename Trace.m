function [result] = Trace(n, d, q, N, theta1, theta2)

d0 = randi([-q,q], 1, d*n);
d1 = randi([-q,q], 1, d*n);

for i = 1:N
    di(i,1:n*d)=mod(d0+i*theta1,q);
end

for i = 1:N
    dj(i,1:n*d)=mod(d1+i*theta2,q);
end


for i=1:N
    if  (di == dj)
        result = 1;
    else
        result = 0;
    end
end
