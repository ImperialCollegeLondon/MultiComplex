function A = arr4matmulti(p,q) %both p and q need to be in the same size
x = p.zn;
y = q.zn;
    for i = 1: length(x)
        for j = 1 : length(x)
            b = j;
            a = i;
            n = length(x);
            while length(x) > 1
                if a <= n/2 & b <= n/2 % upper left
                    x = x(1:n/2);
                    y = y(1:n/2);
                elseif a <= n/2 & b > n/2 % upper right
                    x = -x(n/2 +1 :end);
                    y = -y(n/2 +1 :end);
                    b = b - n/2;
                elseif a > n/2 & b <= n/2 % lower left
                    x = x(n/2 +1 :end);
                    y = y(n/2 +1 :end);
                    a = a - n/2;
                elseif a > n/2 & b > n/2 % lower right
                    x = x(1:n/2);
                    y = y(1:n/2);
                    a = a - n/2;
                    b = b - n/2;
                end
                n = n/2;
            end
            xx(j)=x; % up to this point
            yy(j)=y;
        end 
        z(i) = dot(xx,yy);
    end
A = z;
end
