function data = leg_pol_rec(u, order)

%Recursevely compute legendre polinomial

for n = 0:order
    for m = 0:n

        if n == 0 && m == 0
            data = [n, m, 1];
        else 
            if n == 1 && m == 0
                data = [data; n m sqrt(3)*u];
            else 
                if n == 1 && m == 1
                    data = [data; n m sqrt(3)*sqrt(1-u^2)];
                else
                    if m == n && n > 1
                        i = find(data(:, 1) == n-1 & data(:, 2) == n-1);

                        data = [data; n m sqrt((2*n+1)/(2*n))*sqrt(1-u^2)*data(i, 3)];
                    else 
                        if m == n-1
                            i = find(data(:, 1) == n-1 & data(:, 2) == n-1);
                            data = [data; n m sqrt(2*n+1)*u*data(i, 3)];
                        else
                            if m < (n-1)
                                chi1 = sqrt((2*n+1)*(2*n-1)/((n-m)*(n+m)));
                                chi2 = sqrt((2*(n-1)+1)*(2*(n-1)-1)/(((n-1)-m)*((n-1)+m)));
                                
                                i = find(data(:, 1) == n-1 & data(:, 2) == m);
                                j = find(data(:, 1) == n-2 & data(:, 2) == m);

                                data = [data; n m chi1*u*data(i,3) - chi1/(chi2)*data(j, 3)];
                            end
                        end
                    end     
                end
            end
        end
    end
end
