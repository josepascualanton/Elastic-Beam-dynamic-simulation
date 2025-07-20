function [x, iter] = Gauss_Seidel(A, b);

    n = length(b);
    x = rand(n,1);
    normVal = 1E10;
    tol = 1E-5;
    iter = 0;
    itmax = 250;

    while normVal>tol & iter<itmax
        x_init = x;

        for i = 1:n 

            suma = A(i,1:i-1)*x(1:i-1) + A(i,i+1:n)*x_init(i+1:n);
            x(i) = (b(i) - suma)/A(i, i);
        end
        iter = iter + 1;
        normVal = norm(x_init - x);

    end

end