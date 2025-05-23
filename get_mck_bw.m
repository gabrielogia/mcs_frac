function [M, C, K] = get_mck_bw(m, c, k, a, ndof,y0)
    M = zeros(2 * ndof, 2 * ndof);
    C = eye(2 * ndof);
    K = zeros(2 * ndof, 2 *ndof);
 
    for i=1:ndof
        C(i,i) = c(i);
        K(i,i) = a(i)*k(i);
        K(i, i + ndof) = (1-a(i)) * y0 * k(i);

        if i < ndof
            C(i,i+1) = -c(i+1);
            K(i, i + 1) = -a(i+1)*k(i+1);
            K(i, i + ndof + 1) = -(1 - a(i+1)) * y0 * k(i+1); 
        end

        for j=1:i
            M(i,j) = m(i);
        end
    end
end