function vx = get_variance(A)

    [r,c] = size(A);

    % Loop in the time dimension. Make sure the time dimension correspond
    % to the columns of the array A.
    for i=1:c
        vx(i) = var(A(:,i));
    end
