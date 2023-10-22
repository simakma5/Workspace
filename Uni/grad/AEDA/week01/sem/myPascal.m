function P = myPascal(N)
    if ~isnumeric(N)
        error("Invalid input: Not a number!")
    elseif N < 0
        error("Invalid input: Negative number!")
    else
        P = ones(N);
        for i = 2:N
            for j = 2:N
                P(i,j) = P(i-1,j) + P(i,j-1);
            end
        end
    end
end