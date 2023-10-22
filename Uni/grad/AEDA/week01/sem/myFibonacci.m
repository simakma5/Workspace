function Fn = myFibonacci(N)
    if ~isnumeric(N)
        error("Invalid input: Not a number!")
    elseif N < 2
        if N < 0
            error("Invalid input: Negative number!")
        else
            Fn = N;
        end
    else
        Fn_2 = 0;
        Fn_1 = 1;
        for n = 2:N
            Fn = Fn_1 + Fn_2;
            Fn_2 = Fn_1;
            Fn_1 = Fn;
        end
    end
end