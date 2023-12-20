function output = mapGrayFunc(input)
    %input: 1xN binary vector
    %output: 1xN/2 complex vector

    intermediate = reshape(input,2,[]).';
    output = zeros(length(intermediate),1);
    output(all(intermediate == [0 1], 2)) = -1/sqrt(2) + 1/sqrt(2) *1i;
    output(all(intermediate == [1 1], 2)) = 1/sqrt(2) + 1/sqrt(2) *1i;
    output(all(intermediate == [1 0], 2)) = 1/sqrt(2) - 1/sqrt(2) *1i;
    output(all(intermediate == [0 0], 2)) = -1/sqrt(2) - 1/sqrt(2) *1i;
end