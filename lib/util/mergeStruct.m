    function A = mergeStruct(A, B)
        f = fieldnames(B);
        for i = 1:numel(f), A.(f{i}) = B.(f{i}); end
    end
