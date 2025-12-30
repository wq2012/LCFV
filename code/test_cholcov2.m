% Test suite for cholcov2
function test_cholcov2()
    test_positive_definite();
    test_positive_semi_definite();
    test_non_square();
    fprintf('All cholcov2 tests passed.\n');
end

function test_positive_definite()
    Sigma = [2, 1; 1, 2];
    [T, p] = cholcov2(Sigma);
    assert(p == 0, 'Should be positive definite');
    assert(norm(T'*T - Sigma) < 1e-10, 'Reconstruction error too high');
end

function test_positive_semi_definite()
    Sigma = [1, 1; 1, 1];
    [T, p] = cholcov2(Sigma);
    assert(p == 0, 'Should be positive semi-definite');
    assert(norm(T'*T - Sigma) < 1e-10, 'Reconstruction error too high');
end

function test_non_square()
    Sigma = [1, 2, 3; 4, 5, 6];
    try
        cholcov2(Sigma);
        error('Should have thrown error for non-square matrix');
    catch ME
        assert(strcmp(ME.message, 'Sigma is not square'), 'Wrong error message');
    end
end
