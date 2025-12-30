% Test suite for solve_LCFV1
function test_solve_LCFV1()
    test_synthetic_data();
    fprintf('All solve_LCFV1 tests passed.\n');
end

function test_synthetic_data()
    m = 5;
    n = 10;
    G = rand(m, n);
    C = rand(n, n);
    C = C * C'; % Make symmetric
    alpha = 1.0;

    [M, W] = solve_LCFV1(G, C, alpha);

    assert(all(size(M) == [m, m]), 'M should be m x m');
    assert(all(size(W) == [m, m]), 'W should be m x m');
end
