% Test suite for solve_LCFV2
function test_solve_LCFV2()
    test_synthetic_data();
    fprintf('All solve_LCFV2 tests passed.\n');
end

function test_synthetic_data()
    m = 5;
    n = 10;
    G = rand(m, n);
    C = rand(n, n);
    C = C * C'; % Make symmetric
    alpha = 1.0;

    M = solve_LCFV2(G, C, alpha);

    % Expected size of M is depends on rank of C which is random but likely full rank?
    % The output M is [eye(m); E] where E comes from cholcov(C*alpha) / G
    % cholcov(C) returns A where A'*A approx C.
    % If C is n x n, A is k x n.
    % E is A / G.
    % M = [eye(m); E]. E will have k rows.
    % So M has m+k rows.

    assert(size(M, 2) == m, 'M should have m columns');
    assert(size(M, 1) >= m, 'M should have at least m rows');
end
