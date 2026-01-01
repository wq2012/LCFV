function test_sift()
%TEST_SIFT Unit tests for Dense SIFT extraction

    addpath(genpath(pwd));
    
    test_sift_output_dimension();
    test_sift_rgb_handling();
    test_sift_zero_image();
    test_sift_simple_pattern();
    
    fprintf('All SIFT tests passed.\n');
end

function test_sift_output_dimension()
    fprintf('Testing SIFT output dimension...\n');
    % Create random image 64x64
    img = rand(64, 64);
    
    step_size = 4;
    bin_size = 4;
    
    descs = compute_dense_sift(img, step_size, bin_size);
    
    [D, N] = size(descs);
    assert(D == 128, 'Descriptor dimension should be 128');
    assert(N > 0, 'Should return some descriptors');
    
    % Check N count
    patch_size = 4 * bin_size; % 16
    grid_x = patch_size/2:step_size:64-patch_size/2;
    grid_y = patch_size/2:step_size:64-patch_size/2;
    expected_N = length(grid_x) * length(grid_y);
    
    assert(N == expected_N, sprintf('Expected %d descriptors, got %d', expected_N, N));
end

function test_sift_rgb_handling()
    fprintf('Testing SIFT RGB handling...\n');
    img = rand(32, 32, 3);
    descs = compute_dense_sift(img);
    assert(size(descs, 1) == 128, 'RGB input should produce 128D descriptors');
end

function test_sift_zero_image()
    fprintf('Testing SIFT on zero image...\n');
    img = zeros(32, 32);
    descs = compute_dense_sift(img);
    assert(all(descs(:) == 0), 'Zero image should produce zero descriptors');
end

function test_sift_simple_pattern()
    fprintf('Testing SIFT on simple pattern...\n');
    % Vertical stripes
    img = zeros(32, 32);
    img(:, 1:2:end) = 1;
    
    descs = compute_dense_sift(img, 4, 4);
    
    % Should have some non-zero gradients
    assert(sum(descs(:)) > 0, 'Vertical stripes should produce non-zero descriptors');
    
    % Check normalization (norm should be 1 or 0)
    col_norms = sqrt(sum(descs.^2, 1));
    % We accept 0 norms (features in flat regions) or ~1 norms
    is_valid = (abs(col_norms - 1) < 1e-4) | (col_norms < 1e-4);
    assert(all(is_valid), 'Descriptors should be normalized');
end
