function v = fvecs_read(filename)
%FVECS_READ Read vector files in .fvecs format
%
%   v = FVECS_READ(filename)
%
%   Input:
%       filename - Path to .fvecs file
%
%   Output:
%       v        - D x N matrix
%
%   Copyright (C) 2013 Quan Wang <wangq10@rpi.edu>,
%   Signal Analysis and Machine Perception Laboratory,
%   Department of Electrical, Computer, and Systems Engineering,
%   Rensselaer Polytechnic Institute, Troy, NY 12180, USA
%
%   You are free to use this software for academic purposes if you cite our paper:
%   Quan Wang, Xin Shen, Meng Wang, Kim L. Boyer,
%   Label Consistent Fisher Vectors for Supervised Feature Aggregation,
%   22nd International Conference on Pattern Recognition (ICPR), 2014.
%
%   For commercial use, please contact the authors.

    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    
    % Read dimension of first vector
    d = fread(fid, 1, 'int32');
    
    if isempty(d)
        v = [];
        fclose(fid);
        return;
    end
    
    % Rewind
    fseek(fid, 0, 'bof');
    
    % Each vector has 4 bytes for dim + d * 4 bytes for floats = (d+1)*4 bytes
    fseek(fid, 0, 'eof');
    file_size = ftell(fid);
    fseek(fid, 0, 'bof');
    
    num_vecs = file_size / (4 * (d + 1));
    
    if mod(num_vecs, 1) ~= 0
        fclose(fid);
        error('File size does not match expected format');
    end
    
    % Read all data
    % Data is stored as (dim, v1, v2, ..., vd) repeated
    % We read as integers first to handle the structure
    
    v = zeros(d, num_vecs);
    
    for i = 1:num_vecs
        dim = fread(fid, 1, 'int32');
        if dim ~= d
            fclose(fid);
            error('Dimension mismatch at vector %d', i);
        end
        v(:, i) = fread(fid, d, 'float32');
    end
    
    fclose(fid);
end
