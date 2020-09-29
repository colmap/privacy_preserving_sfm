% Copyright (c) 2020, ETH Zurich.
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%
%     * Neither the name of ETH Zurich nor the names of its contributors may
%       be used to endorse or promote products derivedfrom this software
%       without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
% Authors: Johannes L. Schoenberger (jsch-at-demuc-dot-de)
%          Viktor Larsson (viktor.larsson@inf.ethz.ch)
%          Marcel Geppert (marcel.geppert@inf.ethz.ch)

function generate_gravity_file(model_path)
% For a given model, computes a consistent gravity direction for all
% cameras and stores it in a new file 'gravity.txt' along with the model
% files.
%
% Input arguments:
% model_path    The path to the sparse model folder, this will be handed to
% read_model

[cameras, images, points3D] = read_model(model_path);

gravity = extract_gravity(images);

% Copy the data into a matrix to make the file writing simpler
gravity_keys = gravity.keys;
gravity_mat = zeros(size(gravity_keys,1),4);
for i = 1:numel(gravity_keys)
    key = gravity_keys{i};
    val = gravity(key);
    gravity_mat(i,:) = [double(key), gravity(key)'];
end
file = fopen(fullfile(model_path, 'gravity.txt'), 'w+');

assert(images.Count == size(gravity_mat,1));

image_ids = images.keys();

for i = 1:numel(image_ids)
    image_id = image_ids{i};
    gravity_id = gravity_keys{i};
    assert(gravity_id == image_id)
    fprintf(file, '%u ', gravity_mat(i, 1));
    fprintf(file, '%s ', images(image_id).name);
    fprintf(file, '%.6f %.6f %.6f\n', gravity_mat(i, 2:4));
end

% fprintf(file, '%u %.6f %.6f %.6f\n', gravity_mat');
fclose(file);