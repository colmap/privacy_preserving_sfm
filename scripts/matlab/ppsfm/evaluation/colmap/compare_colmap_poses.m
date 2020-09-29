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

function diffs = compare_colmap_poses(model_path1, model_path2)
% Computes the pose differences between two colmap models.
% output: <image_id> <rotation_diff (axang)> <translation_diff>

[~, images1, ~] = read_model(model_path1);
[~, images2, ~] = read_model(model_path2);

image_ids1 = cell2mat(images1.keys());
image_ids2 = cell2mat(images2.keys());

all_im_ids = unique([image_ids1 image_ids2]);
im_ids_both = intersect(image_ids1, image_ids2);
im_ids_1_not_2 = setdiff(all_im_ids, image_ids2);
im_ids_2_not_1 = setdiff(all_im_ids, image_ids1);

% We assume that model 1 is the reference and we want to check the error of
% model 2. We compute the errors for images that are both contained, set
% the errors of images that are missing in model 1 to 0 and the ones of
% images missing in model 2 to inf.
% We store the differences as [image_id axang_diff pos_diff]
diffs = zeros(numel([im_ids_both im_ids_1_not_2 im_ids_2_not_1]), 8);

diff_idx = 1;
for image_id = im_ids_both
    R_diff = images2(image_id).R * images1(image_id).R';
    pos_diff = R_diff * images1(image_id).t - images2(image_id).t;
    im_diffs = [double(image_id) rotm2axang(R_diff) pos_diff'];
    diffs(diff_idx,:) = im_diffs;
    diff_idx = diff_idx + 1;
end

for image_id = im_ids_2_not_1
    diffs(diff_idx,:) = [image_id zeros(1,7)];
    diff_idx = diff_idx + 1;
end

for image_id = im_ids_1_not_2
    diffs(diff_idx,:) = [image_id inf(1,7)];
    diff_idx = diff_idx + 1;
end