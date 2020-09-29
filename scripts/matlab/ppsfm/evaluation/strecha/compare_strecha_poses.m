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

function diffs = compare_strecha_poses(strecha_groundtruth_folder, model_path)

[~, model_images, ~] = read_model(model_path);

model_image_ids = cell2mat(model_images.keys());
num_model_images = numel(model_image_ids);

gt_images = get_strecha_image_gt(strecha_groundtruth_folder);
num_gt_images = numel(gt_images.keys());
% gt_image_names = gt_images.keys();

% Frist set all errors to inf, this way we don't need to search for issing
% images later
diffs = inf(num_gt_images, 8);

diff_idx = 1;
for image_id = model_image_ids
    % The convention for poses is inverse in the models.
    % We report errors in the strecha convention
    image_name = model_images(image_id).name;
    
    R_diff = quat2rotm(quatmultiply(rotm2quat(gt_images(image_name).R'), rotm2quat(model_images(image_id).R')));
    pos_diff = -R_diff * model_images(image_id).t - gt_images(image_name).R' * gt_images(image_name).pos;
    
    im_diffs = [double(image_id) rotm2axang(R_diff) pos_diff'];
    diffs(diff_idx,:) = im_diffs;
    diff_idx = diff_idx + 1;
end