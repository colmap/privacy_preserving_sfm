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

function plot_model_comparison_impl(colmap_images, privacy_images)

% Compute image positions to simplify the visualization later
colmap_images = update_image_positions(colmap_images);
privacy_images = update_image_positions(privacy_images);

% Check which images are contained in both models and which are unique
colmap_image_ids = cell2mat(colmap_images.keys);
privacy_image_ids = cell2mat(privacy_images.keys);

common_image_ids = intersect(colmap_image_ids, privacy_image_ids);
extra_colmap_image_ids = setdiff(colmap_image_ids, common_image_ids);
extra_privacy_image_ids = setdiff(privacy_image_ids, common_image_ids);

num_common = size(common_image_ids, 2);

fig = figure;
ax = axes(fig);
axis(ax, 'equal');
hold(ax, 'on');
view(ax, [0 -1 0]);

% For common images, plot the colmap model ones and the offset to the
% privacy ones
for im_id = common_image_ids
    image = colmap_images(im_id);
    plotCamera( ...
        'Parent', ax, 'Orientation', image.R, 'Location', image.position, ...
        'Size', 0.1);
end

pos_diffs = zeros([2, num_common, 3]);
for i = 1:num_common
    id = common_image_ids(i);
    pos1 = colmap_images(id).position;
    pos2 = privacy_images(id).position;
    pos_diffs(:,i,:) = permute([pos1 pos2]', [1 3 2]);
end

line(ax, pos_diffs(:,:,1), pos_diffs(:,:,2), pos_diffs(:,:,3), 'Color', 'red');

% Plot the extra colmap images
for im_id = extra_colmap_image_ids
    image = colmap_images(im_id);
    plotCamera( ...
        'Parent', ax, 'Orientation', image.R, 'Location', image.position, ...
        'Size', 0.1, 'Color', 'yellow');
end

% Plot the additional privacy images
for im_id = extra_privacy_image_ids
    image = privacy_images(im_id);
    plotCamera( ...
        'Parent', ax, 'Orientation', image.R, 'Location', image.position, ...
        'Size', 0.1, 'Color', 'blue');
end