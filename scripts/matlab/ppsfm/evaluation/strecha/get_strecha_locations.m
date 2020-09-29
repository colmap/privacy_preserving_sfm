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

function get_strecha_locations2(strecha_groundtruth_dir)

images = get_strecha_image_gt(strecha_groundtruth_dir);

image_names = images.keys();
num_files = numel(image_names);
image_positions = zeros(num_files,3);


for i=1:num_files
    image_name = image_names{i};
    image = images(image_name);
    image_positions(i,:) = image.pos';
end

f = fopen(fullfile(strecha_groundtruth_dir, 'image_locations.txt'), 'w+');
for i=1:num_files
    fprintf(f, '%s %.8f %.8f %.8f\n', image_names{i}, image_positions(i,:));
end
fclose(f);