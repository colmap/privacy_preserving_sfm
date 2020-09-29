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

function run_evaluations(dataset_folder, model_folder_name)

for set=["Castle", "Castle_large", "Entry", "Herzjesu", "Herzjesu_large", "Fountain"]
    gt_folder_path = fullfile(dataset_folder, char(set), 'groundtruth');
    model_path = fullfile(dataset_folder, char(set), model_folder_name);
    errors = compare_strecha_poses(gt_folder_path, model_path);
    
    % Only consider registered images
    errors = errors(~isinf(errors(:,1)), :);
    
    % 5th column is rotation error angle in rad
    mean_rot_err = mean(abs(errors(:,5)));
    rot_err_std_dev = std(abs(errors(:,5)));
    
    pos_err = sqrt(sum(errors(:,6:8).^2,2));
    mean_pos_err = mean(pos_err);
    pos_err_std_dev = std(pos_err);
    
    fprintf('%s mean orientation error (deg): %f (%f)\n', char(set), rad2deg(mean_rot_err), rad2deg(rot_err_std_dev));
    fprintf('%s mean position error (cm): %f (%f)\n', char(set), mean_pos_err * 100, pos_err_std_dev * 100);
end