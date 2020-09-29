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

function err_all = init_eval(inits_folder, num_inits, rel_path, reference_model_path)

% addpath ~/work/thirdparty/colmap/scripts/matlab/
% gt_path = '~/datasets/lund_cath_large/sparse_baseline/0/';
[~, images, ~] = read_model(reference_model_path);

%%

rot_angle = @(R) acos(max(-1.0,min(1.0,(trace(R)-1.0)/2)))*180/pi;
err_all = [];
err_max = [];

pos_err_all = [];
scales_all = [];

for i = 1:num_inits
    
    init_model_path = fullfile(inits_folder, num2str(i), rel_path);
       
    % Read model for init
    [~, init_img, ~] = read_model(init_model_path);
    ids = double(cell2mat(init_img.keys()));
        
    Pgt = {};
    P = {};
    for i = 1:4
        P{i} = [init_img(ids(i)).R init_img(ids(i)).t];        
        Pgt{i} = [images(ids(i)).R images(ids(i)).t];
    end
    
    
    H = [P{1}; 0 0 0 1];
    Hgt = [Pgt{1}; 0 0 0 1];
    for i = 1:4
        P{i} = P{i} * inv(H); 
        Pgt{i} = Pgt{i} * inv(Hgt);
    end
            
    % Compute errors
    err2 = rot_angle(P{2}(:,1:3)'*Pgt{2}(:,1:3));
    err3 = rot_angle(P{3}(:,1:3)'*Pgt{3}(:,1:3));
    err4 = rot_angle(P{4}(:,1:3)'*Pgt{4}(:,1:3));
       
    err_max = [err_max max([err2 err3 err4])];
    err_all = [err_all err2 err3 err4];
    
    % Compute relative positions
    pos = zeros(3,3);
    posgt = zeros(3,3);
    ata = 0;
    atb = 0;
    for i = 1:3
        proj_mat = P{i+1};
        pos(i,:) = (-proj_mat(1:3,1:3) * proj_mat(1:3,4))';
        proj_mat_gt = Pgt{i+1};
        posgt(i,:) = (-proj_mat_gt(1:3,1:3) * proj_mat_gt(1:3,4))';
        
        ata = ata + pos(i,:) * pos(i,:)';
        atb = atb + pos(i,:) * posgt(i,:)';
    end
    
    scale = atb / ata;
    
    scales_all = [scales_all scale];
    
    pos_diffs = zeros(3);
    pos_errs = zeros(1,3);
    pos = pos * scale;
    
    for i = 1:3
        pos_diffs(i,:) = pos(i,:) - posgt(i,:);
        pos_errs(i) = norm(pos_diffs(i,:));
    end
    
    pos_err_all = [pos_err_all pos_errs];
    
%     err_max = [err_max max([err2 err3])];
%     err_all = [err_all err2 err3];    
end


num_rot_bins = 20;
rot_edges = linspace(0,20,num_rot_bins+1);
bin_counts = zeros(1,num_rot_bins);
for i = 1:num_rot_bins
    bin_counts(i) = sum(err_all < rot_edges(i+1) & err_all >= rot_edges(i));
end
bin_counts(num_rot_bins) = size(err_all,2) - sum(bin_counts(1:(num_rot_bins-1)));

set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')

rot_fig = figure;
rot_ax = axes(rot_fig);
plot(rot_ax, rot_edges, [0 bin_counts] / size(err_all,2),'LineWidth', 2);
set(rot_fig, 'Position', [1691 780 321 146]);
xlabel(rot_ax, 'Orientation error (degrees)');

rot_ax2 = axes(figure);
histogram(rot_ax2, err_all, rot_edges);

num_pos_bins = 100;
pos_edges = linspace(0,1,num_pos_bins+1);
pos_bin_counts = zeros(1,num_pos_bins);
for i=1:num_pos_bins
    pos_bin_counts(i) = sum(pos_err_all < pos_edges(i+1) & pos_err_all >= pos_edges(i));
end
pos_bin_counts(num_pos_bins) = size(pos_err_all,2) - sum(pos_bin_counts(1:(num_pos_bins-1)));

pos_fig = figure;
pos_ax = axes(pos_fig);
plot(pos_ax, pos_edges*100, [0 pos_bin_counts] / size(pos_err_all,2), 'LineWidth', 2);
set(pos_fig, 'Position', [1691 780 321 146]);
xlabel(pos_ax,'Position error (cm)');

pos_ax2 = axes(figure);
histogram(pos_ax2, pos_err_all, pos_edges);
% scales_all