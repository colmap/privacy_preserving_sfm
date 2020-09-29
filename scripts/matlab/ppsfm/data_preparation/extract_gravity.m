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

function gravity = extract_gravity(images)
% This function extracts a gravity vector for each image of a COLMAP model.
% This is solely based on the image poses in the model, so the extracted
% gravity direction is consistent for the whole model, but not necessary
% the correct gravity direction in the images.
% The vector for each image is given in the image's coordinate system.
%
% Input arguments:
% images    The image data map as given by 'read_model'
%
% Output arguments:
% gravity   Map with image IDs as keys, gravity vector as value

image_ids = images.keys;
gravity = containers.Map('KeyType','int64','ValueType','any');

for id_cell = image_ids
    image_id = id_cell{1};
    image = images(image_id);
    
    % Gravity direction in COLMAP is the positive Y-axis and the image
    % poses are given as transformation from world to image, so we only
    % need to apply the rotation
    gravity(image_id) = image.R * [0 1 0]';
end