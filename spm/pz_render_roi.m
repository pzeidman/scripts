function pz_render_roi(roi_file, render_file)
    % Renders an ROI on a 3D surface-rendered brain
    %
    % Inputs: 
    % roi_file - the ROI image to render
    % render_file - (optional) the surface rendering to draw on
    %
    % ---------------------------------------------------------------------
    % Copyright (C) 2012 Peter Zeidman
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    % 
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    % 
    % You should have received a copy of the GNU General Public License
    % along with this program.  If not, see <http://www.gnu.org/licenses/>.   
    % ---------------------------------------------------------------------    

    assert(~isempty(roi_file) && exist(roi_file,'file'), 'Could not find ROI file');

    if nargin < 2
        render_file = fullfile(spm('Dir'),'canonical','cortex_20484.surf.gii');
    end

    % Read the ROI
    V = spm_vol(roi_file);
    [t xyz] = spm_read_vols(V);

    % Convert from mm space to voxel space
    M   = inv(V.mat);
    xyz = M(1:3,:)*[xyz ; ones(1,size(xyz,2))];
    
    % Display
    dat    = struct('XYZ',  xyz,...
      't',    t,...
      'mat',  V.mat,...
      'dim',  V.dim(1:3)');    
    spm_render(dat,NaN,render_file);
end