function significant_connections = ...
    pz_bma_significances(BMS, matrix_letter, matrix_number, significance_cutoff)
    % Tests all connections in a BMS for significance against zero. Follows
    % the DCM convention of the most significant p-value being 1. Based 
    % on code by Mohamed Seghier.
    %
    % Inputs:
    % BMS                 - the loaded BMS structure
    % matrix_letter       - 'A' or 'B'
    % matrix_number       - (optional) the B-matrix index (e.g. 2 or 3)
    % significance_cutoff - (optional) 1-pvalue, e.g. 0.99 => p=0.01
    %
    % Outputs:
    % significant_connection_matrix - binary matrix of significant
    % connections (>= significance_cutoff)
    %
    % very_significant_connection_matrix - binary matrix of very
    % significant connections (>= very_significant_cutoff)
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

    % Validate input
    assert(~isempty(BMS) && isstruct(BMS), 'BMS must be a structure');
    assert(~isempty(matrix_letter) && ismember(matrix_letter, ['A','B','C']));    

    if (nargin < 3 || isempty(matrix_number))
        if strcmp(matrix_letter, 'B') 
            error('A matrix number must be provided for the B-matrix');
        else
            matrix_number = 1;
        end
    end      
    
    if (nargin < 4 || isempty(significance_cutoff))
        significance_cutoff = 0.95;
    end
    
    % Get priors from an example subject
    sample_dcm = get_sample_dcm(BMS);
    priors = sample_dcm.DCM.M.pE.(matrix_letter);
    if strcmp(matrix_letter, 'B') 
        priors = priors(:,:,matrix_number);
    end
     
    % Get 10,000 samples of each connection
    connectivity_matrix = get_connectivity_matrix(matrix_letter, matrix_number, BMS);
    nsamp = get_num_samples(BMS);    

    % Compute probability of each connection. This is the proportion of
    % samples greater than the prior.    
    p_connectivity = zeros(size(connectivity_matrix,1), size(connectivity_matrix,2));
    for samp = 1:nsamp
        p_connectivity = p_connectivity + (connectivity_matrix(:,:,samp) > priors );
    end
    p_connectivity = p_connectivity / nsamp;
            
    % Correct p-value for sign
    mean_bms = mean(connectivity_matrix,3);
    p_connectivity(mean_bms < priors) = 1 - p_connectivity(mean_bms < priors);

    % Identify significant connections
    significant_connections = p_connectivity >= significance_cutoff;
    
    % Output some nice facts and figures   
    disp(' ');
    disp('Priors:');
    disp( num2str(priors) );
    disp(' ');    
    
    disp('BMS Mean values:');    
    disp( num2str(mean_bms) );
    disp(' ');

    disp('p-values:');
    disp( num2str(p_connectivity) );
    disp(' ');
    
    fprintf('Significant connections (p > %2.2f):\n', significance_cutoff);
    disp( num2str(mean_bms .* significant_connections) );
    disp(' ');      
    
    %----------------------------------------------------------------------
    function matrix = get_connectivity_matrix(matrix_letter, matrix_number, BMS)
        % Retrieves the requested A,B or C matrix. Dimensions are the
        % original size of the matrix x 10,000 samples from the probability
        % distribution.
        
        % rfx or ffx?
        if isfield(BMS.DCM, 'rfx') == 1
            disp('Using random effects results');
            BMA = BMS.DCM.rfx.bma;
        else
            disp('Using fixed effects results');
            BMA = BMS.DCM.ffx.bma;
        end;
        
        % matrix a or b?
        if strcmp(matrix_letter, 'A') == 1
            matrix = BMA.a;
        elseif strcmp(matrix_letter, 'B') == 1
            matrix = squeeze( BMA.b(:,:,matrix_number,:) );
        elseif strcmp(matrix_letter, 'C') == 1
            matrix = BMA.c;
        else
            error('The matrix letter have value A or B');
        end
    end

    %----------------------------------------------------------------------
    function nsamp = get_num_samples(BMS)
        % Gets the number of samples generated for the probability
        % distribution of each parameter
        if isfield(BMS.DCM, 'rfx') == 1
            nsamp = BMS.DCM.rfx.bma.nsamp;
        else
            nsamp = BMS.DCM.ffx.bma.nsamp;
        end
    end

    %----------------------------------------------------------------------
    function sample_dcm = get_sample_dcm(BMS)
        % Loads an example DCM from the BMS
        if isfield(BMS.DCM, 'rfx') == 1
            model_space = load(BMS.DCM.rfx.data);
        else
            model_space = load(BMS.DCM.ffx.data);
        end
        sample_dcm = load( model_space.subj(1).sess.model(1).fname );
    end
end