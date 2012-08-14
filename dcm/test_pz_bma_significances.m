function passed = test_pz_bma_significances

    nsamp = 10000;        
    a_matrix = zeros(3,3, nsamp);
    
    % The connection from 1->1 is centred on zero
    a_matrix(1,1,1:nsamp) = pickNormalNumbers(1, 1, nsamp, 0, 1);
    
    % The connection from 1->2 is weakly positive
    a_matrix(2,1,1:nsamp) = pickNormalNumbers(1, 1, nsamp, 0.75, 1);
    
    % The connection from 1->3 is weakly negative
    a_matrix(3,1,1:nsamp) = pickNormalNumbers(1, 1, nsamp, -0.75, 1);

    % The connection from 2->2 is strongly positive
    a_matrix(2,2,1:nsamp) = pickNormalNumbers(1, 1, nsamp, 2, 1);    
    
    % The connection from 3->3 is strongly negative
    a_matrix(3,3,1:nsamp) = pickNormalNumbers(1, 1, nsamp, -2, 1);
        
    % Build BMS
    BMS = struct();    
    BMS.DCM.rfx.bma.a = a_matrix;    
    BMS.DCM.rfx.bma.nsamp = nsamp;    
    
    % Run
    significant_connections = ...
        pz_bma_significances(BMS, 'A');
    
    % Check
    assert(significant_connections(1,1) == 0);
    assert(significant_connections(2,1) == 0);
    assert(significant_connections(3,1) == 0);
    assert(significant_connections(2,2) == 1);
    assert(significant_connections(3,3) == 1);
    
    passed = true;

    % ---------------------------------------------------------------------
    function numbers = pickNormalNumbers(n, m, p, mean, sd)
        % m,n,p - dimensions
        % mean, sd - params
        numbers = mean + sd.*randn(n,m,p);
    end
end