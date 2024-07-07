    function [u,epsi,epsi1,epsi3,sigma,sigma1,sigma3] = solid3D(nodes,elements,load1,bcs)

    ks = [];
    n = size(nodes, 1);
    K = zeros(n*3, n*3);
    r = zeros(n*3, 1);


    
    for i = 1:size(elements, 1)

        element = elements(i,:);
        X = nodes(element(2:5), 2:4);

        E = element(6);
        nu = element(7);

        V = volume(X);
        B = calculateBMatrix(X);
        C = calculateCMatrix(E, nu);
        k = V * B * C * B';
        idxs = [];

        % Add to global Stiffness matrix;
        for j = 0:3
            idx = (element(j+2)-1)*3+1;
            idxs = [idxs, idx, idx+1, idx+2];
        end
        K(idxs, idxs) = K(idxs, idxs) + k;
    end

    % assemble force vector
    load_idxs = [];
    for i = 1:size(load1,1)
        f = load1(i,:);
        idx = (f(1)-1)*3+1;
        load_idxs = [load_idxs, idx:idx+2];
        r(idx:idx+2) = f(2:4);
    end
    unknown_load_idxs = setdiff(1:(n*3), load_idxs);

    % add BCs
    bc_idxs = [];
    non_zero_bcs = [];
    zero_bcs = [];
    u = NaN(n*3, 1);
    for i = 1:size(bcs,1)
      bc = bcs(i,:);
      idx = (bc(1)-1)*3+1:(bc(1)-1)*3+3;
      % remove NaN
      not_Nan = find(~isnan(bc(2:4)));
      bc_idxs = [bc_idxs, idx(not_Nan)];
      u(idx(not_Nan)) = bc(not_Nan+1);
    end
    unknown_bc_idxs = setdiff(1:(n*3), bc_idxs);

    % Take into account bcs    
    
    % Create add 1 to diagonal where we have a BC
    m = size(bc_idxs,2);
    K_diag = zeros(m, n*3);
    for i = 1:m
       K_diag(i,bc_idxs(i)) = 1; 
    end
    K(bc_idxs,:) = K_diag;

    %add non-0 BC values to force vector
    r(bc_idxs) = u(bc_idxs);
    known_idxs = [non_zero_bcs ,load_idxs];

    
    % Solve for unknown displacements
    u = K\r;

    % Calculate strains & stressesc

    epsi = zeros(size(elements, 1), 7);
    sigma = zeros(size(elements, 1), 7);
    sigma1 = [];
    sigma3 = [];
    epsi1 = [];
    epsi3 = [];

    for i = 1:size(elements, 1)

        element = elements(i,:);
        X = nodes(element(2:5), 2:4);

        E = element(6);
        nu = element(7);

        V = volume(X);
        B = calculateBMatrix(X);
        C = calculateCMatrix(E, nu);
        idxs = [];
        
        for j = 0:3
            idx = (element(j+2)-1)*3+1;
            idxs = [idxs, idx:idx+2];
        end

        eps = B' * u(idxs);
        sig = C * eps;

        % add final strains and stresses
        
        epsi(i, 1) = element(1);
        epsi(i, 2:7) = eps;

        sigma(i, 1) = element(1);
        sigma(i, 2:7) = sig;



        % calculate principal strains

        epsilon_mat = [eps(1), eps(4)/2, eps(6)/2;...
            eps(4)/2, eps(2), eps(5)/2;...
            eps(6)/2, eps(5)/2, eps(3)];

        sigma_mat = [sig(1), sig(4), sig(6);...
            sig(4), sig(2), sig(5);...
            sig(6), sig(5), sig(3)];

        [V, D_eps] = eig(epsilon_mat);
        [V, D_sig] = eig(sigma_mat);
        principal_strains = sort(diag(D_eps), 'ascend');
        principal_stresses = sort(diag(D_sig), 'ascend');

        sigma1 = [sigma1; principal_stresses(3)];

        sigma3 = [sigma3; principal_stresses(1)];

        epsi1 = [epsi1; principal_strains(3)];

        epsi3 = [epsi3; principal_strains(1)];

    end

 
    
end

function V = volume(X)

    x21 = X(2,:) - X(1,:);
    x31 = X(3,:) - X(1,:);
    x41 = X(4,:) - X(1,:);
    

    V = abs(dot(x21, cross(x31,x41)))/6;

end

function B = calculateBMatrix(X)
    x21 = X(2,:) - X(1,:);
    x31 = X(3,:) - X(1,:);
    x41 = X(4,:) - X(1,:);
    x32 = X(3,:) - X(2,:);
    x42 = X(4,:) - X(2,:);
    x43 = X(4,:) - X(3,:);
    C = num2cell(x21);
    [xx21, y21, z21] = C{:};
    C = num2cell(x31);
    [xx31, y31, z31] = C{:};
    C = num2cell(x41);
    [xx41, y41, z41] = C{:};
    C = num2cell(x32);
    [xx32, y32, z32] = C{:};
    C = num2cell(x42);
    [xx42, y42, z42] = C{:};
    C = num2cell(x43);
    [xx43, y43, z43] = C{:};

    a1 = X(2,2)*z43 - X(3,2)*z42 + X(4,2)*z32;
    a2 = -X(1,2)*z43 + X(3,2)*z41 - X(4,2)*z31;
    a3 = X(1,2)*z42 - X(2,2)*z41 + X(4,2)*z21;
    a4 = -X(1,2)*z32 + X(2,2)*z31 - X(3,2)*z21;
    
    b1 = -X(2,1)*z43 + X(3,1)*z42 - X(4,1)*z32;    
    b2 = X(1,1)*z43 - X(3,1)*z41 + X(4,1)*z31;
    b3 = -X(1,1)*z42 + X(2,1)*z41 - X(4,1)*z21;
    b4 = X(1,1)*z32 - X(2,1)*z31 + X(3,1)*z21;
    
    c1 = X(2,1)*y43 - X(3,1)*y42 + X(4,1)*y32;
    c2 = -X(1,1)*y43 + X(3,1)*y41 - X(4,1)*y31;
    c3 = X(1,1)*y42 - X(2,1)*y41 + X(4,1)*y21;
    c4 = -X(1,1)*y32 + X(2,1)*y31 - X(3,1)*y21;

    V = volume(X);

    B = zeros(12, 6);
    B([1 4 7 10], 1) = [a1 a2 a3 a4];
    B([2 5 8 11], 2) = [b1 b2 b3 b4];
    B([3 6 9 12], 3) = [c1 c2 c3 c4];
    B([2 1 5 4 8 7 11 10], 4) = [a1 b1 a2 b2 a3 b3 a4 b4];
    B([3 2 6 5 9 8 12 11], 5) = [b1 c1 b2 c2 b3 c3 b4 c4];
    B([1 3 4 6 7 9 10 12], 6) = [c1 a1 c2 a2 c3 a3 c4 a4];

    B = B/(6*V);

end

function C = calculateCMatrix(E, nu)
    % E is Young's modulus
    % nu is Poisson's ratio

    % Precompute scale factor
    scaleFactor = E / ((1 + nu) * (1 - 2 * nu));

    % Initialize the C matrix
    C = zeros(6,6);

    % Fill the diagonal blocks
    C(1:3, 1:3) = nu * ones(3, 3);
    C(1,1) = 1 - nu;
    C(2,2) = 1 - nu;
    C(3,3) = 1 - nu;

    % Fill the shear components
    C(4,4) = (1 - 2 * nu) / 2;
    C(5,5) = (1 - 2 * nu) / 2;
    C(6,6) = (1 - 2 * nu) / 2;

    % Scale the entire matrix
    C = scaleFactor * C;
end