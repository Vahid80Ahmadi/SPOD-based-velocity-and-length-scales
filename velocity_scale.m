clc;
clear all;
close all;
%% 

test_condition = {'Phi10U03H000', 'Phi10U07H000', 'Phi10U11H000', 'Phi10U15H000', 'Phi10U19H000', 'Phi10U23H000',...
                  'Phi10U03H040', 'Phi10U07H040', 'Phi10U11H040', 'Phi10U15H040', 'Phi10U19H040', 'Phi10U23H040'};
pixel_size = 81.52e-06;

for k = 1:12

    load(['C:\Users\vahidah\OneDrive - UBC\Codes\OH\SPOD\' ...
    cell2mat(test_condition(k)) '\nfft1024_novlp512_nblks20\spod_energy.mat']);

    u_conv = zeros(1,206);
    for i = 1:206
    
        load(['C:\Users\vahidah\OneDrive - UBC\Codes\OH\SPOD\' ...
            cell2mat(test_condition(k)) '\nfft1024_novlp512_nblks20\' sprintf('spod_f%04d.mat', i)]);

        Psi  = rot90(reshape(Psi,[1024, 1024]), 2);
        
        amp = abs(Psi);
        amp = amp ./ max(amp(:));

        mask = amp > 0.05; 

        ph = angle(Psi);

        mask_u = true(1024,1024);
        mask_u(ph==0) = false;

        ph1 = unwrap2D_masked_LS(ph, mask_u);
       
        nanmask = isnan(ph1);

        ph1(isnan(ph1)) = 0;

        

        [Dx, Dy] = make_differentiation_matrices(1024, 1024);

        grad_I = tv_regularized_gradient(...
                                        ph1, ...
                                        'Dx', Dx, 'Dy', Dy, ...
                                        'tuning_parameter', 1, ...
                                        'splitting_parameter', 1, ...
                                        'iterations', 10, ...
                                        'diagnostics', true);


        k_loc = hypot(squeeze(grad_I(1, :, :)), squeeze(grad_I(2, :, :))) ./ pixel_size;

        Umap = (2*pi*f(i)) ./ k_loc;
        Umap(~mask) = NaN;

        u_conv(i) = mean(Umap(:), 'omitnan');
        disp(i)

    end
    disp(cell2mat(test_condition(k)));
    save(['U_conv_3\' cell2mat(test_condition(k)) '.mat'], "u_conv");
end

function phi = unwrap2D_masked_LS(ph, mask)


    wrapToPi = @(a) mod(a + pi, 2*pi) - pi;

    [nr,nc] = size(ph);

    idx = zeros(nr,nc);
    idx(mask) = 1:nnz(mask);
    nUnknown = nnz(mask);

    maskH = mask(:,1:end-1) & mask(:,2:end);
    [rH,cH] = find(maskH);
    iH = idx(sub2ind([nr,nc], rH, cH));
    jH = idx(sub2ind([nr,nc], rH, cH+1));
    dH = wrapToPi(ph(sub2ind([nr,nc], rH, cH+1)) - ph(sub2ind([nr,nc], rH, cH)));

    maskV = mask(1:end-1,:) & mask(2:end,:);
    [rV,cV] = find(maskV);
    iV = idx(sub2ind([nr,nc], rV, cV));
    jV = idx(sub2ind([nr,nc], rV+1, cV));
    dV = wrapToPi(ph(sub2ind([nr,nc], rV+1, cV)) - ph(sub2ind([nr,nc], rV, cV)));

    nEqH = numel(iH);
    nEqV = numel(iV);
    nEq  = nEqH + nEqV + 1;

    rows = zeros(2*nEqH + 2*nEqV + 1,1);
    cols = zeros(2*nEqH + 2*nEqV + 1,1);
    vals = zeros(2*nEqH + 2*nEqV + 1,1);
    b    = zeros(nEq,1);

    p = 0;

    for k = 1:nEqH
        eq = k;
        p = p + 1; rows(p) = eq; cols(p) = jH(k); vals(p) =  1;
        p = p + 1; rows(p) = eq; cols(p) = iH(k); vals(p) = -1;
        b(eq) = dH(k);
    end

    for k = 1:nEqV
        eq = nEqH + k;
        p = p + 1; rows(p) = eq; cols(p) = jV(k); vals(p) =  1;
        p = p + 1; rows(p) = eq; cols(p) = iV(k); vals(p) = -1;
        b(eq) = dV(k);
    end

    eq = nEq;
    p = p + 1; rows(p) = eq; cols(p) = 1; vals(p) = 1;
    b(eq) = 0;

    A = sparse(rows(1:p), cols(1:p), vals(1:p), nEq, nUnknown);

    phi_vec = A\b;

    phi = nan(nr,nc);
    phi(mask) = phi_vec;
end



function gradient = tv_regularized_gradient(image, varargin)
    p = inputParser;
    addRequired(p, 'image');
    addParameter(p, 'tuning_parameter', 1e0);
    addParameter(p, 'splitting_parameter', []);
    addParameter(p, 'iterations', 10);
    addParameter(p, 'p', 1.0);
    addParameter(p, 'epsilon', 0.0);
    addParameter(p, 'diagnostics', true);
    addParameter(p, 'Dx', []);
    addParameter(p, 'Dy', []);
    addParameter(p, 'dtype', 'single');
    addParameter(p, 'scaling', false);
    addParameter(p, 'p_only', false);
    parse(p, image, varargin{:});
    opts = p.Results;

    if strcmpi(opts.dtype, 'single')
        castfun = @single;
    else
        castfun = @double;
    end

    image0 = castfun(image);
    if opts.scaling
        scale = max(abs(image(:)));
        image0 = image0 / scale;
    else
        scale = 1.0;
    end

    [rows, columns] = size(image);

    if opts.p_only
        Dx = speye(rows * columns);
        Dy = speye(rows * columns);
    else
        Dx = opts.Dx;
        Dy = opts.Dy;
        if isempty(Dx) || isempty(Dy)
            [dx, dy] = make_differentiation_matrices(rows, columns, 1, true, 'periodic');
            if isempty(Dx)
                Dx = dx;
            end
            if isempty(Dy)
                Dy = dy;
            end
        end
    end

    gradient = p_admm(image0, opts.tuning_parameter, opts.splitting_parameter, ...
        opts.iterations, Dx, Dy, opts.p, opts.epsilon, opts.diagnostics, opts.dtype, opts.p_only);

    gradient = scale * gradient;
    gradient = castfun(gradient);
end


function U = p_admm(array, mu, lmbda, iterations, Dx, Dy, p, epsilon, diagnostics, dtype, p_only)
    [rows, columns] = size(array);
    num = rows * columns;

    array = double(array);

    U = zeros(2, rows, columns, 'double');
    W = zeros(4, num, 'double');
    B = zeros(4, num, 'double');

    if isempty(lmbda)
        lmbda = 1 / double(mu);
    end

    c = 1.6;
    [xker, yker, xadjker, yadjker] = make_kernels(rows, columns, mu, lmbda, dtype, p_only);
    [muKxTb, muKyTb] = apply_adjoint(array, xadjker, yadjker, mu, dtype);

    for itr = 1:iterations
        rhs = (Dx' * (W(1, :)' - B(1, :)') + Dy' * (W(2, :)' - B(2, :)')) / lmbda;
        rhs = muKxTb + reshape_rowmajor(rhs, rows, columns);
        out = fft2(rhs);
        U(1, :, :) = real(ifft2(out .* xker));

        rhs = (Dx' * (W(3, :)' - B(3, :)') + Dy' * (W(4, :)' - B(4, :)')) / lmbda;
        rhs = muKyTb + reshape_rowmajor(rhs, rows, columns);
        out = fft2(rhs);
        U(2, :, :) = real(ifft2(out .* yker));

        U1 = squeeze(U(1, :, :));
        U2 = squeeze(U(2, :, :));
        DU = [Dx * vec_rowmajor(U1), ...
              Dy * vec_rowmajor(U1), ...
              Dx * vec_rowmajor(U2), ...
              Dy * vec_rowmajor(U2)]';

        W = p_shrink(DU + B, lmbda, p, epsilon);
        B = B + c * (DU - W);

        if diagnostics
            objective = sum(sum(DU.^2, 1).^(0.5 * p));
            cstr = sqrt(sum((W(:) - DU(:)).^2));
            % fprintf('iteration = %2d, objective = %7.3e, constraint = %7.3e\n', itr - 1, objective, cstr);
        end
    end
end


function [xker, yker, xadjker, yadjker] = make_kernels(rows, columns, mu, lmbda, dtype, p_only)
    if strcmpi(dtype, 'single')
        obj = zeros(rows, columns, 'single');
    else
        obj = zeros(rows, columns, 'double');
    end

    obj(:, :) = 0.0;
    obj(1, 1) = -1.0;
    obj(1, 2) = 1.0;
    xker = fft2(obj);
    xadjker = xker;
    xadjker(:, 2:end) = 1.0 ./ xadjker(:, 2:end);
    xadjker(:, 1) = 0.0;
    xker = abs(xker).^2;

    obj(:, :) = 0.0;
    obj(1, 1) = -1.0;
    obj(2, 1) = 1.0;
    yker = fft2(obj);
    yadjker = yker;
    yadjker(2:end, :) = 1.0 ./ yadjker(2:end, :);
    yadjker(1, :) = 0.0;
    yker = abs(yker).^2;

    if p_only
        laplace_ker = 2.0 / lmbda;
    else
        laplace_ker = (xker + yker) / lmbda;
    end

    yker = yker ./ (yker .* laplace_ker + mu);
    xker = xker ./ (xker .* laplace_ker + mu);
end


function [muKxTd, muKyTd] = apply_adjoint(data, xadjker, yadjker, mu, ~)
    datahat = fft2(data);
    KxTd = real(ifft2(xadjker .* datahat));
    KyTd = real(ifft2(yadjker .* datahat));
    muKxTd = mu * KxTd;
    muKyTd = mu * KyTd;
end


function W = p_shrink(X, lmbda, p, epsilon)
    magnitude = sqrt(sum(X.^2, 1));
    nonzero = magnitude;
    nonzero(magnitude == 0.0) = 1.0;
    magnitude = max(magnitude - lmbda^(2.0 - p) .* (nonzero.^2 + epsilon).^(p / 2.0 - 0.5), 0) ./ nonzero;
    W = X .* magnitude;
end


function [Dx, Dy, Dz] = make_differentiation_matrices(rows, columns, channels, no_z, boundary_conditions)
    if nargin < 3 || isempty(channels)
        channels = 1;
    end
    if nargin < 4 || isempty(no_z)
        no_z = true;
    end
    if nargin < 5 || isempty(boundary_conditions)
        boundary_conditions = 'periodic';
    end

    D = spdiags([-ones(columns, 1), ones(columns, 1)], [0 1], columns, columns);
    switch lower(boundary_conditions)
        case 'neumann'
            D(columns, columns) = 0.0;
        case 'periodic'
            D(columns, 1) = 1.0;
        otherwise
    end
    S = speye(rows);
    Sz = speye(channels);
    Dx = kron(kron(S, D), Sz);

    D = spdiags([-ones(rows, 1), ones(rows, 1)], [0 1], rows, rows);
    switch lower(boundary_conditions)
        case 'neumann'
            D(rows, rows) = 0.0;
        case 'periodic'
            D(rows, 1) = 1.0;
        otherwise
    end
    S = speye(columns);
    Dy = kron(kron(D, S), Sz);

    if no_z
        Dz = [];
        return;
    end

    D = spdiags([-ones(channels, 1), ones(channels, 1)], [0 1], channels, channels);
    switch lower(boundary_conditions)
        case 'neumann'
            D(channels, channels) = 0.0;
        case 'periodic'
            D(channels, 1) = 1.0;
        otherwise
    end
    Sr = speye(rows);
    Sc = speye(columns);
    Dz = kron(kron(Sr, Sc), D);
end


function v = vec_rowmajor(A)
    v = reshape(A.', [], 1);
end


function A = reshape_rowmajor(v, rows, columns)
    A = reshape(v, columns, rows).';
end

