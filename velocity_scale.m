clc;
clear all;
close all;
%% 

test_condition = {'Phi10U03H000', 'Phi10U07H000', 'Phi10U11H000', 'Phi10U15H000', 'Phi10U19H000', 'Phi10U23H000',...
                  'Phi10U03H040', 'Phi10U07H040', 'Phi10U11H040', 'Phi10U15H040', 'Phi10U19H040', 'Phi10U23H040'};
pixel_size = 81.52e-06;

for k = 9:9

    load(['C:\Users\vahidah\OneDrive - UBC\Codes\OH\SPOD\' ...
    cell2mat(test_condition(k)) '\nfft1024_novlp512_nblks20\spod_energy.mat']);
    
    u_conv = zeros(1,206);
    for i = 25:25%1:206
    
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
       
        [ph_y, ph_x] = gradient(ph1, pixel_size, pixel_size);
        k_loc = hypot(ph_x, ph_y);

        Umap = (2*pi*f(i)) ./ k_loc;
        Umap(~mask) = NaN;
        % unwrap_ph = unwrap_phase(ph);
        
        
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