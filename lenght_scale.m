clc;
close all;
clear all;
%%

test_condition = {'Phi10U03H000', 'Phi10U07H000', 'Phi10U11H000', 'Phi10U15H000', 'Phi10U19H000', 'Phi10U23H000',...
                  'Phi10U03H040', 'Phi10U07H040', 'Phi10U11H040', 'Phi10U15H040', 'Phi10U19H040', 'Phi10U23H040'};


for k = 9:9

    load(['C:\Users\vahidah\OneDrive - UBC\Codes\OH\SPOD\' ...
    cell2mat(test_condition(k)) '\nfft1024_novlp512_nblks20\spod_energy.mat']);

    L_C = zeros(99,129);
    final_area = zeros(99,129);
    num = zeros(99,129);

    for i = 25:25

        load(['C:\Users\vahidah\OneDrive - UBC\Codes\OH\SPOD\' ...
            cell2mat(test_condition(k)) '\nfft1024_novlp512_nblks20\' sprintf('spod_f%04d.mat', i)]);

        Psi  = rot90(reshape(Psi,[1024, 1024]), 2);

        freq = f(i);
        dt = 0.01*1/freq;
        nCycles = 1;
        fps = 20;
        T  = nCycles / freq;
        t  = 0:dt:T;
        
        
        [images] = mode_movie(Psi, freq, dt, nCycles, fps);
        der      = second_order(images, t);
        threshold = 0.1;
    
        for j = 1:size(der, 3)
            [L_C(j,i), final_area(j,i),  num(j,i)]  = source_size(squeeze(der(:,:,j)), threshold);
        end
    
    end
    disp(cell2mat(test_condition(k)));
    save(['Results\Source_size\' cell2mat(test_condition(k)) '.mat'], "L_C", "final_area", "num");
end



function [images] = mode_movie(Phi, f, dt, nCycles, fps)
    
    if nargin<3, dt = 1e-4; end
    if nargin<4, nCycles = 2; end
    if nargin<5, fps = 20; end
    
    A  = abs(Phi);
    th = angle(Phi);

    T  = nCycles / f; 
    t  = 0:dt:T; 
    omega = 2*pi*f;
    
    images = zeros(1024,1024, numel(t));

    for k = 1:length(t)
        q = A .* cos(omega*t(k) + th);
        images(:,:,k) = q;

    end
end

function der = second_order(images, time)


    dim=size(images);
    der = (images(:,:,3:dim(3)) - images(:,:,1:dim(3)-2))./(2*(time(2) - time(1)));

end


function [L, final_area,  num_tot] = source_size(der, threshold)
        copy = der;

        sigma = 1.0;
        deriv_img_smoothed = imgaussfilt(copy, sigma);
        
        deriv_img_smoothed = deriv_img_smoothed./max(abs(deriv_img_smoothed(:)));
        
        binary_tot = abs(deriv_img_smoothed) > threshold;
        
        [labeled_tot, num_tot] = bwlabel(binary_tot, 8);

        area = zeros(1,num_tot);
        
        for i = 1:num_tot
            ind = find(labeled_tot == i);
            area(i) = numel(ind);
        end
        
        area_weight = area ./ sum(area);

        final_area = sum(area_weight .* area);
                
        L = sqrt(final_area * 4 / pi);

end