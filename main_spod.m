clc;
clear all;
close all; 

folders = {'2022-09-14','2022-09-15'};
h2_per  = {'040', '000'};
files   = {'03','07','11','15','19','23'};

load Background.mat
load WF.mat

background   = load("Background.mat");
background   = struct2cell(background);
background   = cell2mat(background);

% Loading white field
white_field  = load("WF.mat");
white_field  = struct2cell(white_field);
white_field  = cell2mat(white_field);

median_window_size = 11;

% white_field1 = medfilt2(white_field, [median_window_size median_window_size]);
white_field  = white_field/max(max(white_field));
white_field  = white_field.*(white_field > 0.22);

for index = 1:numel(folders)

    for sub_index = 1:numel(files)

        load(['C:\Users\vahidah\OneDrive - UBC\Data Processing\PremixedMicrophone_September 2022\OH\Rotatioin\' cell2mat(folders(index)) '\Phi10U' cell2mat(files(sub_index)) 'H' cell2mat(h2_per(index)) '.mat'])
        % Clearing useless variables
        clear center centroidX centroidY;
        
        number         = 1000001:1:1010898;
        end_idx        = 5120;
        p              = zeros(end_idx,1024,1024);

        for i = 1:end_idx

            filename_image =  ['D:\UBC-PremixedMicrophone_September 2022\' ...
                              cell2mat(folders(index)) '\Conditions' '\Phi10U' cell2mat(files(sub_index)) 'H' cell2mat(h2_per(index)) '\' ...
                              'OHs\1111\1_C001H001S000' num2str(number(i)) '.tif'];
        
            mg_file = im2double(imread(strcat(filename_image)));
             
            % [1] Subtracting the background
            mg_file = mg_file - background;
            % [2] Subtracting the mean
            mg_file = mg_file - average_image; 
            % [3] Normaalizing with white field
            mg_file = mg_file./white_field;
            % [4] Deleting the fringes
            mg_file  = mg_file .* (mg_file < 1);
            mg_file  = mg_file .* (mg_file > 0 );
            % [5] Removing NaN values
            mg_file(isnan(mg_file)) = 0;
            % [6] Rotating the image
            % copy = rot90(image(:,:,i), 2);
            % [7] Median filter
            mg_file = medfilt2(mg_file, [median_window_size median_window_size]);
            
            p(i,:,:) = mg_file;        
        end
    
        %% SPOD


        opts.savefft   = true;
        opts.deletefft = true;
        opts.nsave     = 3;
        opts.savedir   = ['F:\SPOD-modes-sajjad-test\Phi10U' cell2mat(files(sub_index)) 'H' cell2mat(h2_per(index))];
        
        nDFT    = 1024;
        nOvlp   = nDFT/2;
        dt      = 1/10000;
        window  = hammwin(nDFT);
        [~,~,~] = spod(p,window,ones(1024*1024,1),nOvlp,dt,opts);

        clear p;
        disp(['Phi10U' cell2mat(files(sub_index)) 'H' cell2mat(h2_per(index)) '-done!'])
    end
    
end
