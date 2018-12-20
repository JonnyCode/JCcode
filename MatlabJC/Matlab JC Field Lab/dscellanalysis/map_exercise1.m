
% /Users/sravi/matlab/DS cell analysis/Raster Plots 8s datarun002/DS
% 
% /Users/sravi/matlab/DS cell analysis/Raster Plots 8s datarun002/DS/Likes Both Speeds
% 
% /Users/sravi/matlab/DS cell analysis/Raster Plots 8s datarun002/DS/Likes Both Speeds but 32 > 256
% 
% /Users/sravi/matlab/DS cell analysis/Raster Plots 8s datarun002/DS/Likes both Speeds but 256 > 32

    plot_rf_summaries(datarun, nonzeros(listin(a,:)));

        plot_time_courses(datarun, nonzeros(listin(a,:)), 1, 1);
        

        plot_rf_summaries(datarun, DSCELLS(2,(DSCELLS(2,:))~=0))
            plot_time_courses(datarun, nonzeros(listin(a,:)), 1, 1);

        
cellindic = DSCELLS(2,(DSCELLS(2,:))~=0);
cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
figure(4)
plot_time_courses(datarun00, A, 1, 1)
datarun02 = get_autocorrelations(datarun02, cellids);
plot_autocorrelograms(datarun02, cellids, 'foa', 0);
datarun00 = get_autocorrelations(datarun00, A);
plot_autocorrelograms(datarun00, A, 'foa', 0);

        

%Cell Type: Likes Only 256 Slow Speed
cellindic = [58 72 73 75 387 697 761 786];
cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.8) ;
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
        plot_time_courses(datarun00, A, 1, 1)

        
        

% subplot(2,1,1)
% plot_rf_summaries(datarun02, get_cell_indices(datarun02, 3109))
% subplot(2,1,2)
% plot_rf_summaries(datarun00, get_cell_indices(datarun00, 3111))

cellindic = [58 72 73 75 387 697 761 786];
cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun01, 'master_cell_type', cellids, 'corr_threshold', 0.8);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun01, A)
figure(2)
plot_rf_portraits(datarun01,A)
figure(3)
plot_rfs(datarun01, A)
        

%256 > 32
% ls('/Users/sravi/matlab/DS cell analysis/Raster Plots 8s datarun002/DS/Likes both Speeds but 256 > 32')
%  dirData = dir('/Users/sravi/matlab/DS cell analysis/Raster Plots 8s datarun002/DS/Likes both Speeds but 256 > 32'); 
%  
% [pathstr, name, ext] = fileparts(ans)

cellindic = [5 10 15 20 23 40 44 48 60 103 120 149 163 174 197 232 234 237 242 248 249 250 263 264 265 278 281 290 292 293 294 304 329 330 374 381 385 388 403 429 433 448 480 487 493 508 516 518 545 547 550 561 564 580 584 611 629 632 635 645 651 656 659 670 701 720 741 748 753 757 758 762 787 798 814 815 832 858 863 79 162 203 230 449];
cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
        plot_time_courses(datarun00, A, 1, 1)


cellindic = [5 10 15 20 23 40 44 48 60 103 120 149 163 174 197 232 234 237 242 248 249 250 263 264 265 278 281 290 292 293 294 304 329 330 374 381 385 388 403 429 433 448 480 487 493 508 516 518 545 547 550 561 564 580 584 611 629 632 635 645 651 656 659 670 701 720 741 748 753 757 758 762 787 798 814 815 832 858 863 79 162 203 230 449];
cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun01, 'master_cell_type', cellids, 'corr_threshold', 0.8);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun01, A)
figure(2)
plot_rf_portraits(datarun01,A)
figure(3)
plot_rfs(datarun01, A)

%256 > 32 UP


cellindic = [10 40 44 149 163 174 232 281 292 293 294 374 381 385 480 545 629 748 230];

cellindic = [10 40 44 149 163 174 232 281 292 293 294 374 381 385 480 545 629 748 230 72 75 58 73 786]; %256 > 32 PLUS ONLY 256 UP



cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
plot_time_courses(datarun00, A, 1, 1)

%256 > 32 DOWN

cellindic = [162 60 329 388 448 518 547 611 645 651 701 798 863];

cellindic = [162 60 329 388 448 518 547 611 645 651 701 798 863 697]; %256 > 32 + ONLY 256 DOWN


cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
plot_time_courses(datarun00, A, 1, 1)

        
 %256 > 32 right
cellindic = [5 23 237 242 330 429 508 561 564 584 670 720 741 753 787 815 79 203 449];

cellindic = [5 23 237 242 330 429 508 561 564 584 670 720 741 753 787 815 79 203 449]; %ADDED ONLY 256 RIGHT = NONE


cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
plot_time_courses(datarun00, A, 1, 1)

 %256 > 32 LEFT
cellindic = [15 48 103 120 197 234 278 290 304 403 487 516 550 580 632 635 656 659 814 858];

 %ONLY 256 LEFT = NONE


cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
plot_time_courses(datarun00, A, 1, 1)



%both speeds

cellindic = [6 31 45 55 67 98 112 114 115 121 134 146 150 171 180 181 195 196 215 221 226 233 255 258 295 325 358 363 365 371 379 380 409 430 437 447 455 456 458 482 483 497 502 507 517 523 529 539 544 554 560 579 592 619 623 653 658 688 689 690 703 704 709 710 749 751 763 764 788 24 228];
cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(5)
plot_rf_summaries(datarun00, A)
figure(6)
plot_rf_portraits(datarun00,A)
figure(7)
plot_rfs(datarun00, A)
        plot_time_courses(datarun00, A, 1, 1)



cellindic = [6 31 45 55 67 98 112 114 115 121 134 146 150 171 180 181 195 196 215 221 226 233 255 258 295 325 358 363 365 371 379 380 409 430 437 447 455 456 458 482 483 497 502 507 517 523 529 539 544 554 560 579 592 619 623 653 658 688 689 690 703 704 709 710 749 751 763 764 788 24 228];
cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun01, 'master_cell_type', cellids, 'corr_threshold', 0.8);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun01, A)
figure(2)
plot_rf_portraits(datarun01,A)
figure(3)
plot_rfs(datarun01, A)

%down both speeds
cellindic = [6 45 98 121 146 180 221 325 380 437 517 523 539 619]; 


cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0)); 
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
 plot_time_courses(datarun00, A, 1, 1)

%up both speeds
cellindic = [112 114 115 150 195 196 226 233 258 295 371 379 447 455 456 458 502 544 592 653 688 689 690 749 751];
cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
 plot_time_courses(datarun00, A, 1, 1)
 
 
%right both speeds
cellindic = [67 171 358 430 482 560 658 788 228];

cellindic = [67 171 358 430 482 560 658 788 228 238]; %ADDED 32 >256 RIGHT


cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
 plot_time_courses(datarun00, A, 1, 1)
 
 %left both speeds
cellindic = [24 31 55 134 181 215 255 554 579 623 703 704 709 710 763 764];

cellindic = [24 31 55 134 181 215 255 554 579 623 703 704 709 710 763 764 760]; %ADDED 32 >256 LEFT


cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
 plot_time_courses(datarun00, A, 1, 1)


%Sustained non DS stimuli
cellindic = [59 100 102 184 286 306 488 510 537 634 759 791 820 843];
cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.8);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
        plot_time_courses(datarun00, A, 1, 1)


 

        

cellindic = [59 100 102 184 286 306 488 510 537 634 759 791 820 843];
cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun01, 'master_cell_type', cellids, 'corr_threshold', 0.8);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(2)
plot_rf_summaries(datarun01, A)




%Sustained DS stimuli

cellindic = [2 151 161 210 213 214 289 372 378 432 513 707 711 717];
cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.8);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
        plot_time_courses(datarun00, A, 1, 1)



cellindic = [2 151 161 210 213 214 289 372 378 432 513 707 711 717];
cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun01, 'master_cell_type', cellids, 'corr_threshold', 0.8);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(2)
plot_rf_summaries(datarun01, A)

%Orientation Selective

cellindic = [30 131 137 158 159 167 168 169 216 219 266 285 315 335 347 415 416 438 468 469 470 512 531 556 602 644 667 691 692 732 828 830];
cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun00, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun00, A)
        plot_time_courses(datarun00, A, 1, 1)

figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)

cellindic = [30 131 137 158 159 167 168 169 216 219 266 285 315 335 347 415 416 438 468 469 470 512 531 556 602 644 667 691 692 732 828 830];cellids = datarun02.cell_ids(1,cellindic);
cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun02, datarun01, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun01, A)
figure(2)
plot_rf_portraits(datarun01,A)
figure(3)
plot_rfs(datarun01, A)
figure(4)
plot_time_courses(datarun01, A, 1, 1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT DATA000 DANIEL'S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CLASSIFICATION


%ON type 3 & type 4

A = [271 1066 1546 2492 2956 3182 3571 4396 5206 6798 7113 1711 2551 4247 5191 5476 6841];
A = [271 1066 1546 2492 2956 3182 3571 4396 5206 6798 1711 2551 4247 5476 6841]; %subset
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
plot_time_courses(datarun00, A, 1, 1)
datarun00 = get_autocorrelations(datarun00, A);
plot_autocorrelograms(datarun00, A, 'foa', 0);


%OFF type 1
A = [4591 5344 6934 664 1426 2581 3107 4067 4741 5446 587 961 1592 2611 3376 4246 4816 5671 7517 1231 1816 2911 3436 4381 5086 6616 7591 1261 2147 2957 3511];
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
plot_time_courses(datarun00, A, 1, 1)
datarun00 = get_autocorrelations(datarun00, A);
plot_autocorrelograms(datarun00, A, 'foa', 0);

%OFF type 2

A = [4921 6196 707 1531 2116 2716 3541 4951 6226 886 1727 2146 2761 3707 5418 541 901 1891 2161 3046 3752 5911 543 1141 1936 2476 3181 4066 5912 7156 1277 181 2641 3514 4412];
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
plot_time_courses(datarun00, A, 1, 1)
datarun00 = get_autocorrelations(datarun00, A);
plot_autocorrelograms(datarun00, A, 'foa', 0);

%OFF type 1 & type 2
A = [4591 664 1426 2581 3107 4067 4741 5446 587 961 1592 2611 4816 5671 7517 1231 1816 2911 3436 4381 5086 6616 1261 2147 2957 3511 4921 6196 1531 2116 2716 4951 6226 886 1727 2146 2761 3707 5418 541 1891 2161 3752 5911 543 1936 2476 3181 4066 5912 7156 1277 181 2641 4412];
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
plot_time_courses(datarun00, A, 1, 1)
datarun00 = get_autocorrelations(datarun00, A);
plot_autocorrelograms(datarun00, A, 'foa', 0);



%subsets of OFF type 1 & type 2
A = [961 1261 1592 2581 2611 2911 2957 3436 3511 4381 4591 4816 5086 5446 5671 6616 ]; %subsets of OFF type 1
A = [181 541 886 1531 1891 2116 2146 2161 2476 2641 2716 3181 4066 4921 4951 5911 6226 7156];%subsets of OFF type 2
A = [961 1261 1592 2581 2611 2911 2957 3436 3511 4381 4591 4816 5086 5446 5671 6616 181 541 886 1531 1891 2116 2146 2161 2476 2641 2716 3181 4066 4921 4951 5911 6226 7156];
%subsets of OFF type 1 & type 2
figure(1)
plot_rf_summaries(datarun00, A)
figure(2)
plot_rf_portraits(datarun00,A)
figure(3)
plot_rfs(datarun00, A)
plot_time_courses(datarun00, A, 1, 1)
datarun00 = get_autocorrelations(datarun00, A);
plot_autocorrelograms(datarun00, A, 'foa', 0);

%%%%%%%%%%%%%%%%%%%%%

cellindic = [10 24 40 48 55 72 75 79 112 113 114 115 134 149 195 226 232 237 292 293 295 304 374 379 502 507 554 659 703 709 749 760  761 763]; %ON
cellindic = [15 20 31 44 60 67 103 146 182 203 221 228 242 250 290 358 409 433 493 497 550 632 635 656 757 762 788 815]; %ON OFF
cellindic = [15 31 103 290 550 632 635 656]; %ON OFF LEFT
cellindic = [67 182 203 228 242 358 788 815]; %ON OFF RIGHT
cellindic = [20 60 146 221 250 409 433 493 497 757 762]; %ON OFF DOWN RIGHT

cellindic = [24 79 134 226 379 502 507 709]; %more selective oN

cellids = datarun02.cell_ids(1,cellindic);
[cell_list_map, failed_cells] = map_ei(datarun, datarun02, 'master_cell_type', cellids, 'corr_threshold', 0.95);
emptyCells = cellfun(@isempty,cell_list_map);
cell_list_map(1,(emptyCells==0))
A = cell_list_map(1,(emptyCells==0));
A = cell2mat(A);
figure(1)
plot_rf_summaries(datarun000, A)
figure(2)
plot_rf_portraits(datarun000,A)
figure(3)
plot_rfs(datarun000, A)
plot_time_courses(datarun000, A, 1, 1)
datarun000 = get_autocorrelations(datarun000, A);
plot_autocorrelograms(datarun000, A, 'foa', 0);

