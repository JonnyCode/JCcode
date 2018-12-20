function Tform = map_Ei_to_camera(ArrayImagePath, DataPath)

% Creates appropriate coordinate transform, 'Tform' that will transform eis into monitor space 
% hacked from map_array_from_BW_image
% JC 2017-12-07

%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);
data = load_data(DataPath, opt);

im_array = imread(ArrayImagePath);

figure
imshow(im_array);

%% get array location in display coordinates

%                 EI                               DISPLAY
%
%               195 (1)                         386(5)  264(6)
%                 / \                               ______
%               /     \                            /      \
%   264 (6)    |       |    126 (2)               /        \
%   386 (5)    |       |    4   (3)       455(4)  \        / 195(1)
%               \     /                            \      /
%                 \ /                               ------
%                455 (4)                          4(3)  126(2)
%             +____                  
%     129 (1) |    + 249(2)
%             |    | 
%             |    |
%     512 (4) +____| 392 (3) 
%                  +

disp('select first 4 array corner electrodes (see image in function)') 
array_location_display = ginput;

% get array location in ei coordinates
elec_corner = [195 126 4 455]; % dense array
%elec_corner = [129 249 392 512]; % sparse array
array_location_ei = data.ei.position(elec_corner,:);
Tform = maketform('projective', array_location_ei, array_location_display);

disp('test should be near zero')
test = tformfwd(Tform, array_location_ei)-array_location_display % should be equal or close to zeros
