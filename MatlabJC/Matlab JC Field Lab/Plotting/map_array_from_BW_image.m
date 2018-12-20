function Tform = map_array_from_BW_image(BwImagePath,stixelSize,BwXmlPath,ArrayImagePath, DataPath)

% Creates appropriate coordinate transform, 'Tform' that will transform eis into monitor space 

% 'BwImagePath' (image of BWN frame taken by the camera)
% 'stixelSize' (width of stixel)
% 'BwXmlPath' (xml file for BW noise)
% 'ArrayImagePath'(image of electrodes array, taken in same x,y position as the BW image).
% 'DataPath' (any data run with available ei)

% JC 11/7/2016 created function from XY script,'map_script'

%% get array location in display coordinates

% load stimulus picture taken by camera
im_s = imread(BwImagePath);
%im_s = flipud(im_s) ; % CAUTION: only when cammera image is flipped

% get stimulus frame in display coordinates
mov = get_movie(BwXmlPath, 0, 1);
mov_frame = matrix_scaled_up(squeeze(mov(:,:,1)), stixelSize);

% select control points
disp('select 4 pairs of points')
[movingPoints,fixedPoints] = cpselect(im_s, mov_frame,'Wait', true) ;

%% load data
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);
data = load_data(DataPath, opt);

%% register two images
tform = fitgeotrans(movingPoints, fixedPoints, 'projective');
registered = imwarp(im_s, tform,'OutputView',imref2d(size(mov_frame)));
figure 
imshow(registered);
figure
imshowpair(mov_frame,registered,'blend');

% load array image taken by camera
im_array = imread(ArrayImagePath);

% transform array image into display coordinates
registered_array = imwarp(im_array, tform, 'OutputView', imref2d(size(mov_frame)));
figure
imshow(registered_array);

% get array location in display coordinates

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

disp('select first 4 array corner electrodes (see image in function)') 
array_location_display = ginput;

% get array location in ei coordinates
elec_corner = [195 126 4 455];
array_location_ei = data.ei.position(elec_corner,:);
Tform = maketform('projective', array_location_ei, array_location_display);

disp('test should be near zero')
test = tformfwd(Tform, array_location_ei)-array_location_display % should be equal or close to zeros

