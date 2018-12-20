function PlotConfocalImage(BaseFileName, MaxProjectionFlag, XYFlag, LocalFlag, verbose)

% PlotConfocalImage
%
%   Reads ics file generated by Nikon confocol and generates 4 outputs:
%       (1) movie of fluorescence channel #1
%       (2) combined movie of both fluorescence channels
%       (3) single projected image
%       (4) parameters file for recomputing images
%
%   BaseFileName is name of image data file (without extension).
%   verbose option prompts for start/end frames to compute projected image

green = gray;
green(:, 1) = 0;
green(:, 3) = 0;
colormap(green);

if (LocalFlag)
    RawImageFolderPath = '~/data/images/raw/';
    ImageFolderPath = '~/data/images/';
else
    [status, groupsPath] = unix('df | grep /Volumes/Groups | grep ^afp | awk "{print \$6}"');
    RawImageFolderPath = strcat(groupsPath, '/Lab Staff/Images/raw/confocal/to-process-and-backup/');
    ImageFolderPath = strcat(groupsPath, '/Lab Staff/Images/processed/confocal/');
end

%**************************************************************************
% open ics parameters file and extract key information
%**************************************************************************
ImageFileName = strcat(RawImageFolderPath, strcat(BaseFileName, '.ics'));
icsfid = fopen(ImageFileName, 'r');

% get x, y, z size of zstack
TextLine = fgetl(icsfid);
TargetTextLine1 = 'layout';
TargetTextLine2 = 'sizes';
while (isempty(strfind(TextLine, TargetTextLine1)) | isempty(strfind(TextLine, TargetTextLine2)))
    TextLine = fgetl(icsfid);
    if (TextLine == -1)
        break;
    end
end
ImageParameters = sscanf(TextLine(14:length(TextLine)), '%f %f %f %f %f');
NumChannels = ImageParameters(2);
XSize = ImageParameters(3);
YSize = ImageParameters(4);
ZSize = ImageParameters(5);

% get x, y, z size of zstack
TextLine = fgetl(icsfid);
TargetTextLine1 = 'history';
TargetTextLine2 = 'PassCount';
while (isempty(strfind(TextLine, TargetTextLine1)) | isempty(strfind(TextLine, TargetTextLine2)) | strfind(TextLine, TargetTextLine2) > 20)
    TextLine = fgetl(icsfid);
    if (TextLine == -1)
        break;
    end
end
PassCount = sscanf(TextLine(18:length(TextLine)), '%f');
gain2 = 0;
if (PassCount > 1)
    gain2 = 1;
end

%**************************************************************************
% open stored image file parameters if applicable
%**************************************************************************
cd(ImageFolderPath);
clear ImageParameters
ImageParametersFileName = strcat(BaseFileName, 'parameters.mat');
if (exist(ImageParametersFileName))
    load(ImageParametersFileName)
    fprintf(1, 'Using stored image parameters\n');
end
if (~exist('ImageParameters') || ~isfield(ImageParameters, 'CompressionExponent'))
    ImageParameters.CompressionExponent = 0.2;     % compression exponent for projected image
end
if (~exist('ImageParameters') || ~isfield(ImageParameters, 'MaxValue'))
    ImageParameters.MaxValue = 1;                   % max value for projected image
end
if (~exist('ImageParameters') || ~isfield(ImageParameters, 'MinValue'))
    ImageParameters.MinValue = 0.25;                  % min value for projected image
end
if (~exist('ImageParameters') || ~isfield(ImageParameters, 'ZStackCompressionExponent'))
    ImageParameters.ZStackCompressionExponent = 0.25;    % compression exponent for zstack
end
if (~exist('ImageParameters') || ~isfield(ImageParameters, 'ZStackMaxValue'))
    ImageParameters.ZStackMaxValue = 1.0;               % max for zstack
end
if (~exist('ImageParameters') || ~isfield(ImageParameters, 'ZStackMinValue'))
    ImageParameters.ZStackMinValue = 0.15;            % min for zstack
end

%**************************************************************************
% open z-stack file and read images into matrix
%**************************************************************************
ImageFileName = strcat(RawImageFolderPath, strcat(BaseFileName, '.ids'));
fid = fopen(ImageFileName, 'r');
clear im im2;
fseek(fid, 4, 'bof');
for r = 1:ZSize
    im(:, :, r) = fread(fid, [XSize YSize], 'uint16', 2*(NumChannels-1), 'l');
end
if (gain2 > 0)
    fseek(fid, 2, 'bof');
    for r = 1:ZSize
        im2(:, :, r) = fread(fid, [XSize YSize], 'uint16', 2*(NumChannels-1), 'l');
    end
end
fclose(fid);

%**************************************************************************
% make movie of z-stack
%**************************************************************************
clear ZStackMovie
for r = 1:ZSize
    tempim = im(:, :, r);
    tempim = tempim.^ImageParameters.ZStackCompressionExponent;
    tempim = tempim - min(min(tempim));
    tempim = tempim ./ max(max(tempim));
    indices = find(tempim > ImageParameters.ZStackMaxValue);
    tempim(indices) = ImageParameters.ZStackMaxValue;
    indices = find(tempim < ImageParameters.ZStackMinValue);
    tempim(indices) = ImageParameters.ZStackMinValue;
    tempim = tempim - ImageParameters.ZStackMinValue;
    tempim = (tempim/max(tempim(:))) * 63 + 1;
    scim(:, :, r) = tempim;
    ZStackMovie(r) = im2frame(tempim, green);
end
scrsz = get(0, 'ScreenSize');
figure(1); clf;
set(1, 'Position', [1 scrsz(4)/2 scrsz(3)/2.2 scrsz(4)/1.6]);
movie(ZStackMovie, 1)

if (gain2 > 0)
    clear ZStackMovie2
    for r = 1:ZSize
        tempim = im2(:, :, r);
        tempim = tempim.^ImageParameters.ZStackCompressionExponent;
        tempim = tempim - min(min(tempim));
        tempim = tempim ./ max(max(tempim));
        indices = find(tempim > ImageParameters.ZStackMaxValue);
        tempim(indices) = ImageParameters.ZStackMaxValue;
        indices = find(tempim < ImageParameters.ZStackMinValue);
        tempim(indices) = ImageParameters.ZStackMinValue;
        tempim = tempim - ImageParameters.ZStackMinValue;
        tempim = (tempim/max(tempim(:))) * 63 + 1;
        scim2(:, :, r) = tempim;
        ZStackMovie2(r) = im2frame(tempim, green);
    end
    scrsz = get(0, 'ScreenSize');
    figure(1); clf;
    set(1, 'Position', [1 scrsz(4)/2 scrsz(3)/2.2 scrsz(4)/1.6]);
    movie(ZStackMovie2, 1)
end

%**************************************************************************
% collapse into single frame
%**************************************************************************
clear tempim
fprintf(1, '%d frames in z-stack\n', ZSize);
ImageParameters.StartFrame = 1;
ImageParameters.EndFrame = ZSize;

clear ProjectedImage ScaledProjectedImage ProjectedImage2 ScaledProjectedImage2 CombinedProjectedImage;

% channel 1
if (XYFlag)
    if (MaxProjectionFlag)
        for x = 1:size(scim, 1)
            for y = 1:size(scim, 2)
                ProjectedImage(x, y) = max(im(x, y, :));
            end
        end
    else
        for r = ImageParameters.StartFrame:ImageParameters.EndFrame
            tempim = scim(:, :, r);
            if (r == ImageParameters.StartFrame)
                ProjectedImage = tempim;
            else
                ProjectedImage = ProjectedImage + tempim;
            end
        end
        ProjectedImage = ProjectedImage / (ImageParameters.EndFrame - ImageParameters.StartFrame + 1);
    end
else
    if (MaxProjectionFlag)
        for x = 1:size(scim, 1)
            for z = 1:size(scim, 3)
                ProjectedImage(x, z) = max(im(x, :, z));
            end
        end
    else
        for r = 1:size(im, 2);
            tempim = squeeze(im(r, :, :));
              if (r == 1)
                ProjectedImage = tempim;
            else
                ProjectedImage = ProjectedImage + tempim;
            end
        end
    end
    ProjectedImage = ProjectedImage';
end

% if verbose option chosen, repeat until satisfied with projected image
EndFlag = 0;
while (EndFlag == 0)
    ScaledProjectedImage = (ProjectedImage - min(min(ProjectedImage))) / (max(max(ProjectedImage)) - min(min(ProjectedImage)));
    ScaledProjectedImage = ScaledProjectedImage.^ImageParameters.CompressionExponent - ImageParameters.MinValue;
    ScaledProjectedImage = ScaledProjectedImage / (max(max(ScaledProjectedImage)));
    ScaledProjectedImage = ScaledProjectedImage / ImageParameters.MaxValue;
    indices = find(ScaledProjectedImage(:) < 0);
    ScaledProjectedImage(indices) = 0;
    indices = find(ScaledProjectedImage(:) > 1);
    ScaledProjectedImage(indices) = 1;
    clf;
    imagesc(ScaledProjectedImage, [0 1]);
    colormap(green);
    axis equal
    fprintf(1, '\n**********CHANNEL 1**********\n');
    if (verbose == 0)
        EndFlag = 1;
    else
        EndFlag = input('\nok (enter 1) or change parameters (enter 0)?');
    end
    if (~isnumeric(EndFlag))
        EndFlag = 0;
    end
    if (EndFlag == 0 && verbose)
        fprintf(1, 'CURRENT: compression exponent = %d min = %d max = %d\n\n', ImageParameters.CompressionExponent, ImageParameters.MinValue, ImageParameters.MaxValue);
        ImageParameters.CompressionExponent = input('Enter new compression exponent ');
        ImageParameters.MinValue = input('Enter new minimum value ');
        ImageParameters.MaxValue = input('Enter new maximum value ');       
    end    
end

% channel 2
clear ProjectedImage;
if (gain2 > 0)

    if (XYFlag)
        if (MaxProjectionFlag)
            for x = 1:size(scim2, 1)
                for y = 1:size(scim2, 2)
                    ProjectedImage(x, y) = max(im2(x, y, :));
                end
            end
        else
            for r = ImageParameters.StartFrame:ImageParameters.EndFrame
                tempim = scim2(:, :, r);
                if (r == ImageParameters.StartFrame)
                    ProjectedImage = tempim;
                else
                    ProjectedImage = ProjectedImage + tempim;
                end
            end
            ProjectedImage = ProjectedImage / (ImageParameters.EndFrame - ImageParameters.StartFrame + 1);
        end
    else
        if (MaxProjectionFlag)
            for x = 1:size(scim2, 1)
                for z = 1:size(scim2, 3)
                    ProjectedImage(x, z) = max(im2(x, :, z));
                end
            end
        else
            for r = 1:size(im2, 2);
                tempim = squeeze(im2(r, :, :));
                  if (r == 1)
                    ProjectedImage = tempim;
                else
                    ProjectedImage = ProjectedImage + tempim;
                end
            end
        end
        ProjectedImage = ProjectedImage';
    end
    
    % if verbose option chosen, repeat until satisfied with projected image
    EndFlag = 0;
    while (EndFlag == 0)
        ScaledProjectedImage2 = (ProjectedImage - min(min(ProjectedImage))) / (max(max(ProjectedImage)) - min(min(ProjectedImage)));
        ScaledProjectedImage2 = ScaledProjectedImage2.^ImageParameters.CompressionExponent - ImageParameters.MinValue;
        ScaledProjectedImage2 = ScaledProjectedImage2 / (max(max(ScaledProjectedImage2)));
        ScaledProjectedImage2 = ScaledProjectedImage2 / ImageParameters.MaxValue;
        indices = find(ScaledProjectedImage2(:) < 0);
        ScaledProjectedImage2(indices) = 0;
        indices = find(ScaledProjectedImage2(:) > 1);
        ScaledProjectedImage2(indices) = 1;
        clf;
        imagesc(ScaledProjectedImage2, [0 1]);
        colormap(green);
        axis square
        fprintf(1, '\n**********CHANNEL 2**********\n');
        if (verbose == 0)
            EndFlag = 1;
        else
            EndFlag = input('\nok (enter 1) or change parameters (enter 0)?');
        end
        if (EndFlag == 0 && verbose)
            fprintf(1, 'CURRENT: compression exponent = %d min = %d max = %d\n\n', ImageParameters.CompressionExponent, ImageParameters.MinValue, ImageParameters.MaxValue);
            ImageParameters.CompressionExponent = input('Enter new compression exponent ');
            ImageParameters.MinValue = input('Enter new minimum value ');
            ImageParameters.MaxValue = input('Enter new maximum value ');       
        end    
    end

    CombinedProjectedImage(:, :, 1) = ScaledProjectedImage;
    CombinedProjectedImage(:, :, 2) = ScaledProjectedImage;
    CombinedProjectedImage(:, :, 3) = ScaledProjectedImage2;
    imagesc(CombinedProjectedImage, [0 1]);

    if (isfield(ImageParameters, 'BFBrightness'))
        BFBrightness = ImageParameters.BFBrightness;
    else
        BFBrightness = 6 * mean(im(:))/mean(im2(:));
    end
end
ImageFileName = strcat(ImageFolderPath, BaseFileName);
imwrite(ScaledProjectedImage, strcat(ImageFileName, '.jpg'), 'jpeg', 'Quality', 100);    
movie2avi(ZStackMovie, ImageFileName);
clear ZStackMovie
if (gain2 > 0)
    ImageFileName = strcat(ImageFolderPath, strcat(BaseFileName, 'ch2'));
    imwrite(ScaledProjectedImage2, strcat(ImageFileName, '.jpg'), 'jpeg', 'Quality', 100);    
    movie2avi(ZStackMovie2, ImageFileName);
    clear ZStackMovie2;
end

%**************************************************************************
% make combined movie 
%**************************************************************************
if (gain2 > 0)
    FloBrightness = 1;
    EndFlag = 0;
    while (EndFlag == 0)
        clear CombinedZStackMovie colorimage;
        for r = 1:ZSize
            tempim = im(:, :, r);
            tempim = tempim.^ImageParameters.ZStackCompressionExponent;
            tempim = tempim - min(min(tempim));
            tempim = tempim ./ max(max(tempim));
            indices = find(tempim > ImageParameters.ZStackMaxValue);
            tempim(indices) = ImageParameters.ZStackMaxValue;
            indices = find(tempim < ImageParameters.ZStackMinValue);
            tempim(indices) = ImageParameters.ZStackMinValue;
            tempim = tempim - ImageParameters.ZStackMinValue;
            tempim = tempim * FloBrightness ./ max(tempim(:));
            indices = find(tempim > 1);
            tempim(indices) = 1;
            tempim = tempim * 63 + 1;
            tempim2 = im2(:, :, r);
            tempim2 = tempim2.^ImageParameters.ZStackCompressionExponent;
            tempim2 = tempim2 - min(min(tempim2));
            tempim2 = tempim2 ./ max(max(tempim2));
            indices = find(tempim2 > ImageParameters.ZStackMaxValue);
            tempim2(indices) = ImageParameters.ZStackMaxValue;
            indices = find(tempim2 < ImageParameters.ZStackMinValue);
            tempim2(indices) = ImageParameters.ZStackMinValue;
            tempim2 = tempim2 - ImageParameters.ZStackMinValue;
            tempim2 = tempim2 * BFBrightness ./ max(tempim2(:));
            indices = find(tempim2 > 1);
            tempim2(indices) = 1;
            tempim2 = tempim2 * 63 + 1;
            tempim = tempim ./ max(max(tempim));
            tempim2 = tempim2 ./ max(max(tempim2));
            colorimage(:, :, 1) = tempim2;
            colorimage(:, :, 2) = tempim;
            colorimage(:, :, 3) = tempim2;
            CombinedZStackMovie(r) = im2frame(colorimage);
        end
        figure(2); clf;
        set(2, 'Position', [scrsz(3)/2 scrsz(4)/2 scrsz(3)/2.2 scrsz(4)/1.6]);
        movie(CombinedZStackMovie, 1)
        fprintf(1, '\n**********COMBINED CHANNELS**********\n');
        if (verbose == 0)
            EndFlag = 1;
        else
            EndFlag = input('\nok (enter 1) or change relative brightness (enter 0)?');
        end
        if (EndFlag == 0 && verbose)
            fprintf(1, 'CURRENT: channel 1 brightness (cell) = %d channel 2 brightness (background) = %d\n\n', BFBrightness, FloBrightness);
            BFBrightness = input('Enter new channel 1 brightness ');
            FloBrightness = input('Enter new channel 2 brightness ');
        end    
    end    
    ImageParameters.BFBrightness = BFBrightness;
    ImageParameters.FloBrightness = FloBrightness;
end

%**************************************************************************
% save projection and avi movie
%**************************************************************************
if (gain2 > 0)
    ImageFileName = strcat(ImageFolderPath, strcat(BaseFileName, 'combined'));
    imwrite(CombinedProjectedImage, strcat(ImageFileName, '.jpg'), 'jpeg', 'Quality', 100);    
    movie2avi(CombinedZStackMovie, ImageFileName);
end
ImageParametersFileName = strcat(ImageFolderPath, strcat(BaseFileName, 'parameters'));
save(ImageParametersFileName, 'ImageParameters');

