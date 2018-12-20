function PlotConfocalImageG230(BaseFileName, verbose)

% PlotConfocalImageG230
%
%   Reads ics file generated by Andor camera and generates 4 outputs:
%       (1) movie of fluorescence channel
%       (2) combined movie of fluorescence and bright field
%       (3) single projected image
%       (4) parameters file for recomputing images
%
%   BaseFileName is name of image data file (without extension).
%   verbose option prompts for start/end frames to compute projected image

green = gray;
green(:, 1) = 0;
green(:, 3) = 0;
colormap(green);

%RawImageFolderPath = '/users/maura/data/images/raw/';
%ImageFolderPath = '/users/maura/data/images/';
RawImageFolderPath = '~/Data/images/';
ImageFolderPath = '~/Data/images/';

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
ImageParameters = sscanf(TextLine(14:length(TextLine)), '%f %f %f %f');
XSize = ImageParameters(2);
YSize = ImageParameters(3);
ZSize = ImageParameters(4);

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
    ImageParameters.CompressionExponent = 0.25;     % compression exponent for projected image
end
if (~exist('ImageParameters') || ~isfield(ImageParameters, 'MaxSTD'))
    ImageParameters.MaxSTD = 3.5;                   % max value for projected image
end
if (~exist('ImageParameters') || ~isfield(ImageParameters, 'MinSTD'))
    ImageParameters.MinSTD = 0.75;                  % min value for projected image
end
if (~exist('ImageParameters') || ~isfield(ImageParameters, 'ZStackCompressionExponent'))
    ImageParameters.ZStackCompressionExponent = 0.25;    % compression exponent for zstack
end
if (~exist('ImageParameters') || ~isfield(ImageParameters, 'BFCompressionExponent'))
    ImageParameters.BFCompressionExponent = 2;    % compression exponent for bf image
end
if (~exist('ImageParameters') || ~isfield(ImageParameters, 'ZStackMaxSTD'))
    ImageParameters.ZStackMaxSTD = 3.5;               % max for zstack
end
if (~exist('ImageParameters') || ~isfield(ImageParameters, 'ZStackMinSTD'))
    ImageParameters.ZStackMinSTD = 0.75;            % min for zstack
end

%**************************************************************************
% open z-stack file and read images into matrix
%**************************************************************************
ImageFileName = strcat(RawImageFolderPath, strcat(BaseFileName, '.ids'));
fid = fopen(ImageFileName, 'r');
BrightFieldFileName = strcat(RawImageFolderPath, strcat(BaseFileName, 'brightfield.ids'));
[bffid, errmsg] = fopen(BrightFieldFileName, 'r');
clear im;
for r = 1:ZSize
    im(:, :, r) = fread(fid, [XSize YSize], 'uint16', 0, 'l');
end
fclose(fid);

if (isempty(errmsg))
    clear bfim;
    for r = 1:ZSize
        bfim(:, :, r) = fread(bffid, [XSize YSize], 'uint16', 0, 'l');
    end
    fclose(bffid);
end

%**************************************************************************
% make movie of z-stack
%**************************************************************************
clear ZStackMovie
for r = 1:ZSize
    tempim = im(:, :, r);
    tempim = tempim.^ImageParameters.ZStackCompressionExponent;
    tempim = tempim - min(min(tempim));
    tempim = tempim ./ max(max(tempim));
    ThresholdValue = max([0 (mean(tempim(:)) - ImageParameters.ZStackMinSTD * std(tempim(:)))]);
    MaxValue = min([1 (mean(tempim(:)) + ImageParameters.ZStackMaxSTD * std(tempim(:)))]);
    indices = find(tempim > MaxValue);
    tempim(indices) = MaxValue;
    tempim = tempim / MaxValue;
    thresh = mean(tempim(:, 1)) * ThresholdValue;
    indices = find(tempim < thresh);
    tempim(indices) = thresh;
    tempim = tempim - thresh;
    tempim = tempim ./ max(tempim(:));
    tempim = tempim * 63 + 1;
    ZStackMovie(r) = im2frame(tempim, green);
end
scrsz = get(0, 'ScreenSize');
figure(1); clf;
set(1, 'Position', [1 scrsz(4)/2 scrsz(3)/2.2 scrsz(4)/1.6]);
movie(ZStackMovie, 1)

%**************************************************************************
% collapse into single frame
%**************************************************************************
clear tempim
fprintf(1, '%d frames in z-stack\n', ZSize);
if (isfield(ImageParameters, 'StartFrame'))
    fprintf(1, 'using frames %d to %d for projected image\n', ImageParameters.StartFrame, ImageParameters.EndFrame);
else
    ImageParameters.StartFrame = 1;
    ImageParameters.EndFrame = ZSize;
end
if (~isfield(ImageParameters, 'StartFrame') && verbose == 0)
    ImageParameters.StartFrame = input('Start frame for projected image ');
    ImageParameters.EndFrame = input('End frame for projected image ');
end
EndFlag = 0;

% if verbose option chosen, repeat this look until satisfied with projected
% image
while (EndFlag == 0)
    for r = ImageParameters.StartFrame:ImageParameters.EndFrame
        tempim = im(:, :, r);
        tempim = tempim.^ImageParameters.CompressionExponent;
        tempim = tempim - min(min(tempim));
        tempim = tempim ./ max(max(tempim));
        ThresholdValue = max([0 (mean(tempim(:)) - ImageParameters.MinSTD * std(tempim(:)))]);
        MaxValue = min([1 (mean(tempim(:)) + ImageParameters.MaxSTD * std(tempim(:)))]);
        indices = find(tempim > MaxValue);
        tempim(indices) = MaxValue;
        tempim = tempim / MaxValue;
        thresh = mean(tempim(:, 1)) * ThresholdValue;
        indices = find(tempim < thresh);
        tempim(indices) = thresh;
        tempim = tempim - thresh;
        if (r == ImageParameters.StartFrame)
            ProjectedImage = tempim;
        else
            ProjectedImage = ProjectedImage + tempim;
        end
    end
    ProjectedImage = ProjectedImage / (ImageParameters.EndFrame - ImageParameters.StartFrame + 1);
    clf;
    imagesc(ProjectedImage, [0 1]);
    colormap(green);
    axis square
    if (verbose == 0)
        EndFlag = 1;
    else
        EndFlag = input('ok (enter 1) or change start/end (enter 0)?');
    end
    if (EndFlag == 0 && verbose)
        fprintf(1, 'using frames %d to %d compression exponent = %d\n', ImageParameters.StartFrame, ImageParameters.EndFrame, ImageParameters.CompressionExponent);
        fprintf(1, 'min = %d max = %d\n', ImageParameters.MinSTD, ImageParameters.MaxSTD);
        ImageParameters.StartFrame = input('Start frame for projected image ');
        ImageParameters.EndFrame = input('End frame for projected image ');
        ImageParameters.CompressionExponent = input('Compression exponent ');
        ImageParameters.MinSTD = input('How many st devs to go below mean ');
        ImageParameters.MaxSTD = input('How many st devs to go above mean ');       
    end    
end

if (isempty(errmsg))
    clear tempim
    for r = 1:ZSize
        tempim = bfim(:, :, r);
        tempim = tempim.^ImageParameters.BFCompressionExponent;
        tempim = tempim - min(min(tempim));
        tempim = tempim ./ max(max(tempim));
        if (r == 1)
            ProjectedImage2 = tempim;
        else
            ProjectedImage2 = ProjectedImage2 + tempim;
        end
    end
    ProjectedImage2 = ProjectedImage2 / ZSize;
    CombinedProjectedImage(:, :, 1) = ProjectedImage;
    CombinedProjectedImage(:, :, 2) = ProjectedImage;
    CombinedProjectedImage(:, :, 3) = ProjectedImage2;
    imagesc(CombinedProjectedImage, [0 1]);
    if (isfield(ImageParameters, 'BFBrightness'))
        BFBrightness = ImageParameters.BFBrightness;
    else
        BFBrightness = 6 * mean(im(:))/mean(bfim(:));
    end
end

%**************************************************************************
% make combined movie 
%**************************************************************************
FloBrightness = 1;
if (isempty(errmsg))
    EndFlag = 0;
    while (EndFlag == 0)
        clear BFZStackMovie CombinedZStackMovie colorimage;
        for r = 1:ZSize
            tempim = im(:, :, r);
            tempim = tempim.^ImageParameters.ZStackCompressionExponent;
            tempim = tempim - min(min(tempim));
            tempim = tempim ./ max(max(tempim));
            tempim = tempim * FloBrightness;
            ThresholdValue = max([0 (mean(tempim(:)) - ImageParameters.ZStackMinSTD * std(tempim(:)))]);
            MaxValue = min([1 (mean(tempim(:)) + ImageParameters.ZStackMaxSTD * std(tempim(:)))]);
            indices = find(tempim > MaxValue);
            tempim(indices) = MaxValue;
            tempim = tempim / MaxValue;
            thresh = mean(tempim(:, 1)) * ThresholdValue;
            indices = find(tempim < thresh);
            tempim(indices) = thresh;
            tempim = tempim - thresh;
            tempim = tempim ./ max(tempim(:));
            tempim = tempim * 63 + 1;
            tempim2 = bfim(:, :, r);
            tempim2 = tempim2.^ImageParameters.BFCompressionExponent;
            tempim2 = tempim2 - min(min(tempim2));
            tempim2 = tempim2 * BFBrightness ./ max(max(tempim2));
            indices = find(tempim2 > 1);
            tempim2(indices) = 1;
            tempim2 = tempim2 * 63 + 1;
            BFZStackMovie(r) = im2frame(tempim2, gray);
            tempim = tempim ./ max(max(tempim));
            tempim2 = tempim2 ./ max(max(tempim2));
            colorimage(:, :, 1) = tempim;
            colorimage(:, :, 2) = tempim;
            colorimage(:, :, 3) = tempim2;
            CombinedZStackMovie(r) = im2frame(colorimage);
        end
        figure(2); clf;
        set(2, 'Position', [scrsz(3)/2 scrsz(4)/2 scrsz(3)/2.2 scrsz(4)/1.6]);
        movie(BFZStackMovie, 1)
        movie(CombinedZStackMovie, 1)
        if (verbose == 0)
            EndFlag = 1;
        else
            EndFlag = input('ok (enter 1) or change start/end (enter 0)?');
        end
        if (EndFlag == 0 && verbose)
            fprintf(1, 'bf brightness = %d flo brightness = %d\nflo compression exponent = %d bf compression exponent = %d min = %d max = %d\n', BFBrightness, FloBrightness, ImageParameters.ZStackCompressionExponent, ImageParameters.BFCompressionExponent, ImageParameters.ZStackMinSTD, ImageParameters.ZStackMaxSTD);
            BFBrightness = input('Bright field brightness ');
            FloBrightness = input('fluorescence image brightness ');
            ImageParameters.ZStackCompressionExponent = input('Flo Compression exponent ');
            ImageParameters.BFCompressionExponent = input('BF Compression exponent ');
            ImageParameters.ZStackMinSTD = input('How many st devs to go below mean ');
            ImageParameters.ZStackMaxSTD = input('How many st devs to go above mean ');       
        end    
    end    
end

%**************************************************************************
% save projection, avi movie and parameters
%**************************************************************************
ImageFileName = strcat(ImageFolderPath, BaseFileName);
imwrite(ProjectedImage, strcat(ImageFileName, '.jpg'), 'jpeg', 'Quality', 100);    
movie2avi(ZStackMovie, ImageFileName);
if (isempty(errmsg))
    ImageParameters.BFBrightness = BFBrightness;
    ImageFileName = strcat(ImageFolderPath, strcat(BaseFileName, 'combined'));
    imwrite(CombinedProjectedImage, strcat(ImageFileName, '.jpg'), 'jpeg', 'Quality', 100);    
    movie2avi(CombinedZStackMovie, ImageFileName);
    ImageFileName = strcat(ImageFolderPath, strcat(BaseFileName, 'bf'));
    movie2avi(BFZStackMovie, ImageFileName);    
end
ImageParametersFileName = strcat(ImageFolderPath, strcat(BaseFileName, 'parameters'));
save(ImageParametersFileName, 'ImageParameters');
