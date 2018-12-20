function ReturnedCellInfo = PlotConfocalImageFile(CellInfo, ImageFolderPath, RawImageFolderPath)

green = gray;
green(:, 1) = 0;
green(:, 3) = 0;
colormap(green);

if (~isfield(CellInfo, 'ProjectedImage'))

    %**************************************************************************
    % open ics parameters file and extract key information
    %**************************************************************************
    ImageFileName = input('ics file name (with extension): ', 's');
    ImageFileName = strcat(RawImageFolderPath, ImageFileName);
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
    XSize = ImageParameters(3);
    YSize = ImageParameters(4);
    ZSize = ImageParameters(5);

    %**************************************************************************
    % open z-stack file and read images into matrix
    %**************************************************************************
    ImageFileName = input('ids file name (with extension): ', 's');
    ImageFileName = strcat(RawImageFolderPath, ImageFileName);
    fid = fopen(ImageFileName, 'r');
    clear im;
    for r = 1:ZSize
        im(:, :, r) = fread(fid, [XSize YSize], 'uint16', 6, 'l');
    end
    fclose(fid);

    %**************************************************************************
    % collapse into single frame
    %**************************************************************************
    clear tempim
    for r = 1:ZSize
        if (r == 1)
            ProjectedImage = im(:, :, r);
        else
            ProjectedImage = ProjectedImage + im(:, :, r);
        end
    end
    ProjectedImage = ProjectedImage / ZSize;

    %**************************************************************************
    % make movie of z-stack
    %**************************************************************************
    for r = 1:ZSize
        tempim = im(:, :, r);
        tempim = tempim.^0.4;
        tempim = tempim - min(min(tempim));
        tempim = tempim ./ max(max(tempim));
        thresh = 0.8;
        indices = find(tempim > thresh);
        tempim(indices) = thresh;
        tempim = tempim / thresh;
        thresh = mean(tempim(:, 1)) * 0.9;
        indices = find(tempim < thresh);
        tempim(indices) = thresh;
        tempim = tempim - thresh;
        tempim = tempim * 63 + 1;
        ZStackMovie(r) = im2frame(tempim, green);
    end


    %**************************************************************************
    % update CellInfo
    %**************************************************************************
    Indices = find(CellInfo.CellFile == '/');
    CellFileName = CellInfo.CellFile(max(Indices)+1:length(CellInfo.CellFile));
    ImageFileName = strcat(CellFileName, 'zstack');
    ImageFileName = strcat(ImageFolderPath, ImageFileName);
    save(ImageFileName, 'ZStackMovie');
    CellInfo.ProjectedImage = ProjectedImage;
    CellInfo.ZStack = ImageFileName;
end

%**************************************************************************
% plot stored zstack and projection
%**************************************************************************
load(CellInfo.ZStack);
scrsz = get(0, 'ScreenSize');
figure(1); clf;
set(1, 'Position', [scrsz(3)/2 scrsz(4)/2 scrsz(3)/2.2 scrsz(4)/1.6]);
movie(ZStackMovie, 1)
clf;
imagesc(CellInfo.ProjectedImage.^0.35, [4 124].^0.3);
colormap(green);
axis square

ReturnedCellInfo = CellInfo;