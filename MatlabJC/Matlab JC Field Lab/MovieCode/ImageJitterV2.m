function Output = ImageJitterV2(Image, AxisShiftStd, varargin) 

% This function will take a mat file "Image" and make a movie where each
% (NumFrames)frame is a shifted version of the original Image.  "AxisShiftStd" (Image
% stixels) std of vert-horizontal image location distribution (or delta - see DeltaFlag).

% adapted from 'ImageJitter.m' to add rotation and used imtranslate instead
% of ad hoc matrix shift

% JC 12/12/2016

p = inputParser;
p.addParamValue('ImageRotateStd', 0, @isnumeric) ; % (degrees) std of rotation (default is no rotation)
p.addParamValue('Seed', [], @isnumeric) ; % seed to determine random x,y jitter
p.addParamValue('AxisShiftMax', [], @isnumeric) ; % maximum Axis shift
p.addParamValue('MaskFlag',false, @islogical) ; % this flag provides a mask (caution: makes Output.mov larger)
p.addParamValue('MaskSize', [], @isnumeric) ; % number of stixels on each side (defaults to as large as necessary)
p.addParamValue('MinMovSize', 300, @isnumeric) ; % (pix) minimum size of one movie edge before code stops
p.addParamValue('SmoothNumFrames', 1, @isnumeric) ; % average deta over this number of frames (slows jitter)
p.addParamValue('NumFrames', 60*60*1, @isnumeric) ; % Max number of frames to jitter images (default ~1 minutes) 
p.addParamValue('DeltaFlag', false, @islogical) ; % this flag causes the 'AxisShiftStd' and smoothNumFrames 
% to specify the change in position not the position

p.parse(varargin{:});
params = p.Results;

% image size        
ImageNumX = size(Image,1) ;
ImageNumY = size(Image,2) ;

% construct image jitter vectors
if ~isempty(params.Seed) ;
    rng(params.Seed) ; % set random seed so it is consistent
end

if params.DeltaFlag==0 ; % if you Axis params are defining the position 
    Xdelta_SumTemp = normrnd(0,1,1,params.NumFrames) ; % position (relative to original start)
    Ydelta_SumTemp = normrnd(0,1,1,params.NumFrames) ; 
    Rdelta_SumTemp = normrnd(0,1,1,params.NumFrames) ; % rotation 
    
    if params.SmoothNumFrames>1 ; % smooth jitter vectors
        Xdelta_SumTemp = smooth(Xdelta_SumTemp,params.SmoothNumFrames) ;
        Ydelta_SumTemp = smooth(Ydelta_SumTemp,params.SmoothNumFrames) ;
        Rdelta_SumTemp = smooth(Rdelta_SumTemp,params.SmoothNumFrames) ;
    end
    
    Xdelta_Sum = (Xdelta_SumTemp - mean(Xdelta_SumTemp))*AxisShiftStd/std(Xdelta_SumTemp) ; % sets std as specified
    Ydelta_Sum = (Ydelta_SumTemp - mean(Ydelta_SumTemp))*AxisShiftStd/std(Ydelta_SumTemp) ; 
    Rdelta_Sum = (Rdelta_SumTemp - mean(Rdelta_SumTemp))*params.ImageRotateStd/std(Rdelta_SumTemp) ;

    Xdelta_Sum = round(Xdelta_Sum) ;
    Ydelta_Sum = round(Ydelta_Sum) ;
    
    if ~isempty(params.AxisShiftMax) ;
        Xdelta_Sum(Xdelta_Sum>params.AxisShiftMax) = params.AxisShiftMax ;
        Ydelta_Sum(Ydelta_Sum>params.AxisShiftMax) = params.AxisShiftMax ;
    end

    Xdelta = [Xdelta_Sum(1),diff(Xdelta_Sum)'] ;
    Ydelta = [Ydelta_Sum(1),diff(Ydelta_Sum)'] ;
    Rdelta = [Rdelta_Sum(1),diff(Rdelta_Sum)'] ;
    
elseif params.DeltaFlag==1 ; % if Axis params are defining the change position
    XdeltaTemp = normrnd(0,1,1,params.NumFrames) ; % change in position
    YdeltaTemp = normrnd(0,1,1,params.NumFrames) ;
    RdeltaTemp = normrnd(0,1,1,params.NumFrames) ;

    if params.SmoothNumFrames>1 ; % smooth jitter vectors
        XdeltaTemp = smooth(XdeltaTemp,params.SmoothNumFrames) ;
        YdeltaTemp = smooth(YdeltaTemp,params.SmoothNumFrames) ;
        RdeltaTemp = smooth(RdeltaTemp,params.SmoothNumFrames) ;
    end

    Xdelta = (XdeltaTemp - mean(XdeltaTemp))*AxisShiftStd/std(XdeltaTemp) ; % sets std as specified
    Ydelta = (YdeltaTemp - mean(YdeltaTemp))*AxisShiftStd/std(YdeltaTemp) ;
    Rdelta = (RdeltaTemp - mean(RdeltaTemp))*params.ImageRotateStd/std(RdeltaTemp) ;

    Xdelta = round(Xdelta) ; % make indicies
    Ydelta = round(Ydelta) ;

    if ~isempty(params.AxisShiftMax) ;
        Xdelta(Xdelta>params.AxisShiftMax) = params.AxisShiftMax ;
        Ydelta(Ydelta>params.AxisShiftMax) = params.AxisShiftMax ;
    end
    
    Xdelta_Sum = cumsum(Xdelta) ; % position relative to start
    Ydelta_Sum = cumsum(Ydelta) ;
    Rdelta_Sum = cumsum(Rdelta) ;
end
    
% number of frames that you can actually get before image shifts too much
for f=1:params.NumFrames ;
    if ImageNumX-range(Xdelta_Sum(1:f))<params.MinMovSize | ...
            ImageNumY-range(Ydelta_Sum(1:f))<params.MinMovSize ;
        params.NumFrames = min(f-1,params.NumFrames) ;
    end
end

% contstrain jitter vectors to time before image jitters too much
Xdelta = Xdelta(1:params.NumFrames) ; 
Ydelta = Ydelta(1:params.NumFrames) ;
Rdelta = Rdelta(1:params.NumFrames) ;
Xdelta_Sum = Xdelta_Sum(1:params.NumFrames) ;
Ydelta_Sum = Ydelta_Sum(1:params.NumFrames) ;
Rdelta_Sum = Rdelta_Sum(1:params.NumFrames) ;

% mask
mask = zeros(ImageNumX,ImageNumY) ;
maskPrep = zeros(ImageNumX,ImageNumY) ;
if isempty(params.MaskSize) ; % if no mask is specified
    mask(1:max(Xdelta_Sum),:) = mean(Image(:)) ;
    mask(end+min(Xdelta_Sum):end,:) = mean(Image(:)) ;
    mask(:,1:max(Ydelta_Sum)) = mean(Image(:)) ;
    mask(:,end+min(Ydelta_Sum):end) = mean(Image(:)) ;
    
    maskPrep(max(Xdelta_Sum)+1:end+min(Xdelta_Sum)-1,...
        max(Ydelta_Sum)+1:end+min(Ydelta_Sum)-1) = 1 ;
else
    mask(1:params.MaskSize,:) = mean(Image(:)) ;
    mask(end-params.MaskSize:end,:) = mean(Image(:)) ;
    mask(:,1:params.MaskSize) = mean(Image(:)) ;
    mask(:,end-params.MaskSize:end) = mean(Image(:)) ;
    maskPrep(params.MaskSize+1:end-params.MaskSize-1,...
        params.MaskSize+1:end-params.MaskSize-1) = 1 ;
end

% final mask (to center shifted image)
[mr,mc] = find(maskPrep==1) ;
maskXnum = min(range(mr),600) ;
maskYnum = min(range(mc),800) ;
FinalMaskTemplate = ones(600,800)*mean(Image(:)) ;
maskXstart = ceil((600-maskXnum)/2)+1 ;
maskYstart = ceil((800-maskYnum)/2)+1 ;

maskedMovBoundsX = [min(mr):min(min(mr)+600-1,max(mr)-1)] ;
maskedMovBoundsY = [min(mc):min(min(mc)+800-1,max(mc)-1)] ;

% prep output
Output.params = params ;
Output.Image = Image ;
Output.Xdelta = Xdelta ;
Output.Ydelta = Ydelta ;
Output.Rdelta = Rdelta ;
if params.MaskFlag ; % if there is a mask 
    Output.mov = nan(600,800,params.NumFrames) ;
else
    Output.mov = nan(maskXnum,maskYnum,params.NumFrames) ;
end


for f=1:params.NumFrames ;
    movFrame = Image ;
    movFrame = imtranslate(movFrame,[Xdelta_Sum,Ydelta_Sum]) ; % option ('FillValues',mean(Image(:))
    movFrame = imrotate(movFrame,Rdelta_Sum(f),'crop') ; % rotation,translation order changes shift
    
    if params.MaskFlag ;
        movFrame = movFrame.*maskPrep + mask ; % movie with mask
        Output.mov(:,:,f) = FinalMaskTemplate ;
        Output.mov(maskXstart:maskXstart+maskXnum-1,maskYstart:maskYstart+maskYnum-1,f) = ...
            movFrame(maskedMovBoundsX,maskedMovBoundsY) ; % movie centered
    else
        Output.mov(:,:,f) = movFrame(maskedMovBoundsX,maskedMovBoundsY) ; 
    end
end
    

    
    
    
