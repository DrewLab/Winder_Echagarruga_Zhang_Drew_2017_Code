function [ImageStack] = GetCBVFrameSubset(filename, image_height, image_width, FrameInds)
%   function [ImageStack] = GetCBVFrameSubset(filename, image_height, image_width, FrameInds)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Reads in camera images from a binary file with defined
%   width and height and at defined positions in time.
%   
%_______________________________________________________________
%   PARAMETERS: 
%               filename - [string] the name of the binary file to read
%               from
%
%               image_height - [double] the height of the image in pixels
%
%               image_width - [double] the width of the image in pixels
%
%               FrameInds - [array] the positions of the desired frames in
%               time, as an index                       
%_______________________________________________________________
%   RETURN:                     
%               ImageStack - [array] a 3D array of images                   
%_______________________________________________________________

pixels_per_frame=image_width*image_height; % number to give to fread along with 16-bit pixel depth flag
skipped_pixels = pixels_per_frame*2; % Multiply by two because there are 16 bits (2 bytes) per pixel
NumFrames = length(FrameInds);

%open the file , get file size , back to the begining
fid=fopen(filename);
fseek(fid,0, 'eof')
thefile_size=ftell(fid);
fseek(fid,0, 'bof');

nframesInFile=floor(thefile_size/(2*pixels_per_frame))

% Handle Error
if fid == -1
    error(wraptext(['Error. ReadDalsaBinary_Matrix: fopen.m cannot open '...
        filename]))
end

%preallocate memory
ImageStack = NaN*ones(image_height, image_width, NumFrames);

% Loop over each frame
t1 = tic;
for n=1:NumFrames
    fseek(fid,(FrameInds(n)-1)*skipped_pixels,'bof');
    z=fread(fid, pixels_per_frame,'*int16','b');
    % Convert linear array into a 256x256x1 frame
    img=reshape(z,image_height,image_width);
    
    % Orient the frame so that rostral is up
    ImageStack(:,:,n) = rot90(img',2);
end
elapsed = toc(t1);
wraptext(['GetCBVFrameSubset: Frames acquired in ' num2str(elapsed) ' seconds.'])
fclose('all');