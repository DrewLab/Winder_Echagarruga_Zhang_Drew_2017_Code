function [ImageStack]=ReadDalsaBinary_Matrix(filename, image_height, image_width)
%   function [theimage]=ReadDalsaBinary_Matrix(thefile, image_height, image_width)
%
%   Author: Aaron Winder, adapted from ReadDalsaBinary2 by Patrick Drew
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Converts a binary file containing images from a Dalsa 1M60
%   pantera camera into an 256x256xn stack of images with proper orientation.
%   
%_______________________________________________________________
%   PARAMETERS:             
%                   filename - [string] the name of the saved binary file, 
%                   including the .bin extension
%
%                   image_height - [int] the height of the image in pixels
%
%                   image_width - [int] the width of the image in pixels
%_______________________________________________________________
%   RETURN:                     
%                   ImageStack - [array] 256x256xn array containing 
%                   measured absolute intensity values
%_______________________________________________________________


pixels_per_frame=image_width*image_height;

%open the file , get file size , back to the begining
fid=fopen(filename);

% Handle Error
if fid == -1
    error(wraptext(['Error. ReadDalsaBinary_Matrix: fopen.m cannot open '...
        filename]))
end

fseek(fid,0, 'eof');
thefile_size=ftell(fid);
fseek(fid,0, 'bof');

% identify the number of frames to read. Each frame has a previously
% defined width and height (as inputs), along with a grayscale "depth" of
% 2"
NumFrames=floor(thefile_size/(2*pixels_per_frame));
display(['ReadDalsaBinary_Matrix: ' filename '; '  num2str(NumFrames) ' Frames']);

%preallocate memory
ImageStack = NaN*ones(image_height, image_width, NumFrames);

% Loop over each frame
for n=1:NumFrames
    z=fread(fid, pixels_per_frame,'*int16','b');
    
    % Convert linear array into a 256x256x1 frame
    img=reshape(z(1:pixels_per_frame),image_height,image_width);
    
    % Orient the frame so that rostral is up
    ImageStack(:,:,n) = rot90(img',2);
end

fclose('all');