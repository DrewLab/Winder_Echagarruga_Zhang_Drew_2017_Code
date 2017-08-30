function [image_grad]=ReadBinFile_U8MatrixGradient(filename, height, width)
%%  [image_grad]=ReadBinFile_U8MatrixGradient(filename, height, width)
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   SUMMARY: Reads binary file with little endian ordering and 8-bit
%   grayscale formatting, calculates the gradient of the image to
%   emphasize the edges and outputs a cell array of images.
%_______________________________________________________________
%   PARAMETERS:
%               filename - [string] complete file name including extension
%
%               height - [integer] the height of the recorded image in
%               pixels
%
%               width - [integer]  the width of the recorded image in
%               pixels
%
%   RETURN:
%               imageout - [array, u8 integer] the measured intensity
%               values organized as [image width, image height, frame
%               number]
%_______________________________________________________________

% Calculate pixels per frame for fread
pixels_per_frame=width*height;

% open the file , get file size , back to the begining
fid=fopen(filename);
fseek(fid,0, 'eof');
thefile_size=ftell(fid);
fseek(fid,0, 'bof');

% Identify the number of frames to read. Each frame has a previously
% defined width and height (as inputs), U8 has a depth of 1.
nframes_to_read=floor(thefile_size/(pixels_per_frame));
display(['ReadBinFile_U8MatrixGradient: ' num2str(nframes_to_read)...
    ' frames to read.'])

% PreAllocate
image_grad = int8(zeros(width,height,nframes_to_read));
for n=1:nframes_to_read
    z=fread(fid,pixels_per_frame,'*uint8',0,'l');
    ind_img=reshape(z(1:pixels_per_frame),width,height);
    image_grad(:,:,n) = int8(gradient(double(ind_img)));
end
fclose(fid);
