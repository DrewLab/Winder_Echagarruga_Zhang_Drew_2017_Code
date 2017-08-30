function [theimage]=ReadDalsaBinary(thefile, image_height, image_width)

%   Adapted from ReadDalsaBinary2 by Aaron Winder, Drew Lab, ESM, 
%   Penn State University, Nov 2013
%   Version 1

%   SUMMARY: Converts a binary file containing images from a Dalsa 1M60
%   pantera camera into an image with proper orientation.
%_______________________________________________________________
%   INPUTS:
%                       thefile - the name of the binary file including the
%                       extension

%                       image_height - the height of the image in pixels

%                       image_width - the width of the image in pixels
%_______________________________________________________________
%   OUTPUTS:
%                       theimage - a cell array containing a stack of
%                       images
%_______________________________________________________________
%   REQUIRED SCRIPTS:
%_______________________________________________________________
%   CALLED BY:
%_______________________________________________________________
%   FUTURE VERSIONS:
%_______________________________________________________________
%   CHANGES FROM PREV VERS:
%_______________________________________________________________

%imagebasics
pixels_per_frame=image_width*image_height;
image_type='int16';
%open the file , get file size , back to the begining
fid=fopen(thefile)
fseek(fid,0, 'eof')
thefile_size=ftell(fid);
fseek(fid,0, 'bof');

% identify the number of frames to read. Each frame has a previously
% defined width and height (as inputs), along with a grayscale "depth" of
% 2"

nframes_to_read=floor(thefile_size/(2*pixels_per_frame))
%preallocate memory
% blocksize=nframes_to_read;
% maxrow=image_height;
% maxcol=image_width;
% A=zeros(maxrow,maxcol,'int16');
theimage=cell(1,nframes_to_read);
for n=1:nframes_to_read
    z=fread(fid, pixels_per_frame,'*int16','b');
    img=reshape(z(1:pixels_per_frame),image_height,image_width);
    theimage{n} = rot90(img',2);
%     theimage{n}=A;
end

fclose('all');