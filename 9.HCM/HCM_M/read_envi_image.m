%+
% NAME:
%   read_envi_image
%
% PURPOSE:
%   This function reads an ENVI image directly into MATLAB.
%
% CALLING SEQUENCE:
%   image = read_envi_image( filename )
%
% INPUT:
%   filename
%      The ENVI image filename to be read
%
% RETURN VALUE:
%   An array containing the image data. This array will be of the correct
%   data type for the provided image.
%
% MODIFICATION HISTORY:
%   Written by       Sarah Paul and Carl Salvaggio
%   December, 2009   Original Code
%   February, 2011   Code restructuring
%
% DISCLAIMER:
%   This source code is provided "as is" and without warranties as to
%   performance or merchantability. The author and/or distributors of
%   this source code may have made statements about this source code.
%   Any such statements do not constitute warranties and shall not be
%   relied on by the user in deciding whether to use this source code.
%
%   This source code is provided without any express or implied warranties 
%   whatsoever. Because of the diversity of conditions and hardware under
%   which this source code may be used, no warranty of fitness for a
%   particular purpose is offered. The user is advised to test the source
%   code thoroughly before relying on it. The user must assume the entire
%   risk of using the source code.
%-

function image = read_envi_image( filename )

   headerFilename = strcat( filename, '.hdr' );
   [samples, lines, bands, dataType, interleave, byteOrder] = ...
      read_envi_header( headerFilename );

   fid = fopen( filename );
   if fid == -1
      error( 'Error opening image file, bad file ID' );
   end

   if byteOrder == 0
      image = fread( fid, dataType, 'ieee-le' );
   else
      image = fread( fid, dataType, 'ieee-be' );
   end

   interleave( isspace( interleave ) ) = [];
   switch interleave
      case 'bip'
         image = reshape( image, bands, samples, lines );
      case 'bsq'
         image = reshape( image, samples, lines, bands );
      case 'bil'
         image = reshape( image, samples, bands, lines );
      otherwise
         error( 'Error reshaping the image array, invalid interleave' );
   end
   
   image = permute( image, [2 1 3] );

   image = squeeze( image );

   switch dataType
      case 'uint8'
         image = uint8( image );
      case 'int16'
         image = int16( image );
      case 'int32'
         image = int32( image );
      case 'float32'
         image = single( image );
      case 'float64'
         image = double( image );
      case 'uint16'
         image = uint16( image );
      case 'uint32'
         image = uint32( image );
      case 'int64'
         image = int64( image );
      case 'uint64'
         image = uint64( image );
   end

   fclose( fid );
