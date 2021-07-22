%+
% NAME:
%   read_envi_header
%
% PURPOSE:
%   This function reads an ENVI header file and parses out the 
%   number of samples, lines, bands, data type, interleave type and byte
%   order. These parameters are returned in an array.
%
% CALLING SEQUENCE:
%   [samples, lines, bands, dataType, interleave, byteOrder] = ...
%                                   read_envi_header( headerFilename )
%
% INPUT:
%   headerFilename
%      A string containing the full file path of the ENVI header file
%
% OUTPUT:
%   samples
%      A variable containing the number of samples in the image associated
%      with the ENVI header file
%   lines
%      A variable containing the number of lines in the image associated
%      with the ENVI header file
%   bands 
%      A variable containing the number of bands in the image
%      associated with the ENVI header file
%   dataType 
%      A variable containing the data type of the image associated
%      with the ENVI header file
%   interleave
%      A variable containing the files band interleave type; either
%      BIP, BSQ, or BIL are possible
%   byteOrder
%      A variable containing the byte order (0 is little endian [least
%      significant byte first], 1 is big endian [most significant byte
%      first])
%
% MODIFICATION HISTORY:
%   Written by       Sarah Paul and Carl Salvaggio
%   December, 2009   Original Code
%   February, 2011   Modified to parse byte order
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

function[samples, lines, bands, dataType, interleave, byteOrder] = ...
   read_envi_header( headerFilename)

   headerFID = fopen( headerFilename );
   if headerFID == -1
      error( 'Error opening header file, bad file ID' );
   end

   while 1

      tline = fgetl( headerFID );

      if ~ischar( tline )
         break;
      end

      [parameter, value] = strtok( tline, '=' );

      if strncmpi( parameter, 'samples', length( 'samples' ) )
         [delimeter, value] = strtok( value );
         samples = str2num( value );
      end

      if strncmpi( parameter, 'lines', length( 'lines' ) )
         [delimeter, value] = strtok( value );
         lines = str2num( value );
      end

      if strncmpi( parameter, 'bands', length( 'bands' ) )
         [delimeter, value] = strtok( value );
         bands = str2num( value );
      end

      if strncmpi( parameter, 'data type', length( 'data type' ) )
         [delimeter, value] = strtok( value );
         dataType = str2num( value );

         switch dataType
            case 1
               dataType = 'uint8';
            case 2
               dataType = 'int16';
            case 3
               dataType = 'int32';
            case 4
               dataType = 'float32';
            case 5
               dataType = 'float64';
            case 12
               dataType = 'uint16';
            case 13
               dataType = 'uint32';
            case 14
               dataType = 'int64';
            case 15
               dataType = 'uint64';
            otherwise
               error( 'Unknown data type' );
         end

      end

      if strncmpi( parameter, 'interleave', length( 'interleave' ) )
         [delimeter, value] = strtok( value );
         interleave = strtrim( char( value ) );
      end

      if strncmpi( parameter, 'byte order', length( 'byte order' ) )
         [delimeter, value] = strtok( value );
         byteOrder = str2num( value );
      end

   end

   fclose( headerFID );