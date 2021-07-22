function [ output_args ] = WriteGridASCII( filename,ncols,nrows,data)
fid = fopen([filename],'w');
for i=1:nrows
     for j=1:ncols
        fprintf(fid, '%f ', data(i,j));
     end
   fprintf(fid, '\n');
end
fclose(fid);
clear temp;