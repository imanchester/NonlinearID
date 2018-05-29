function [im, dims] = datread_jack(filename, nrows, ncols, nsets, datatype)
% function [im dims] = datread(filename, nrows, ncols, nsets, datatype)
% Requires: Image processing toolbox.
%           Currently: reads only 1 channel images
%           If datatype is specified, it must be a legal Matlab datatype.
% Effects: Reads a datfile, returning the result in im, 
%          and the associated dimensions in dims.
%          nrows, ncols, nsets, and datatype are OPTIONAL and will be read
%          from the header unless specified by the user.
% 
%          To display, must normalize all values between 0 and 1.
%          Example:
%          imagesc(imriv) OR
%          montage(imrivseq/255, [115 170 5])
%          truesize


error(nargchk(1, 5, nargin));

fid = fopen([filename,'/data'], 'r');
if fid == -1
  error(['Could not open file: ' filename])
end

if nargin < 5
  [status, DATtypeff] = unix(['getkey ' filename ' _data_type']);
  % discard newline at end of char
  DATtype = sscanf(DATtypeff, '%s');
  if status ~= 0, error(['Could not read datatype of DATfile ' filename]); end
matlab_types = [ 
 'uchar         '
 'schar         '
 'short         '
 'ushort        '
 'int           '
 'uint          '
 'float         '
 'double        '];

dat_types = [
    'unsigned_1'  
    'signed_1  '  
    'signed_2  '  
    'unsigned_2'  
    'signed_4  '  
    'unsigned_4'  
    'float     '  
    'double    '];

  datatype = translate_str(dat_types, matlab_types, DATtype);
end

if nargin < 4
  [status, str] = unix(['getkey ' filename ' _data_sets']);
  if status ~= 0
     nsets = 1;
     disp('Warning: could not read _data_sets from file; assuming 1 data set')
  else
    nsets = sscanf(str, '%d');
    if prod(size(nsets)) ~= 1, error('getkey returned more than 1 number'); end
  end
end

if nargin < 3
  [status, str] = unix(['getkey ' filename ' _dimensions']);
  if status ~= 0, error(['Could not read dimensions of DAT file ' filename]); end
  datdims = sscanf(str, '%d');
  [datdimsrows, datdimscol] = size(datdims);
  if datdimsrows > 2
    error(['Datfile ' file ' can have at most 2 dimensions'])
  end
  if datdimscol ~= 1,  error('getkey returned something weird'); end
  nrows = datdims(1);  
  if datdimsrows == 1
    ncols = 1;
  else
    ncols = datdims(2);
  end
end

[status, str] = unix(['getkey ' filename ' _channels']);
if status == 0 
  nchannels = sscanf(str, '%d');
  if prod(size(nchannels)) ~= 1, error('getkey returned more than 1 number'); end
  if nchannels ~= 1, error('Datfile can have at most one channel!'); end
end    
  
% disp(['#rows = ' int2str(nrows) '  #cols = ' int2str(ncols) '  #frames = ' ... 
%     int2str(nsets) '  datatype: ' datatype ' (' DATtype ')']);

% Old version
% -------------------------------------------------------------------------
% % Preallocate memory
% im = zeros(nrows, ncols*nsets);
% 
% [tmp_im, noRead] = fread(fid, [ncols, nrows], datatype);
% 
% if noRead ~= ncols*nrows 
% %     error(['While reading frame ', int2str(i), ' could only read ', int2str(noRead), ' bytes out of ', int2str(ncols*nrows)]);
% end
% 
% im(imslice([nrows ncols nsets], 1)) = tmp_im';
% 
% for i = 2:nsets
%   [tmp_im, noRead] = fread(fid, [ncols, nrows], datatype);
%   if noRead ~= ncols*nrows 
%     error(['While reading frame ', int2str(i), ' could only read ', int2str(noRead), ' bytes out of ', int2str(ncols*nrows)]);
%   end
%   im(imslice([nrows ncols nsets], i)) = tmp_im';
% end  


% New version
% -------------------------------------------------------------------------
% Preallocate memory: 3d array
im = zeros(nrows, ncols, nsets);


for i = 1:nsets
    
    [tmp_im, noRead] = fread(fid, [ncols, nrows], datatype);
    
    im(:,:,i) = tmp_im';
    
end



% if noRead ~= ncols*nrows 
% %     error(['While reading frame ', int2str(i), ' could only read ', int2str(noRead), ' bytes out of ', int2str(ncols*nrows)]);
% end
% 
% im(imslice([nrows ncols nsets], 1)) = tmp_im';

% for i = 2:nsets
%   [tmp_im, noRead] = fread(fid, [ncols, nrows], datatype);
%   if noRead ~= ncols*nrows 
%     error(['While reading frame ', int2str(i), ' could only read ', int2str(noRead), ' bytes out of ', int2str(ncols*nrows)]);
%   end
%   im(imslice([nrows ncols nsets], i)) = tmp_im';
% end 

fclose(fid);  
dims = [nrows ncols nsets];

