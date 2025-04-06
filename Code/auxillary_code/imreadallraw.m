function ImgSeq = imreadallraw(filename,x,y,nFrames,precision)

    % '*uint8' 8 bit imaging, raw data, behavioral camera
    % '*uint16' 16 bit imaging, raw data, VSD camera
    % '*float32' 32 bit, filtered data, VSD camera
    
    if isempty(nFrames)
        nFrames = 1000000;
    end
    if isempty(precision)
        precision = '*float32';
    end
    fid0 = fopen(filename, 'r', 'b'); 
%     fid0 = fopen(filename, 'r'); %#RG: removed "b" because valid file identifier was not created on the laptop
 
    ImgSeq = fread(fid0,[x*y nFrames],precision);
    fclose(fid0);
 
    nFrames = size(ImgSeq,2);
    ImgSeq = reshape(ImgSeq,x,y,nFrames);
    ImgSeq = permute(ImgSeq, [2, 1, 3]);
end
