
function Accepted = MAGNET(RefSeq, ReadSeq, ErrorThreshold)
%function [Accepted, ErrorCount, MagnetMask1,mask1,mask2,mask3,mask4,mask5,mask6,mask7] = MAGNET(RefSeq, ReadSeq, ErrorThreshold)

% ErrorThreshold = 5;
% RefSeq   = 'AAAAAAACGTATATCCTCTTTATTTGGGGTGGAGAGTTCTGTAGATGTCTATTAGGTCCACTTGGTGCAGAGCTGAGTTCAATTCCTGGGTATCCTTGTT';
% ReadSeq  = 'AAAAAAATGTATATCCTCTTAATTTGGGGTGGACAGTTCTGTAGATGTCTATTATGTCCACTTGGTGCAGAGATGAGTTCAATTCCTGGGTATCCTTTTT';
% 
% RefSeq   = 'AAAAAAATGTATATCCTCTTTATTTGGGGTGGAGAGTTCTGTAGATGTCTATTAGGTCCACTTGGTGCAGAGCTGAGTTCAATTCCTGGGTATCCTTGTT';
% ReadSeq  = 'AAAAAAATGTATATTCTGTTGATTTGGGGTGGAGAGTTCTGTAGATGTCTATTAGGTCTGCTTGGTGCAGAGCTGAGTTCAATTCCTGGGTATCCTTGTT';
% 
% RefSeq   = 'AAAAAAATGTATATCCTCTTTATTTGGGGTGGAGAGTTCTGTAGATGTCTATTAGGTCCACTTGGTGCAGAGCTGAGTTCAATTCCTGGGTATCCTTGTT';
% ReadSeq  = 'AAAAAAATGTATATTCTGTTGATTTGGGGTGGAGAGTTCTGTAGATGTCTATTAGGTCTGCTTGGTGCAGAGCTGAGTTCAATTCCTGGGTATCCTTGTT';
% AlignStruct = localalign(ReadSeq, RefSeq);

Accepted=0;
HammingMask=zeros((2*ErrorThreshold)+1,length(ReadSeq));
% Original Hamming Mask
ReadSeq1 = ReadSeq;
for n=1:length(ReadSeq)
    HammingMask(1,n)= not(strcmp(ReadSeq1(n), RefSeq(n)));
end
% Check Hamming Distance
ErrorCount = length(find(HammingMask(1,:)));
if ErrorCount<=ErrorThreshold
    Accepted=1;
end


if Accepted==0
    % Incremental-step right shifted Hamming Masks (Deletion)
    for e=1:ErrorThreshold
        ReadSeq1 =  [RefSeq(1:e), ReadSeq(1:end-e)];
        for n=1:length(ReadSeq)
            HammingMask(e+1,n)= not(strcmp(ReadSeq1(n), RefSeq(n)));
        end
        % Check Hamming Distance
        
        ErrorCount = length(find(HammingMask(e+1,:)));
        if ErrorCount<=ErrorThreshold
            Accepted=1;
        end
        
        %if HammingMask(e+1,e+1)==0
        HammingMask(e+1,1:e) = zeros(1,e);
        %end
    end
end

if Accepted==0
    % Incremental-step left shifted Hamming Masks (Insertion)
    for e=1:ErrorThreshold
        ReadSeq1 =  [ReadSeq(e+1:end), RefSeq(end-e+1:end)];
        for n=1:length(ReadSeq)
            HammingMask(ErrorThreshold+e+1,n)= not(strcmp(ReadSeq1(n), RefSeq(n)));
        end
        % Check Hamming Distance
        ErrorCount = length(find(HammingMask(ErrorThreshold+e+1,:)));
        if ErrorCount<=ErrorThreshold
            Accepted=1;
        end
        
        %if HammingMask(ErrorThreshold+e+1,end-e)==0
        HammingMask(ErrorThreshold+e+1,end-e+1:end) = zeros(1,e);
        % end
    end
end
% RunsMask1 = '';
if Accepted ==0
    % Search for Longest Consecutive Zeros
    MagnetMask = ones(1,length(HammingMask(1,:)));
    GlobalMaxZerosStartIndex=1;
    GlobalMaxZerosEndIndex=length(ReadSeq);
    
    MagnetMask = Extraction_Encapsulation(HammingMask,GlobalMaxZerosStartIndex, GlobalMaxZerosEndIndex, MagnetMask, ErrorThreshold);
    ErrorCount = length(find(MagnetMask));
    if ErrorCount<=ErrorThreshold
        Accepted=1;
    else
        Accepted=0;
    end
    % RunsMask1 = regexprep(mat2str(double(MagnetMask(1,:))),'[^\w'']','');
end



% mask1=regexprep(mat2str(double(HammingMask(1,:))),'[^\w'']','');
% mask2=regexprep(mat2str(double(HammingMask(2,:))),'[^\w'']','');
% mask3=regexprep(mat2str(double(HammingMask(3,:))),'[^\w'']','');
% mask4=regexprep(mat2str(double(HammingMask(4,:))),'[^\w'']','');
% mask5=regexprep(mat2str(double(HammingMask(5,:))),'[^\w'']','');
% mask6=regexprep(mat2str(double(HammingMask(6,:))),'[^\w'']','');
%
% mask7=regexprep(mat2str(double(HammingMask(7,:))),'[^\w'']','');
% % mask8=regexprep(mat2str(double(HammingMask(8,:))),'[^\w'']','');
% % mask9=regexprep(mat2str(double(HammingMask(9,:))),'[^\w'']','');

end