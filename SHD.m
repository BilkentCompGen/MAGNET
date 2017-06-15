
function SHDAccepted = SHD(RefSeq, ReadSeq, ErrorThreshold)
%  ErrorThreshold = 0;
% RefSeq  = 'TCCATTGACATTCGTGAGCTGCTCCTTCTCTCCCACCCCTTTGCCC';
% ReadSeq = 'TCCATTGACATTCGTGACTCTCCTTCTCTCCCACCCCTTTGCCCCC';

% RefSeq  = 'GATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAGTTGGATGTGTGCATTTTCCTGAGAGGAA';
% GATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAGTTGGATGTGTGCATTTTCCTGAGAGGAA - 0
% GATCCTTGAAGCGCCCCCACGGGCTTCTTCTCAAAGTTGGATGTGTGCATTTTCCTGAGAGGAA - 2
% GATCCTTGAAGCGCCCCCAAGGGGTTCTTCTCAAAGTTGGATGTGTGCATTTTCCTGAGAGGAA - 2
% GATCCTTGAAGCGCCCCCAAGGCCATCTTCTCAAAGATGGATGTGTGCATTTTCCTGAGAGGAA - 2
% GATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAGTTGGATGTGTGCAAACCTTCCTGAGAGG - 4
% GATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAGTAATGGATGTGTGCATTTTCCTGAGAGG - 2
% GATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAGTTGGATGTGTGCATTTTCCTTTTTGGAA - 4
% GATCCTTGAAGCGCCAAGGGCATCTTCTCAAAGTTGGATGTGTGCATTTTCCTGAGAGGAAAGC - 3
% GATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAGTTGGATGTGTACGTTTTTTCCTGAGAGG - 4
% GATCCTTGAAGCGCCCCCAAAAAAATCTTCTCAAAGTTGGATGTGTGCATTTTCCTGAGAGGAA - 4
% GATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAGCCGGATGTGTGCATTTTCCTGAGAGGAA - 2
% GATCCTTGAAGCGCCGGCAAGGGCATCTTCTCAAAGTTGGATGTGTCGCTTTTCCTGAGAGGAA - 4
% GATCCTTGAAGCGCGGGCAAGGGCATCTTCTCAAAGAAGGATGTGTGCATTTTCCTGAGAGGAA - 5
% GATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAGTGGGATGTGTGCATTTTCCTGAGAGGAA - 1
% GATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAACCTGGATGTGTGCATTTTCCTGAGAGGAA - 2
% GATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAAAATGGATGTGTGCATTTTCCTGAGAGGAA - 2
% GATCCTTGAAGCGCCCCCAAGGGCATCTTCTCAAACCTGGATGTGTGCAGGTTCCTGAGAGGAA - 4
% GATCCTTGAAGCGCCCCCAATTGCATCTTCTCAAAGTTGGATGTGTGCACCTTCCTGAGAGGAA - 4
% GATCCTTGAAGCGGGGACAAGGGCATCTTCTCAAAGTTGGATGTGTGCATTTTCCTGAGAGGAA - 4
% ReadSeq = 'GATCCTTGAAGCGCCAAGGGCATCTTCTCAAAGTTGGATGTGTGCATTTTCCTGAGAGGAAAGC';% - 16

HammingMask=zeros((2*ErrorThreshold)+1,length(ReadSeq));
% Original Hamming Mask
ReadSeq1 = ReadSeq;
for n=1:length(ReadSeq)
    HammingMask(1,n)= not(strcmp(ReadSeq1(n), RefSeq(n)));
end
HammingMask(1,:) = SRSoSHD(HammingMask(1,:));
% Incremental-step right shifted Hamming Masks
for e=1:ErrorThreshold
    ReadSeq1 =  [RefSeq(1:e), ReadSeq(1:end-e)];
    
    for n=1:length(ReadSeq)
        HammingMask(e+1,n)= not(strcmp(ReadSeq1(n), RefSeq(n)));
    end
    HammingMask(e+1,:) = SRSoSHD(HammingMask(e+1,:));
end


% Incremental-step left shifted Hamming Masks
for e=1:ErrorThreshold
    ReadSeq1 =  [ReadSeq(e+1:end), RefSeq(end-e+1:end)];
    
    for n=1:length(ReadSeq)
        HammingMask(ErrorThreshold+e+1,n)= not(strcmp(ReadSeq1(n), RefSeq(n)));
    end
    HammingMask(ErrorThreshold+e+1,:) = SRSoSHD(HammingMask(ErrorThreshold+e+1,:));
end
% Generating the AND Mask
ANDMask = logical(HammingMask(1,:));
for n=2:((2*ErrorThreshold)+1)
    ANDMask = and(ANDMask, HammingMask(n,:));
end

% Speculative removal of short-matches (SRS) Count
MinErrors=0;
i=1;
%ANDMaskStr=regexprep(mat2str(double(ANDMask)),'[^\w'']','');
while i < length(ANDMask)
    if (i<=(length(ANDMask)-3))
        if (ANDMask(1,i)==0  && ANDMask(1,i+1)==0 && ANDMask(1,i+2)==0 && ANDMask(1,i+3)==0) %case '0000'
            MinErrors=MinErrors+0;
        elseif (ANDMask(1,i)==0  && ANDMask(1,i+1)==1 && ANDMask(1,i+2)==0 && ANDMask(1,i+3)==1) % {'0101','0110','1001','1010','1011','1101' }
            MinErrors=MinErrors+2;
        elseif (ANDMask(1,i)==0  && ANDMask(1,i+1)==1 && ANDMask(1,i+2)==1 && ANDMask(1,i+3)==0) % {'0101','0110','1001','1010','1011','1101' }
            MinErrors=MinErrors+2;
        elseif (ANDMask(1,i)==1  && ANDMask(1,i+1)==0 && ANDMask(1,i+2)==0 && ANDMask(1,i+3)==1) % {'0101','0110','1001','1010','1011','1101' }
            MinErrors=MinErrors+2;
        elseif (ANDMask(1,i)==1  && ANDMask(1,i+1)==0 && ANDMask(1,i+2)==1 && ANDMask(1,i+3)==0) % {'0101','0110','1001','1010','1011','1101' }
            MinErrors=MinErrors+2;
        elseif (ANDMask(1,i)==1  && ANDMask(1,i+1)==0 && ANDMask(1,i+2)==1 && ANDMask(1,i+3)==1) % {'0101','0110','1001','1010','1011','1101' }
            MinErrors=MinErrors+2;
        elseif (ANDMask(1,i)==1  && ANDMask(1,i+1)==1 && ANDMask(1,i+2)==0 && ANDMask(1,i+3)==1) % {'0101','0110','1001','1010','1011','1101' }
            MinErrors=MinErrors+2;
        else
            MinErrors=MinErrors+1;
        end
    end
    i = i + 4;
end

if MinErrors<=ErrorThreshold
    SHDAccepted=1;
else
    SHDAccepted=0;
end
end