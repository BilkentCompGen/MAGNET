function MagnetMask = Extraction_Encapsulation(HammingMask,GlobalMaxZerosStartIndex, GlobalMaxZerosEndIndex, MagnetMask, ErrorThreshold)

% GlobalMaxZerosStartIndex=1;
% GlobalMaxZerosEndIndex=length(ReadSeq);

if (GlobalMaxZerosStartIndex <= GlobalMaxZerosEndIndex)
    % Search for Longest Consecutive Zeros in the Middle
    ConsecutiveZerosInfo = ConsecutiveZeros(HammingMask(1,GlobalMaxZerosStartIndex:GlobalMaxZerosEndIndex),1);
    for i=2:((2*ErrorThreshold)+1)
        ConsecutiveZerosInfo = [ConsecutiveZerosInfo, ConsecutiveZeros(HammingMask(i,GlobalMaxZerosStartIndex:GlobalMaxZerosEndIndex),i)];
    end
    [MaxZeros, MaxZerosIndex] = max(ConsecutiveZerosInfo(1,:));
    if GlobalMaxZerosStartIndex==1 && GlobalMaxZerosEndIndex==length(MagnetMask)
        if(MaxZeros<((length(MagnetMask)/(ErrorThreshold+1))))
            ErrorCount = length(MagnetMask)-MaxZeros;
            Accepted=0;
            return;
        end
    end
    MaxZerosStartIndex = GlobalMaxZerosStartIndex-1+ConsecutiveZerosInfo(2, MaxZerosIndex);
    MaxZerosEndIndex = GlobalMaxZerosStartIndex-1+ConsecutiveZerosInfo(3, MaxZerosIndex);
    % MaxZerosMaskID = ConsecutiveZerosInfo(4, MaxZerosIndex);
    MagnetMask(1,MaxZerosStartIndex:MaxZerosEndIndex) = zeros(1,MaxZeros);
    
    % Search for Longest Consecutive Zeros in the Left
    MagnetMask = Extraction_Encapsulation(HammingMask,GlobalMaxZerosStartIndex, MaxZerosStartIndex-2, MagnetMask, ErrorThreshold);
    
    % Search for Longest Consecutive Zeros in the Right
    MagnetMask = Extraction_Encapsulation(HammingMask,MaxZerosEndIndex+2, GlobalMaxZerosEndIndex, MagnetMask, ErrorThreshold);
end
% %%
% % Search for Longest Consecutive Zeros in the Right
% while(MaxZerosStartIndex>GlobalMaxZerosStartIndex+1)
%     ConsecutiveZerosInfo = ConsecutiveZeros(HammingMask(GlobalMaxZerosStartIndex,1:MaxZerosStartIndex-2),1);
%     for i=2:((2*ErrorThreshold)+1)
%         ConsecutiveZerosInfo = [ConsecutiveZerosInfo, ConsecutiveZeros(HammingMask(i,GlobalMaxZerosStartIndex:MaxZerosStartIndex-2),i)];
%     end
%     [MaxZeros, MaxZerosIndex] = max(ConsecutiveZerosInfo(1,:));
%     MaxZerosStartIndex = GlobalMaxZerosStartIndex-1+ConsecutiveZerosInfo(2, MaxZerosIndex);
%     MaxZerosEndIndex1 = GlobalMaxZerosStartIndex-1+ConsecutiveZerosInfo(3, MaxZerosIndex);
%     %MaxZerosMaskID = ConsecutiveZerosInfo(4, MaxZerosIndex);
%     MagnetMask(1,MaxZerosStartIndex:MaxZerosEndIndex1) = zeros(1,MaxZeros);
% end
%
% % Search for Longest Consecutive Zeros in the Left
% while(MaxZerosEndIndex<GlobalMaxZerosEndIndex-1)
%     ConsecutiveZerosInfo = ConsecutiveZeros(HammingMask(1,MaxZerosEndIndex+2:length(ReadSeq)),1);
%     for i=2:((2*ErrorThreshold)+1)
%         ConsecutiveZerosInfo = [ConsecutiveZerosInfo, ConsecutiveZeros(HammingMask(i,MaxZerosEndIndex+2:length(ReadSeq)),i)];
%     end
%     [MaxZeros, MaxZerosIndex] = max(ConsecutiveZerosInfo(1,:));
%     MaxZerosStartIndex = MaxZerosEndIndex+1+ConsecutiveZerosInfo(2, MaxZerosIndex);
%     MaxZerosEndIndex = MaxZerosEndIndex+1+ConsecutiveZerosInfo(3, MaxZerosIndex);
%     %MaxZerosMaskID = ConsecutiveZerosInfo(4, MaxZerosIndex);
%     MagnetMask(1,MaxZerosStartIndex:MaxZerosEndIndex) = zeros(1,MaxZeros);
% end

% MagnetMask1 = regexprep(mat2str(double(MagnetMask(1,:))),'[^\w'']','');
% mask1=regexprep(mat2str(double(HammingMask(1,:))),'[^\w'']','');
% mask2=regexprep(mat2str(double(HammingMask(2,:))),'[^\w'']','');
% mask3=regexprep(mat2str(double(HammingMask(3,:))),'[^\w'']','');
% mask4=regexprep(mat2str(double(HammingMask(4,:))),'[^\w'']','');
% mask5=regexprep(mat2str(double(HammingMask(5,:))),'[^\w'']','');
% mask6=regexprep(mat2str(double(HammingMask(6,:))),'[^\w'']','');
% mask7=regexprep(mat2str(double(HammingMask(7,:))),'[^\w'']','');
% mask8=regexprep(mat2str(double(HammingMask(8,:))),'[^\w'']','');
% mask9=regexprep(mat2str(double(HammingMask(9,:))),'[^\w'']','');

end