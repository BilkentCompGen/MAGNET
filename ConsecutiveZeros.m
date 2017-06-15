function ConsecutiveZerosInfo = ConsecutiveZeros(InputSequence,MaskID)

% Inputs:
% InputSequence
% Value

% Output:
% returns the length of the longest repeated subsequence of Value's

%InputSequence = [ 0 0 1 1 0 0 0 1 1 0 0 0 0 1 1 1 0 0 0 0 1 0 1 ]; % data


w = [ 1 InputSequence 1 ]; % change the ones to zeros to find consecutives ones
EndIndex = find(diff(w)==1)-1;
StartIndex = find(diff(w)==-1);
ConsecutiveZerosVector = (EndIndex+1)-StartIndex; % lenghts of runs of 0's
%MaskIDVector = MaskID*ones(1,length(ConsecutiveZerosVector));

%ConsecutiveZerosInfo = [ConsecutiveZerosVector; StartIndex; EndIndex; MaskIDVector];
ConsecutiveZerosInfo = [ConsecutiveZerosVector; StartIndex; EndIndex];

% Value=0; % we want longest consecutive subsequence of Zeros
% y = diff([0 (InputSequence==Value) 0]);
% u=find(y==1)-find(y==-1);
% ConsecutiveZerosVector = -u;
%
% ZerosLengthMax = max(ConsecutiveZerosVector); % the length of the longest repeated subsequence
% ConsecutiveZerosVector(ConsecutiveZerosVector==max(ConsecutiveZerosVector))=0;
% ZerosLength2Max = max(ConsecutiveZerosVector);
% ConsecutiveZerosVector(ConsecutiveZerosVector==max(ConsecutiveZerosVector))=0;
% ZerosLength3Max = max(ConsecutiveZerosVector);
% %StartIndex = v(p); EndIndex = ConsecutiveZerosVector(p)-1; % this for the start and end indices

end