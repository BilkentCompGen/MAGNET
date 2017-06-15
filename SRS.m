function HammingMask = SRS(HammingMask)


%HammingMask=[1 0 0 0 0 0 0 0 0 0 0 1 0 1 1 1 0 1 0 1 0 0 1 1 1 1 0 1 1 0 1 1 1 0 1 1 1 0 1 1 1 0 0 1 0 0 0];

% Speculative removal of short-matches (SRS)
for i = 1:length(HammingMask)
    if (i<=length(HammingMask)-3)
        
        % pattern = '101'
        if(HammingMask(i) && (~HammingMask(i+1)) && HammingMask(i+2))
            HammingMask(i+1) = 1;
        end
        
        if (i<=length(HammingMask)-4)
            
            % pattern = '1001'
            if(HammingMask(i) && (~HammingMask(i+1)) && (~HammingMask(i+2)) && HammingMask(i+3))
                
                HammingMask(i+1) = 1;
                HammingMask(i+2) = 1;
            end
        end
    end
end

end