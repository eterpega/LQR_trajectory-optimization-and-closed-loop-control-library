function serialCallback(obj, ~)
    global messageData messageReceived
    persistent byteCount buffer
    
    if(isempty(byteCount)) 
        byteCount = 0;
        buffer = zeros(1, 34);
    end
    
%     str = fscanf(obj);
%     disp(str);
	byteIn = fscanf(obj);
    buffer(byteCount) = byteIn;
    byteCount = byteCount + 1;
    if(byteCount == 34)
        messageReceived = true;
    
end