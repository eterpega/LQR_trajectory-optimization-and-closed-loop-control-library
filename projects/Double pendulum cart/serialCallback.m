function serialCallback(obj, ~)
    global messageData messageReceived
    buffer = uint8(fread(obj, 34));
    data = buffer(1:2)
    typecast(buffer(3:end), 'single')
end