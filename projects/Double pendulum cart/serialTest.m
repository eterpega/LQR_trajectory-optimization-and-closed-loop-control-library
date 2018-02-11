s = serial('COM5');
set(s, 'BaudRate', 921600);
s.Terminator = 'CR';
s.BytesAvailableFcnMode = 'terminator'
s.BytesAvailableFcn = @serialCallback;
fopen(s);
for i=1:20
%     i
    pause(0.1);
end
fclose(s);
delete(s);
clear s;