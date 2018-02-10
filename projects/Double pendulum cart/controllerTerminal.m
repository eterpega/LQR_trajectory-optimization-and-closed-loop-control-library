%State vector
% CMD


s = serial('COM5');
set(s, 'BaudRate', 921600);%115200);
s.Terminator = 'LF';%'CR';
s.BytesAvailableFcnMode = 'terminator';
s.BytesAvailableFcn = @serialCallback;
fopen(s);

x = 100*rand(1, 1000);%[10.5 pi 67 89];
data = typecast(single(x), 'uint8');
% data = 'abcdefgh12345';
% fwrite(s, data);
startTime = tic;
for j = 1:length(data)
%     disp(['Writing ' num2str(data(j))]);
    fwrite(s, data(j));
%     out = fscanf(s)
%     fwrite(s, data);
%     sprintf('%d', data(j));
%     pause(0.1);
end
toc(startTime)
disp('Transfer complete')
x(end)
pause(2);
% out = fscanf(s)
% out = fscanf(s)
% out = fscanf(s)
s.BytesAvailableFcn = '';
fclose(s);
delete(s);
clear s;