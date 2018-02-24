% State vector
% CMD
clear serialCallback
global messageData messageReceived

messageData = [];
messageReceived = true;

s = serial('COM5');
set(s, 'BaudRate', 921600);

% s.Terminator = 'LF';%'CR';
% s.BytesAvailableFcnMode = 'terminator';
% s.BytesAvailableFcn = @serialCallback;

% s.Terminator = 'LF';%'CR';
s.BytesAvailableFcnMode = 'byte';
s.BytesAvailableFcnCount = 34;
s.BytesAvailableFcn = @serialCallback;


fopen(s);

% x = 100*rand(1, 1000);%[10.5 pi 67 89];
% x = [uint8(3) uint32(5)];
% data = typecast(single(x), 'uint8');

% traj = 1:400;%100*rand(1, 10);
% trajBytes = typecast(single(traj), 'uint8');
% data = [uint8(3) typecast(uint32(length(trajBytes)+5), 'uint8') trajBytes];

data = [uint8(0) typecast(uint32(5), 'uint8')];

startTime = tic;
for j = 1:length(data)
%     disp(['Writing ' num2str(data(j))]);
%     fwrite(s, data(j));
end
toc(startTime)
disp('Transfer complete')
% x(end)
pause(2);
% out = fscanf(s)
% out = fscanf(s)
% out = fscanf(s)
s.BytesAvailableFcn = '';
fclose(s);
delete(s);
clear s;

% try
%     fprintf(s1, 'Vroem');
% end
