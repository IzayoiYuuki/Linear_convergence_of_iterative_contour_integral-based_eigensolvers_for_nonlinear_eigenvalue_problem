clear;
rng(0);

% data
TestMatrices = ["spring", "aw2", "butterfly", "ls", "photonic", "r2r", "hadeler", "gun", "pdde"];
FileNames = "./figure/9examples/" + TestMatrices + ".txt";

mytable = [];
for fileNo = 1 : 9

    % read from files
    FullFile = readlines(FileNames(fileNo));
    timenow = str2num(FullFile(2));

    % write to mytable
    mytable = [mytable; timenow];
end

writelines(string(datetime('now')), './figure/9examples/timetable.txt', WriteMode="append");
dlmwrite('./figure/9examples/timetable.txt', mytable, '-append', 'delimiter', ' ', 'precision', 4);