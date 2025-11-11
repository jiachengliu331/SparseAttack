function maxNumber = findMaxNumber(folderPath)

    files = dir(fullfile(folderPath,'*test*'));

    numbers = [];

    for i=1:length(files)

        fileName = files(i).name;

        token = regexp(fileName,'test(\d+)','tokens');

        if ~isempty(token)
            num = str2double(token{1}{1});
            numbers = [numbers,num];
        end

    end

    if ~isempty(numbers)
        maxNumber = max(numbers) + 1;
    else
        maxNumber = 1;
    end

end

