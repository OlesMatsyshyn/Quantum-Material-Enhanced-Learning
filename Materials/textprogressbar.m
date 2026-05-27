function textprogressbar(c)
%TEXTPROGRESSBAR Minimal MATLAB-compatible console progress bar.

persistent strCR;

if isempty(strCR) && ~ischar(c)
    error('The text progress bar must be initialized with a string.');
elseif isempty(strCR) && ischar(c)
    fprintf('%s', c);
    strCR = -1;
elseif ~isempty(strCR) && ischar(c)
    strCR = [];
    fprintf([c '\n']);
elseif isnumeric(c)
    c = max(0, min(100, floor(c)));
    percentageOut = [num2str(c) '%%'];
    percentageOut = [percentageOut repmat(' ', 1, 5 - length(percentageOut))];
    strOut = percentageOut;

    if strCR == -1
        fprintf(strOut);
    else
        fprintf([strCR strOut]);
    end

    strCR = repmat('\b', 1, length(strOut) - 1);
else
    error('Unsupported argument type.');
end
end