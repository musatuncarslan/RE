function [div, divNum] = partitionPeriods(Simparams)

    divNo = find(diff(Simparams.simPeriods) ~= 1); % index of division
    if (isempty(divNo) ~= 1) % divide trajectory into its piecewise partitions
        divNo = [0 divNo length(Simparams.simPeriods)];
        divNum = length(divNo)-1; % number of divions
        div = [];
        for k=1:divNum
            div = [div divNo(k)+1 divNo(k+1)]; % start and end index of each division
        end
    else % trajectory is continuous, it can be divided into partitions from anywhere
        divNum = 1;
        div = [1 length(Simparams.simPeriods)];
    end

    % check if any of the divisions have more than some set amount of periods,
    % if so partition them further
    divTemp = div; % assign temporary division vector
    newDiv = [];
    divSize = 1000;
    for k=1:divNum
            divIdx = [divTemp((k-1)*2+1) divTemp(k*2)];
            divL = divIdx(2) - divIdx(1); % length of divison
            if divL > divSize % some number divSize, so that each division has at most divSize periods for example
                pDiv = divIdx(1):divIdx(2);
                newDiv = [newDiv, divIdx(1) sort([pDiv(divSize:divSize:end) pDiv(divSize:divSize:end)+1]) divIdx(2)];
            else
                newDiv = [newDiv, divIdx(1), divIdx(2)];
            end
    end
    div = sort(newDiv);
    divNum = length(div)/2;

end