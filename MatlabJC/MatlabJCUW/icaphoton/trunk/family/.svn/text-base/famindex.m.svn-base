function indexnumbs = famindex(condition,epochnumbs)
% function indexnumbs = famindex(condition,epochnumbs)
% Given epochnumbers finds the index numbers in the FamilyCondition that
% corresponds to these numbers.

%  allepochs is a list of all epochs in the familycondition (real epoch #s)
%  mkmk  March, 2002

allepochs = condition.EpochNumbers;
epochlength = length(epochnumbs);
indexnumbs = zeros(1,epochlength);
for i = 1:epochlength
  indexnumbs(i) = find(allepochs == epochnumbs(i));
end

