cd ~/Desktop % change to needed ovation file
load try2.mat % load name of ovation file

import auimodel.* % imports needed functions
import vuidocument.*

list = EpochList(OvationExport) ; % creates an object called 'list' : unorganized  epochs

tree = EpochTree(list,{'protocolSettings.acquirino:cellBasename',...
    'protocolSettings.notes:Mean',...
    'protocolSettings.notes:Amp',...
    'protocolSettings.notes:PreSynapticHold'}) % create an object tree structure called 'tree': organized epochs

tree.visualize




% 
% figure(1); clf;
% for leaf = 1:length(tree.leafNodes) % for each leaf (unique combo of mean, amp, holdsignal as specied above)
%     PreEpochData = tree.leafNodes{leaf}.epochList.responsesByStreamName('Amp_1') ; % this is a method which gets all the epochs from a particular leaf from Amp1
%     meanLeaf(leaf,:) = mean(PreEpochData) ; % the mean of these leaf epochs
% %     plot([1:size(PreEpochData, 2)], meanLeaf(leaf,:)) ; % plot these leaf means
% %     pause(1);
%     
%     traitStructure{leaf} = getLeafTraits(tree.leafNodes{leaf}); % get leaf traits
%     hp(leaf) = traitStructure{leaf}.pshp ; %  holding  potential
% % 
% %     if hp>10
% end 









