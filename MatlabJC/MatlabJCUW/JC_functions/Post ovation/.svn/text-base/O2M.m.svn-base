% this script will import deal with ovation segregated data

cd ~/ % change to needed ovation file
load % load name of ovation file

import auimodel.* % imports needed functions
import vuidocument.*

list = EpochList(OvationExport) ; % creates an object called 'list'

tree = EpochTree(list,{'protocolSettings.acquirino:cellBasename'}) % create an object tree structure called 'tree'

tree.visualize % visualize the tree structure you made

plot(tree.leafNodes{1}.epochList.elements{1}.responses.Amp_1.data) % example of plotting


