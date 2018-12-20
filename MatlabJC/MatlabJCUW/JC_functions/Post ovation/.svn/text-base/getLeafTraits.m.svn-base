function traitStructure = getLeafTraits(leaf)
% GETLEAFTRAITS
%   traitStructure = getLeafTraits(leaf)
%   pulls out parameters from an Epoch object, typically for analysis on
%   an leaf containing data of a similar type.
%   traitStructure is a struct with the following fields:
%       flfd - a basename suitable for a structure field name eg c2_100708
%       dvec - date in vector form, needed for getting the calibration
%       dura - duration
%       prep - prepnts
%       taip - tail points
%       stmp - stimulation points
%       ampl - amplitude in DACC counts of the traces
%       samp - sampling interval in seconds
%
% TA 10/20/08

% edited JC 12/29?08
% Made for input tree.leafNodes{leaf}
% added pshp - presynaptic holding potential


%--- Get leaf parameters
epli = leaf.epochList;
resn = epli.responseStreamNames;
stin = epli.stimuliStreamNames;

elex = epli.elements{1};
dvec = elex.startDate;
prep = elex.protocolSettings.(strcat('stimuli:',stin{1},':prepts'));
taip = elex.protocolSettings.(strcat('stimuli:',stin{1},':tailpts'));
stmp = elex.protocolSettings.(strcat('stimuli:',stin{1},':stmpts'));
samp = elex.protocolSettings.('acquirino:samplingInterval')/1e6; %sec
dura = elex.protocolSettings.('acquirino:dataPoints')*samp; %sec;
ampl = elex.protocolSettings.(strcat('stimuli:',stin{1},':amp'));
file = elex.protocolSettings.('acquirino:cellBasename');
pshp = elex.protocolSettings.('notes:PreSynapticHold');

% file field name for this leaf, eg c12_082408_327_2010
flfd = sprintf('%s_%s_%g_%g',file(7:end),file(1:end-2),ampl,dura);

traitStructure = struct(...
    'resn',resn,...
    'stin',stin,...
    'dvec',dvec,...
    'prep',prep,...
    'taip',taip,...
    'stmp',stmp,...
    'samp',samp,...
    'dura',dura,...
    'ampl',ampl,...
    'file',file,...
    'flfd',flfd,...
    'pshp',pshp);
