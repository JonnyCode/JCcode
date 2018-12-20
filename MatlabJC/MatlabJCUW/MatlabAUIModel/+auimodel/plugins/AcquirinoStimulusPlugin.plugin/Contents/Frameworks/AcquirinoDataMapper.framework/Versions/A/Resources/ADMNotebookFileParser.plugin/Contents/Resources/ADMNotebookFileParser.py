#!/usr/bin/env python
#
#  ADMNotebookFileParser.py
#  AcquirinoDataMapper
#
#  Created by Daniel Carleton on 12/20/07.
#  Copyright (c) 2007 __MyCompanyName__. All rights reserved.
#

from Foundation import *
from enum import Enum
import logging as log
import objc
import reflex
import re

# allow parser to be run independent from the data mapper during development
if ( __name__ == '__main__' ):
  objc.loadBundle('AcquirinoDataMapper', globals(),
    bundle_path='/Applications/AcquirinoDataMapper.app/Contents/Frameworks/AcquirinoDataMapper.framework')
  log.basicConfig(level = log.DEBUG)
else:
  log.basicConfig(level = log.WARN)

AbstractParser = objc.lookUpClass('ADMAbstractNotebookFileParser')
Cell = objc.lookUpClass('ADMAcquirinoCell')
EpochFamily = objc.lookUpClass('ADMAcquirinoEpochFamily')
Note = objc.lookUpClass('ADMAcquirinoNote')

def iscontinuous(a_list):
  l = a_list[0] - 1
  for i in a_list:
    if not i == l + 1:
      return False
    l = i
  return True

class ADMNotebookFileParser(AbstractParser):
  # the tokens this parser is capable of extracting from a notebook file
  Token = Enum('Purpose_Line', 'Cell_Header', 'Comment_Line',
    'Comment_Timestamp', 'Tail_Comment', 'Protocol_Comment',
    'Epoch', 'Parameters_Line', 'Family_Start')
  
  # matches timestamp format used with many tokens e.g. "3:28:37 PM"
  TIMESTAMP_RE = '\d{1,2}:\d{1,2}:?\d{1,2}? [AMP]{2}'
  
  EXPERIMENT_TIMEZONE = 'PST'
  
  def init(self):
    self = super(ADMNotebookFileParser, self).init()
    
    # fields accessed by objc via properties of the same name minus
    # underscores
    self.experimentPurpose_ = None
    self.experimentStartDate_ = None
    self.notes_ = None
    self.acquirinoCells_ = None
    
    self.dateFormatter = NSDateFormatter.alloc()
    self.currentExperimentDate = ""
    
    # this date format is based on concatenating a trimmed cell
    # data file name and a timestamp. e.g. "091807 12:02:47 PM"
    # timezone is appended
    self.dateFormatter.initWithDateFormat_allowNaturalLanguage_\
      ("%m%d%y %I:%M:%S %p %Z", None)
    
    self.initScanner()
    return self

  def parse(self):
    file = open(self.filePath(), 'Ur')
    iterator = self.scanner(file, self)
    
    Token = ADMNotebookFileParser.Token
    states = self.scanner.states
    
    self.notes_ = []
    self.cells_ = []

    current_cell = None
    cell_basenames = set()
    last_epoch_number = -1
    token = iterator.next()
    
    # this loop processes the token stream and constructs the object graph.
    # the file is assumed to be well-formed at this stage, because action
    # functions are monitoring the stream, and so only business logic
    # appears below.
    while token != None:
      log.debug("token: " + str(token))

      if token.id == Token.Purpose_Line:
        purpose = ''
        
        while token.id == Token.Purpose_Line:
          purpose += token.value
          token = iterator.next()
        
        self.experimentPurpose_ = purpose.strip()
        self.experimentPurpose_ += "\n\n"
        
        # we're on a new token that follows the purpose lines now -- go
        # back to the top and process it.
        continue
      
      # a tail comment outside of an epoch family is simply a note
      if token.id == Token.Tail_Comment and\
        iterator.state == states['in-cell']:
        
        comment, timestamp = self.parseTailComment(token.value)
        self.appendNote(comment, timestamp)
      
      if token.id == Token.Comment_Timestamp:
        log.debug("time: " + str(token))
        timestamp = self.timestampFromTime(token.value)

        comment = ""
        token = iterator.next()        
        while token != None and token.id == Token.Comment_Line:
          comment += token.value.strip() + "\n"
          token = self.safeNext(iterator)
        
        self.appendNote(comment, timestamp)
        continue
      
      if token.id == Token.Cell_Header:
        basename, timestamp = self.parseCellHeader(token.value)

        # experiment start date is defined as that of the first cell
        if self.experimentStartDate_ == None:
          self.experimentStartDate_ = timestamp

        if current_cell != None:
          self.labelCell(current_cell)
          log.debug(current_cell)
          self.cells_.append(current_cell)
          cell_basenames.add(current_cell.dataFileBaseName())

        if cell_basenames.issuperset(set([basename])):
          raise Exception, 'Cell data file %s appears twice in notes.'\
            % basename

        current_cell = Cell.new()
        current_cell.setDataFileBaseName_(basename)
        current_cell.setStartDate_(timestamp)
        current_cell.setEpochFamilies_([])

      if token.id == Token.Family_Start:
        start_date = self.parseFamilyStart(token.value)
        
        token = iterator.next()
        if token.id == Token.Family_Start:
          continue
        elif not token.id == Token.Protocol_Comment:
          log.debug(token)
          raise Exception, 'Protocol comment must follow family starting at %s'\
            % start_date

        family = EpochFamily.new()
        family.setStartDate_(start_date)
        family.setProtocolComment_(token.value.strip("*\n"))
        family.setEpochSettings_(
          NSMutableDictionary.dictionaryWithCapacity_(10))

        epochs = dict()
        parameters = dict()
        
        token = iterator.next()
        while token != None and token.id == Token.Parameters_Line:
          self.parseFamilyParameters(token.value, parameters)
          token = self.safeNext(iterator)

        family.setParameters_(parameters)

        epoch_number = None
        while token != None and token.id in [Token.Epoch, Token.Tail_Comment]:
          if token.id == Token.Epoch:
            epoch_number, epoch_settings = self.parseEpoch(token.value)
            epochs[epoch_number] = epoch_settings
          elif token.id == Token.Tail_Comment:
            comment, stamp = self.parseTailComment(token.value)

            note = Note.new()
            note.setComment_(comment.strip())
            note.setTimestamp_(stamp)

            family.setTailNote_(note)
 
          token = self.safeNext(iterator)

        # ignore empty families.  assume no epochs were collected for them.
        # (see assertion below)
        if epoch_number == None:
          continue

        # this family should start up where the last one left off
        start_epoch_number = last_epoch_number + 1
        last_epoch_number = epoch_number

        epoch_numbers = list(epochs.keys())
        epoch_numbers.sort()

        missing_epochs = list(set(range(start_epoch_number, last_epoch_number))\
          .difference(set(epoch_numbers)))
        missing_epochs.sort()

        log.debug('start: %d, last: %d, epochs: %s, missing: %s',
          start_epoch_number, last_epoch_number, str(epoch_numbers),
          str(missing_epochs))

        # workaround for Igor bug described in ticket:7 (sometimes first
        # epoch is left off of a family)
        if (len(missing_epochs) == 1 and
            missing_epochs[0] == start_epoch_number):
          log.warn("Missing single epoch (#%d) from start of group " +
            "in notes file (Igor bug) -- including anyway.", missing_epochs[0])
          
          epochs[start_epoch_number] = dict()

        # epochs are missing from the middle of this family.  this is a failure
        # case and the file must be repaired.
        elif not iscontinuous(epoch_numbers):
          log.error("Epochs missing from middle of group (found #s %s, cell " +
            " %s) in notes file.", str(epoch_numbers),
            current_cell.dataFileBaseName())

          raise Exception, 'Malformed notes file -- please repair.'
  
        # a contiguous block of epochs is missing.  this can be caused by a
        # calibration step's epochs not being recorded.  note that blocks
        # missing from the end of a cell can only be detected by the caller.
        elif (len(missing_epochs) > 1 and
              epoch_numbers[0] == missing_epochs[-1] + 1):
          log.warn("Epoch group (#s %s) missing from cell %s in notes file " +
            "-- excluding.", str(missing_epochs),
            current_cell.dataFileBaseName())
            
          start_epoch_number = epoch_numbers[0]

        family.setEpochNumberRange_\
          (self.makeRangeFrom_to_(start_epoch_number,
          last_epoch_number))
        
        family.epochSettings().addEntriesFromDictionary_(epochs)

        current_cell.epochFamilies().append(family)
        
        if token != None:
          continue

      token = self.safeNext(iterator)

    self.labelCell(current_cell)
    self.cells_.append(current_cell)

  def parseFamilyParameters(self, parameters_line, parameters):
    # parses lines of the the following forms:
    # Amp = 0.1 V   Dur = 500 ms   Mean = 0 V  Pre= 500 ms  Tail=  500 ms
    # PreSynaptic: Hold= 0 mV;  PostSynaptic: Hold= 0 mV
    # PostSynaptic: Hold= -60 mV; Pre/Tail = -60 mV; Dur= 100 ms;  Pre= 100 ms;
    # PreSynaptic: Hold= -60 mV  Hold During Step= -60 mV
    pattern = re.compile('([^:=]+ *= *[-0-9.]+ ?\w*|[^: ]+:)')

    group = ''
    for match in pattern.findall(parameters_line):
      if match[-1] == ':':
        group = match.strip('\t:')
        continue

      key, value = match.split('=')
      key = key.strip('; ')
      value = value.strip(' ')
      
      key = key.replace('/', '') # for e.g. "Pre/Tail"
      key = key.replace(' ', '') # for e.g. "
      key = group + key
      
      parameters[key] = value

  def parseEpochParameters(self, line, params):
    # used for epoch settings lines like:
    # seed = 1  Noise SD = 0.15 Temp = 31.8 C

    # FIXME: units are lost from last setting before a flag.  assuming low
    # importance as only unit seen in wild is "C" for celcius.
    pattern = re.compile('([^:=]+ *= *[-0-9.]+ *(?:\w{1}(?: |$))?|[^: ]+:)')

    group = ''
    for match in pattern.findall(line):
      key, value = match.split('=')
      key = key.strip(' ')
      value = value.strip(' ')

      params[key] = value
    
    # used to grab trailing flags e.g. for gclamp experiments "OffStepInhG1"
    # you can have one flag at the end of an epoch line and it will go into
    # settings key "Flag".  note that confusion could result from units of
    # more than one character in a different param. (usuall "C" for celcius)
    pattern = re.compile('.+[^=]+\s+(\S\S+)$')
    m = pattern.match(line)

    if m:
      params['Flag'] = m.groups()[0]

  def labelCell(self, cell):
      fams = cell.epochFamilies()
      notes = self.notes_
      
      log.debug("cell->label fams: %s, notes: %s" % (fams, notes))

      if len(fams) > 0 and len(notes) > 0\
         and fams[-1].startDate().compare_(notes[-1].timestamp()) < 0:
        cell.setLabel_(notes[-1].comment().strip())
      else:
        cell.setLabel_(cell.dataFileBaseName())

  def parseEpoch(self, epoch_line):
    pattern = re.compile("^\s+Epoch # (\d+)\s*(.*)")
    match = pattern.match(epoch_line)
    
    settings = dict()
    if match.group(2):
      self.parseEpochParameters(match.group(2), settings)

    return int(match.group(1)), settings

  def safeNext(self, iterator):
    try:
      token = iterator.next()
    except StopIteration:
      token = None
    
    return token
  
  def appendNote(self, comment, timestamp):
    note = Note.new()
    note.setComment_(comment)
    note.setTimestamp_(timestamp)

    self.notes_.append(note)
  
  def timestampFromTime(self, time):
    zone = ADMNotebookFileParser.EXPERIMENT_TIMEZONE
    
    log.debug("time from: |%s|" % time)
    
    # time does not include seconds
    if ( len(time.split(':')) == 2 ):
      parts = time.split(' ')
      time = parts[0] + ':00' + ' ' + parts[1]
    
    datestr = self.currentExperimentDate + " " + time.strip() + " " + zone
    log.debug("date->parse string: |%s|" % datestr)
    
    return self.dateFormatter.dateFromString_(datestr)
    
  def parseFamilyStart(self, start_line):
    timestamp_re = ADMNotebookFileParser.TIMESTAMP_RE

    log.debug("famStart->parse line: |%s|" % start_line)

    pattern = re.compile('.+?(' + timestamp_re + ')')
    match = pattern.match(start_line)

    time = match.group(1)
    
    return self.timestampFromTime(time)
  
  def parseTailComment(self, tail_comment):
    timestamp_re = ADMNotebookFileParser.TIMESTAMP_RE

    pattern = re.compile('^\s*\*(.+)\s+(' + timestamp_re + ')$')
    match = pattern.match(tail_comment)

    comment = match.group(1).strip()
    time = match.group(2)
    
    return comment, self.timestampFromTime(time)
  
  def parseCellHeader(self, cell_header):
    # extract a time + date stamp from cell_header, which is
    # e.g. "Macintosh HD:Users:fred:acquisition:091807Ac1  12:02:47 PM"

    basename = cell_header.split("\t")[0].rpartition(":")[2]
    date = basename[0 : 6]
    time = cell_header.split("\t")[1]
   
    formatter = NSDateFormatter.alloc()\
      .initWithDateFormat_allowNaturalLanguage_('%m%d%y', None)

    parsed = formatter.dateFromString_(date)

    # all this footwork is needed because dates get wrapped around.  e.g.
    # having a "month" of 13 yields jan.
    if parsed == None or not int(date[0:2]) == parsed.monthOfYear()\
       or not int(date[2:4]) == parsed.dayOfMonth()\
       or not date[4:6] == str(parsed.yearOfCommonEra())[2:4]:

      raise Exception, "Data file filenames must start with MMDDYY "\
        + '(%s invalid)' % basename

    # CODE DEBT: this is implemented as a side-effect -- should be explicit
    self.currentExperimentDate = date
    
    stamp = self.dateFormatter.dateFromString_\
      (date + " " + time + " " + ADMNotebookFileParser.EXPERIMENT_TIMEZONE)
    
    if stamp == None:
      raise Exception, "Unable to extract timestamp from cell header: "\
        + cell_header
    
    return basename, stamp

  def initScanner(self):
    # lexical tokenizer for notebook files -- configured below.
    self.scanner = reflex.scanner('in-purpose')
    scanner = self.scanner

    Token = ADMNotebookFileParser.Token
    timestamp_re = ADMNotebookFileParser.TIMESTAMP_RE

    # matches the line written when a new cell file is started.
    # e.g. "Macintosh HD:Users:fred:acquisition:091807Ac2  2:33:38 PM"
    cell_header_re = '^\S[^:]+:.+\s+' + timestamp_re + '\\n'

    scanner.state('in-purpose')
    scanner.rule(cell_header_re, token = Token.Cell_Header,
      tostate = 'in-cell')
    scanner.rule('.*\\n', token = Token.Purpose_Line)

    # a string of asterisk denote the beginning and end of comments,
    # which can appear any time outside of epoch families.  most have
    # 50, but we are lenient.
    comment_delimit_re = '^' + ''.join('\\*' for i in range(20)) + '.+\\n'

    # these setting change related comments can appear alone in a cell or
    # at the tail end of an epoch family.
    # e.g. "*Green LED Mean  = 0 V   1:55:48 PM"
    tail_comment_re = '^\s*\\*.+\s+' + timestamp_re + '\\n'

    # these single asterisk delimited comments appear before epochs are
    # listed in a family. e.g. "*LED Pulse using Green LED*" or, in the
    # case of gclamp, "*Conductance clamp step, exc amp = 0 exc mean = 0
    # inh amp = 25"
    protocol_comment_re = '^\\*[^*]+\\*?\\n'

    # epoch families start this way. e.g. "  USER  1:43:51 PM"
    family_start_re = '^\s+[A-Z]+.+' + timestamp_re + '.*\\n'
    empty_line_re = '^\s*$'

    log.debug("family start re: " + family_start_re)

    scanner.state('in-cell')
    scanner.rule(cell_header_re, token = Token.Cell_Header)
    scanner.rule(tail_comment_re, token = Token.Tail_Comment)
    scanner.rule(protocol_comment_re, token = Token.Protocol_Comment)
    scanner.rule(comment_delimit_re, tostate = 'in-cell-comment')
    scanner.rule(family_start_re, tostate = 'in-epoch-family',
      token = Token.Family_Start)
    scanner.rule(empty_line_re)

    # all comment blocks start with a single line containing a timestamp
    comment_timestamp_re = '^' + timestamp_re + '\\n'

    scanner.state('in-cell-comment')
    scanner.rule(comment_delimit_re, tostate = 'in-cell')
    scanner.rule(comment_timestamp_re, token = Token.Comment_Timestamp)
    scanner.rule('.*\\n', token = Token.Comment_Line)

    epoch_re = '^\s+Epoch # \d+.*\\n'

    # if we see something being set equal to something on a line that's not an
    # epoch, treat it as one or more parameters.
    parameters_re = '\s*\w*:?\s*[^\t =]+\s*=.+\\n'

    scanner.state('in-epoch-family')
    scanner.rule(protocol_comment_re, token = Token.Protocol_Comment)
    scanner.rule(epoch_re, token = Token.Epoch)
    scanner.rule(tail_comment_re, token = Token.Tail_Comment)
    scanner.rule(parameters_re, token = Token.Parameters_Line)
    scanner.rule(empty_line_re, tostate = 'in-cell')
    scanner.rule('.+\\n')

# more support for running the parser independently for dev purposes
if ( __name__ == '__main__' ):
  log.basicConfig(level = log.INFO)

  parser = ADMNotebookFileParser.parserForFile_\
    ('/Users/dacc/Desktop/101808-notes.2.txt')

  parser.parse()

#  for family in parser.cells()[0].epochFamilies():
#    print family

#  print len(parser.cells()[0].epochFamilies())

#  for note in parser.notes():
#    print str(note.timestamp()) + ": " + note.comment()

#  for family in parser.cells()[0].epochFamilies():
#    print family.startDate()

#  print parser.cells()[0]

#  print parser.cells()[0].epochFamilies()[0].parameters()
#  print parser.cells()[0].epochFamilies()[1].parameters()

