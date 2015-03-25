#-*- coding: UTF-8 -*- 

import re

############################################################
# Utils functions:
#			Useful common bioinformatic functions
############################################################

def ReverseComplete(seq):
	rSeq = seq[::-1]
	rcSeq = ''
	for p in rSeq :
		if(p == 'A'):
			rcSeq += 'T'
		elif(p == 'C'):
			rcSeq += 'G'
		elif(p == 'G'):
			rcSeq += 'C'
		elif(p == 'T'):
			rcSeq += 'A'
		else:
			rcSeq += p

	return rcSeq

############################################################
# Constants:
#			Constants about regular expressions for sam 
#			files
############################################################

# precompiled regular expression

blankLineRe = re.compile('^\s*$')
numRe = re.compile('^[-+]?[0-9]+$')
floatRe = re.compile('^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$')
printRe = re.compile('^[ !-~]+$')
hexRe = re.compile('^[0-9A-F]+$')
arrayRe = re.compile('^[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+$')

samHeaderRe = re.compile('(^@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$)|(^@CO\t.*)$')

qnameRe = re.compile('^[ -?A-~]{1,255}$')
intRe = re.compile('^[0-9]+$')
rnameRe = re.compile('^\*|[!-()+-<>-~][!-~]*$')
cigarRe = re.compile('^\*|([0-9]+[MIDNSHPX=])+$')
rnextRe = re.compile('^\*|=|[!-()+-<>-~][!-~]*$')
tlenRe = re.compile('^-?[0-9]+$')
seqRe = re.compile('^\*|[A-Za-z=.]+$')
qualRe = re.compile('^[!-~]+$')
optionRe = re.compile('^[A-Za-z][A-Za-z0-9]:[AifZHB]:[!-~]+$')

mdRe = re.compile('^[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*$')
tagNumRe = re.compile('^[0-9]+$')
mdTagRe = re.compile('^[0-9]+|\^[A-Z]+|[A-Z]+$')
cigarTagRe = re.compile('^[MIDNSHPX=]$')

############################################################
# Class： 	AlignmentOptional
#       	Alignment Optional Field Object
############################################################

class SamAlignmentTag(object):
	def __init__(self):
		self.type = ''
		self.value = None

	def parse(self, tagStr):
		fields = tagStr.split(':')
		if(len(fields) != 3):
			return False
		self.type = fields[1]
		if(self.type == 'A'):
			self.value = fields[2]
		elif(self.type == 'i'):
			if(not numRe.match(fields[2])):
				return False
			self.value = int(fields[2])
		elif(self.type == 'f'):
			if(not floatRe.match(fields[2])):
				return False
			self.value = float(fields[2])
		elif(self.type == 'Z'):
			if(not printRe.match(fields[2])):
				return False
			self.value = fields[2]
		elif(self.type == 'H'):
			if(not hexRe.match(fields[2])):
				return False
			self.value = int(fields[2], 16)
		elif(self.type == 'B'):
			if(not arrayRe.match(fields[2])):
				return False
			# parse B type values later
			self.value = fields[2]
		else:
			return False
		return True

############################################################
# Class:	Alignment
#			Alignment record for sam file 
############################################################

class SamAlignment(object):
	def __init__(self):
		# Qname: Query template name, string, [!-?A-~]{1,255}
		self.qname = ''
		
		# Flag: Bitwise flag, int, [0,2^16 -1]
		self.flag = -1
		
		# Rname: Reference sequence name, string, \*|[!-()+-<>-~][!-~]*
		self.rname = ''
		
		# Pos: 1-based leftmost mapping position, int, [0,2^31 -1]
		self.pos = -1
		
		# MapQ: Mapping Quality, int, [0,2^8 -1]
		self.mapq = 255
		
		# CIGAR: CIGAR string, string, \*|([0-9]+[MIDNSHPX=])+
		self.cigar = ''
		
		# Rnext: Ref. name of the mate/next read, string, \*|=|[!-()+-<>-~][!-~]*
		self.rnext = ''
		
		# Pnext: Position of the mate/next read, int, [0,2^31 -1]
		self.pnext = -1
		
		# Tlen: Observed template length, int, [-2^31 +1,2^31 -1]
		self.tlen = 0
		
		# Seq: Segment Sequence, string, \*|[A-Za-z=.]+
		self.seq = ''
		
		# Qual: ASCII of Phred-scaled base Quality+33, string, [!-~]+
		self.qual = ''

		# Optional: optional fields

		self.tags = {}

	def parse(self, astr):
		fields = astr.split('\t')

		# there are 11 mandatory fields in a valid aligment record 

		if(len(fields) < 11):
			return False

		# validate qname field

		if(not qnameRe.match(fields[0])):
			self.qname = ''
			return False
		self.qname = fields[0].strip()

		# validate flag field

		if(not intRe.match(fields[1])):
			self.flag = -1
			return False
		self.flag = int(fields[1])
		if(self.flag >= 0xFFFFFF):
			self.flag = -1
			return False
		
		# validate rname field

		if(not rnameRe.match(fields[2])):
			self.rname = ''
			return False
		self.rname = fields[2].strip()

		# validate pos field
		
		if(not intRe.match(fields[3])):
			self.pos = -1
			return False
		self.pos = int(fields[3])
		if(self.pos >= 0x7FFFFFFF):
			self.pos = -1
			return False

		# validate mapq field
		
		if(not intRe.match(fields[4])):
			self.mapq = 255
			return False
		self.mapq = int(fields[4])
		if(self.mapq > 0xFF):
			self.mapq = 255
			return False

		# validate cigar field

		if(not cigarRe.match(fields[5])):
			self.cigar = ''
			return False
		self.cigar = fields[5].strip()

		# validate rnext field

		if(not rnextRe.match(fields[6])):
			self.rnext = ''
			return False
		self.rnext = fields[6].strip()

		# validate pnext field

		if(not intRe.match(fields[7])):
			self.pos = -1
			return False
		self.pnext = int(fields[7])
		if(self.pnext >= 0x7FFFFFFF):
			self.pos = -1
			return False

		# validate tlen field

		if(not tlenRe.match(fields[8])):
			self.tlen = 0
			return False
		self.tlen = int(fields[8])

		# validate seq field

		if(not seqRe.match(fields[9])):
			self.seq = ''
			return False
		self.seq = fields[9].strip()

		# validate qual field

		if(not qualRe.match(fields[10])):
			self.qual = ''
			return False
		self.qual = fields[10].strip()

		# Optional fields

		i = 11
		while i < len(fields):
			if(not optionRe.match(fields[i])):
				return False
			else:
				# extract tag
				tag = fields[i][:2]
				tagContent = SamAlignmentTag()
				if(tagContent.parse(fields[i])):
					self.tags[tag] = tagContent
				else:
					return False
			i = i + 1

		return True

	# Adjusted sequence by CIGAR value

	def seqByCigar(self):
		if(self.cigar == ''):
			return self.seq
		
		nums = tagNumRe.findall(self.cigar)
		tags = cigarTagRe.findall(self.cigar)
		numLen = len(nums)
		tagLen = len(tags)
		if(numLen != tagLen):
			return ''

		orgSeq = list(self.seq)
		orgSeq.reverse()
		adjSeq = ''
		try:
			for i in range(tagLen):
				if(tags[i] == 'M'):
					for j in range(int(nums[i])):
						adjSeq += orgSeq.pop()
				elif(tags[i] == 'I'):
					for j in range(int(nums[i])):
						orgSeq.pop()
				elif(tags[i] == 'D'):
					for j in range(int(nums[i])):
						adjSeq += 'N'
				elif(tags[i] == 'N'):
					for j in range(int(nums[i])):
						adjSeq += 'N'
				elif(tags[i] == 'S'):
					for j in range(int(nums[i])):
						orgSeq.pop()
				elif(tags[i] == '='):
					for j in range(int(nums[i])):
						adjSeq += orgSeq.pop()
				elif(tags[i] == 'X'):
					for j in range(int(nums[i])):
						adjSeq += orgSeq.pop()
		except IndexError:
			return ''

		return adjSeq

	# Adjusted sequence by MD value

	def seqByMD(self):
		if(not ('MD' in self.tags)):
			return self.seq

		mdstr = self.tags['MD'].value
		tags = mdTagRe.findall(mdstr)

		orgSeq = list(self.seq)
		orgSeq.reverse()
		adjSeq = ''
		try:
			for tag in tags:
				if(tag.isdigit()):
					for j in range(int(tag)):
						adjSeq += orgSeq.pop()
				elif(tag.isalpha()):
					for j in range(len(tag)):
						adjSeq += orgSeq.pop()
				elif(tag[0] == '^'):
					for j in range(len(tag) - 1):
						adjSeq += 'N'
		except IndexError:
			return ''

		return adjSeq

	# return alignment

	def str(self):
		alignmentStr = (self.qname + '\t' + str(self.flag) + '\t' + self.rname + 
						'\t' + str(self.pos) + '\t' + str(self.mapq) + '\t' + self.cigar +
						'\t' + self.rnext + '\t' + str(self.pnext) + '\t' + str(self.tlen) +
						'\t' + self.seq + '\t' + self.qual + '\t')
		for tag in self.tags :
			alignmentTag = self.tags[tag]
			alignmentStr = alignmentStr + tag + ':' + alignmentTag.type + ':' + str(alignmentTag.value)
			alignmentStr = alignmentStr + '\t'
		return alignmentStr


############################################################
# Class:	Sam Header
#			Sam header class
############################################################

class SamHeader(object):
	def __init__(self):
		self.tag = ''
		self.value = ''
		# self.value = {}
	
	def parse(self, tagStr):
		if(not samHeaderRe.match(tagStr)):
			return False
		self.tag = tagStr[:3]
		self.value = tagStr[3:].strip()
		# fields = tagStr[3:].split()
		# i = 0
		# while (i < len(fields)):
		# 	value = fields[i].split(':')
		# 	try:
		# 		self.value[value[0]] = value[1]
		# 	except:
		# 		self.value[value[0]] = ''
		# 	i = i + 1
		return True

	def str(self):
		headerStr = self.tag + '\t' + self.value
		return headerStr


############################################################
# Class:	Sam File
#			Sam File structure
############################################################

class Sam(object):
	def __init__(self):
		self.header = []
		self.aligment = []

	# Parse a line and store the parse result in Sam object

	def parseLine(self, line):
		if(blankLineRe.math(line)):

			# Ignor blank lines

			return True
		elif(line[0] == '@'):

			# Try to parse the line as header

			header = SamHeader()
			if(header.parse(line)):
				self.header.append(header)
				return True
			else:
				return False
		else:

			# Try to parse the line as alignment

			alignment = SamAlignment()
			if(alignment.parse(line)):
				self.alignment.append(alignment)
				return True
			
		return False

	# Parse a line and return the parse object, for external call

	def parseObject(self, line):
		if(blankLineRe.match(line)):

			# Ignor blank lines

			return None
		elif(line[0] == '@'):

			# Try to parse the line as header

			header = SamHeader()
			if(header.parse(line)):
				return header
		else:

			# Try to parse the line as alignment

			alignment = SamAlignment()
			if(alignment.parse(line)):
				return alignment

		return None