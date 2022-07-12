
"""
#####################################################################
Class to read the usual outputs from SmallDataIO C++ GRChombo class
#####################################################################

The expected format is one or more rows of headers (starting with #)
and then subsequent rows of coordinates and data

You can do a bunch of things with this File class. You can add files, subtract
files, complexify columns and a bunch of other stuff you'll probably end up
doing manually because it takes less time than understand how to use this.

Nomenclature:
- FileSet is a set of files (e.g. all 'coords' files of the AH)
- File is well, a .dat file
- Block: each file is composed of blocks. Most files have only 1 block. Except
the WeylExtractionOut files, which have multiple blocks. These files have
one block per extraction radii (each block being all the 2D spherical points).
Again, most files have only 1 block.

For these "Extraction" dump type of files, the data structure of blocks is
("header1|data1 \n\n header2|data2 ...")
Here assuming that a single empty line is what is dividing blocks
(in reality files typically have 2 of dividing empty lines between blocks)

Usage File class:

filename = "data/Weyl_integral_64.dat"
reader = File(filename)
print(reader.getHeaders())
print(reader.getData()[0:2])
print(len(reader.blocks))
print(len(reader.getHeaders()))
print(len(reader.getData()))

Usage FileSet class:

file_prefix = "data/Weyl_integral_*.dat"
reader = FileSet(file_prefix)
print(reader.getFile(0).getName())
print(reader.getFile(0).wasRead())
print(reader.getFile(0).getHeaders())

"""

import glob # to read directories based on a regex expression
import operator # do things like operator.iadd(a,b)
import copy # use deepcopy
import numpy # turn lists into numpy arrays
import types # to get types.FunctionType

class Block:
    """
    Contains headers and data of a block of data of SmallDataIO files
    not to use directly. Reading to a 'Block' is managed by 'read()' function of 'File'
    """
    
    ### PUBLIC ###
    
    def __init__(self): 
        """constructor"""
        self._headers = []
        self._data = []
        self._num_labels = None # undefined
        self._headerMap = {}

    def getData(self):    return self._data
    def getHeaders(self): return self._headers
    def getHeaderMap(self): return self._headerMap

    def numRows(self):    return len(self._data)
    def numColumns(self): return len(self._data[0])
    def numLabels(self):  return self._num_labels
    def numValues(self):  return self.numColumns() - self.numLabels()

    def complexifyColumns(self, col_start : int = -1, col_end : int = -1):
        """Complexify columns in pairs, useful for Weyl4 (by default do from the 1st 'value' - non label - column til the end)"""
        if col_start < 0: col_start = self.numLabels()
        if col_end   < 0: col_end   = self.numColumns() - 1

        assert isinstance(col_end, int) and isinstance(col_start, int)
        assert col_end >= 0 and col_end < self.numColumns()
        assert col_start >= 0 and col_start < self.numColumns()
        assert (col_end - col_start + 1) % 2 == 0 # we complexify every pair
        columns = self.transposeData(self._data)
        newColumns = []
        for col in range(col_start):
            newColumns.append(columns[col])
        for col in range(col_start, col_end+1, 2):
            newColumns.append([(columns[col][row] + 1j * columns[col+1][row]) for row in range(len(columns[col]))])
        for col in range(col_end+1, self.numColumns()):
            newColumns.append(columns[col])

        # transpose back to row-like list
        self._data = self.transposeData(newColumns)
        # self.toNumpy() # numpy converts everything to complex...
        return self.getData()

    def addColumn(self, formula : types.FunctionType):
        rows = []
        for row in self:
            if isinstance(row, numpy.ndarray):
                row = row.tolist()
            row.append(formula(row, self._headerMap))
            rows.append(row)
        self._data = rows

    def removeColumn(self, col_start : int, col_end : int = -1):
        assert col_start >= self.numLabels() and col_start < self.numColumns()
        if col_end < 0: col_end = col_start
        assert col_end >= self.numLabels() and col_end < self.numColumns()

        rows = []
        for row in self:
            if isinstance(row, numpy.ndarray):
                row = row.tolist()
            for i in range(col_start, col_end+1):
                del row[col_start] # it is always at the same index, because we delete it
            rows.append(row)
        self._data = rows

    def transposeData(self, data):
        return [list(col) for col in zip(*data)]

    def findHeader(self, label : str):
        """Return (row, column) index of a given header that contains 'label'"""
        for key, val in self._headerMap.items():
            if label in key:
                return val
        return None

    def getHeaderValue(self, label : str):
        """If a header is written as 'label = value', retrieve value"""
        index = self.findHeader(label)
        if index is None:
            return None
        header = self._headers[index[0]][index[1]]
        if '=' not in header:
            return None
        return float(header.split('=')[1].strip())

    def __len__(self):               return self.numRows()  # ask for 'len(Block)'
    def __getitem__(self, idx : int):      return self._data[idx] # ask for row 'i' as 'block[i]'
    def __setitem__(self, idx: int, val : float): self._data[idx] = val  # 'block[idx] = val'

    # allow to use operators +-*/ between blocks (does not change the 'labels')
    # 'iadd' refers to the '+=' operator (that change the object itself), likewise for the other
    def __add__(self, other):      return self.__operator(operator.add,      other, assign = False)
    def __iadd__(self, other):     return self.__operator(operator.iadd,     other, assign = True)
    def __radd__(self, other):     return self.__operator(operator.add,     other, assign = False)
    def __sub__(self, other):      return self.__operator(operator.sub,      other, assign = False)
    def __isub__(self, other):     return self.__operator(operator.isub,     other, assign = True)
    def __mul__(self, other):      return self.__operator(operator.mul,      other, assign = False)
    def __imul__(self, other):     return self.__operator(operator.imul,     other, assign = True)
    def __rmul__(self, other):     return self.__operator(operator.mul,     other, assign = False)
    def __truediv__(self, other):  return self.__operator(operator.truediv,  other, assign = False)
    def __itruediv__(self, other): return self.__operator(operator.itruediv, other, assign = True)

    def addLine(self, line : str):
        """Used in 'File' class below"""
        line = line.strip() # just to make sure
        self.__add_header_line(line) if line.startswith('#') else self.__add_data_line(line)

    def toNumpy(self):
        """
        Used in 'File' class below
        Convert to numpy at the end of writing all lines
        It's easier to convert to numpy array only at the end
        because otherwise we would need to know in advance how many rows and columns the block has
        """
        self._data = numpy.array(self._data)

    ### PRIVATE ###

    def __add_header_line(self, line):
        line = line[1:] # remove '#'
        # split with >=2 whitespaces such that pieces as 'r = ' or 'time =' stay together
        # (there is a risk of heaving 2 very big headers, change separator to 'None' if you have that problem)
        self._headers.append( [s.strip().replace(',','') for s in line.split('  ') if s] )
        row = len(self._headers) - 1
        for c, col in enumerate(self._headers[row]):
            self._headerMap[col] = (row, c)

    def __add_data_line(self, line):
        if self._num_labels == None: # count number of labels if not yet done
            self._num_labels = self.__get_num_labels_from_line(line)
        self._data.append( [float(s.strip()) for s in line.split() if s] )

    def __get_num_labels_from_line(self, line): # retrieve number of labels from a line
        all_lengths = [(len(s.strip(' -'))) for s in line.split() if s] # width of each number, remove minus signs
        lengths = list(set(all_lengths)) # count how many of each
        if len(lengths) == 2: # labels have only length, data another length
            return all_lengths.count(min(lengths)) # assume coordinates have the minimum length
        else:
            print("Labels could not be detected from width of columns (maybe first line has nans)")
            return 0

    def __operator(self, oper, obj, assign): # generic operator
        if isinstance(obj, Block):                    return self.__operator__Block(oper, obj, assign)
        elif hasattr(obj, "__len__"):                 return self.__operator__ListOfConstants(oper, obj, assign)
        elif isinstance(obj, (int, float, complex)):  return self.__operator__Constant(oper, obj, assign)
        else:                                         raise TypeError

    def __operator__Constant(self, oper, constant, assign): # generic operator that receives a constant
        return self.__operator__ListOfConstants(oper, [constant]*self.numValues(), assign)

    def __operator__ListOfConstants(self, oper, listOfConstants, assign):
        """Generic operator that receives a list / array"""
        assert len(listOfConstants) == self.numValues()
        num_labels = self.numLabels()
        num_columns = self.numColumns()

        out = self if assign else copy.deepcopy(self)
        for row in out:
            for col in range(num_labels, num_columns):
                row[col] = oper(row[col], listOfConstants[col - num_labels])

        return out

    def __operator__Block(self, oper, other_block, assign):
        """Generic operator that receives another Block"""
        num_labels = self.numLabels()
        num_columns = self.numColumns()
        assert other_block.numLabels() == num_labels
        assert other_block.numColumns() == num_columns
        assert other_block.numRows() == self.numRows()

        out = self if assign else copy.deepcopy(self)
        for row, other_row in zip(out._data, other_block):
            for col in range(num_labels, num_columns):
                row[col] = oper(row[col], other_row[col])

        out.headers = [] # doesn't make sense to mix headers
        return out

class File:
    """Contains all the data of a SmallDataIO file, as a list of blocks"""

    ### PUBLIC ###

    def __init__(self, filename : str, read : bool = True):
        self.filename = filename
        self.blocks = []
        self._was_read = False
        if read: self.read()

    def numBlocks(self): return len(self.blocks)
    def numLabels(self): return self.getBlock(0).numLabels() # assume same in all blocks
    def numValues(self): return self.getBlock(0).numValues() # assume same in all blocks
    def wasRead(self):   return self._was_read
    def getName(self):   return self.filename
    def getBlock(self, idx = 0):
        assert self._was_read
        return self.blocks[idx]

    # methods assuming there is only 1 block
    def getData(self):      return self.getBlock(0).getData()
    def getHeaders(self):   return self.getBlock(0).getHeaders()
    def getHeaderMap(self): return self.getBlock(0).getHeaderMap()
    def numRows(self):      return self.getBlock(0).numRows()
    def numColumns(self):   return self.getBlock(0).numColumns()

    def transposeData(self, data):         return self.getBlock(0).transposeData(data)
    def findHeader(self, label : str):     return self.getBlock(0).findHeader(label)
    def getHeaderValue(self, label : str): return self.getBlock(0).getHeaderValue(label)

    def complexifyColumns(self, col_start : int = -1, col_end : int = -1):
        """Complexify columns in pairs, useful for Weyl4 (by default do from the 1st 'value' - non label - column til the end)"""
        for block in self.blocks:
            block.complexifyColumns(col_start, col_end)

    def addColumn(self, formula : types.FunctionType):
        for block in self.blocks:
            block.addColumn(formula)

    def removeColumn(self, col_start : int, col_end : int = -1):
        for block in self.blocks:
            block.removeColumn(col_start, col_end)

    def removeBlock(self, col_start : int, col_end : int = -1):
        assert col_start >= 0 and col_start < self.numBlocks()
        if col_end < 0: col_end = col_start
        assert col_end >= 0 and col_end < self.numBlocks()

        for i in range(col_start, col_end+1):
            del self.blocks[col_start] # it is always at the same index, because we delete it

    def __len__(self):               return self.numBlocks()   # ask for 'len(file)'
    def __getitem__(self, idx):      return self.getBlock(idx) # ask for block 'i' as 'file[i]'
    def __setitem__(self, idx, val):
        assert self._was_read
        self.blocks[idx] = val

    # allow to use operators +-*/ between files (does not change the 'labels')
    def __add__(self, other):      return self.__operator(operator.add,      other, assign = False)
    def __iadd__(self, other):     return self.__operator(operator.iadd,     other, assign = True)
    def __radd__(self, other):     return self.__operator(operator.radd,     other, assign = False)
    def __sub__(self, other):      return self.__operator(operator.sub,      other, assign = False)
    def __isub__(self, other):     return self.__operator(operator.isub,     other, assign = True)
    def __mul__(self, other):      return self.__operator(operator.mul,      other, assign = False)
    def __imul__(self, other):     return self.__operator(operator.imul,     other, assign = True)
    def __rmul__(self, other):     return self.__operator(operator.mul,     other, assign = False)
    def __truediv__(self, other):  return self.__operator(operator.truediv,  other, assign = False)
    def __itruediv__(self, other): return self.__operator(operator.itruediv, other, assign = True)
    
    def read(self):
        assert not self._was_read
        self.blocks = [] # reset
        found_new_block = True # start by adding 1st Block
        with open(self.filename) as file:
            for line in file:
                line = line.strip()
                if line: # if line is not empty
                    if found_new_block:
                        self.blocks.append(Block())
                        found_new_block = False
                    self.blocks[-1].addLine(line)
                else:
                    found_new_block = True
        # convert everything to numpy at the end
        for b in range(self.numBlocks()):
            self.blocks[b].toNumpy()
        self._was_read = True

    ### PRIVATE ###

    # def __lines_per_block(self): # it works, just not needed
    #     count = 0
    #     with open(self.filename) as file:
    #         for line in file:
    #             line = line.strip()
    #             if not line:                 break
    #             if not line.startswith('#'): count += 1
    #     return count

    def __operator(self, oper, obj, assign): # generic operator
        if isinstance(obj, File):                    return self.__operator__File(oper, obj, assign)
        if isinstance(obj, Block):                    return self.__operator__Block(oper, obj, assign)
        elif hasattr(obj, "__len__"):                 return self.__operator__ListOfConstants(oper, obj, assign)
        elif isinstance(obj, (int, float, complex)):  return self.__operator__Constant(oper, obj, assign)
        else:                                         raise TypeError

    def __operator__Constant(self, oper, constant, assign): # generic operator that receives a constant
        return self.__operator__ListOfConstants(oper, [constant]*self.numValues(), assign)

    def __operator__ListOfConstants(self, oper, listOfConstants, assign):
        """Generic operator that receives a list / array"""
        out = self if assign else copy.deepcopy(self)
        for b, block in enumerate(out.blocks):
            out.blocks[b] = oper(block, listOfConstants)
        return out

    def __operator__Block(self, oper, other_block, assign):
        """Generic operator that receives another Block"""
        out = self if assign else copy.deepcopy(self)
        for b, block in enumerate(out.blocks):
            out.blocks[b] = oper(block, other_block)
        return out

    def __operator__File(self, oper, other_file, assign):
        """Generic operator that receives another File"""
        out = self if assign else copy.deepcopy(self)
        for b, block, other_block in zip(range(self.numBlocks()), out.blocks, other_file):
            out.blocks[b] = oper(block, other_block)
        out.filename = "MixedFile" # doesn't make sense to keep either filename
        return out

class FileSet:
    """List of Files, reads a set of SmallDataIO files based on a prefix (can be just 1 as well)"""

    ### PUBLIC ###

    def __init__(self, file_prefix = None, read : bool = True, skip : int = 1, files=None):
        assert file_prefix is not None or files is not None
        if files is None:
            file_list = glob.glob(file_prefix)
            file_list.sort()
            self.files = [File(file_list[index], read = read) for index in range(0, len(file_list), skip)]
        else:
            self.files = files

    def numFiles(self):         return len(self.files)
    def getFile(self, idx = 0): return self.files[idx]

    def __len__(self):          return self.numFiles()   # ask for 'len(reader)'
    def __getitem__(self, idx): return self.getFile(idx) # ask for file 'i' as 'reader[i]'
    def __setitem__(self, idx, file):
        self.files[idx] = file

    # def complexifyColumns(self, col_start : int = -1, col_end : int = -1):
    #     """Complexify columns in pairs, useful for Weyl4 (by default do from the 1st 'value' - non label - column til the end)"""
    #     for file in self.files:
    #         file.complexifyColumns(col_start, col_end)

    # def addColumn(self, formula : types.FunctionType):
    #     for file in self.files:
    #         file.addColumn(formula)

    # def removeColumn(self, col_start : int, col_end : int = -1):
    #     for file in self.files:
    #         file.removeColumn(col_start, col_end)

    # def removeBlock(self, col_start : int, col_end : int = -1):
    #     for file in self.files:
    #         file.removeBlock(col_start, col_end)
