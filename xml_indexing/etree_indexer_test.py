import xml.etree.ElementTree as etree

class IterParseIterator(etree._IterParseIterator):
    """Implements the etree iterparse with indexing
    
    This works by reading one byte at a time, feeding 
    that byte to the parser, and finally checking if
    the parser has registered any events via the 
    read_events() call.
    
    This only generates xml tag beginnings and ends, 
    the resulting data structure must still be produced
    manually.
    """
    
    def __init__(self, *args, **kwargs):
        self.tag = ""   # addition
        self._old_position = 0 # addition
        self.position = 0   # addition
        etree._IterParseIterator.__init__(self, *args, **kwargs)
        
    def __next__(self):
        self._old_position = self.position # addition
        currentposition = self.position # addition
        while 1:
            for event in self._parser.read_events():
                self.position = currentposition
                return event
            if self._parser._parser is None:
                self.root = self._root
                if self._close_file:
                    self._file.close()
                raise StopIteration
            # load event buffer
            data = self._file.readline()  # alteration
            #data = self._file.read(1)
            currentposition += len(data) # addition
            if data:
                self._parser.feed(data)
            else:
                self._root = self._parser._close_and_return_root()

    next = __next__

    def read_prev_line(self):
        readlen = self.position - self._old_position
        self._file.seek(-1*readlen, 1)
        line = self._file.read(readlen)
        assert self._file.tell() == self.position
        return line


if __name__ == "__main__":
    #initialize xml file
    xmlfile = open("simple.xml", 'rb')
    a = IterParseIterator(xmlfile, events=('start', 'end'), parser=None, close_source=True)
    print("initial position = {}".format(a._file.tell()))
    
    while True:
        out = next(a)
        line = a.read_prev_line()
        print("foundline;  {}".format(repr(line)))
        print("current event {}; elem {} @ {}\n".format(out[0], out[1].tag, a._file.tell()))
        b = input("press any key to end iteration: ")
        if b:
            break