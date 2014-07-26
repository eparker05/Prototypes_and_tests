from xml.parsers.expat import ParserCreate

class Element(object):
    """A simple to use element class for indexing"""
    def __init__(self, name, begin=None, end=None):
        self.parent = None
        self.name = name
        self.text = ""
        self.attributes = None
        self.children = []
        self.indexbegin = begin
        self.indexend = end
        
    def add_child(self, child):
        child.parent = self
        self.children.append(child)
    
    def get_all_children_by_name(self,name):
        validchildren = []
        for child in self.children:
            if child.name == name:
                validchildren.append(child)
            validchildren.extend(child.get_all_children_by_name(name))
        return validchildren
            
    def first_child(self):
        if not self.children:
            return None
        else:
            return self.children[0]
            
    def __repr__(self):
        return "< Element name={}, position=({},{}) >". \
                format(self.name, self.indexbegin, self.indexend)
                
    def depth(self):
        """this is not efficient and is mostly for debugging"""
        if self.parent is None:
            return 0
        else:
            return self.parent.depth() + 1

class ExpatHandler(object):
    """ExpatHandler class will return an indexed Element tree
    
    The  targetfield attribute is the fundamental unit of indexing,
    these xml tags must be a sequential list with no tags in between
    having a different identity. The namestoparse list contains tags
    that will be extracted into the element tree.
    
    ExpatHandler assumes compilant well formmated XML, several types
    of formatting errors will result in difficult to decipher 
    errors while other sorts of errors will not be detected. Best
    practice is to use a secondary parser for actual parsing and
    validation of data.
    """
    
    def __init__(self, handle, parser_class=ParserCreate):
        #set up parser
        self._handle = handle
        self._parser_class = parser_class
        
        self.targetfield = "PLANT"
        self.namestoparse = ["PLANT", "NAMES", "N", "ZONE", "PRICE"] 
    
    def parse_from_position(self, position=0):
        handle = self._handle
        #make the parser
        parser = self._parser = self._parser_class()
        parser.StartElementHandler = self.start_element
        parser.EndElementHandler = self.end_element
        parser.CharacterDataHandler = self.char_data
        handle.seek(position)
        self.baseposition = position
        
        rootelem = Element(name="ROOT", begin=position)
        self.rootelem = rootelem
        self.currentelem = rootelem
        self.target_tag_met = False
        
        self.tags = {}
        self.tagcounts = {}
        try:
            parser.ParseFile(handle)
        except StopIteration:
            return rootelem
        
        #A return should have happened at this point
        raise ValueError("Check that file contains target element")
        
    def start_element(self, name, attrs):
        print("{}new name {}".format("-"*(self.currentelem.depth()+1), name))
        if self.currentelem.indexend is True:
            self._finish_element()
            
        if name in self.namestoparse:
            self.savetext = True
            byteindex = self._parser.CurrentByteIndex + self.baseposition
            newelement = Element(name, begin=byteindex)
            newelement.attributes = attrs
            self.currentelem.add_child(newelement)
            self.currentelem = newelement
        else:
            self.savetext = False
        
    def end_element(self, name):
        if name == self.targetfield:
            end = self._parser.CurrentByteIndex + len(self.targetfield) + 3 \
                  + self.baseposition
            self.currentelem.indexend = end
            self.rootelem.indexend = end
            raise StopIteration()
        print("{}end name {}".format("-"*(self.currentelem.depth()+1), name))
        if self.currentelem.indexend is True:
            self._finish_element()
        if name == self.currentelem.name:
            self.currentelem.indexend = True        

    def char_data(self, data):
        if data.strip() and self.savetext:
            self.currentelem.text += data.strip()
            
    def _finish_element(self):
        """ any element eligible for finishing is saved here
        
        An element has ended; fix the end byte index and fetch the parent node.
        This will possibly also stop iteration if the end of the target node is found"""
        assert self.currentelem.indexend is True
        self.currentelem.indexend = self._parser.CurrentByteIndex + self.baseposition
        self.currentelem = self.currentelem.parent
    

#open file in binary mode for robuster byte offsets.
xmlFile = open("simple.xml", 'rb')

h = ExpatHandler(xmlFile)
a = h.parse_from_position()
print(repr(a))
for c in a.children:
    print("child {} exists".format(c.name))
    print("child {} has {} children and 1 parents".format(c.name, len(c.children)))

print("\n\nnext section\n\n")
a = h.parse_from_position(c.indexend)
print(repr(a))
for c in a.children:
    print("child {} exists".format(c.name))
    print("child {} has {} children and 1 parents".format(c.name, len(c.children)))