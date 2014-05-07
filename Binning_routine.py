# There's a good chance that quite a bit of this code will eventually become
# part of my biopython feature branch, so I haven't been shy about commenting,
# formatting, and unittests. 
#
# Copyright Evan Parker 2014
# this code is released under the the Biopython license 
# see http://www.biopython.org/DIST/LICENSE for the complete license


from math import ceil, floor, log
import sys

#set up integer_types variable for interpreter neutral code
#that accomodates integers larger than 2**31
if sys.version[0] == '2':
    integer_types = (int, long)
else:
    integer_types = (int,)


class FeatureBinCollection(object):
    """this class manages the creation and maintenance of feature indices

       This class is used to organize feature data in a qucikly retrievable
       data structure. The feature data must be added as a tuple containing
       at least two indicies: first anotated residue and the last as a half
       open half closed interval [first, last). The indices are assumed to be
       the first two elements of the stored tuple, but they may be re-asigned
       on instantiation via the beginindex and endindex kwarks.

       EXAMPLE
       -------
       defined below is a 3-tuple format of (beginindex, endindex, fileidx)
       three features are added to a newly initialized featurebin 

       >>> ft0 = (5574, 5613, 2300) 
       >>> ft1 = (0, 18141, 1300 )
       >>> ft2 = (5298, 6416, 3540)
       >>> featurebin = FeatureBinCollection()
       >>> featurebin.insert( ft0 )
       >>> featurebin.insert( ft1 )
       >>> featurebin.insert( ft2 )
       >>> len(featurebin)
       3
       
       Now that the 'featurebin' instance has some features, they can be
       retrieved with a standard getter using single integer indices or
       slice notation.

       >>> featurebin[1]
       [(0, 18141, 1300)]
       >>> sliceresult = featurebin[5200:5300]
       >>> sliceresult.sort()
       >>> sliceresult
       [(0, 18141, 1300), (5298, 6416, 3540)]


       BACKGROUND:
       -----------
       The basic idea of using feature bins is to group features into 
       bins organized by their span and sequence location. These bins then allow
       only likely candiate features to be queried rather than all features. The 
       example below illustrated with Figure 1 shows a similar scheme where feature1 
       is stored in bin-0, feature2 in bin-4 and feature3 in bin-2. Each sequence is
       stored in the smallest bin that will fully contain the sequence. A query of 
       all features in the region denoted by query1 could be quickly performed by 
       only searching through bins 0, 2, 5, and 6. Were this data structure many 
       levels deep, the performance savings would be large

       ___Figure 1_________________________________________________
       |                                                           |
       |    feature1  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~              |
       |    feature2  |       ~~~~                  |              |
       |    feature3  |       |  |               ~~~~~~~~~~~~~     |
       |              |       |  |               |  |        |     |
       | bins:        |       |  |               |  |        |     |
       |    0_________|_______|__|_______._______|__|____.___|_    |
       |    1_________________|__|__   2_:_______|_______:___|_    |
       |    3__________  4____|__|__   5_:________  6____:_____    |
       |                                 :               :         |
       |                                 :               :         |
       |    query1                       [ ? ? ? ? ? ? ? ]         |
       |...........................................................|

       Further reading on the math behind the idea can be found in: 
           Journal:  Bioinformatics Vol 27 no. 5 2011, pages 718-719
           Article:  "Tabix: fast retrieval of sequence features from generic
                      tab delimited files"
           Author:   Heng Li

       The implementation by Li has, as its largest bin ~500 million (2^29) and its smallest
       bin ~16000 (2^14). Each level of binning is separated by a factor of 8 (2^3).
       The implementation herein abandons a static binning scheme and instead
       starts with the smallest and largest bins as 256 and 8 million respectively. 
       These bins can then be dynamically expanded increasing by a factor of 8
       every time new data is found to be larger than the largest bin. As a practical
       matter of sanity checking, bin sizes are capped at 200 billion basepairs (2^41).

       Under some circumstances the exact size of a sequence and all related annotations
       is known beforehand. If this is the case the length kwarg allows the binning object
       to be solidified on instantiation at the correct length.
       
       This structure knows nothing about the global sequence index and is indexed 
       at zero. Any index transformation must be done at a higher level. It is important
       that all sequences and features stored here are indexed to zero.
       """
    
    def __init__(self, length = None, beginindex=0, endindex=1):
        """ initialize the class and set standard attributes

        kwargs:

        length:
            when length == None, the bins are dynamically sized.
            when length is a positive integer, the appropriate bin 
            size is selected and locked. Exceeding this value will
            cause exceptions when the max bin size is locked

        beginindex:
            the index of the first residue within the tuple that will
            be stored with the FeatureBinCollection.

        endindex:
            the index of the last residue (as a open interval) inside 
            the tuple that will be stored with FeatureBinCollection
        """
        
        #these should not be changed
        self._bin_level_count = 6
        self._bins = [[] for i in range(37449)]

        # this defines the indices of the begin and end sequence info
        # in the tuple structures stored in the bins
        self._beginindex = beginindex
        self._endindex = endindex

        #default action: start small (8M) and allow expansion
        self._dynamic_size = True
        if length is None:
            self._set_max_bin_power(23)
            
        #alternate action if a sequence length is provided
        # set to smallest power able to fully contain
        elif isinstance(length, integer_types) and length > 0:
            default_powers = [23,26,29,32,35,38,41]
            for power in default_powers:
                if length <= 2**power:
                    self._set_max_bin_power(power)
                    self._dynamic_size = False
                    break
            if self._dynamic_size: #this should have been set to False
                error_string = "Sequence length is {}: must be less than 2^41".format(length)
                raise ValueError(error_string)


        
    def _increase_bin_sizes(self):
        """increase max bin size 8x (2**3) and re-organize existing binned data
        
        In order to increase the total maximum bin size, the lowest set
        of bins must be merged up one level (first step) then the entire set
        must be moved down one level without disturbing the organization scheme.
        
        An assertion in this routine blocks sequences larger than 2**41 from
        being created.
        """
        
        oldsizepower = self._max_bin_power
        newsizepower = oldsizepower + 3
        assert newsizepower <= 41
        self._set_max_bin_power(newsizepower)
        
        # first, remove the lowest level
        # by merging it up to the previous level
        level = 5
        oL = int((2**(3*level) - 1)/7.0)
        new_level = 4
        new_oL = int((2**(3*new_level) - 1)/7.0)
        old_size = 2**(oldsizepower - 3*level)
        new_size = 2**(oldsizepower - 3*new_level)
        for k in range(4681, 37449):    
            bin_begin_old = (k-oL)*old_size
            k_new = int(floor(new_oL + (bin_begin_old/new_size)))
            #extend required to save existing data
            self._bins[k_new].extend(self._bins[k])
            self._bins[k] = []
        
        #then, move everything down.
        for k_inverse in range(4681):
            k = 4680 - k_inverse
            level = int( floor( log((7*k + 1),2)/3.0 ) )
            new_level = level + 1
            oL = int((2**(3*level) - 1)/7.0)
            new_oL = int((2**(3*new_level) - 1)/7.0)
            
            new_index = k - oL + new_oL 
            self._bins[new_index] = self._bins[k]
            self._bins[k] = []
            
        
    def _set_max_bin_power(self, power):
        """sets the maximum bin power and fixes other necessary attributes"""
        
        self._max_bin_power = power
        self._min_bin_power = self._max_bin_power - 3*(self._bin_level_count-1)
        self._size_list = [2**(self._min_bin_power+3*n) for n in range(self._bin_level_count)]
        self.maxSeqLength = self._size_list[-1]
    
    def insert(self, testTuple):
        """inserts a tuple with a sequence range into the feature bins
        
        data is assumed to be somewhat scrubbed, coming from a parser
        or a parser consumer."""
        
        beginindex = self._beginindex
        endindex = self._endindex


        begin = testTuple[beginindex]
        end = testTuple[endindex]
        assert isinstance(begin, integer_types)
        assert isinstance(end, integer_types)
        assert begin < end
        span = end-begin
        
        bin_index = self._calculate_bin_index(begin, span)
        self._bins[bin_index].append(testTuple)
        
    def __len__(self):
        return sum(len(bin) for bin in self._bins)    

    def __getitem__(self, key):

        #set some locals
        beginindex = self._beginindex
        endindex = self._endindex

        #any integers are just converted to a 'len() == 1' slice
        if isinstance(key, integer_types):
            key = slice(key, key+1)

        #check that it is a slice and it has no step property (or step==1)
        if not isinstance(key, slice):
            raise TypeError("lookups in the feature bin must use slice or int keys")
        if key.step is not None and key.step != 1:
            raise KeyError("lookups in the feature bin may not use slice stepping")

        #check that the key is within boundaries:
        if key.start < 0 or key.stop > self._size_list[-1]:
            raise IndexError("key out of bounds")
        if key.start > key.stop:
            raise IndexError("key not valid, slice.start > slice.stop")

        #code taken from self._calculate_bin_index(), comments removed
        possible_entries = []
        bin_level_count = self._bin_level_count
        max_bin_power = self._max_bin_power 
        for l_inverse in range(bin_level_count):
            L = bin_level_count - 1 - l_inverse
            oL = (2**(3*L) - 1)/7
            sL = float(2**(max_bin_power-3*L))
            k1 = int(floor(oL + (key.start/sL)))
            #k2 is incremented since range is used
            k2 = int(ceil(oL - 1 + (key.stop)/sL)) + 1
            for bin in range(k1,k2):
                possible_entries.extend(self._bins[bin])
        
        #check for bound and overlapping sequences
        return_entries = []
        for feature in possible_entries:
            #this covers fully bound sequence and left overlap   ssssssFsFsFsFFFFF
            if key.start <= feature[beginindex] < key.stop:
                return_entries.append(feature)
            #this covers left sequence right sequence overlap of (F)  FFFFFFFsFsFsFsssss
            elif key.start < feature[endindex] <= key.stop:
                return_entries.append(feature)
            #this covers features fully bound by a sequence    ssssssFsFsFsFsFssssss
            elif key.start <= feature[beginindex] and key.stop >= feature[endindex]:
                return_entries.append(feature)
            #this covers seqyebces fully bound by a feature      FFFFFFFsFsFsFsFsFFFFFFF
            elif key.start >= feature[beginindex] and key.stop <= feature[endindex]:
                return_entries.append(feature)
        return return_entries


    def _calculate_bin_index(self, begin,span):
        """ This function returns a bin index given a (begin, span) interval
        
        The equations for determination of bin index derived from this publication: 
           Journal:  Bioinformatics Vol 27 no. 5 2011, pages 718-719
           Article:  "Tabix: fast retrieval of sequence features from generic
                      tab delimited files"
           Author:   Heng Li
        
        This function should only be used privately in the context of having no
        easier relationship to assign bin index. Placing this in a loop for any
        task other than arbitrary assignments is a bad idea since many of the
        common tasks can provide bin index through other relationships.
        _increase_bin_sizes is an example of a routine that would suffer
        performance penalties were it to use this routine yet can be run
        efficiently using alternate bin index relationships.
        """
        
        #take care of bad data first
        assert span > 0   # no function may pass a span <= 0
        assert begin >= 0  # all indices must be positive
        if self._dynamic_size:
            assert begin+span <= 2**41   # len(seq) > 219 billion is not reasonable
        elif not self._dynamic_size and begin+span > 2**self._max_bin_power:
            error_string = "feature index at {}: must be less than 2^{}".format \
                                                (begin+span, self._max_bin_power)
            raise ValueError(error_string)

        #run the assignment loop
        while True:
            bin_level_count = self._bin_level_count
            max_bin_power = self._max_bin_power 
            for l_inverse in range(bin_level_count):
                L = bin_level_count - 1 - l_inverse
                # calculate offset (oL) of the list at level L
                oL = (2**(3*L) - 1)/7
                #calculate size (sL) of the list: the number of residues in width
                sL = float(2**(max_bin_power-3*L))
                # interval[0] >= (k - oL)*sL
                # rearrange to form
                # k =< oL + (interval[0])/sL
                k1 = int(floor(oL + (begin/sL)))
                # interval[1] < (k - oL + 1)*sL
                # rearrange to form
                #k > 1 + oL + (interval[1])/sL
                #k2 = oL - 1 + (begin+span)/sL  
                k2 = int(ceil(oL - 1 + (begin+span)/sL))
                if k1 == k2:
                    return k1
            # no suitable bin has been found, expand bins
            # this will not be invoked if the data satisfies
            # criteria checked prior to the while loop
            self._increase_bin_sizes()


           
if __name__ ==  "__main__":
    """ the following unit tests will eventually be used outside of this"""

    import unittest
    
    class TestFeatureBinCollection(unittest.TestCase):
        def setUp(self):
            self.bins = FeatureBinCollection()
            
        def test_binning_default_level(self):
            k_0_1 = self.bins._calculate_bin_index(0,1)
            k_0_256 = self.bins._calculate_bin_index(0,256)
            k_1_256 = self.bins._calculate_bin_index(1,256)
            k_256_256 = self.bins._calculate_bin_index(256,256)
            largeBin = 15000
            kLarge = self.bins._calculate_bin_index(largeBin*256,256)
            
            self.assertEqual(k_0_1, 4681) #this is the lowest lv5 offset
            self.assertEqual(k_1_256, 585) #lowest for lv4
            self.assertEqual(k_0_1, k_0_256)
            self.assertEqual(k_256_256-1, k_0_256)
            self.assertEqual(k_256_256-1, k_0_256)
            self.assertEqual(kLarge, 4681+15000)
            
        def test_binning_size_changer(self):
            #test largest of the default size
            k_full = self.bins._calculate_bin_index(0,8388608)
            self.assertEqual(k_full, 0)
            self.assertEqual(self.bins.maxSeqLength, 8388608) 
            self.assertEqual(self.bins._max_bin_power, 23)
            #trigger a size rearrangement 
            k_overFull = self.bins._calculate_bin_index(1,8388608)
            self.assertEqual(self.bins.maxSeqLength, 8388608*8)
            k_full = self.bins._calculate_bin_index(0,8388608)
            self.assertEqual(k_overFull, 0)
            self.assertEqual(k_full, 1)
            self.assertEqual(self.bins._max_bin_power, 26)
            #trigger a size rearrangement twice
            k_overFull = self.bins._calculate_bin_index(1,8388608*8)
            self.assertEqual(self.bins.maxSeqLength, 8388608*8*8)
            k_full = self.bins._calculate_bin_index(0,8388608*8)
            self.assertEqual(k_overFull, 0)
            self.assertEqual(k_full, 1)
            self.assertEqual(self.bins._max_bin_power, 29)
            
        def test_insertion_simple(self):
            testTuple1 = (0, 256)
            testTuple2 = (1, 257)
            testTuple3 = (256, 256+256)
            testTuple4 = (4194305, 5194305)
            self.bins.insert(testTuple1)
            self.bins.insert(testTuple2)
            self.bins.insert(testTuple3)
            self.bins.insert(testTuple4)
            self.assertIn(testTuple1, self.bins._bins[4681])
            self.assertIn(testTuple2, self.bins._bins[585])
            self.assertIn(testTuple3, self.bins._bins[4682])
            self.assertIn(testTuple4, self.bins._bins[5])
            
        def test_insertion_and_rearrangement(self):
            testTuple1 = (0, 256)
            testTuple2 = (1, 257)
            testTuple3 = (256, 256+256)
            testTuple4 = (4194305, 5194305)
            testTuple6 = (0, 8388608)
            testTuple0 = (20000, 50000)
            self.bins.insert(testTuple1)
            self.bins.insert(testTuple2)
            self.bins.insert(testTuple3)
            self.bins.insert(testTuple4)
            self.bins.insert(testTuple6)
            self.bins.insert(testTuple0)
            self.assertIn(testTuple1, self.bins._bins[4681])
            self.assertIn(testTuple2, self.bins._bins[585])
            self.assertIn(testTuple3, self.bins._bins[4682])
            self.assertIn(testTuple4, self.bins._bins[5])
            self.assertIn(testTuple6, self.bins._bins[0])
            #trigger a size rearrangement 
            testTuple5 = (0, 9000000)
            self.bins.insert(testTuple5)
            testTuple7 = (9000000, 9002500)
            testTuple8 = (0, 67108864)
            self.bins.insert(testTuple7)
            self.bins.insert(testTuple8)
            # recalculate bin positions (note that _calculate_bin_index is not
            # used during the rearrangement procedure) this step allows me to check
            # through two methods that the rearranged bins are still organized
            a = lambda tt: tt[0]
            b = lambda tt: tt[1] - tt[0]
            bin2 = self.bins._calculate_bin_index(a(testTuple2),b(testTuple2))
            bin3 = self.bins._calculate_bin_index(a(testTuple3),b(testTuple3))
            bin4 = self.bins._calculate_bin_index(a(testTuple4),b(testTuple4))
            bin6 = self.bins._calculate_bin_index(a(testTuple6),b(testTuple6))
            bin1 = self.bins._calculate_bin_index(a(testTuple1),b(testTuple1))
            bin7 = self.bins._calculate_bin_index(a(testTuple7),b(testTuple7))
            bin8 = self.bins._calculate_bin_index(a(testTuple8),b(testTuple8))
            bin5 = self.bins._calculate_bin_index(a(testTuple5),b(testTuple5))
            bin0 = self.bins._calculate_bin_index(a(testTuple0),b(testTuple0))
            self.assertIn(testTuple1, self.bins._bins[4681])
            self.assertIn(testTuple2, self.bins._bins[4681])
            self.assertIn(testTuple3, self.bins._bins[4681])
            self.assertIn(testTuple4, self.bins._bins[5+8])
            self.assertIn(testTuple5, self.bins._bins[0])
            self.assertIn(testTuple6, self.bins._bins[1])
            self.assertIn(testTuple1, self.bins._bins[bin1])
            self.assertIn(testTuple2, self.bins._bins[bin2])
            self.assertIn(testTuple3, self.bins._bins[bin3])
            self.assertIn(testTuple4, self.bins._bins[bin4])
            self.assertIn(testTuple5, self.bins._bins[bin5])
            self.assertIn(testTuple6, self.bins._bins[bin6])
            self.assertIn(testTuple7, self.bins._bins[bin7])
            self.assertIn(testTuple0, self.bins._bins[bin0])
            correctedsize_list = [2048, 16384, 131072, 1048576, 8388608, 67108864]
            self.assertEqual(self.bins._size_list, correctedsize_list)
            #border case; last small bin
            testTuple8_1 = (67108863, 67108864)
            self.bins.insert(testTuple8_1)
            bin8_1 = self.bins._calculate_bin_index(a(testTuple8_1),b(testTuple8_1))
            self.assertIn(testTuple8_1, self.bins._bins[bin8_1])
            #border case; one size up
            testTuple8_2 = (0, 67108865)
            self.assertNotIn(67108864*8, self.bins._size_list)
            self.bins.insert(testTuple8_2)
            self.assertIn(67108864*8, self.bins._size_list)
            bin8_1 = self.bins._calculate_bin_index(a(testTuple8_1),b(testTuple8_1))
            self.assertIn(testTuple8_1, self.bins._bins[bin8_1])

        def test_static_lists(self):
            #testing the maximum length of
            staticbins = FeatureBinCollection(length=67108864)
            a = lambda tt: tt[0]
            b = lambda tt: tt[1] - tt[0]
            testTuple0 = (0, 2048)
            testTuple1 = (1, 2049)
            testTuple2 = (20000, 50000)
            testTuple3 = (67108864-3, 67108864)
            testTuple4 = (4194305, 5194305)
            testTuple5 = (0, 8388608)
            testTuple6 = (0, 67108864)
            staticbins.insert(testTuple0)
            staticbins.insert(testTuple1)
            staticbins.insert(testTuple2)
            staticbins.insert(testTuple3)
            staticbins.insert(testTuple4)
            staticbins.insert(testTuple5)
            staticbins.insert(testTuple6)
            bin0 = staticbins._calculate_bin_index(a(testTuple0),b(testTuple0))
            bin1 = staticbins._calculate_bin_index(a(testTuple1),b(testTuple1))
            bin2 = staticbins._calculate_bin_index(a(testTuple2),b(testTuple2))
            bin3 = staticbins._calculate_bin_index(a(testTuple3),b(testTuple3))
            bin4 = staticbins._calculate_bin_index(a(testTuple4),b(testTuple4))
            bin5 = staticbins._calculate_bin_index(a(testTuple5),b(testTuple5))
            bin6 = staticbins._calculate_bin_index(a(testTuple6),b(testTuple6))
            self.assertIn(testTuple0, staticbins._bins[4681])
            self.assertIn(testTuple0, staticbins._bins[bin0])
            self.assertIn(testTuple1, staticbins._bins[585])
            self.assertIn(testTuple1, staticbins._bins[bin1])
            self.assertIn(testTuple2, staticbins._bins[bin2])
            self.assertIn(testTuple3, staticbins._bins[bin3])
            self.assertIn(testTuple4, staticbins._bins[bin4])
            self.assertIn(testTuple5, staticbins._bins[bin5])
            self.assertIn(testTuple6, staticbins._bins[bin6])
            #test that exceptions are raised if a over-sized bin is used
            overSizedTuple1 = (0, 67108865)
            overSizedTuple2 = (67108864, 67108865)
            reallyReallyoversizedEnd = 1+2**41
            self.assertRaises(ValueError, staticbins.insert, overSizedTuple1)
            self.assertRaises(ValueError, staticbins.insert, overSizedTuple2)
            self.assertRaises(ValueError, FeatureBinCollection, reallyReallyoversizedEnd)

        def test_getter_base_cases(self):
            testTuple0 = (0, 2048)
            testTuple1 = (1, 2049)
            testTuple2 = (0, 2**41)
            self.bins.insert(testTuple0)
            self.bins.insert(testTuple1)
            self.assertRaises(TypeError, self.bins.__getitem__, "hello")
            self.assertRaises(TypeError, self.bins.__getitem__, 5.45)
            self.assertRaises(KeyError, self.bins.__getitem__, slice(0,25,2))
            self.assertRaises(IndexError, self.bins.__getitem__, 2**27)
            self.assertRaises(IndexError, self.bins.__getitem__, slice(55,25))
            #Resize the array to test large numbers
            #  this is required to check that long type integers (python2) do not interfere
            #  with the getter routine, also that py3 is not failing when testing integers
            self.bins.insert(testTuple2)
            resultsForTuple2 = self.bins[2**40]
            self.assertIn(testTuple2, resultsForTuple2)


        def test_getter_functionality(self):
            testTuple0 = (0, 2048)
            testTuple1 = (1, 2049)
            testTuple2 = (20000, 50000)
            testTuple3 = (0, 8388608)
            testTuple4 = (0, 67108864)
            self.bins.insert(testTuple0)
            self.bins.insert(testTuple1)
            self.bins.insert(testTuple2)
            self.bins.insert(testTuple3)
            #   test length=1 retrievals 
            #int key
            self.assertNotIn(testTuple0, self.bins[2049])
            self.assertNotIn(testTuple0, self.bins[2048])
            self.assertNotIn(testTuple1, self.bins[0])
            self.assertNotIn(testTuple2, self.bins[50000])
            self.assertNotIn(testTuple2, self.bins[19999])
            #slice key
            self.assertNotIn(testTuple0, self.bins[2049:2050])
            self.assertNotIn(testTuple0, self.bins[2048:2049])
            self.assertNotIn(testTuple1, self.bins[0:1])
            self.assertNotIn(testTuple2, self.bins[50000:50001])
            self.assertNotIn(testTuple2, self.bins[19999:20000])
            #fully contained
            self.assertIn(testTuple3, self.bins[1000:1500])
            self.assertNotIn(testTuple2, self.bins[1000:1500])
            self.assertIn(testTuple1, self.bins[1000:1500])
            #feature overlap of sequence right side
            self.assertIn(testTuple2, self.bins[19000:21000])
            self.assertIn(testTuple2, self.bins[19000:20001])
            self.assertNotIn(testTuple2, self.bins[19000:20000])
            #feature overlap of sequence left side
            self.assertIn(testTuple2, self.bins[29000:55000])
            self.assertIn(testTuple2, self.bins[49999:55000])
            self.assertNotIn(testTuple2, self.bins[50000:55000])
            #
            # trigger a dynamic expansion event and re-run tests
            #
            self.bins.insert(testTuple4)
            self.assertNotIn(testTuple0, self.bins[2049])
            self.assertNotIn(testTuple0, self.bins[2048])
            self.assertNotIn(testTuple1, self.bins[0])
            self.assertNotIn(testTuple2, self.bins[50000])
            self.assertNotIn(testTuple2, self.bins[19999])
            #slice key
            self.assertNotIn(testTuple0, self.bins[2049:2050])
            self.assertNotIn(testTuple0, self.bins[2048:2049])
            self.assertNotIn(testTuple1, self.bins[0:1])
            self.assertNotIn(testTuple2, self.bins[50000:50001])
            self.assertNotIn(testTuple2, self.bins[19999:20000])
            #fully contained
            self.assertIn(testTuple3, self.bins[1000:1500])
            self.assertNotIn(testTuple2, self.bins[1000:1500])
            self.assertIn(testTuple1, self.bins[1000:1500])
            #feature overlap of sequence right side
            self.assertIn(testTuple2, self.bins[19000:21000])
            self.assertIn(testTuple2, self.bins[19000:20001])
            self.assertNotIn(testTuple2, self.bins[19000:20000])
            #feature overlap of sequence left side
            self.assertIn(testTuple2, self.bins[29000:55000])
            self.assertIn(testTuple2, self.bins[49999:55000])
            self.assertNotIn(testTuple2, self.bins[50000:55000])


    unittest.main( exit=False )

    import doctest
    doctest.testmod()
    