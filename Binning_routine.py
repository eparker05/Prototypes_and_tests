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

class StupidFeatureBinCollection(object):
    """this class manages a flat list of features and retrieves them

    Interaction with this class should be the same as FeatureBinCollection
    but this class lacks the binning strategy that improves performance.
    
    no stability checks and base cases are managed, this is unstable
    and is only useful for performance comparison."""
    
    def __init__(self, length = None, beginindex=0, endindex=1):
        self._bins = []
        self._beginindex = beginindex 
        self._endindex = endindex
        self._is_sorted = False
        
    def insert(self, feature_tuple):
        beginindex = self._beginindex
        endindex = self._endindex

        begin = feature_tuple[beginindex]
        end = feature_tuple[endindex]
        assert isinstance(begin, integer_types)
        assert isinstance(end, integer_types)
        assert begin <= end
        span = end-begin
        
        self._is_sorted = False
        self._bins.append(feature_tuple)
    
    def sort(self):
        self._bins.sort()
        self._is_sorted = True 
        
    def __getitem__(self, key):

        if not self._is_sorted:
            self.sort()
        
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
        
        #fix begin or end index for slicing of forms: bins[50:] or bins[:50] or even bins[:]
        if key.stop is None:
            key = slice(key.start, self._max_sequence_length)
        if key.start is None:
            key = slice(0, key.stop)
        
        #check that the key is within boundaries:
        if key.start < 0:
            raise IndexError("key out of bounds")
        if key.start > key.stop:
            raise IndexError("key not valid, slice.start > slice.stop")

        #check for bound and overlapping sequences
        possible_entries = self._bins
        return_entries = []
        for feature in possible_entries:
            #this covers fully bound sequence and left overlap   ssssssFsFsFsFFFFF
            if key.start <= feature[beginindex] < key.stop:
                return_entries.append(feature)
            #this covers left sequence right sequence overlap of (F)  FFFFFFFsFsFsFsssss
            elif key.start < feature[endindex] <= key.stop:
                return_entries.append(feature)
            #this covers seqyebces fully bound by a feature      FFFFFFFsFsFsFsFsFFFFFFF
            elif key.start >= feature[beginindex] and key.stop <= feature[endindex]:
                return_entries.append(feature)
            #ends the iteration once no more values are possible
            if key.stop < feature[beginindex]:
                break
        return return_entries

        
class FeatureBinCollection(object):
    """this class manages the creation and maintenance of feature indices

       This class is used to organize feature data in a quickly retrievable
       data structure. The feature data must be added as a tuple containing
       at least two indices: first annotated residue and the last as a half
       open half closed interval [first, last). The indices are assumed to be
       the first two elements of the stored tuple, but they may be re-assigned
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
       only likely candidate features to be queried rather than all features. The 
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
       matter of sanity checking, bin sizes are capped at 2.2 trillion residues (2^41).

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
        self._sorted = False
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
        self._max_sequence_length = self._size_list[-1]
    
    def insert(self, feature_tuple):
        """inserts a tuple with a sequence range into the feature bins
        
        data is assumed to be somewhat scrubbed, coming from a parser
        or a parser consumer."""
        
        beginindex = self._beginindex
        endindex = self._endindex

        #reset sorted quality
        self._sorted = False

        begin = feature_tuple[beginindex]
        end = feature_tuple[endindex]
        assert isinstance(begin, integer_types)
        assert isinstance(end, integer_types)
        assert begin <= end
        span = end-begin
        
        bin_index = self._calculate_bin_index(begin, span)
        self._bins[bin_index].append(feature_tuple)
        
    def __len__(self):
        return sum(len(bin) for bin in self._bins) 

    def sort(self):
        """this performs bin-centric sorting, necessary for faster retrieval"""
        #bins must be sorted by the begin index, this is fastest
        if self._beginindex == 0:
            for i in range(len(self._bins)):
                self._bins[i].sort()
        #this is a bit slower but accomodates diverse data structures
        else:
            for i in range(len(self._bins)):
                self.bins[i].sort(key = lambda tup: tup[self._beginindex])
        #reset sorted quality
        self._sorted = True
            

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
        
        #fix begin or end index for slicing of forms: bins[50:] or bins[:50] or even bins[:]
        if key.stop is None:
            key = slice(key.start, self._max_sequence_length)
        if key.start is None:
            key = slice(0, key.stop)
        
        #check that the key is within boundaries:
        if key.start < 0:
            raise IndexError("key out of bounds")
        if key.start > key.stop:
            raise IndexError("key not valid, slice.start > slice.stop")
        if key.stop >= self._max_sequence_length:
            key = slice(key.start, self._max_sequence_length-1)

        #pre-sort if necessary
        if not self._sorted:
            self.sort()
            
        #code taken from self._calculate_bin_index(), comments removed
        return_entries = []
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
            if k2-k1 > 2:
                for bin in range(k1+1, k2-1):
                    return_entries.extend( self._bins[bin]) 
            for binn in set([k1,k2-1]):
                #for binn in range(k1,k2):
                for feature in self._bins[binn]:
                    #this covers fully bound sequence and left overlap   ssssssFsFsFsFFFFF
                    if key.start <= feature[beginindex] < key.stop:
                        return_entries.append(feature)
                    #this covers left sequence right sequence overlap of (F)  FFFFFFFsFsFsFsssss
                    elif key.start < feature[endindex] <= key.stop:
                        return_entries.append(feature)
                    #this covers seqyebces fully bound by a feature      FFFFFFFsFsFsFsFsFFFFFFF
                    elif key.start >= feature[beginindex] and key.stop <= feature[endindex]:
                        return_entries.append(feature)
                    if key.stop < feature[beginindex]:
                        break
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
        
        #take care base cases with the span parameter
        assert span >= 0   # no function may pass a span < 0
        #this is required for zero length determination of bin location
        span = max(1, span)
        
        # all indices must be positive
        assert begin >= 0  
        if self._dynamic_size:
            assert begin+span <= 2**41   # len(seq) > 2.19 trillion is not reasonable
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
                group_length = (2**(3*(L+1)) - 1)/7
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
                if k1 == k2 and k1 < group_length:
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
            
        def test_bin_index_finder_smallest_bin_and_leftmost(self):
            k_0_256 = self.bins._calculate_bin_index(0,256)
            self.assertEqual(k_0_256, 4681) #this is the leftmost level 5 bin
        
        def test_bin_index_finder_smallest_bin_and_2ndleftmost(self):
            k_256_256 = self.bins._calculate_bin_index(256,256)
            self.assertEqual(k_256_256, 4682) #this is the leftmost level 5 bin
            
        def test_bin_index_finder_smallest_bin_and_2ndleftmost_zero_length(self):
            k_256_0 = self.bins._calculate_bin_index(256,0)
            self.assertEqual(k_256_0, 4682) 
        
        def test_bin_index_finder_smallest_bin_last_index(self):
            k_1_256 = self.bins._calculate_bin_index(1,256)
            self.assertEqual(k_1_256, 585) 
         
        def test_bin_index_finder_smallest_bin_last_index(self):
            k_big_256 = self.bins._calculate_bin_index(8388352,256)
            self.assertEqual(k_big_256, 37448) # this is the rightmost lv5 bin
        
        def test_bin_index_finder_largest_bin_and_ensure_no_recalculation(self):
            k_small_big = self.bins._calculate_bin_index(0, 8388608)
            self.assertEqual(k_small_big, 0)
            #value below should not have changed
            k_big_256 = self.bins._calculate_bin_index(8388352,256)
            self.assertEqual(k_big_256, 37448)
          
        def test_binning_size_changer_once(self):
            #test largest of the default size
            k_full = self.bins._calculate_bin_index(0,8388608)
            self.assertEqual(k_full, 0)
            #trigger a size rearrangement 
            k_overFull = self.bins._calculate_bin_index(1,8388608)
            k_full = self.bins._calculate_bin_index(0,8388608)
            self.assertEqual(self.bins._max_sequence_length, 8388608*8)
            self.assertNotIn(256, self.bins._size_list)
            self.assertIn(2048, self.bins._size_list)
            self.assertEqual(k_overFull, 0)
            self.assertEqual(k_full, 1)
            self.assertEqual(self.bins._max_bin_power, 26)
        
        def test_binning_size_changer_multiple_steps(self):
            #this sets some inital values up 
            self.test_binning_size_changer_once()
            #trigger a size rearrangement the second time
            k_overFull = self.bins._calculate_bin_index(1,8388608*8)
            self.assertEqual(self.bins._max_sequence_length, 8388608*8*8)
            k_full = self.bins._calculate_bin_index(0,8388608*8)
            self.assertEqual(k_overFull, 0)
            self.assertEqual(k_full, 1)
            self.assertEqual(self.bins._max_bin_power, 29)
            self.assertNotIn(2048, self.bins._size_list)
        
        def test_insertion_bad_negative_val(self):
            testTuple1 = (-1, 56)
            self.assertRaises(AssertionError,self.bins.insert,testTuple1)
            
        def test_insertion_once_smallest(self):
            testTuple1 = (0, 256)
            self.bins.insert(testTuple1)
            self.assertIn(testTuple1, self.bins._bins[4681])
            
        def test_insertion_once_smallest_but_overlaps(self):
            testTuple2 = (1, 257)
            self.bins.insert(testTuple2)
            self.assertIn(testTuple2, self.bins._bins[585])
            
        def test_insertion_once_smallestbin_rightmost(self):
            testTuple3 = (8388608-256, 8388608)
            self.bins.insert(testTuple3)
            self.assertIn(testTuple3, self.bins._bins[37448])

        def test_insertion_zero_length_smallestbin_rightmost(self):
            testTuple3 = (8388607, 8388607)
            self.bins.insert(testTuple3)
            self.assertIn(testTuple3, self.bins._bins[37448])
            
        
        def test_insertion_once_smallestbin_rightmost_lv4(self):
            testTuple4 = (8388608-257, 8388608)
            self.bins.insert(testTuple4)
            self.assertIn(testTuple4, self.bins._bins[4680])    

        def test_insertion_with_rearrangement(self):
            testTuple3 = (256, 256+256)
            self.bins.insert(testTuple3)
            self.assertIn(testTuple3, self.bins._bins[4682])
            #trigger a size rearrangement with an insertion
            testTuple5 = (0, 9000000)
            self.bins.insert(testTuple5)
            self.assertIn(testTuple3, self.bins._bins[4681])
            self.assertIn(testTuple5, self.bins._bins[0])
            self.assertEqual(self.bins._max_sequence_length, 8388608*8)
            self.assertNotIn(256, self.bins._size_list)
            self.assertIn(2048, self.bins._size_list)
            self.assertEqual(self.bins._max_bin_power, 26)
        
        def test_insertion_level3_and_level5_then_resize_to_fit_in_same_bin(self):
            test_tuple_lv3 = (16384 , 16384+16384)
            test_tuple_lv5 = (16384+16384-256 , 16384+16384)
            self.bins.insert(test_tuple_lv3)
            self.bins.insert(test_tuple_lv5)
            self.assertIn(test_tuple_lv3, self.bins._bins[74])
            self.assertIn(test_tuple_lv5, self.bins._bins[4681+127])
            #trigger rearrangement by 2 levels
            k_overFull = self.bins._calculate_bin_index(1,8388608*8)
            self.assertIn(test_tuple_lv3, self.bins._bins[4682])
            self.assertIn(test_tuple_lv5, self.bins._bins[4682])
            self.assertEqual(self.bins._max_bin_power, 29)
            
            
        def test_overflows_of_static_defined_lists(self):
            staticbins = FeatureBinCollection(length=67108864)
            #insert a chunk of data
            testTuple1 = (1, 2049)
            staticbins.insert(testTuple1)
            self.assertIn(testTuple1, staticbins._bins[585])
            #test that exceptions are raised if a over-sized bin is used
            overSizedTuple1 = (0, 67108865)
            overSizedTuple2 = (67108864, 67108865)
            self.assertRaises(ValueError, staticbins.insert, overSizedTuple1)
            self.assertRaises(ValueError, staticbins.insert, overSizedTuple2)
            
        def test_oversized_definition_of_the_collection(self):    
            reallyReallyoversizedEnd = 1+2**41
            self.assertRaises(ValueError, FeatureBinCollection, reallyReallyoversizedEnd)

        def test_getter_get_values_from_empty_set(self):
            resultsEmpty = self.bins[2**20]
            self.assertEqual([], resultsEmpty)
            
        def test_getter_get_values_edge_cases(self):
            resultsEdgeRight = self.bins[-1+2**23]
            resultsEdgeLeft = self.bins[0]
            self.assertEqual([], resultsEdgeRight)
            self.assertEqual([], resultsEdgeLeft)
            
        def test_getter_get_values_from_out_of_bounds(self):
            self.assertRaises(IndexError, self.bins.__getitem__, -1)
            self.assertRaises(IndexError, self.bins.__getitem__, 1+2**23)
            
        def test_getter_typeError_string(self):
            self.assertRaises(TypeError, self.bins.__getitem__, "hello")
            
        def test_getter_typeError_float(self):
            self.assertRaises(TypeError, self.bins.__getitem__, 5.45)          
            
        def test_getter_typeError_stepped_slice(self):
            self.assertRaises(KeyError, self.bins.__getitem__, slice(0,25,2))
            
        def test_getter_half_slices(self):
            emptyL = self.bins[:50]
            emptyR = self.bins[50:]
            emptyM = self.bins[:]
            self.assertEqual([], emptyL)
            self.assertEqual([], emptyR)
            self.assertEqual([], emptyM)

        def test_getter_half_slices_with_result_right_border(self):
            resultR = (20000,30000)
            self.bins.insert(resultR)
            emptyL = self.bins[:20000]
            emptyR = self.bins[20000:]
            emptyM = self.bins[:]
            self.assertEqual([], emptyL)
            self.assertIn(resultR, emptyR)
            self.assertIn(resultR, emptyM)
        
        def test_getter_half_slices_with_result_left_border(self):
            resultL = (10000,20000)
            self.bins.insert(resultL)
            emptyL = self.bins[:20000]
            emptyR = self.bins[20000:]
            emptyM = self.bins[:]
            self.assertEqual([], emptyR)
            self.assertIn(resultL, emptyL)
            self.assertIn(resultL, emptyM)
            
        def test_insertion_where_medium_sized_bin_is_out_of_bounds(self):
            resultL = (8388608, 8388608+2047)
            self.bins.insert(resultL)
            self.assertIn(resultL, self.bins[838840:8388610])
            self.assertEqual(self.bins._max_bin_power, 26)

        def test_getter_overlap_left(self):
            feature = (100000,200000)
            self.bins.insert(feature)
            result = self.bins[99000:101000]
            self.assertIn(feature, result)
        
        def test_getter_overlap_right(self):
            feature = (100000,200000)
            self.bins.insert(feature)
            result = self.bins[199000:201000]
            self.assertIn(feature, result)

        def test_getter_feature_inside_region(self):
            feature = (100000,200000)
            self.bins.insert(feature)
            result = self.bins[99000:201000]
            self.assertIn(feature, result) 
            
        def test_getter_region_inside_feature(self):
            feature = (100000,200000)
            self.bins.insert(feature)
            result = self.bins[101000:199000]
            self.assertIn(feature, result) 
        
        def test_getter_zero_length_inside(self):
            testTuple3 = (8388605, 8388605)
            self.bins.insert(testTuple3)
            self.assertIn(testTuple3,self.bins[8388604:8388606])
            self.assertIn(testTuple3,self.bins[8388604:8388605])
            self.assertIn(testTuple3,self.bins[8388605:8388606])
            self.assertEqual([],self.bins[8388603:8388604])
            
    unittest.main( exit=False )

    print("now running doctests")
    
    import doctest
    if doctest.testmod():
        print("DOCTESTS: work")
    
    #
    #
    # the following was hacked togeather to be used for basic performance testing
    # 
    #
    #
    sbins = StupidFeatureBinCollection()
    bins = FeatureBinCollection()
    
    import random
    import time
    
    #measuring insertion time
    numberToInsert = 500000
    count = numberToInsert
    tinitialize = time.clock()
    sbins_insertion_time = 0
    bins_insertion_time = 0
    size = 8*67108864
    while count > 0:
        count -= 1
        #models a distribution where 60% are small features (one or several AAs)
        # 20% are medium features
        # 20% are large features
        if count > 0.8*numberToInsert:
            num1 = random.randint(0,size)
            num2 = random.randint(0,size)
        elif count > 0.6*numberToInsert:
            num1 = random.randint(0,size/16)
            num2 = random.randint(0,size/16)
            randomAddition = random.randint(0,15*size/16)
            num1 += randomAddition
            num2 += randomAddition
        else:
            num1 = random.randint(0,size/128)
            num2 = random.randint(0,size/128)
            randomAddition = random.randint(0,127*size/128)
            num1 += randomAddition
            num2 += randomAddition
        idx1 = min(num1, num2)
        idx2 = max(num1, num2)
        newFeature = (idx1, idx2,)
        sbin0 = time.clock()
        sbins.insert(newFeature)
        sbins_insertion_time += time.clock()- sbin0
        bin0 = time.clock()
        bins.insert(newFeature)
        bins_insertion_time += time.clock()- bin0
    #the sbins sort is actually a part of the insertion process since 
    #data retrieval requires this to be sorted
    bin0 = time.clock()
    sbins.sort()    
    sbins_insertion_time += time.clock()- bin0
    bin0 = time.clock()
    bins.sort()    
    bins_insertion_time += time.clock()- bin0
    print("insertiontime bins {}, stupidbins {}".format(bins_insertion_time,sbins_insertion_time))
    
    #measuring retrieval time
    number_to_retrieve = 100
    count = number_to_retrieve
    tinitialize = time.clock()
    sbins_insertion_time = 0
    bins_insertion_time = 0
    while count > 0:
        count -= 1
        if count > 0.5* number_to_retrieve:
            num1 = random.randint(0,size)
            num2 = random.randint(0,size)
        else:
            num1 = random.randint(0,size/128)
            num2 = random.randint(0,size/128)
            randomAddition = random.randint(0,127*size/128)
            num1 += randomAddition
            num2 += randomAddition
        idx1 = min(num1, num2)
        idx2 = max(num1, num2)
        sbin0 = time.clock()
        sbinsResult = sbins[idx1:idx2]
        sbins_insertion_time += time.clock()- sbin0
        bin0 = time.clock()
        binsResult = bins[idx1:idx2]
        bins_insertion_time += time.clock()- bin0
        binsResult.sort()
        sbinsResult.sort()
        assert binsResult == sbinsResult
        
    sbins.sort()    
    print("retrieval time: bins {}, stupidbins {}".format(bins_insertion_time,sbins_insertion_time))
    