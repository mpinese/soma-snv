#!/usr/bin/env python3
import collections
import intervaltree
import logging
import math
import operator
import re
import time


_COMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
_MAXIT = 100
_FPMIN = 1e-30
_EPS = 3e-7
_MpileupField = collections.namedtuple('MpileupField', ['chrom', 'pos', 'ref', 'alt_depth', 'depth', 'alt_counts'])


def _complement(seq):
    return ''.join([_COMPLEMENTS[base] for base in seq])


def _normalise_triplet(triplet):
    # Normalises the triplet of bases in triplet (supplied as a 3-member list)
    # by reverse complementing if necessary, to ensure that the middle element 
    # is in {C, T, N}.  Returns the normalised triplet as well as a flag to 
    # indicate whether normalisation took place.
    if triplet[1] == 'C' or triplet[1] == 'T' or triplet[1] == 'N':
        return triplet, False
    return _complement(triplet)[::-1], True


def _tally_sequence_string(seq_string, bq_string, mq_string, min_bq, min_mq):
    # Format of seq_string, from the samtools mpileup docs:
    # "a dot stands for a match to the reference base on the forward strand, a comma for 
    # a match on the reverse strand, a '>' or '<' for a reference skip, `ACGTN' for a 
    # mismatch on the forward strand and `acgtn' for a mismatch on the reverse strand. 
    # A pattern `\\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this 
    # reference position and the next reference position. The length of the insertion 
    # is given by the integer in the pattern, followed by the inserted sequence. 
    # Similarly, a pattern `-[0-9]+[ACGTNacgtn]+' represents a deletion from the 
    # reference. The deleted bases will be presented as `*' in the following lines. 
    # Also at the read base column, a symbol `^' marks the start of a read. The ASCII 
    # of the character following `^' minus 33 gives the mapping quality. A symbol `$' 
    # marks the end of a read segment."
    string_uppercase = seq_string.strip().upper()
    i = 0
    j = 0
    ref_count = 0
    nonref_counter = collections.Counter()

    while i < len(string_uppercase):
        if string_uppercase[i] == 'A' or string_uppercase[i] == 'C' or string_uppercase[i] == 'G' or string_uppercase[i] == 'T':
            if ord(bq_string[j]) - 33 >= min_bq and ord(mq_string[j]) - 33 >= min_mq:
                nonref_counter.update(string_uppercase[i])
            j += 1
        elif string_uppercase[i] == ',' or string_uppercase[i] == '.':
            if ord(bq_string[j]) - 33 >= min_bq and ord(mq_string[j]) - 33 >= min_mq:
                ref_count += 1
            j += 1
        elif string_uppercase[i] == '$': # End of a read
            pass
        elif string_uppercase[i] == 'N': # Unknown base
            j += 1
            pass
        elif string_uppercase[i] == '*': # Placeholder for deletion described at a previous locus
            pass
        elif string_uppercase[i] == '^': # Start of a read
            i += 1                       # Skip the subsequent mapping quality character
        elif string_uppercase[i] == '+' or string_uppercase[i] == '-':  # Indel
            indel_length_field = re.match('[0-9]+', string_uppercase[i+1:]).group()
            indel_length = int(indel_length_field)
            i += len(indel_length_field) + indel_length
        else:
            warnings.warn('Unexpected symbol {} in read string: {}'.format(read_string[i], read_string))
        i += 1
    return ref_count, sum(nonref_counter.values()), nonref_counter


def _parse_mpileup_line(line, min_bq, min_mq):
    # samtools mpileup -C 50 -a -E -I -q 30 -Q 30 -s -f <ref.fa> <data.bam/cram>
    chrom, pos_string, ref, _, seq_string, bq_string, mq_string = line.rstrip().split('\t')

    pos = int(pos_string)
    ref_depth, alt_depth, alt_counts = _tally_sequence_string(seq_string, bq_string, mq_string, min_bq, min_mq)
    depth = ref_depth + alt_depth

    return _MpileupField(chrom, pos, ref, alt_depth, depth, alt_counts)


def _most_common_alt(alt_counts):
    """
    Returns the most common alt allele in the collections.Counter object
    alt_counts.  In the event of a tie for most common, breaks ties by 
    base alphabetical order.

    >>> _most_common_alt(collections.Counter({'A': 10}))
    ('A', 10)
    >>> _most_common_alt(collections.Counter({'A': 10, 'C': 5}))
    ('A', 10)
    >>> _most_common_alt(collections.Counter({'A': 10, 'C': 20}))
    ('C', 20)
    >>> _most_common_alt(collections.Counter({'A': 10, 'C': 10}))
    ('A', 10)
    >>> _most_common_alt(collections.Counter({'T': 10, 'C': 10}))
    ('C', 10)
    >>> _most_common_alt(collections.Counter({'T': 10, 'C': 10, 'A': 20}))
    ('A', 20)
    """
    most_common_pairs = alt_counts.most_common()
    most_common_count = most_common_pairs[0][1]
    most_common_alts = [base for base, count in most_common_pairs if count == most_common_count]
    return (sorted(most_common_alts)[0], most_common_count)


def _betacf(a, b, x):
    """
    >>> _betacf(1, 2, 0.1)
    1.1728395061728394
    >>> _betacf(0, 2, 0.1)
    1.2345679012345678
    >>> _betacf(3, 5, 0.9)
    3918.5714285777376
    """
    # Accessory function for incomplete Beta function from NRC 6.4
    qab=a+b
    qap=a+1.0
    qam=a-1.0
    c=1.0
    d=1.0-qab*x/qap
    if abs(d) < _FPMIN:
        d=_FPMIN
    d=1.0/d
    h=d
    for m in range(1, _MAXIT+1):
        m2=2.0*m
        aa=m*(b-m)*x/((qam+m2)*(a+m2))
        d=1.0+aa*d
        if math.fabs(d) < _FPMIN:
            d=_FPMIN
        c=1.0+aa/c
        if math.fabs(c) < _FPMIN: 
            c=_FPMIN
        d=1.0/d
        h = h*d*c
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
        d=1.0+aa*d
        if math.fabs(d) < _FPMIN: 
            d=_FPMIN
        c=1.0+aa/c
        if math.fabs(c) < _FPMIN:
            c=_FPMIN
        d=1.0/d
        d_el=d*c
        h = h*d_el
        if math.fabs(d_el-1.0) < _EPS:
            break

    if m > _MAXIT:
        raise ValueError("a or b too big, or MAXIT too small in betacf")
    
    return h


def _betai(a, b, x):
    # Incomplete Beta function from NRC 6.4
    assert x >= 0.0 and x <= 1.0

    if x == 0.0 or x == 1.0:
        bt = 0.0
    else:
        bt = math.exp(math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)+a*math.log(x)+b*math.log(1.0-x))

    if x < (a+1.0)/(a+b+2.0):
        return bt*_betacf(a,b,x)/a
    else:
        return 1.0-bt*_betacf(b,a,1.0-x)/b


def _pbinom(q, size, prob, lower_tail=True):
    """
    >>> _pbinom(2, 4, 0.5)
    0.6875000000000001
    >>> _pbinom(2, 10, 0.5)
    0.05468749999999989
    >>> _pbinom(0, 10, 0.9)
    1.000000082740371e-10
    >>> _pbinom(10, 10, 0.9)
    1.0
    >>> _pbinom(0, 10, 0.1)
    0.3486784400015138
    >>> _pbinom(3, 2, 0.5)
    Traceback (most recent call last):
      File "/apps/python3/3.4.3/lib/python3.4/doctest.py", line 1318, in __run
        compileflags, 1), test.globs)
      File "<doctest __main__._pbinom[5]>", line 1, in <module>
        _pbinom(3, 2, 0.5)
      File "soma-snv.py", line 197, in _pbinom
        assert q >= 0 and q <= size
    AssertionError
    """
    assert prob >= 0.0 and prob <= 1.0
    assert size > 0
    assert q >= 0 and q <= size

    if q == size:
        upper_p = 0.0
    else:
        upper_p = _betai(q+1, size-q, prob)

    if lower_tail:
        return 1.0 - upper_p
    else:
        return upper_p



class Blacklist:
    def __init__(self):
        self.trees = {}


    def load_bed(self, instream):
        self.trees = {}
        for line in instream:
            chrom, start0, end1 = line.split('\t', 3)[:3]
            if chrom not in self.trees:
                self.trees[chrom] = intervaltree.IntervalTree()
            self.trees[chrom].add(intervaltree.Interval(int(start0), int(end1)))


    def intersects(self, chrom, pos):
        if chrom not in self.trees:
            return False

        return len(self.trees[chrom][pos]) > 0



class LssDetector:
    TestResult = collections.namedtuple('TestResult', ['was_powered', 'was_detection', 'critval_E', 'critval_H', 'sensitivity', 'snr'])


    def __init__(self, blacklist = Blacklist(), rate_het=1.0e-3, rate_sc=5.0e-7, vaf=0.1, min_snr=10.0, epsilon=1.0e-2, min_bq=30, min_mq=30, max_dp=100):
        self.blacklist = blacklist
        self.rate_het = rate_het
        self.rate_sc = rate_sc
        self.epsilon = epsilon
        self.vaf = vaf
        self.min_snr = min_snr
        self.min_mq = min_mq
        self.min_bq = min_bq
        self.max_dp = max_dp
        threshold_cache_array_size = 101
        if max_dp + 1 < threshold_cache_array_size:
            threshold_cache_array_size = max_dp + 1
        self.threshold_cache_dense = [None] * threshold_cache_array_size
        self.threshold_cache_sparse = {}
#        self._precompute_dense_thresholds()


    def _precompute_dense_thresholds(self):
        for depth in range(len(self.threshold_cache_dense)):
            self.threshold_cache_dense[depth] = self._calc_optimal_thresholds(depth)


    def _get_optimal_thresholds(self, depth):
        if depth < len(self.threshold_cache_dense):
            if self.threshold_cache_dense[depth] == None:
                self.threshold_cache_dense[depth] = self._calc_optimal_thresholds(depth)
            return self.threshold_cache_dense[depth]

        if depth not in self.threshold_cache_sparse:
            self.threshold_cache_sparse[depth] = self._calc_optimal_thresholds(depth)
        return self.threshold_cache_sparse[depth]


    @staticmethod
    def _calc_pn_pr(ce, ch, depth, epsilon, rc, rh, rr, pa):
        ae = _pbinom(ce, depth, epsilon, lower_tail=False)
        ah = _pbinom(ch + 1, depth, 0.5+epsilon/3.0, lower_tail=True)
        pn = rr*(_pbinom(ch, depth, pa) - _pbinom(ce - 1, depth, pa))
        pr = (rc/(1-rc)) * (pn / (rr*ae + rh*ah))
        return (pn, pr)


    def _calc_optimal_thresholds(self, depth):
        """
        Find ce and ch that yield optimal pn, subject to:
          1 <= ce <= ch <= n-1, and
          pr >= gr.
        Use an exhaustive search over the grid of valid ce, ch pairs.

        >>> cl = LssDetector(rate_het=1.0e-3, rate_sc=5.0e-7, vaf=0.1, min_snr=10.0, epsilon=1.0e-2)
        >>> cl._calc_optimal_thresholds(40)
        (None, None, 0, None)
        >>> cl._calc_optimal_thresholds(100)
        (10, 28, 0.6601016080599246, 18.0580053085918)
        >>> cl._calc_optimal_thresholds(135)
        (11, 43, 0.8803333325220946, 11.584716396883863)
        """
        rr = 1.0-self.rate_het-self.rate_het*self.rate_het
        pa = (1-4/3*self.epsilon)*self.vaf + self.epsilon
        rc = self.rate_sc
        rh = self.rate_het

        # 1. Identify a starting point for the optimization that is 
        #    within the feasible region defined by pr >= self.min_snr.
        #    We search along the diagonal ce = ch.
        for ch in range(1, depth):
            pn, pr = self._calc_pn_pr(ch, ch, depth, self.epsilon, rc, rh, rr, pa)
            if pr >= self.min_snr:
                break
        else:
            # No feasible region
            return (None, None, 0, None)

        # 2. Starting from the feasible point above, optimise pn by
        #    coordinate ascent.
        best_ce, best_ch, best_pr, best_pn = ch, ch, pr, pn

        while True:
            for ce, ch in ((best_ce-1,best_ch), (best_ce,best_ch+1), (best_ce-1,best_ch+1)):
                if ce <= 0 or ch <= 0 or ce >= depth or ch >= depth:
                    continue
                pn, pr = self._calc_pn_pr(ce, ch, depth, self.epsilon, rc, rh, rr, pa)
                if pr >= self.min_snr and pn > best_pn:
                    best_ce, best_ch, best_pr, best_pn = ce, ch, pr, pn
                    break
            else:
                break

        return (best_ce, best_ch, best_pn, best_pr)


    def _test_for_lsv(self, depth_alt, depth):
        if depth == 0:
            return self.TestResult(False, None, None, None, None, None)
        cE, cH, sens, snr = self._get_optimal_thresholds(depth)
        if cE == None:
            return self.TestResult(False, None, None, None, None, None)     # Cannot achieve target SNR; don't proceed with test.

        # This locus passes test S1 if depth_alt >= cE, and test S2 if depth_alt <= cH.
        # See LSS-5 for test definitions.
        return self.TestResult(True, depth_alt >= cE and depth_alt <= cH, cE, cH, sens, snr)


    @staticmethod
    def _emit_detection_line(fields, result, ref_triplet, outstream):
        # Using multiple calls to outstream.write because I believe that's fastest,
        # but I've not tested this.  Alternatives would be %/format/f-strings, or
        # simple concatenation.  TODO.

        # Preparations: what was the alt allele, is normalisation needed?
        most_common_alt, most_common_alt_dp = _most_common_alt(fields.alt_counts)
        norm_ref, norm_applied = _normalise_triplet(ref_triplet)
        if norm_applied:
            most_common_alt = _complement(most_common_alt)

        # Format: CHROM     POS     REF     ALT     Strand  DP      RD      AD      sensitivity     snr
        # tab-separated
        outstream.write(fields.chrom)
        outstream.write('\t')
        outstream.write(str(fields.pos))
        outstream.write('\t')
        outstream.write(norm_ref[0])
        outstream.write(norm_ref[1])
        outstream.write(norm_ref[2])
        outstream.write('\t')
        outstream.write(most_common_alt)
        outstream.write('\t')
        outstream.write('-1' if norm_applied else '+1')
        outstream.write('\t')
        outstream.write(str(fields.depth))
        outstream.write('\t')
        outstream.write(str(fields.depth - fields.alt_depth))
        outstream.write('\t')
        outstream.write(str(most_common_alt_dp))
        outstream.write('\t')
        outstream.write(str(result.sensitivity))
        outstream.write('\t')
        outstream.write(str(result.snr))                
        outstream.write('\n')


    @staticmethod
    def _emit_sensitivity_tally(sens_tally, outstream):
        # FORMAT: REF    sens
        for base1 in ('A', 'C', 'G', 'T'):
            for base2 in ('C', 'T'):
                for base3 in ('A', 'C', 'G', 'T'):
                    outstream.write('{}{}{}\t{}\n'.format(base1, base2, base3, sens_tally[base1 + base2 + base3]))


    def process_mpileup_stream(self, instream, outstream_vars, outstream_sens):
        # Stream the mpileup output and emit variants as they are detected.
        # In parallel, keep a running tally of sensitivity for each ref
        # triplet at non-blacklisted loci.
        #
        # There is complexity in this in that the ref triplet for a line 
        # contains the *next* base, which has not yet been processed.
        # Therefore variant writing lags behind reading by one base.
        #
        # The use of a blacklist also adds logic -- loci in the blacklist
        # should not have any variants emitted nor contribute to the 
        # sensitivity tally, but should be used in the ref triplets.

        sensitivity_tally = collections.Counter()

        ref_triplet = ['N'] * 3
        last_fields = None
        last_result = self.TestResult(False, False, None, None, None, None)

        line_count = 0
        start_time = time.clock()
        last_print_time = start_time
        last_print_count = line_count
        print_interval = 1000000
        parsing_time = 0
        test_time = 0
        blacklist_time = 0

        for line in instream:
            line_count += 1
            if line_count % print_interval == 0:
                this_print_time = time.clock()
                delta_print_time = this_print_time - last_print_time
                delta_print_count = line_count - last_print_count
                # TotalLoci  TotalTime  Rate  ParseTotal  TestTotal  BlacklistTotal  LociSinceLastStatus  TimeSinceLastStatus  RateSinceLastStatus
                logging.info('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
                    line_count, this_print_time - start_time, line_count / (this_print_time - start_time),
                    parsing_time, test_time, blacklist_time, 
                    delta_print_count, delta_print_time, delta_print_count / delta_print_time))
                last_print_time = this_print_time
                last_print_count = line_count

            temp = time.clock()
            fields = _parse_mpileup_line(line, self.min_bq, self.min_mq)    # type(fields) = MpileupField (chrom, pos, ref, alt_depth, depth, alt_counts)
            parsing_time += (time.clock() - temp)

            # Update the tally of the reference triplet, centred at the previous base
            del ref_triplet[0]
            ref_triplet.append(fields.ref)

            # Emit the previous line if a variant was found
            if last_result.was_detection == True:
                self._emit_detection_line(last_fields, last_result, ref_triplet, outstream_vars)
                # Clear last_result to ensure no tripping on future iterations.
                # This is necessary as the last_x = x code below is skipped for
                # unpowered and blacklisted sites, so last_x is not always updated.
                last_result = self.TestResult(False, False, None, None, None, None)

            # Skip further processing for high depth sites
            if fields.depth > self.max_dp:
                continue

            # Test for LSV
            temp = time.clock()
            result = self._test_for_lsv(fields.alt_depth, fields.depth)     # type(result) = self.TestResult (was_powered, was_detection, critval_E, critval_H, sensitivity, snr)
            test_time += (time.clock() - temp)

            # Skip futher processing for unpowered sites
            if result.was_powered == False:
                continue

            # Skip further processing for blacklisted sites
            # Note that this check is after the LSV test, although it could be 
            # placed before.  The reason for the current placement is that the 
            # interval intersection code is very slow at the moment, and it's 
            # actually faster to calculate LSV first and check overlap later.
            temp = time.clock()
            blacklist_intersection = self.blacklist.intersects(fields.chrom, fields.pos - 1)
            blacklist_time += (time.clock() - temp)
            if blacklist_intersection:
                continue

            # Update the sensitivity tally for the current line.
            sensitivity_tally[''.join(_normalise_triplet(ref_triplet)[0])] += result.sensitivity

            # Copy results to buffer for possible emission in the next iteration
            last_fields = fields
            last_result = result

        # Tail processing of the last variant
        if last_result.was_detection == True:
            self._emit_detection_line(last_fields, last_result, ref_triplet, outstream_vars)

        # Emit the sensitivities
        if outstream_sens != None:
            self._emit_sensitivity_tally(sensitivity_tally, outstream_sens)



if __name__ == '__main__':
    import argparse
    import doctest
    import sys

    doctest.testmod()

    parser = argparse.ArgumentParser(description='Find subclonal SNVs.')

    parser.add_argument('--mpileup', help='Input mpileup result (- for stdin) [-]', type=str, default='-')
    parser.add_argument('--blacklist', help='BED of blacklist regions [not used]', type=open)
    parser.add_argument('--variants', help='Destination for detected variants (- for stdout) [-]', type=str, default='-')
    parser.add_argument('--background', help='Destination for background coverage tally [not written]', type=str)
    parser.add_argument('--snr', help='Target minimum SNR for detections [10.0]', type=float, default=10.0)
    parser.add_argument('--vaf', help='Target subclonal variant AF [0.1]', type=float, default=0.1)
    parser.add_argument('--error', help='Per-read error rate [1e-2]', type=float, default=1.0e-2)
    parser.add_argument('--bq', help='Minimum base quality threshold (Phred-scale) [30]', type=int, default=30)
    parser.add_argument('--mq', help='Minimum mapping quality threshold (Phred-scale) [30]', type=int, default=30)
    parser.add_argument('--het', help='Heterozygous locus rate [1e-3]', type=float, default=1.0e-3)
    parser.add_argument('--scv', help='Subclonal variation rate [5e-7]', type=float, default=5.0e-7)
    parser.add_argument('--maxdp', help="Ignore loci with depth exceeding this value [100]", type = int, default=100)

    args = parser.parse_args()

    if args.mpileup == '-':
        instream = sys.stdin
    else:
        instream = open(args.mpileup, 'rt')

    blacklist = Blacklist()
    if args.blacklist:
        blacklist.load_bed(args.blacklist)

    if args.variants == '-':
        outstream = sys.stdout
    else:
        outstream = open(args.variants, 'wt')

    if args.background:
        background = open(args.background, 'wt')
    else:
        background = None

    logging.basicConfig(format='%(asctime)s\t%(relativeCreated)d\t%(levelname)s\t%(message)s', level=logging.INFO)

    detector = LssDetector(blacklist, args.het, args.scv, args.vaf, args.snr, args.error, args.bq, args.mq, args.maxdp)
    detector.process_mpileup_stream(instream, outstream, background)

    logging.shutdown()
