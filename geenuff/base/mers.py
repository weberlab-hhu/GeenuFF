import itertools
import copy

from . import helpers

class MerCounter(object):
    amb = 'ambiguous_mers'

    def __init__(self, k):
        self.k = k
        self.counts = {}
        # calculate all possible mers of this length, and set counter to 0
        for mer in itertools.product('ATCG', repeat=k):
            mer = ''.join(mer)
            self.counts[mer] = 0
        self.counts[MerCounter.amb] = 0
        # most recent base pairs of up to length k
        self.sliding_mer = []

    def export(self):
        out = copy.deepcopy(self.counts)
        for key in self.counts:
            if key != MerCounter.amb:
                # collapse to cannonical kmers
                rc_key = helpers.reverse_complement(key)
                if key != min(key, rc_key):
                    out[rc_key] += out[key]
                    out.pop(key)
        return out

    def add_sequence(self, sequence):
        for i in range(len(sequence) -  self.k + 1):
            mer = sequence[i:i+self.k]
            try:
                self.counts[mer] += 1
            except KeyError:
                self.counts[MerCounter.amb] += 1
