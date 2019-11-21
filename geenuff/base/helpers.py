import logging
import copy
import hashlib


##### General #####

def sequence_hash(sequence):
    sha1 = hashlib.sha1(sequence.encode())
    return sha1.hexdigest()


def in_enum_values(x, enum):
    return x in [item.value for item in enum]


def none_to_list(x):
    if x is None:
        return []
    else:
        assert isinstance(x, list)
        return x


def convert2list(obj):
    if isinstance(obj, list):
        out = obj
    elif isinstance(obj, set) or isinstance(obj, GeneratorType) or isinstance(obj, tuple):
        out = list(obj)
    else:
        out = [obj]
    return out


def full_db_path(path):
    sqlite_prefix = 'sqlite:///'
    if not path.startswith(sqlite_prefix):
        return sqlite_prefix + path
    return path


def db_attr_as_dict(orm_obj):
    """removes known, shared non-db entry attributes from a copy of orm_obj.__dict__
    can be used, e.g. to setup new db entry matching the here-by filtered one
    """
    exclude = ["_sa_instance_state", "handler"]
    out = copy.copy(orm_obj.__dict__)
    for item in exclude:
        if item in out:
            del out[item]
    return out


def chunk_str(string, length):
    for i in range(0, len(string), length):
        yield string[i:(i+length)]


def get_repr(class_name, params, addition=''):
    param_str = ', '.join('{}:{}'.format(k, v) for k, v in params.items())
    if addition:
        return class_name + '[' +  param_str + ', ' + addition + ']'
    else:
        return class_name + '[' +  param_str + ']'


##### Mapper #####

class NonMatchableIDs(Exception):
    pass


class Mapper(object):
    def __call__(self, key, *args, **kwargs):
        return key


class CheckMapper(Mapper):
    def __init__(self, keys):
        self.keys = set(keys)

    def __call__(self, key, *args, **kwargs):
        if key in self.keys:
            return key
        else:
            raise KeyError("{} not in known keys {}".format(key, self.keys))


class DictMapper(Mapper):
    def __init__(self, key_vals):
        if not isinstance(key_vals, dict):
            raise ValueError('key_vals must be a dict instance')
        self.key_vals = key_vals

    def __call__(self, key, *args, **kwargs):
        return self.key_vals[key]


def make_key_mapper(known_keys, other_keys):
    to_match_keys = set(known_keys)
    other_keys = set(other_keys)
    # if we have a perfect match
    if to_match_keys == other_keys or to_match_keys.issuperset(other_keys):
        return CheckMapper(to_match_keys)
    elif to_match_keys.issubset(other_keys):
        raise NonMatchableIDs('known is a subset of other keys')

    logging.info("attempting to match up, non-identical IDs")
    oth2known = {}
    # for each tree key, does it have exactly one match?
    for key in other_keys:
        matches = [x for x in to_match_keys if key in x]
        if len(matches) == 1:
            # setup dict[old_key] = new_key
            oth2known[key] = matches[0]
        elif len(matches) == 0:
            # pretending no match is ok, (a warning will be logged by are_keys_compatible) but seriously NCBI?
            logging.debug('no matches found for {} in known keys, e.g. {}'.format(
                key, list(known_keys)[:min(4, len(known_keys))]
            ))
        else:
            raise NonMatchableIDs('could not identify unique match for {}, but instead got {}'.format(key,
                                                                                                      matches))
    if not oth2known:  # bc we can still get an empty dict for complete unrelated gibberish at this point
        raise NonMatchableIDs('could not uniquely match up known: {} and other: {} keys'.format(known_keys,
                                                                                                other_keys))

    # check we matched trees to _unique_ fasta keys
    if len(oth2known.values()) != len(set(oth2known.values())):
        raise NonMatchableIDs('could not uniquely match up known: {} and other: {} keys'.format(known_keys,
                                                                                                other_keys))
    return DictMapper(oth2known)


def two_way_key_match(known_keys, other_keys):
    forward = True
    try:
        mapper = make_key_mapper(known_keys, other_keys)
    except NonMatchableIDs as e:
        try:
            mapper = make_key_mapper(other_keys, known_keys)
            forward = False
        except NonMatchableIDs:
            raise e
    return mapper, forward


def get_seqids_from_gff(gfffile):
    seqids = set()
    with open(gfffile) as f:
        for line in f:
            if not line.startswith('#'):
                seqids.add(line.split('\t')[0])
    return seqids

##### GFF start/end to GeenuFF #####

def get_strand_direction(gffentry):
    if gffentry.strand == '+':
        return True
    elif gffentry.strand == '-':
        return False
    else:
        raise ValueError('cannot interpret strand "{}"'.format(gffentry.strand))


def get_geenuff_start_end(gff_start, gff_end, is_plus_strand):
    gff_start, gff_end = to_count_from_0(gff_start), to_count_from_0(gff_end)

    if is_plus_strand:
        start = gff_start
        end = to_exclusive_end(gff_end, is_plus_strand)
    else:
        start = gff_end
        end = to_exclusive_end(gff_start, is_plus_strand)
    return start, end


def to_count_from_0(coord):
    return coord - 1


def to_exclusive_end(end, is_plus_strand):
    if is_plus_strand:
        return end + 1
    else:
        return end - 1


##### Reverse complement #####

def mk_rc_key():
    fw = "ACGTMRWSYKVHDBN"
    rv = "TGCAKYWSRMBDHVN"
    key = {}
    for f, r in zip(fw, rv):
        key[f] = r
    return key


# so one doesn't recalculate it for every call of revers_complement
REV_COMPLEMENT_KEY = mk_rc_key()


def reverse_complement(seq):
    key = REV_COMPLEMENT_KEY
    rc_seq = []
    for base in reversed(seq):
        try:
            rc_seq.append(key[base])
        except KeyError as e:
            raise KeyError('{} caused by non DNA character {}'.format(e, base))
    return rc_seq


##### Start/Stop codon detection #####

START_CODON = 'ATG'
START_CODON_COMP = ''.join(reverse_complement(START_CODON))
STOP_CODONS = ['TAG', 'TGA', 'TAA']
STOP_CODONS_COMP = [''.join(reverse_complement(c)) for c in STOP_CODONS]


def substr_seq(seq, start, end, is_plus_strand):
    """returns a substring of sequence according to geenuff coordinates and strand type"""
    if is_plus_strand:
        return seq[start:end]
    else:
        return seq[end+1:start+1]


def has_start_codon(seq, start, is_plus_strand):
    if is_plus_strand:
        return substr_seq(seq, start, start + 3, is_plus_strand) == START_CODON
    else:
        return substr_seq(seq, start, start - 3, is_plus_strand) == START_CODON_COMP


def has_stop_codon(seq, end, is_plus_strand):
    if is_plus_strand:
        return substr_seq(seq, end - 3, end, is_plus_strand) in STOP_CODONS
    else:
        return substr_seq(seq, end + 3, end, is_plus_strand) in STOP_CODONS_COMP

##### SQL alchemy core queue control #####

class Counter(object):
    def __init__(self, cl=None, at=0):
        self.cl = cl
        self._at = at

    def sync_with_db(self, session):
        from sqlalchemy import func
        new_at = session.query(func.max(self.cl.id)).one()[0]
        if new_at == None:
            new_at = 0
        self._at = new_at

    def __call__(self, *args, **kwargs):
        self._at += 1
        return self._at


class QueueController(object):
    def __init__(self, session, engine):
        self.session = session
        self.engine = engine
        self.ordered_queues = []

    def execute_so_far(self):
        conn = self.engine.connect()
        for queue in self.ordered_queues:
            if queue.queue:
                conn.execute(queue.action, queue.queue)
                del queue.queue[:]
        self.session.commit()
        logging.info('All core queues executed')


class CoreQueue(object):
    def __init__(self, action):
        self.queue = []
        self.action = action
