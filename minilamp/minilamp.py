import numpy as np
from scipy import stats
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

__author__ = 'Natalia Quinones-Olvera'
__email__ = "nquinones@g.harvard.edu"

# ..............................FUNCTIONS....................................


def get_primers(region_coord, full_region):
    '''
    Builds LAMP primers, as specified from LAMP regions (f3, f2, etc.)
    Saves an additional 'witness' sequence to keep track of the positions
    of the regions.

    In:
    region_coord =
    full_region =
    '''

    # FIP
    f1_seq = full_region[region_coord['f1'][0]:region_coord['f1'][1]]
    f2_seq = full_region[region_coord['f2'][0]:region_coord['f2'][1]]

    fip = (f1_seq.reverse_complement() +
           'TTTT' +
           f2_seq)

    fip_witn = ('1' * len(f1_seq) +
                '0' * 4 +
                '1' * len(f2_seq))

    # BIP
    b1c_seq = full_region[region_coord['b1c'][0]:region_coord['b1c'][1]]
    b2c_seq = full_region[region_coord['b2c'][0]:region_coord['b2c'][1]]

    bip = (b1c_seq +
           'TTTT' +
           b2c_seq.reverse_complement())

    bip_witn = ('1' * len(b1c_seq) +
                '0' * 4 +
                '1' * len(b2c_seq))

    # FOP
    f3_seq = full_region[region_coord['f3'][0]:region_coord['f3'][1]]
    fop = f3_seq

    fop_witn = '1' * len(fop)

    # BOP
    b3c_seq = full_region[region_coord['b3c'][0]:region_coord['b3c'][1]]
    bop = b3c_seq.reverse_complement()

    bop_witn = '1' * len(bop)

    primers_dict = {'fip': fip,
                    'bip': bip,
                    'f3': fop,
                    'b3': bop}

    witn_dict = {'fip': fip_witn,
                 'bip': bip_witn,
                 'f3': fop_witn,
                 'b3': bop_witn}

    return primers_dict, witn_dict


def get_dumbbell(region_coord, full_region):
    '''

    Requires function: 'get_primers()'
    '''

    primers_ds = get_primers(region_coord, full_region)
    primers = primers_ds[0]
    witness = primers_ds[1]

    end_f2 = region_coord['f2'][1]
    start_b2c = region_coord['b2c'][0]

    dumbell_seq = (primers['bip'] +
                   full_region[end_f2:start_b2c].reverse_complement() +
                   primers['fip'].reverse_complement())

    witness_middle = ('0' * (region_coord['f1'][0] - region_coord['f2'][1]) +   # space between f2-f1
                      '1' * (region_coord['f1'][1] - region_coord['f1'][0]) +   # length of f1
                      '0' * (region_coord['b1c'][0] - region_coord['f1'][1]) +  # space between f1-bc1
                      '1' * (region_coord['b1c'][1] - region_coord['b1c'][0]) + # length of bc1
                      '0' * (region_coord['b2c'][0] - region_coord['b1c'][1]))  # space between bc1-bc2

    dumbell_witness = (witness['bip'] +
                       witness_middle[::-1] +
                       witness['fip'][::-1])

    return dumbell_seq, dumbell_witness


def find_loop(witness_seq):
    '''
    Finds the position in a 0's and 1's sequence where
    you find the n=3 transition 0-1 from back to front.
    i.e.
    0011100011100011100
         ^

    Used on the witness string to find the position that
    is not forming the loop (to find what will be appended)
    '''

    # invert
    witness_seq_i = witness_seq[::-1]
    n = 3

    start = witness_seq_i.find('10')
    while start >= 0 and n > 1:
        start = witness_seq_i.find('10', start + 2)
        n -= 1

    position = len(witness_seq_i) - start - 1

    # return (position, witness_seq[:position], witness_seq[position:])
    return position


def amplify(dumbbell_tup, lamp_cycl):
    '''
    '''

    for i in range(lamp_cycl):
        # find the part that corresponds to the loop
        cut_coord = find_loop(dumbbell_tup[1])
        # take the sequence and add the reverse
        # transcribed of what came before the loop
        seq = dumbbell_tup[0] + dumbbell_tup[0][:cut_coord].reverse_complement()
        # same for whitness, but inverted instead
        # of reverse transcribed
        witn = dumbbell_tup[1] + (dumbbell_tup[1][:cut_coord])[::-1]

        # rewrite the dumbbell_tup structure
        dumbbell_tup = (seq, witn)

    # return dumbbell_tup
    return dumbbell_tup[0]


def main_amplification(region_coord, full_region, lamp_cycl):
    '''
    '''
    dumbbell_tup = get_dumbbell(region_coord, full_region)

    amplicon = amplify(dumbbell_tup, lamp_cycl)

    return amplicon


def make_read_v1(amplicon, model):
    '''
    '''
    # read model parameters
    ins_p = model['ins_p']
    del_p = model['del_p']
    sub_p = model['sub_p']

    len_params = model['len_fit']

    # choose length from longnorm distribution with
    # model parameters
    choose_len = 0
    while (choose_len is 0 or
           choose_len > len(amplicon)):
        choose_len = int(stats.lognorm.rvs(len_params[0],
                                           len_params[1],
                                           len_params[2],
                                           1))

    # choose a random start place for the read region
    # (where the read fits)
    start = np.random.choice(range(1, len(amplicon) - choose_len), 1)[0]
    end = start + choose_len

    select_seq = amplicon[start:end]

    # add mutations
    new_seq = []
    wit_seq = []

    for base in select_seq:
        # deletion
        if np.random.uniform() < del_p:
            new_base = ''
            witness = 'D'
        else:
            # substitution
            if np.random.uniform() < sub_p:
                nts = ['A', 'C', 'G', 'T']
                del nts[nts.index(base)]
                new_base = np.random.choice(nts)
                witness = 'M'
            else:
                new_base = base
                witness = '-'

        new_seq.append(new_base)
        wit_seq.append(witness)

        # insertion
        if np.random.uniform() < ins_p:
            nts = ['A', 'C', 'G', 'T']
            ins_base = np.random.choice(nts)
            ins_witness = 'I'

            new_seq.append(ins_base)
            wit_seq.append(ins_witness)

    new_seq = ''.join(new_seq)
    wit_seq = ''.join(wit_seq)

    # return str(select_seq), new_seq, wit_seq
    return new_seq


def make_read_v2(amplicon, model):
    '''
    '''
    # read model parameters
    ins_params = model['ins_fit']
    del_params = model['del_fit']
    sub_params = model['sub_fit']

    len_params = model['len_fit']

    # choose length from longnorm distribution with
    # model parameters
    choose_len = 0
    while (choose_len is 0 or
           choose_len > len(amplicon)):
        choose_len = int(stats.lognorm.rvs(len_params[0],
                                           len_params[1],
                                           len_params[2],
                                           1))

    # choose a random start place for the read region
    # (where the read fits)
    start = np.random.choice(range(1, len(amplicon) - choose_len), 1)[0]
    end = start + choose_len

    select_seq = amplicon[start:end]

    # select a deletion probability
    del_p = 0
    # make sure you don't go bellow 0.
    while del_p <= 0:
        del_p = stats.norm.rvs(loc=del_params[0],
                               scale=del_params[1],
                               size=1)[0]

    sub_p = stats.norm.rvs(loc=sub_params[0],
                           scale=sub_params[1],
                           size=1)[0]

    ins_p = stats.expon.rvs(loc=ins_params[0],
                            scale=ins_params[1],
                            size=1)[0]

    # add mutations
    new_seq = []
    wit_seq = []

    for base in select_seq:
        # deletion
        if np.random.uniform() < del_p:
            new_base = ''
            witness = 'D'
        else:
            # substitution
            if np.random.uniform() < sub_p:
                nts = ['A', 'C', 'G', 'T']
                del nts[nts.index(base)]
                new_base = np.random.choice(nts)
                witness = 'M'
            else:
                new_base = base
                witness = '-'

        new_seq.append(new_base)
        wit_seq.append(witness)

        # insertion
        if np.random.uniform() < ins_p:
            nts = ['A', 'C', 'G', 'T']
            ins_base = np.random.choice(nts)
            ins_witness = 'I'

            new_seq.append(ins_base)
            wit_seq.append(ins_witness)

    new_seq = ''.join(new_seq)
    wit_seq = ''.join(wit_seq)

    # return str(select_seq), new_seq, wit_seq
    return new_seq


def simulate_reads_v1(region_coord, full_region, lamp_cycl, model, number, fasta_file):
    '''
    Generates sequencing reads for a
    LAMP amplification - MinION sequencing process

    Requires functions: 'main_amp()', 'make_read()'

    In:
    region_coord: dictionary with lamp regions' coords
    full_sequence: the full sequence from where the coords come
    lamp_cycl: number of lamp iterations
    model: dictionary with in/del/sub parameters
    number: number of reads to generate
    fasta_file: fasta file to write output
    '''

    # generate amplicon
    amplicon = main_amplification(region_coord, full_region, lamp_cycl)

    # generate read by read
    all_reads = []

    for i in range(number):

        # function make_read makes read by selecting a length,
        # a random position of the amplicon and adding
        # in/del/sub according to model
        seq_obj = Seq(make_read_v1(amplicon, model))

        # names just assigned secuentially
        seq_id = 'seq{}'.format(i)

        # make into SeqRecord object and append to list
        seq_rec = SeqRecord(seq_obj, id=seq_id, description='')
        all_reads.append(seq_rec)

    # write fasta
    SeqIO.write(all_reads, fasta_file, "fasta")

    # return all_reads


def simulate_reads_v2(region_coord, full_region, lamp_cycl, model, number, fasta_file):
    '''
    Generates sequencing reads for a
    LAMP amplification - MinION sequencing process

    Requires functions: 'main_amp()', 'make_read()'

    In:
    region_coord: dictionary with lamp regions' coords
    full_sequence: the full sequence from where the coords come
    lamp_cycl: number of lamp iterations
    model: dictionary with in/del/sub parameters
    number: number of reads to generate
    fasta_file: fasta file to write output
    '''

    # generate amplicon
    amplicon = main_amplification(region_coord, full_region, lamp_cycl)

    # generate read by read
    all_reads = []

    for i in range(number):

        # function make_read makes read by selecting a length,
        # a random position of the amplicon and adding
        # in/del/sub according to model
        seq_obj = Seq(make_read_v2(amplicon, model))

        # names just assigned secuentially
        seq_id = 'seq{}'.format(i)

        # make into SeqRecord object and append to list
        seq_rec = SeqRecord(seq_obj, id=seq_id, description='')
        all_reads.append(seq_rec)

    # write fasta
    SeqIO.write(all_reads, fasta_file, "fasta")

    # return all_reads

# ..........................................................................

if __name__ == '__main__':
    pass
