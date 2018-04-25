import re
import pandas as pd
import subprocess
from Bio import SeqIO

__author__ = 'Natalia Quinones-Olvera'
__email__ = "nquinones@g.harvard.edu"

MAF_CONVERT_PATH = 'maf-convert'
SAMTOOLS_PATH = 'samtools'
LASTAL_PATH = 'lastal'


# ..............................FUNCTIONS....................................


def build_lastal_cmd(db, fastq, out, fmt='fastq'):
    '''
    '''
    if fmt is 'fastq':
        q = 1
    elif fmt is 'fasta':
        q = 0
    else:
        print('invalid format, will read as fastq')

    lastal_cmd = ('{} '
                  '-r6 '
                  '-q12 '
                  '-a15 '
                  '-b3 '
                  '-e150 '
                  '-m100 '
                  '-Q{} '
                  '-j4 '
                  '{} '
                  '{} '
                  '| last-split '
                  '-m1 '
                  '> '
                  '{}'.format(LASTAL_PATH,
                              q,
                              db,
                              fastq,
                              out))

    return lastal_cmd


def build_mafconvert_cmd(maf_file, sam_file):
    '''
    '''
    mafconv_cmd = ('{} '
                   '-d '
                   'sam '
                   '{} '
                   '> '
                   '{}').format(MAF_CONVERT_PATH,
                                maf_file,
                                sam_file)

    return mafconv_cmd


def build_samtobam_cmd(sam_file, bam_file):
    '''
    '''
    samtobam_cmd = ('{} '
                    'view '
                    '-S '
                    '-b '
                    '{} '
                    '| '
                    '{} '
                    'sort '
                    '- > '
                    '{}').format(SAMTOOLS_PATH,
                                 sam_file,
                                 SAMTOOLS_PATH,
                                 bam_file)

    return samtobam_cmd


def build_bamindex_cmd(bam_file):
    '''
    '''
    bamindex_cmd = ('{} '
                    'index '
                    '{}').format(SAMTOOLS_PATH,
                                 bam_file)

    return bamindex_cmd


def master_cmd(db, read_file, fmt, out_files):
    '''
    '''
    out_maf = '{}.maf'.format(out_files)
    out_sam = '{}.sam'.format(out_files)
    out_bam = '{}.bam'.format(out_files)

    cmd1 = build_lastal_cmd(db, read_file, out_maf, fmt=fmt)
    cmd2 = build_mafconvert_cmd(out_maf, out_sam)
    cmd3 = build_samtobam_cmd(out_sam, out_bam)
    cmd4 = build_bamindex_cmd(out_bam)

    subprocess.call(cmd1, shell=True)
    subprocess.call(cmd2, shell=True)
    subprocess.call(cmd3, shell=True)
    subprocess.call(cmd4, shell=True)


def parse_cigar(cigar):
    '''
    '''
    cigar_pattern = re.compile(r'\d+[^\d]')

    cigar_split = re.findall(cigar_pattern, cigar)

    mismatches = 0
    insertions = 0
    deletions = 0
    hardclip = 0

    for element in cigar_split:
        # Mismatches, X
        if 'X' in element:
            splited_x = element.split('X')[0]
            mismatches = mismatches + int(splited_x)

        # Insertions, I
        if 'I' in element:
            splited_i = element.split('I')[0]
            insertions = insertions + int(splited_i)

        # Deletions, D
        if 'D' in element:
            splited_d = element.split('D')[0]
            deletions = deletions + int(splited_d)

        # Clipping, H
        if 'H' in element:
            splited_h = element.split('H')[0]
            hardclip = hardclip + int(splited_h)

    return [insertions, deletions, mismatches, hardclip]


def parse_sam(sam_file):
    '''
    '''
    df = pd.read_table(sam_file,
                       skiprows=[0, 1],
                       header=None,
                       sep='\t')
    df.columns = ['query',
                  'flag',
                  'target',
                  'from',
                  'map_q',
                  'cigar',
                  'rnext',
                  'pnext',
                  'tlen',
                  'seq',
                  'q',
                  'NM',
                  'AS']

    df['len'] = df['seq'].str.len()

    df['NM'] = df['NM'].str.replace('NM:i:', '')
    df['AS'] = df['AS'].str.replace('AS:i:', '')

    df[['I', 'D', 'M', 'H']] = pd.DataFrame(df['cigar'].apply(parse_cigar).values.tolist())

    df['to'] = df['from'] + df['len']

    df = df[['query', 'target', 'AS', 'len', 'I', 'D', 'M', 'H']]
    df['AS'] = pd.to_numeric(df['AS'])
    df['len'] = pd.to_numeric(df['len'])

    return df


def keep_highscore(df):
    '''
    '''
    df = df.sort_values(by='AS', ascending=False)
    df = df.drop_duplicates('query', keep='first')

    return df


def compete_sams(sam_list):
    '''
    '''
    df_merge = pd.concat(sam_list)

    df_merge = df_merge.sort_values(by='AS', ascending=False)
    df_merge = df_merge.drop_duplicates('query', keep='first')

    return df_merge


def fastq_len(file):
    '''
    '''
    ids_l = []
    len_l = []
    status_l = []

    for record in SeqIO.parse(file, "fastq"):
        if '_fail_' in record.description:
            status = 'fail'
        elif '_pass_' in record.description:
            status = 'pass'
        else:
            status = '?'

        ids_l.append(record.id)
        len_l.append(len(record.seq))
        status_l.append(status)

    df = pd.DataFrame()
    df['id'] = ids_l
    df['len'] = len_l
    df['status'] = status_l

    return df

# ..........................................................................

if __name__ == '__main__':
    pass
