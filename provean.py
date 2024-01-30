# Packages
import argparse
import itertools
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree
from pathlib import Path

import Bio.Align
import Bio.SearchIO
import Bio.Seq
import Bio.SeqIO
import Bio.SeqRecord
from Bio.Blast.NCBIWWW import qblast as NCBI_API
from joblib import Parallel, delayed

# Executables paths. If an executable is not available in your PATH, then provide the full path to the executable (e.g: /home/user/Downloads/CDHIT/bin/cdhit)
CDHIT_EXE = 'cdhit'


AA_ALPHABET = list('GAVLMISTCPNQKRHDEFYW')

def replace_extended_IUPAC_residues_with_X(seq, translation_table=str.maketrans({residue:'X' for residue in 'BZJUO'})):
    str_seq = str(seq)
    str_seq = str_seq.translate(translation_table)
    return Bio.Seq.Seq(str_seq)

def get_seq_from_fasta_file(file):
    return replace_extended_IUPAC_residues_with_X(Bio.SeqIO.read(file, "fasta").seq)

def get_min_hit_seq_len(query_seq, min_hit_coverage):
    return int( len(query_seq) * min_hit_coverage )

def run_BLASTP(temp_dir, output_dir, query_seq, e_val_thr, open_penalty, extend_penalty, BLASTP_hitlist_size, save_intermediate_files):
    """
    """
    result = NCBI_API(
        program='blastp', 
        database='nr', 
        sequence=query_seq,
        expect=e_val_thr,
        gapcosts=f'{abs(open_penalty)} {abs(extend_penalty)}', # e.g '10 1'
        hitlist_size=BLASTP_hitlist_size
    )

    blastp_output_file_path = temp_dir / 'blastp_output.xml'
    with open(blastp_output_file_path, 'wt') as file_handle:
        file_handle.write(result.read())

    if save_intermediate_files:
        shutil.copy(src=blastp_output_file_path, dst=output_dir)

    return blastp_output_file_path

def extract_homologous_sequences_from_blastp_output_file(temp_dir, output_dir, min_hit_seq_len, blastp_output_file_path, save_intermediate_files):
    """
    """
    xml_tree = xml.etree.ElementTree.parse(blastp_output_file_path).getroot()

    homologous_sequences = []
    for hit in xml_tree.iter('Hit'):
        hit_id = hit.find('Hit_id').text
        hit_seq = hit.find('Hit_hsps').find('Hsp').find('Hsp_hseq').text
        hit_seq = Bio.Seq.Seq(hit_seq).replace('-','') # Remove gaps

        record = Bio.SeqRecord.SeqRecord(seq=hit_seq, id=hit_id, name='', description='')
        if len(record) >= min_hit_seq_len:
            homologous_sequences.append(record)

    homologous_sequences_fasta_file = temp_dir / 'homologous_sequences.fasta'
    Bio.SeqIO.write(sequences=homologous_sequences, handle=homologous_sequences_fasta_file, format='fasta')

    if save_intermediate_files:
        shutil.copy(src=homologous_sequences_fasta_file, dst=output_dir)

    return homologous_sequences_fasta_file

def get_hitID_sequence_map(homologous_sequences_fasta_file):
    """
    From the given fasta file, returns a dictionary where keys are sequence IDs and values are the corresponding amino acid sequence.
    """
    return Bio.SeqIO.to_dict(Bio.SeqIO.parse(homologous_sequences_fasta_file, "fasta"))

def run_cdhit(temp_dir, output_dir, homologous_sequences_fasta_file, clust_seq_id_thr, save_intermediate_files):
    """
    """
    cdhit_output_file_path = temp_dir / 'cdhit_clusters.clstr'
    command = [
        CDHIT_EXE,
        '-i', str(homologous_sequences_fasta_file),
        '-o', str(cdhit_output_file_path.with_suffix('')),
        '-c', str(clust_seq_id_thr),
        '-d', '0' # Unlimited ID length (by default cdhit truncates IDs longer than 20 characters)
    ]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL)

    if save_intermediate_files:
        shutil.copy(src=cdhit_output_file_path, dst=output_dir)

    return cdhit_output_file_path

def add_cluster_sequences_to_cluster_map(clusters_map, lines, hitID_sequence_map):
    """
    """
    cluster_reference_sequence = None
    cluster_sequences = []
    for line in lines[1:]:
        hitID = line.split(', >')[1].split('...')[0] # e.g: 0	182aa, >UniRef90_A0A7X8D1Y4... at 82.42% => UniRef90_A0A7X8D1Y4... at 82.42% => UniRef90_A0A7X8D1Y4
        hit_sequence_SeqRecord = hitID_sequence_map[hitID]

        if line.strip().endswith('*'):
            cluster_reference_sequence = replace_extended_IUPAC_residues_with_X(hit_sequence_SeqRecord.seq)
        else:
            cluster_sequences.append( replace_extended_IUPAC_residues_with_X(hit_sequence_SeqRecord.seq) )

    assert cluster_reference_sequence is not None, f'Failed to find the reference sequence (*) from the following cluster.\n {lines} '
    clusters_map[cluster_reference_sequence] = cluster_sequences
    return

def parse_cdhit_clusters_file(cdhit_output_file_path, hitID_sequence_map, min_n_clusters):
    """
    """
    # Example of the format of a cdhit .clstr file:
    #   >Cluster 0
    #   0	329aa, >UniRef90_R6E5V7... *
    #   >Cluster 1
    #   0	182aa, >UniRef90_A0A7X8D1Y4... at 82.42%
    #   1	116aa, >UniRef90_A0A7X8D1Y4... at 83.62%
    #   2	327aa, >UniRef90_A0A7C7EQS6... *
    #   etcetc
    with open(cdhit_output_file_path) as file_handle:
        cdhit_clusters_file_lines = file_handle.readlines()
        
    cluster_idexes = [
        index 
        for index, line in enumerate(cdhit_clusters_file_lines)
        if line.startswith('>Cluster')
    ]
    if len(cluster_idexes) < min_n_clusters:
        raise ValueError(f'CDHIT only returned {len(cluster_idexes)} clusters, but {min_n_clusters=}. To solve this problem, adapt some of the parameters, such as min_n_clusters or clust_seq_id_thr.')

    # Its simpler to extract the clusters in reverse order
    clusters_map = {}
    for cluster_index in reversed(cluster_idexes):
        add_cluster_sequences_to_cluster_map(
            clusters_map, 
            cdhit_clusters_file_lines[cluster_index:], 
            hitID_sequence_map
        )

        del cdhit_clusters_file_lines[cluster_index:]

    return clusters_map

def get_semi_global_aligner(open_penalty, extend_penalty, substitution_matrix):
    """
    Semi global alignment as defined in the PROVEAN paper, which the authors define as 'no penalty on end gaps in 
    global alignment'.
    """
    return Bio.Align.PairwiseAligner(
        open_gap_score = open_penalty, 
        extend_gap_score = extend_penalty,
        substitution_matrix = substitution_matrix,
        mode = 'global', query_end_gap_score = 0, target_end_gap_score = 0,
    )

def select_top_clusters(clusters_map, max_n_clusters, semi_global_aligner, query_seq):
    """
    """
    # No need to remove any cluster if there aren't too many
    if len(clusters_map) <= max_n_clusters:
        return

    scores = []
    for cluster_reference_seq in clusters_map.keys():
        best_alignment = semi_global_aligner.align(query_seq, cluster_reference_seq)[0]

        scores.append(best_alignment.score)

        setattr(cluster_reference_seq, 'aligned_query', best_alignment[0])
        setattr(cluster_reference_seq, 'aligned_target', best_alignment[1])
        setattr(cluster_reference_seq, 'unmutated_score', best_alignment.score)


    min_score_of_selected_clusters = sorted(scores, reverse=True)[max_n_clusters-1] # This is the score of the 45th highest scoring cluster reference sequence
    cluster_reference_sequences_to_remove = [
        cluster_reference_seq
        for cluster_reference_seq in clusters_map.keys()
        if cluster_reference_seq.unmutated_score < min_score_of_selected_clusters
    ]

    for cluster_reference_seq in cluster_reference_sequences_to_remove:
        del clusters_map[cluster_reference_seq]

    return
    
def calculate_all_unmutated_alignment_scores(clusters_map, semi_global_aligner, query_seq):
    """
    """
    for seq in clusters_map.keys():
        best_alignment = semi_global_aligner.align(query_seq, seq)[0]

        setattr(seq, 'aligned_query', best_alignment[0])
        setattr(seq, 'aligned_target', best_alignment[1])
        setattr(seq, 'unmutated_score', best_alignment.score)

    for cluster_of_sequences in clusters_map.values():
        for seq in cluster_of_sequences:
            best_alignment = semi_global_aligner.align(query_seq, seq)[0]

            setattr(seq, 'aligned_query', best_alignment[0])
            setattr(seq, 'aligned_target', best_alignment[1])
            setattr(seq, 'unmutated_score', best_alignment.score)

    return

def variant_seq_generator(query_seq, substitutions_only):
    """
    Generates all posible substitutions, insertions and deletions of the given query sequence.
    Returns the mutated sequence and the corresponding variant annotation (ie G2L, G2del, G2insL)
    """
    sequence_length = len(query_seq)
    mutable_seq = Bio.Seq.MutableSeq(query_seq)
    for index, reference_AA in enumerate(mutable_seq):
        assert mutable_seq == query_seq # Make sure we have indeed restored the original sequence before the next iteration

        # Substitutions
        for AA in AA_ALPHABET:
            mutable_seq[index] = AA
            yield Bio.Seq.Seq(mutable_seq), f'{reference_AA}{index+1}{AA}'
        mutable_seq[index] = reference_AA # Reset to original sequence
        
        if substitutions_only:
            continue

        # Insertions
        is_last_residue = index-1 == sequence_length
        if not is_last_residue: # Inserting a residue after the last residue doesn't make any sense
            for AA in AA_ALPHABET:
                mutable_seq.insert(index+1, AA)
                yield Bio.Seq.Seq(mutable_seq), f'{reference_AA}{index+1}ins{AA}'
                del mutable_seq[index+1] # Reset to original sequence

        # Deletion
        del mutable_seq[index]
        yield Bio.Seq.Seq(mutable_seq), f'{reference_AA}{index+1}del'
        mutable_seq.insert(index, reference_AA) # Reset to original sequence

    return

def get_aligned_position(cluster_seq, query_residue, query_position):
    """
    Variant annotation positions (e.g G2L) correspond to the position in the original sequence, but here we are dealing with aligned
    sequences, and therefore positions need to be recalculated to account for possible alignment gaps.
    
    Example
    -------
    Query:  A - N K <- In the alignment, residue N is at position 3, while it's at position 2 in the original sequence
    Target: A L N K
    """
    sequence_position = -1
    for aligned_position, AA in enumerate(cluster_seq.aligned_query):
        if AA != '-':
            sequence_position += 1

        if query_position == sequence_position:
            break
    
    assert query_residue == AA
    return aligned_position

def substitution_rescoring(cluster_seq, variant_annotation, substitution_matrix):
    """
    """
    query_residue, query_position, variant_residue = variant_annotation[0], int(variant_annotation[1:-1]), variant_annotation[-1] # e.g: 'G123A -> 'G', '123', 'A'.
    query_position -= 1 #  Positions in variant annotations are 1-based, unlike Python which is 0-based.

    aligned_position = get_aligned_position(cluster_seq, query_residue, query_position)
    target_residue = cluster_seq.aligned_target[aligned_position] # e.g: 'A'

    # Variants on an aligned gap don't change the score given gap penalties are the same for all residues
    if target_residue == '-':
        return cluster_seq.unmutated_score

    original_position_score = substitution_matrix[query_residue][target_residue]
    variant_position_score = substitution_matrix[variant_residue][target_residue]

    return cluster_seq.unmutated_score - original_position_score + variant_position_score

def calculate_average(values):
    return sum(values) / len(values)

def calculate_variant_seq_PROVEAN_score(variant_seq, variant_annotation, semi_global_aligner, substitution_matrix, clusters_map):
    """
    """
    is_substitution = 'del' not in variant_annotation and 'ins' not in variant_annotation

    clusters_average_delta_scores = []
    for cluster_reference_seq, cluster_of_sequences in clusters_map.items():
        # For substitutions we simply look up the score in the substitution_matrix to recalculate the score
        # For indels we have to redo the alignment from scratch, which is much slower.
        if is_substitution:
            cluster_delta_scores = [
                substitution_rescoring(cluster_seq, variant_annotation, substitution_matrix) - cluster_seq.unmutated_score
                for cluster_seq in itertools.chain([cluster_reference_seq], cluster_of_sequences)
            ]
        else:
            cluster_delta_scores = [
                semi_global_aligner.align(variant_seq, cluster_seq).score - cluster_seq.unmutated_score
                for cluster_seq in itertools.chain([cluster_reference_seq], cluster_of_sequences)
            ]

        clusters_average_delta_scores.append( calculate_average(cluster_delta_scores) )        

    PROVEAN_score = calculate_average(clusters_average_delta_scores)
    
    return round(PROVEAN_score, 3), variant_annotation

def calculate_provean_scores(output_dir, query_seq, substitutions_only, clusters_map, max_n_clusters, substitution_matrix, open_penalty, extend_penalty, n_cores):
    """
    """
    substitution_matrix = Bio.Align.substitution_matrices.load(substitution_matrix)
    semi_global_aligner = get_semi_global_aligner(open_penalty, extend_penalty, substitution_matrix)
    
    select_top_clusters(clusters_map, max_n_clusters, semi_global_aligner, query_seq)
    calculate_all_unmutated_alignment_scores(clusters_map, semi_global_aligner, query_seq)

    with open(output_dir / 'PROVEAN_scores.csv', 'wt') as output_file_handle, Parallel(n_cores, return_as='generator') as parallel:
        output_file_handle.write('variant,PROVEAN_score\n') # Column names
        
        calculate_variant_seq_PROVEAN_score_delayed = delayed(calculate_variant_seq_PROVEAN_score)
        results = parallel(
            calculate_variant_seq_PROVEAN_score_delayed(variant_seq, variant_annotation, semi_global_aligner, substitution_matrix, clusters_map)
            for variant_seq, variant_annotation in variant_seq_generator(query_seq, substitutions_only)
        )
        for PROVEAN_score, variant_annotation in results:
            output_file_handle.write(f'{variant_annotation},{str(PROVEAN_score)}\n')

    return

def provean(
        query_file,
        output_dir,
        substitutions_only=True,
        save_intermediate_files=False,
        homologous_sequences_fasta_file=None,
        BLASTP_hitlist_size=5000,
        e_val_thr=0.1, 
        min_hit_coverage=0.3,
        clust_seq_id_thr=0.8,
        min_n_clusters=15, max_n_clusters=45,
        substitution_matrix='BLOSUM62', 
        open_penalty=-10, extend_penalty=-1,
        n_cores=2
    ):
    """
    Calculates the PROVEAN scores of all possible single residue substitutions of the given query protein. To also calculate
    the scores of all indels, set substitutions_only to False (note that indels are much slower to compute).
    If a FASTA file with homologous sequences is already available, it can be provided through the homologous_sequences_fasta_file
    parameter, which will avoid having to run BLASTP.
    """
    # All intermediate files are kept in a temporary directory. Users can chose to copy these files to 
    # their output directory, such as the BLASTP xml file or the FASTA file with the homologous sequences.
    output_dir.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='PROVEAN_') as temp_dir:
        temp_dir = Path(temp_dir)

        query_seq = get_seq_from_fasta_file(file=query_file)

        if not homologous_sequences_fasta_file:
            print('Running BLASTP via the NCBI REST API ...')
            blastp_output_file_path = run_BLASTP(temp_dir, output_dir, query_seq, e_val_thr, open_penalty, extend_penalty, BLASTP_hitlist_size, save_intermediate_files)
            homologous_sequences_fasta_file = extract_homologous_sequences_from_blastp_output_file(
                temp_dir, output_dir,
                min_hit_seq_len=get_min_hit_seq_len(query_seq, min_hit_coverage),
                blastp_output_file_path=blastp_output_file_path,
                save_intermediate_files=save_intermediate_files
            )
            
        hitID_sequence_map = get_hitID_sequence_map(homologous_sequences_fasta_file)

        cdhit_output_file_path = run_cdhit(temp_dir, output_dir, homologous_sequences_fasta_file, clust_seq_id_thr, save_intermediate_files)
        clusters_map = parse_cdhit_clusters_file(cdhit_output_file_path, hitID_sequence_map, min_n_clusters)
        
        print('Calculating PROVEAN scores ...')
        calculate_provean_scores(output_dir, query_seq, substitutions_only, clusters_map, max_n_clusters, substitution_matrix, open_penalty, extend_penalty, n_cores)
    
    return


def argument_parser():
    """
    """
    parser = argparse.ArgumentParser(
        description = "Description:\n  Python implementation of the PROVEAN pipeline to predict variant deleteriousness.",
        usage = 'provean <query_file> <output_dir> [parameters]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        'query_file', type = str,
        help = 'Full path of a FASTA file containing the query sequence (1 chain only)'
    )

    parser.add_argument(
        'output_dir', type = str,
        help = "Full path of the desired output directory (will be created if it doesn't exist)"
    )
    

    parser.add_argument(
        '-substitutions_only', type = bool, default = True, metavar = '',
        help = "If set to False, PROVEAN scores of indels are also calculated. Note that this is slow"
    )

    parser.add_argument(
        '-save_intermediate_files', type = bool, default = False, metavar = '',
        help = "Save intermediate files generated during the PROVEAN pipeline, such as the BLASTP result or the CDHIT clusters file"
    )

    parser.add_argument(
        '-homologous_sequences_fasta_file', default=None, metavar = '',
        help = 'If already calculated, you can provide a FASTA file containing homologous sequences, which will avoid running BLASTP'
    )

    parser.add_argument(
        '-BLASTP_hitlist_size', type = int, default = 5000, metavar = '',
        help = 'Maximum number of sequences returned by BLASTP'
    )

    parser.add_argument(
        '-e_val_thr', type = float, default = 0.1, metavar = '',
        help = 'E value threshold to use for BLASTP'
    )
    
    parser.add_argument(
        '-min_hit_coverage', type = float, default = 0.3, metavar = '',
        help = 'Minimum relative length of hit sequences with respect to the query sequence. For example, if the query sequence has 100 residues and min_hit_coverage=0.3, then only hit sequences with at least 30 residues are taken from BLASTP'
    )

    parser.add_argument(
        '-clust_seq_id_thr', type = float, default = 0.8, metavar = '',
        help = "Sequence identity threshold used by CDHIT to cluster the homologous sequences"
    )
    
    parser.add_argument(
        '-min_n_clusters', type = int, default = 15, metavar = '',
        help = "Minimum number of clusters. If CDHIT returns fewer clusters, an error will be raised and the program will stop"
    )

    parser.add_argument(
        '-max_n_clusters', type = int, default = 45, metavar = '',
        help = "Maximum number of clusters. If CDHIT returns more clusters, then the top 45 will be selected"
    )

    parser.add_argument(
        '-substitution_matrix', type = str, default = 'BLOSUM62', metavar = '',
        help = "Name of the substitution matrix used to calculate the PROVEAN scores. There is no reason to change it"
    )

    parser.add_argument(
        '-open_penalty', type = int, default = -10, metavar = '',
        help = "Open gap penalty used for BLASTP and to calculate the PROVEAN scores"
    )

    parser.add_argument(
        '-extend_penalty', type = int, default = -1, metavar = '',
        help = "Extension gap penalty used for BLASTP and to calculate the PROVEAN scores"
    )

    parser.add_argument(
        '-n_cores', type = int, default = 2, metavar = '',
        help = "Number of parallele cores to use"
    )
    
    return parser


if __name__ == '__main__':
    parser = argument_parser()
    args = parser.parse_args()

    provean(
        query_file=Path(args.query_file),
        output_dir=Path(args.output_dir),
        substitutions_only=args.substitutions_only,
        save_intermediate_files=args.save_intermediate_files,
        homologous_sequences_fasta_file=args.homologous_sequences_fasta_file,
        BLASTP_hitlist_size=args.BLASTP_hitlist_size,
        e_val_thr=args.e_val_thr, 
        min_hit_coverage=args.min_hit_coverage,
        clust_seq_id_thr=args.clust_seq_id_thr,
        min_n_clusters=args.min_n_clusters, max_n_clusters=args.max_n_clusters,
        substitution_matrix=args.substitution_matrix, 
        open_penalty=args.open_penalty, extend_penalty=args.extend_penalty,
        n_cores=args.n_cores
    )