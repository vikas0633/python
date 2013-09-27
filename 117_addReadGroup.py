import os
import sys
import glob
import pysam
import argparse
import multiprocessing
 
def get_args():
    '''Parse sys.argv'''
    parser = argparse.ArgumentParser()
    parser.add_argument('--cores', type=int,
                        help='the number of avalible processor cores')
    parser.add_argument('-i','--input-dir', required=True,
                        help='The input directory containing the bam files.')
    parser.add_argument('-CN',type=str,
                        help="Name of sequencing center producing the read. \
                        GATK not required.")
    parser.add_argument('-DS',type=str,
                        help="Description. GATK Not Required. ")
    parser.add_argument('-DT',type=str, required=True,
                        help="Date the run was produced (ISO8601 date or date/time). \
                        GATK Not Required. ")
    parser.add_argument('-PI',type=int, required=True,
                        help="Predicted median insert size. GATK Not Required.")
    parser.add_argument('-PL',type=str, required=True,
                        choices = ['CAPILLARY', 'LS454', 'ILLUMINA', 'SOLID', 
                                    'HELICOS', 'IONTORRENT', 'PACBIO'],
                        help="Platform/technology used to produce the reads.")
 
    args = parser.parse_args()
    return args
 
def addRG2Header(filename, sam_info, args):
    """Add read group info to a header."""
    # CREATE TEMPLATE
    # Read group. Unordered multiple @RG lines are allowed.
    RG_template = { 'ID': '',           # Read group identifier. e.g., Illumina flowcell + lane name and number
                    'CN': '',           # GATK Not Required. Name of sequencing center producing the read.
                    'DS': '',           # GATK Not Required. Description
                    'DT': '',           # GATK Not Required. Date the run was produced (ISO8601 date YYYY-MM-DD or YYYYMMDD)
                    'PI': '',           # GATK Not Required. Predicted median insert size.
                    'PU': '',           # GATK Not Required. Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for SOLiD).
                    'SM': '',           # Sample. Use pool name where a pool is being sequenced.
                    'PL': 'ILLUMINA'}   # Platform/technology used to produce the reads.
 
    samfile = pysam.Samfile(filename, 'r')
    new_header = samfile.header.copy()
    samfile.close()
 
    # ADD INFO TO TEMPLATE
    RG_template = RG_template.copy()
    RG_template['ID'] = sam_info["sample_name"]
    if args.CN: RG_template['CN'] = args.CN.upper()
    if args.DS: RG_template['DS'] = args.DS
    RG_template['DT'] = args.DT
    RG_template['LB'] = sam_info["sample_name"]
    RG_template['SM'] = sam_info["sample_name"]
    RG_template['DS'] = "{0}.{1}".format(sam_info['sample_name'], sam_info['locality'])
    RG_template['PI'] = args.PI
    RG_template['PU'] = '{0}.{1}'.format(sam_info['flowcell_id'], sam_info['lane'])
    new_header['RG'] = [RG_template]
    return new_header
 
 
def add_RGs_2_BAMs_runner(data):
    """Generates the correct @RG header and adds a RG field to a bam file."""
    # Make Bam of Sam if it doesn't exist
    filename, new_RG_header = data
    if filename.endswith('sam'):
        sam_handle = pysam.Samfile(filename)
        bam_name = os.path.splitext(filename)[0]+".bam"
        bam_handle = pysam.Samfile( bam_name, "wb", template = sam_handle )
        for s in sam_handle:
            bam_handle.write(s)
        filename = bam_name
 
    # Massage paths and make outputfiles
    pysam.sort(filename, os.path.splitext(filename)[0]+".sorted")
    pysam.index(os.path.splitext(filename)[0]+".sorted.bam")
    filename = os.path.splitext(filename)[0]+".sorted.bam"
    path, filename = os.path.split(filename)
    name, ext = os.path.splitext(filename)
    new_name = name + '.wRG.' + 'bam'
    outfile_name =  os.path.join(path,new_name)
    outfile = pysam.Samfile( outfile_name, 'wb', header = new_RG_header )
 
    # Step 2: Process Samfile adding Read Group to Each Read
    samfile = pysam.Samfile(os.path.join(path, filename))
    samfile.fetch()
    for count, read in enumerate(samfile.fetch()):
        name = read.qname
        read_group = os.path.split(filename)[1].split(".")[0]
        new_tags = read.tags
        new_tags.append(('RG', read_group))
        read.tags = new_tags
        outfile.write(read)
    outfile.close()
 
    # Step 3: Make index of read group enabled samfile
    pysam.index(outfile_name)
    sys.stdout.write(".")
    sys.stdout.flush()
    return
 
def parseFileName(filepath):
    path, filename = os.path.split(filepath)
    sample_name = filename.split(".")[0]
    locality = filename.split(".")[0]
    inline_tag = filename.split(".")[0]
    third_read_tag = filename.split(".")[0]
    sam_handle = pysam.Samfile(filepath,'r')
    print filepath
    #sam_line = sam_handle.next()
    #print filename
    #read_info = sam_line.qname
    #print read_info
    #sam_handle.close()
    #instrument, run_id, flowcell_id, lane = read_info.split(":")[:4]
    instrument, run_id, flowcell_id, lane = "Illumina", "1", "1", "1"
    info = {'sample_name': sample_name,
            'locality': locality,
            'inline_tag': inline_tag,
            'third_read_tag': third_read_tag,
            'instrument':instrument,  
            'run_id': run_id,
            'flowcell_id': flowcell_id,
            'lane': lane}
    return info
 
def add_RGs_2_BAMs(pool, args):
    data_for_map = []
    print 'Making RG headers.'
    for count, filename in enumerate(glob.glob(os.path.join(args.input_dir,'*'))):
        if filename.endswith('sam') or filename.endswith('bam'):
            sam_info = parseFileName(filename)
            new_RG_header = addRG2Header(filename, sam_info, args)
            data_for_map.append([filename, new_RG_header])
 
    sys.stdout.write("\nAdding RGs and making BAMs")
    sys.stdout.flush()
    pool.map(add_RGs_2_BAMs_runner, data_for_map)
    return
 
def main():
    args = get_args()
    cores = args.cores
    pool = multiprocessing.Pool(args.cores)
    add_RGs_2_BAMs(pool, args)
 
if __name__ == '__main__':
    main()