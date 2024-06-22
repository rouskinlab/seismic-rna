import multiprocessing
import os
import pickle
import time
from pathlib import Path

import fastqsplitter
import pandas as pd

from ..align.fqops import FastqUnit
from ..core.arg import (opt_barcode_end,
                        opt_barcode_start,
                        opt_parallel_demultiplexing,
                        opt_clipped_demultiplexing,
                        opt_mismatch_tolerence,
                        opt_index_tolerence,
                        opt_demulti_overwrite,
                        arg_fasta,
                        opt_mask_sections_file,
                        opt_fastqx,
                        opt_refs_meta,
                        opt_tmp_pfx,
                        opt_keep_tmp, )

# from scipy import signal

# secondary barcode names: secondary_signiture
params = [
    # Inputs
    arg_fasta,
    opt_fastqx,
    opt_mask_sections_file,
    opt_barcode_start,
    opt_barcode_end,

    # options
    opt_parallel_demultiplexing,
    opt_clipped_demultiplexing,
    opt_mismatch_tolerence,
    opt_index_tolerence,
    opt_demulti_overwrite,
    opt_refs_meta,
    opt_tmp_pfx,
    opt_keep_tmp,

]


def makes_dict_from_fastq(fpath):
    # do I make a set of the cordinates? ideally that's the only part thats different..?

    f = open(fpath)
    lines = f.readlines()
    f.close()
    reads = {}
    read_lines = []
    for x in range(len(lines)):

        l = lines[x]
        read_lines.append(l)

        if x % 4 == 0:
            spl = l.split()
            identifier_pre = spl[0]  # 5/6
            values = identifier_pre.split(":")
            identifier = (values[4] + ":" + values[5] + ":" + values[6])

        if x % 4 == 3:
            reads[identifier] = read_lines
            read_lines = []
    return reads


def reverse_compliment(sequence):
    rev_seq = sequence[::-1]
    new_str = ""
    for x in range(len(rev_seq)):
        if rev_seq[x] == "A":
            new_str += "T"
        if rev_seq[x] == "G":
            new_str += "C"
        if rev_seq[x] == "C":
            new_str += "G"
        if rev_seq[x] == "T":
            new_str += "A"
    return new_str


class Sequence_Obj():
    def __init__(self,
                 sequence,
                 name,
                 fastq1_path,
                 fastq2_path,
                 workspace,
                 paired=True,
                 fwd_primer="",
                 rev_primer="",
                 secondary_signature="", secondary_signature_start="", secondary_signature_end="",
                 rev_secondary_signature="", rev_secondary_signature_start="", rev_secondary_signature_end="",
                 barcode_start=-1, barcode_end=-1, barcode="",
                 rev_barcode="", rev_barcode_start="", rev_barcode_end=""):
        """
        This object will be hold the reference objects:
            including a dictionary of the ranges and the bv/mh associated with each of those
        """

        # input_df=pd.read_csv()
        self.demutliplex_fasta = ""
        self.sequence = sequence
        self.sample_name = name
        self.paired = paired
        self.sequence = sequence
        self.fwd_primer = fwd_primer
        self.rev_primer = rev_primer
        self.paired = paired
        self.workspace = workspace

        self.fastq1_path = fastq1_path
        self.fastq2_path = fastq2_path
        self.fastq_paths = {1: fastq1_path, 2: fastq2_path}
        # barcode info
        self.barcode_start = barcode_start
        self.barcode_end = barcode_end
        self.barcode = barcode

        self.rev_barcode = rev_barcode
        self.rev_barcode_start = rev_barcode_start
        self.rev_barcode_end = rev_barcode_end

        # secondary signiture info
        self.secondary_signature = secondary_signature
        self.secondary_signature_start = secondary_signature_start
        self.secondary_signature_end = secondary_signature_end

        self.rev_secondary_signature = rev_secondary_signature
        self.rev_secondary_signature_start = rev_secondary_signature_start
        self.rev_secondary_signature_end = rev_secondary_signature_end

        # self.sequence_folder=workspace
        self.fastqs = {1: "", 2: ""}

        self.fastq1_init_bc = set()
        self.fastq2_init_bc = set()

        clip_grepped_fq1 = {}
        clip_grepped_fq2 = {}

        # reads grepped by rev barcode search
        self.fastq1_rev_bc = set()
        self.fastq2_rev_bc = set()

        rev_clip_grepped_fq1 = {}
        rev_clip_grepped_fq2 = {}

        # ids of reads filtered out due to the secondary signiture
        self.fastq1_sec_filtered = set()
        self.fastq2_sec_filtered = set()

        # ids of reads filtered out due to reverse compliment of the secondary signiture
        self.fastq1_sec_rev_filtered = set()
        self.fastq2_sec_rev_filtered = set()

        # this is the total before secondary signiture filtering
        self.fastq1_unfiltered = set()
        self.fastq2_unfiltered = set()

        # this is post filtering or if no sec_sign filtering is going to happen
        self.fastq1_filtered = set()
        self.fastq2_filtered = set()

        # these too will only really need one since theyre the same reads either way
        # set of distnct reads from either fq1 or fq2
        self.fastq_union = set()

        # set of reads found in both fastq
        self.fastq_intersection = set()

        # organized by the amount of NT clipped from barcode and then searched for
        clip_grepped_fq1 = {}
        clip_grepped_fq2 = {}

        # grepped_done_file=


class super_fastq():

    def __init__(self, fpath: str, split_count=10, fastq_name: str = None, super_dir="") -> None:
        # print(super_dir)
        if (fastq_name == None):
            fq = fpath.split("/")[-1]
            dir_name = super_dir + fq.split("_R")[0] + "_" + fq.split("_R")[-1][0] + "_pickle_split/"
        else:
            dir_name = super_dir + fastq_name
        # dir_name+=super_dir
        os.makedirs(dir_name, exist_ok=True)
        try:

            fq_size = os.path.getsize(fpath)
        except:
            raise Exception("File could not be opened, please check if the path you gave is correct!")

        split_list = []
        pickle_file_list = []
        self.original_fastq = fpath
        self.id_to_pickle_dict = {}
        if fq_size < 100000:
            self.split_count = 1
        else:
            self.split_count = split_count

        for x in range(self.split_count):
            sp = dir_name + "split_" + str(x) + ".fastq"
            sp_pickle = dir_name + "split_" + str(x) + ".p"
            self.id_to_pickle_dict[x] = dict()

            pickle_file_list.append(sp_pickle)
            split_list.append(sp)

        self.split_list = split_list
        self.pickle_file_list = pickle_file_list

        self.id_to_pickle_dict_pickle_file = dir_name + "read_ids_per_pickle.p"

        self.split_bool = False
        self.split_fastqs_present = False
        # self.full_id_set_pickle=dir_name+"read_ids_per_pickle.p"

    def check_exists(self):

        for p in self.pickle_file_list:
            path_v = Path(p)
            b = os.path.exists(path_v)
            if (b == False):
                return b
        return True

    def check_set(self):
        return set(make_dict_from_fasta(self.original_fastq))

    def fastq_to_dict(self, fpath):
        f = open(fpath)
        lines = f.readlines()
        f.close()
        reads = {}
        read_lines = []
        for x in range(len(lines)):

            l = lines[x]
            read_lines.append(l)

            if x % 4 == 0:
                spl = l.split()
                identifier_pre = spl[0]  # 5/6
                values = identifier_pre.split(":")
                identifier = (values[4] + ":" + values[5] + ":" + values[6])

            if x % 4 == 3:
                reads[identifier] = read_lines
                read_lines = []
        return reads

    def split_fastq(self, delete_text_fastqs: bool, temp_delete_idSets_to_pickle_dict: bool):
        print("splitting: ", self.original_fastq)
        big_set = set()

        if self.check_exists():

            print("already exists")
            self.split_bool = True
            id_dict = pickle.load(open(self.id_to_pickle_dict_pickle_file, "rb"))

            for k in id_dict.keys():
                # print(k," : ",type(id_dict[k]))

                big_set = big_set.union(id_dict[k])
                # print(len(big_set))
            # print("big set: ",len(big_set))
            return big_set

        fastqsplitter.split_fastqs(self.original_fastq, self.split_list)

        for i in range(len(self.split_list)):
            # print(type(self.split_list[i]))
            temp_dict = self.fastq_to_dict(self.split_list[i])
            big_set = big_set.union(set(temp_dict))
            self.id_to_pickle_dict[i] = set(temp_dict.keys())
            pickle.dump(temp_dict, open(self.pickle_file_list[i], "wb"))

        if (delete_text_fastqs):

            for fq in self.split_list:
                os.remove(fq)
            self.split_fastqs_present = False

        else:

            self.split_fastqs_present = True
        self.split_bool = True

        """
        this set could be pretty big, could be worth temp_deleting
        """

        pickle.dump(self.id_to_pickle_dict, open(self.id_to_pickle_dict_pickle_file, "wb"))
        # self.id_to_pickle_dict.clear()

        # print("in split: ",big_set)
        return big_set

    """
    takes a dict of reads 
    this method is for efficently organizing the reads based on which pickle
    they are in so pickles only need to be opened once

    union dict a diction representing the sequence to its 

    they will write to the "directory_to_write_to" which is just a string 
    names file as construct name (key in union dict) +"_R"+fastq_id+".fastq

    also clears the union dict to save space

    """

    def super_write_fastqs(self, union_dict: dict, directoy_to_write_to: str, fastq_id: int, sequence_objects: dict,
                           sample_name: str):

        """
        organizes reads into sets per k,
        based on which pickle the read is in 
        """
        sample_fq_dir = directoy_to_write_to + sample_name + "/"
        os.makedirs(sample_fq_dir, exist_ok=True)

        organized_joint_dict = {}
        for k in union_dict.keys():
            # print(k)
            joint_set = union_dict[k]

            temp_dictionary_pickle_loc = {}

            self.id_to_pickle_dict = pickle.load(open(self.id_to_pickle_dict_pickle_file, "rb"))
            # print(self.id_to_pickle_dict.keys())

            for i in range(len(self.pickle_file_list)):
                temp_dictionary_pickle_loc[i] = set()

            for ident in joint_set:

                for num in self.id_to_pickle_dict.keys():

                    if (ident in self.id_to_pickle_dict[num]):
                        temp_dictionary_pickle_loc[num].add(ident)

            organized_joint_dict[k] = temp_dictionary_pickle_loc

        # union_dict.clear()
        """
        now that the reads are organized we have to open each pickle file write fastqs 
        """
        # fqs=[]
        for k in sequence_objects.keys():

            if (len(union_dict[k]) > 1):
                fq = sample_fq_dir + k + "_R" + str(fastq_id) + ".fastq"
                fq_buff = open(fq, "wt")
                fq_buff.close()

        for x in range(len(self.pickle_file_list)):

            p_file_object = open(self.pickle_file_list[x], "rb")
            pickled_fastq_dict = pickle.load(p_file_object)
            p_file_object.close()

            for k in organized_joint_dict.keys():

                if (len(union_dict[k]) > 1):

                    # print("writing filtered fastqs"+str(k))
                    fq = sample_fq_dir + k + "_R" + str(fastq_id) + ".fastq"
                    # fqs.append(fq)

                    sequence_objects[k].fastqs[fastq_id] = fq

                    specific_reads = organized_joint_dict[k][x]

                    fq_buff = open(fq, "a")

                    for r in specific_reads:
                        lines = pickled_fastq_dict[r]
                        fq_buff.writelines(lines)
                    fq_buff.close()
        # fqs.sort()

    """
    this methods removes the pickles that contain split fastq information
    """

    def destroy_temp_data(self):
        for x in range(len(self.pickle_file_list)):
            os.remove(self.pickle_file_list[x])


"""
tolerence is added to both ends of search range

if tolerence==-1 then disregards indexes 
    this will increase computation time
no index ranges and high mismatch threshhold will result in extremely high comp time 
"""


def run_seqkit_grep_function(pattern: str,
                             search_start_ind: int,
                             search_end_index: int,
                             fastq_to_search: str,
                             fastq_to_write: str,
                             threads: int = 20,
                             mismatch_threshhold: int = 0,
                             append_bool: bool = False,
                             tolerance: int = 0,
                             delete_fq: bool = False
                             ):
    """
    1 indexed? 
    
    """
    append_char = ">>" if append_bool else ">"

    if (search_start_ind - tolerance < 0):
        search_start_ind = 0
    else:
        search_start_ind -= tolerance
    search_end_index += tolerance
    if (search_start_ind < 1):
        search_end_index += abs(search_start_ind) + 1
        search_start_ind = 1
    cmd = f'seqkit --threads {threads} grep -s -R {search_start_ind}:{search_end_index + tolerance} -p "{pattern}" -m {mismatch_threshhold} -P {fastq_to_search} > {fastq_to_write}'
    print(cmd)  # debug-
    return_code = os.system(cmd)  # debug
    if (delete_fq):
        os.remove(fastq_to_write)
    return (set(makes_dict_from_fastq(fastq_to_write).keys()))


"""
simple fastq dict for retrieveing sequences from fastq
"""


def make_dict_from_fasta(fasta_path) -> dict:
    fa = open(fasta_path, "rt").readlines()
    temp_dict = {}

    for i in range(0, len(fa), 2):
        temp_dict[fa[i][1:].strip()] = fa[i + 1].strip()

    return temp_dict


"""
input csv, mask_sections_file that represents each sequence to be dumultiplexed with many different coloumns 

workspace directory that demultiplexing is being done in 

fq1/fq2 path for big fastq which will be demultiplexed? maybe should be removed


this whole method could be replaced with a dataframe that organizes all of these attributes


"""


def make_sequence_objects_from_csv(input_csv, barcode_start, barcode_end, fasta, fastq1_path, fastq2_path, paired,
                                   workspace) -> dict:
    
    sequence_object_dict = {}
    fasta_dict = make_dict_from_fasta(fasta)
    if (input_csv == ""):
        print("rip")
        for name in fasta_dict.keys():
            seq = fasta_dict[name]
            rev_seq = reverse_compliment(seq)
            bc = seq[barcode_start:barcode_end]

            rev_barcode = reverse_compliment(bc)
            rev_bc_start = rev_seq.index(rev_barcode)
            rev_bc_end = rev_bc_start + len(rev_barcode)

            sequence_object_dict[name] = Sequence_Obj(
                sequence=fasta_dict[name],
                fastq1_path=fastq1_path,
                fastq2_path=fastq2_path,
                name=name,
                paired=paired,
                barcode_start=barcode_start,
                barcode_end=barcode_end,
                barcode=bc,
                rev_barcode=rev_barcode,
                rev_barcode_start=rev_bc_start,
                rev_barcode_end=rev_bc_end,
                secondary_signature_start=-1,
                secondary_signature_end=-1,
                secondary_signature=-1,
                rev_secondary_signature_start=-1,
                rev_secondary_signature_end=-1,
                rev_secondary_signature=-1,
                workspace=workspace + name + "/"
            )

    else:
        print(f"input csv: {input_csv}")
        df = pd.read_csv(input_csv)

        # fasta_dict=make_dict_from_fasta(fasta)

        # sequence_object_dict={}
        cols = set(df.columns)

        if ("Barcode5" not in cols and barcode_start==-1):
            raise Exception("no barcode info given")

        for x in df.index:

            name = df.at[x, "Reference"]
            seq = fasta_dict[name]
            rev_seq = reverse_compliment(seq)

            barcode_start = df.at[x, "Barcode5"]
            barcode_end = df.at[x, "Barcode3"]
            bc = seq[barcode_start:barcode_end]
            rev_barcode = reverse_compliment(bc)
            rev_bc_start = rev_seq.index(rev_barcode)
            rev_bc_end = rev_bc_start + len(rev_barcode)

            if ("Secondary_Signature_Start" in cols):
                secondary_sign_start = df.at[x, "Secondary_Signature_Start"]
                secondary_sign_end = secondary_sign_start + df.at[x, "Secondary_Signature_Length"]
                secondary_sign = fasta_dict[name][secondary_sign_start:secondary_sign_end]

                rev_sec_sign = reverse_compliment(secondary_sign)
                rev_sec_sign_start = rev_seq.index(rev_sec_sign)
                rev_sec_sign_end = rev_sec_sign_start + len(rev_sec_sign)


            else:
                secondary_sign_start = -1
                secondary_sign_end = -1
                secondary_sign = -1

                rev_sec_sign = -1
                rev_sec_sign_start = -1
                rev_sec_sign_end = -1

            sequence_object_dict[name] = Sequence_Obj(
                sequence=fasta_dict[name],
                fastq1_path=fastq1_path,
                fastq2_path=fastq2_path,
                name=name,
                paired=paired,
                barcode_start=barcode_start,
                barcode_end=barcode_end,
                barcode=bc,
                rev_barcode=rev_barcode,
                rev_barcode_start=rev_bc_start,
                rev_barcode_end=rev_bc_end,
                secondary_signature_start=secondary_sign_start,
                secondary_signature_end=secondary_sign_end,
                secondary_signature=secondary_sign,
                rev_secondary_signature_start=rev_sec_sign_start,
                rev_secondary_signature_end=rev_sec_sign_end,
                rev_secondary_signature=rev_sec_sign,
                workspace=workspace + name + "/"
            )
    return sequence_object_dict


"""
seqkit grep is 1 indexed and inclusive of the final value of its range 



clipped int that represents how much 
fastq_id fastq 1 or fastq 2
"""

"""
runs grep and accepts a clipped argument and appends the set to the main dictionary

index tolerence can only apply to the initial but there are cases where that could false 
"""


def run_multi_greps(read_id_dict: dict, clipped: int, index_tolerence: int, delete_fastqs: bool,
                    mismatches_allowed: int, pattern_type: str, pattern: str, pattern_start: int, pattern_end: int,
                    fastq: str, seq_folder: str, front: bool = False):
    # debug#print("running multi_greps for ", pattern_type)
    new_reads_dict = {}
    folder = seq_folder
    # file_name_start=
    new_reads_dict["init_" + pattern_type] = run_seqkit_grep_function(pattern=pattern, search_start_ind=pattern_start,
                                                                      search_end_index=pattern_end,
                                                                      fastq_to_search=fastq,
                                                                      fastq_to_write=seq_folder + "init_" + pattern_type + ".fastq",
                                                                      mismatch_threshhold=mismatches_allowed,
                                                                      tolerance=index_tolerence)
    # if(len(pattern)-pattern_end)<clipped:
    append_these = [seq_folder + "init_" + pattern_type + ".fastq"]
    clipped += 1
    # print("initial grep run")#debug
    if not front:
        for x in range(1, clipped):
            append_these.append(seq_folder + pattern_type + "_clipped_" + str(x) + ".fastq")
            new_reads_dict[pattern_type + "_clipped_" + str(x)] = run_seqkit_grep_function(
                pattern=pattern[:len(pattern) - x], search_start_ind=pattern_start + x, search_end_index=pattern_end,
                fastq_to_search=fastq, fastq_to_write=seq_folder + pattern_type + "_clipped_" + str(x) + ".fastq",
                mismatch_threshhold=0, tolerance=0)
    else:
        for x in range(1, clipped):
            append_these.append(seq_folder + pattern_type + "_clipped_" + str(x) + ".fastq")
            new_reads_dict[pattern_type + "_clipped_" + str(x)] = run_seqkit_grep_function(pattern=pattern[x:],
                                                                                           search_start_ind=pattern_start,
                                                                                           search_end_index=pattern_end - x,
                                                                                           fastq_to_search=fastq,
                                                                                           fastq_to_write=seq_folder + pattern_type + "_clipped_" + str(
                                                                                               x) + ".fastq",
                                                                                           mismatch_threshhold=0,
                                                                                           tolerance=0)
    # deug#print("clipped run")
    # for k in new_reads_dict.keys():
    # keys=list(new_reads_dict.keys())

    for x in range(1, clipped):

        new_reads_dict[pattern_type + "_clipped_" + str(x)] -= new_reads_dict["init_" + pattern_type]
        for i in range(x + 1, clipped):
            new_reads_dict[pattern_type + "_clipped_" + str(x)] -= new_reads_dict[pattern_type + "_clipped_" + str(x)]
    # debug#print("organizing reads")
    for k in new_reads_dict.keys():
        read_id_dict[k] = new_reads_dict[k]
    # debug#print("appending reads")
    cmd = "cat "
    return_fastq = seq_folder + pattern_type + "_appended.fastq"
    # print(append_these  )
    for x in range(len(append_these)):
        cmd += (append_these[x] + " ")
    cmd += "> " + seq_folder + pattern_type + "_appended.fastq"

    os.system(cmd)
    # debug#print("appended!")
    return return_fastq


def append_files(files, new_file_name):
    cmd = "cat "
    for x in range(len(files)):
        cmd += (files[x] + " ")
    cmd += "> " + new_file_name
    os.system(cmd)


def run_seqkit_grep(sequence_object: Sequence_Obj, clipped: int, rev_clipped: int, index_tolerence: int,
                    delete_fastqs: bool, fastq_id: int, mismatches_allowed: int):
    # f=open(str(sequence_object.sample_name)+"_test.txt","wt")
    # f.write(str(vars(sequence_object)))
    fastq = sequence_object.fastq_paths[fastq_id]
    # makes folder for demultiplex
    os.makedirs(sequence_object.workspace, exist_ok=True)

    fastq_folder = sequence_object.workspace + "fq" + str(fastq_id) + "/"
    os.makedirs(fastq_folder, exist_ok=True)

    # fastqs

    fastq_unfiltered = fastq_folder + "unfiltered.fastq"

    fastq_sec_filtered = fastq_folder + "sec_fitlered_filtered.fastq"
    fastq_rev_sec_filtered = fastq_folder + "rev_sec_filtered.fastq"

    fastq_filtered = fastq_folder + "filtered_out.fastq"

    read_ids = {}
    bc_fastq = run_multi_greps(read_id_dict=read_ids, clipped=clipped, index_tolerence=index_tolerence,
                               delete_fastqs=delete_fastqs, mismatches_allowed=mismatches_allowed,
                               pattern_type="barcode", pattern=sequence_object.barcode,
                               pattern_start=sequence_object.barcode_start, pattern_end=sequence_object.barcode_end,
                               fastq=fastq, seq_folder=fastq_folder)
    # debug#print("first_multi_run")
    rev_bc_fastq = run_multi_greps(read_id_dict=read_ids, clipped=0, index_tolerence=index_tolerence,
                                   delete_fastqs=delete_fastqs, mismatches_allowed=mismatches_allowed,
                                   pattern_type="rev_barcode", pattern=sequence_object.rev_barcode,
                                   pattern_start=sequence_object.rev_barcode_start,
                                   pattern_end=sequence_object.rev_barcode_end, fastq=fastq, seq_folder=fastq_folder)

    append_files([bc_fastq, rev_bc_fastq], fastq_unfiltered)
    read_ids["unfiltered"] = set(makes_dict_from_fastq(fastq_unfiltered).keys())

    complete_set = set()

    if (sequence_object.secondary_signature != -1):
        v_threshold = 0
        v_section = pattern = sequence_object.secondary_signature
        if (len(v_section) > 20):
            v_threshold = 4
        elif (len(v_section) <= 20 and len(v_section) > 15):
            v_threshold = 3
        elif (len(v_section) <= 15 and len(v_section) >= 5):
            v_threshold = 2
        elif (len(v_section) < 5):
            v_threshold = 1

        secondary_sign_mismatches = v_threshold  # TODO hardcoded :(

        read_ids["sec_filter"] = run_seqkit_grep_function(pattern=sequence_object.secondary_signature,
                                                          search_start_ind=sequence_object.secondary_signature_start,
                                                          search_end_index=sequence_object.secondary_signature_end,
                                                          fastq_to_search=fastq_unfiltered,
                                                          fastq_to_write=fastq_sec_filtered,
                                                          mismatch_threshhold=secondary_sign_mismatches, tolerance=2)
        read_ids["rev_sec_filter"] = run_seqkit_grep_function(pattern=sequence_object.rev_secondary_signature,
                                                              search_start_ind=sequence_object.rev_secondary_signature_start,
                                                              search_end_index=sequence_object.rev_secondary_signature_end,
                                                              fastq_to_search=fastq_unfiltered,
                                                              fastq_to_write=fastq_rev_sec_filtered,
                                                              mismatch_threshhold=secondary_sign_mismatches,
                                                              tolerance=2)
        complete_set = read_ids["sec_filter"].union(read_ids["rev_sec_filter"])
        read_ids["reads_lost_to_filter"] = read_ids["unfiltered"] - complete_set
        # complete_set=read_ids["sec_filter"].union(read_ids["rev_sec_filter"])
    else:
        complete_set = read_ids["unfiltered"]
    # here the sets need to be compared in order to identify which ids were unique to which search
    # TODO

    # print(vars(sequence_object))
    # complete_set=read_ids["sec_filter"].union(read_ids["rev_sec_filter"])

    complete_set_pickle = fastq_folder + "complete_set_of_reads.p"
    pickle.dump(complete_set, open(complete_set_pickle, "wb"))
    read_id_data_pickle = fastq_folder + "read_id_data.p"
    pickle.dump(read_ids, open(read_id_data_pickle, "wb"))
    # print(read_ids)

    open(sequence_object.workspace + "grepped.txt", "wt").write("done")


"""
latest and greatest
"""


def check_done(sequence_folder: str) -> bool:
    # print("looking here: ",sequence_folder+"grepped.txt")

    comp_set_fq1 = "fq1/complete_set_of_reads.p"

    if (os.path.exists(sequence_folder + "fq1/complete_set_of_reads.p")) \
            and (os.path.exists(sequence_folder + "fq2/complete_set_of_reads.p")) \
 \
            and (os.path.exists(sequence_folder + "fq1/read_id_data.p")) \
            and (os.path.exists(sequence_folder + "fq2/read_id_data.p")) \
 \
            and (os.path.exists(sequence_folder + "grepped.txt")):

        return True
    else:
        return False


def check_all_done(seq_objects: dict()):
    return_dict = {}

    for k in seq_objects.keys():
        print("checking: " + k, end="\r")
        if not (check_done(seq_objects[k].workspace)):
            print("check done not done! ", check_done(seq_objects[k].workspace))
            return_dict[k] = seq_objects[k]
    return return_dict


def grep_both_fastq(sequence_object: Sequence_Obj, clipped: int, rev_clipped: int, index_tolerence: int,
                    delete_fastqs: bool, mismatches_allowed: int):
    run_seqkit_grep(sequence_object=sequence_object, clipped=clipped, rev_clipped=rev_clipped,
                    index_tolerence=index_tolerence, delete_fastqs=delete_fastqs, mismatches_allowed=mismatches_allowed,
                    fastq_id=1)
    run_seqkit_grep(sequence_object=sequence_object, clipped=clipped, rev_clipped=rev_clipped,
                    index_tolerence=index_tolerence, delete_fastqs=delete_fastqs, mismatches_allowed=mismatches_allowed,
                    fastq_id=2)


def parallel_grepping(sequence_objects: dict, fwd_clips: int, rev_clips: int, index_tolerence: int, delete_fastq: bool,
                      paired: bool = True, mismatches: int = 0, threads=10, iteration: int = 0, overwrite: bool = True):
    """
    runs grep in parallel 
    """
    itr_val = iteration
    print("iteration value:XXX ", itr_val)
    countx = 0
    procs = []
    seq_keys = list(sequence_objects.keys())
    seq_count = 0
    THREADS = threads
    # print(seq_keys)
    run_all = True

    while (seq_count < len(seq_keys) + 6):

        print(seq_count)
        procs = []

        seq_index = seq_count
        print(range(seq_index, seq_index + THREADS))
        for i in range(seq_index, (seq_index + THREADS)):
            # seq_index+i
            # print(seq_index)
            if (i < len(seq_keys)):
                seq_keys[i]
                if (not overwrite):
                    previously_run = check_done(sequence_objects[seq_keys[i]].workspace)
                else:
                    previously_run = False
                # print(previously_run)
                if (previously_run and not overwrite):
                    print("grepped")
                else:
                    # (sequence_object:Sequence_Obj,clipped:int,rev_clipped:int,index_tolerence:int,delete_fastqs:bool,mismatches_allowed:int)
                    x = multiprocessing.Process(target=grep_both_fastq, args=(
                        sequence_objects[seq_keys[i]], fwd_clips, rev_clips, index_tolerence, delete_fastq, mismatches))
                    x.start()
                    procs.append(x)
            # seq_index=+1
            else:
                # info#print("no more to run")
                pass

        for p in procs:
            p.join()

        # time.sleep(1000)
        start = time.time()
        seq_count += THREADS
        itr_val += 1

    print("checking")
    not_done_dict = check_all_done(sequence_objects)
    # print(not_done_list)

    while len(not_done_dict) > 0 and itr_val > 4:
        print("rip")
        itr_val = iteration
        # itr_val=iteration
        # not_done_dict=check_all_done(sequence_objects)
        parallel_grepping(sequence_objects=not_done_dict, fwd_clips=fwd_clips, rev_clips=rev_clips,
                          index_tolerence=index_tolerence, delete_fastq=delete_fastq, iteration=itr_val + 1,
                          mismatches=mismatches)
    if (itr_val > 4):
        print("could not finish some of the reads: ", list(not_done_dict.keys()))

        # warning could not finish one of the reads
    return not_done_dict


def regular_grepping(sequence_objects: dict, fwd_clips: int, rev_clips: int, index_tolerence: int, delete_fastq: bool,
                     paired: bool = True, mismatches: int = 0, iteration: int = 0, overwrite: bool = False):
    """
    runs grep in parallel 
    """
    print("regular grepping")
    itr_val = iteration

    for k in sequence_objects.keys():
        previously_run = check_done(sequence_objects[k].workspace)
        if (previously_run and (not overwrite)):
            print("grepped")
        else:
            grep_both_fastq(sequence_objects[k], fwd_clips, rev_clips, index_tolerence, delete_fastq, mismatches)

    not_done_dict = check_all_done(sequence_objects)

    while len(not_done_dict) > 0 and itr_val > 4:
        itr_val = iteration
        """if itr_val>0:
            not_done_dict=check_all_done(sequence_objects)"""

        # not_done_dict=check_all_done(sequence_objects)
        regular_grepping(sequence_objects=not_done_dict, fwd_clips=0, rev_clips=0, index_tolerence=0,
                         delete_fastq=False, iteration=itr_val + 1, mismatches=mismatches)
    if (itr_val > 4):
        print("could not finish one of the reads: ", list(not_done_dict.keys()))


"""
checks each sequence for a grepped.txt and returns true if found 
"""


def finds_multigrepped_reads(sequence_objects: dict, remove: bool = True, resolve: bool = False,
                             print_multi_grep_dict: bool = True, demultiplex_workspace: str = None) -> dict:
    """
    filters reads based on weather or not they map to multiple constructs 
    returns a dictionary mapping each read to a list of the reads it mapped to
    """
    union_dictionary = {}
    fq1_sets = {}
    fq2_sets = {}

    for k in sequence_objects.keys():
        # print(k)
        seq_object_fastq_folder = sequence_objects[k].workspace
        fq1_complete_set_pickle = seq_object_fastq_folder + "fq1/complete_set_of_reads.p"
        fq2_complete_set_pickle = seq_object_fastq_folder + "fq2/complete_set_of_reads.p"

        fq1_complete_set = pickle.load(open(fq1_complete_set_pickle, "rb"))
        fq2_complete_set = pickle.load(open(fq2_complete_set_pickle, "rb"))

        fq1_sets[k] = fq1_complete_set
        fq2_sets[k] = fq2_complete_set
        # joint set of ids from complete sets of fqs
        union_dictionary[k] = fq1_complete_set.union(fq2_complete_set)
        """
        now that all reads are collected into a big dictionary of sets 
        we can see which ones are mapped to multiple reads
        """
    used_reads = {}

    for k in sequence_objects.keys():

        for read in union_dictionary[k]:

            if read not in used_reads.keys():
                used_reads[read] = [k]
            else:
                used_reads[read].append(k)

    multi_grep_dict = {}

    for read in used_reads.keys():
        used_by_list = used_reads[read]

        if len(used_by_list) > 1:
            temp_list = []
            # multi_grep_dict[read]=used_by_list

            for name in used_by_list:

                if (read in fq1_sets[name] and read in fq2_sets[name]):
                    temp_list.append([0, 0, 1])
                elif (read in fq1_sets[name] and read not in fq2_sets[name]):
                    temp_list.append([1, 0, 0])
                elif (read not in fq1_sets[name] and read in fq2_sets[name]):
                    temp_list.append([0, 1, 0])
                else:
                    temp_list.append([0, 0, 0])

                if (remove):
                    union_dictionary[name].remove(read)
            multi_grep_dict[read] = [used_by_list, temp_list]

    if (print_multi_grep_dict):
        pickle.dump(multi_grep_dict, open(demultiplex_workspace + "multigrepped_reads.p", "wb"))

    return union_dictionary


def resolve_or_analyze_multigrepped_reads(union_sets: dict, remove: bool = True, resolve: bool = False):
    pass


def create_report(sequence_objects: dict, fq1: str, fq2: str, working_directory: str, unioned_sets: dict,sample_name:str):
    orginal_len = {}
    mixed_total_dict = {}

    fastqs = {}
    dict_of_lists = {}
    FQS = ["fq1", "fq2"]
    # basic_info=["construct_name_","final_count_","lost_to_multigrep_",]
    temp_read_data_pickle = sequence_objects[list(sequence_objects.keys())[0]].workspace + "fq1/read_id_data.p"
    read_data_keys = list(pickle.load(open(temp_read_data_pickle, "rb")).keys())

    for fq in FQS:
        if (fq == "fq1"):
            fastqs[fq] = fq1
        else:
            fastqs[fq] = fq2

        mixed_total_dict[fq] = set(makes_dict_from_fastq(fastqs[fq]))
        orginal_len[fq] = len(mixed_total_dict[fq])
        print("orginal len: ", orginal_len[fq])
        dict_of_lists[fq + "_construct_name_"] = []
        dict_of_lists[fq + "_final_count_"] = []
        dict_of_lists[fq + "_lost_to_multigrep_"] = []
        for data_key in read_data_keys:
            dict_of_lists[fq + "_" + data_key] = []

    for k in sequence_objects.keys():
        # sequence_objects
        seq_obj = sequence_objects[k]
        final_reads = unioned_sets[k]

        print("final_counts: ", len(final_reads))

        for fq in FQS:
            mixed_total_dict[fq] -= final_reads
            fq_complete_set = pickle.load(open(seq_obj.workspace + fq + "/complete_set_of_reads.p", "rb"))
            fq_read_data = pickle.load(open(seq_obj.workspace + fq + "/read_id_data.p", "rb"))

            dict_of_lists[fq + "_construct_name_"].append(k)
            dict_of_lists[fq + "_final_count_"].append(len(final_reads))
            dict_of_lists[fq + "_lost_to_multigrep_"].append(len(fq_complete_set) - len(final_reads))

            for data_key in read_data_keys:
                dict_of_lists[fq + "_" + data_key].append(len(fq_read_data[data_key]))
    df = pd.DataFrame()

    for k in dict_of_lists.keys():
        df[k] = dict_of_lists[k]
    print("len: mixed: ", len(mixed_total_dict[FQS[0]]))
    print("len: mixed: ", len(mixed_total_dict[FQS[0]]), " / ", orginal_len[FQS[0]], " = ",
          len(mixed_total_dict[FQS[0]]) / orginal_len[FQS[0]])
    len_col = [" "] * len(sequence_objects.keys())
    len_col[0] = len(mixed_total_dict[FQS[0]])
    df["percent_reads_used"] = [(1 - (len(mixed_total_dict[FQS[0]]) / orginal_len[FQS[0]])) * 100] * len(
        sequence_objects.keys())
    df["total_read_count"] = len_col

    print(working_directory + "demultiplex_info.csv")
    df.to_csv(working_directory + sample_name+"_demultiplex_info.csv", index=False)


"""
split is default to 10. disregarding extremes, the higher the split the lighter the memeory load
mask_sections_file csv 
    each construct must have a secondary signiture start index and len in order to process, 
    barcode given in main arguements 
"""


"""



"""

def demultiplex_run(refs_file_csv, demulti_workspace, report_folder, fq_unit: FastqUnit, fasta, barcode_start=0,
                    barcode_end=0, split: int = 10, clipped: int = 0, rev_clipped: int = 0, index_tolerance: int = 0,
                    parallel: bool = False, mismatch_tolerence: int = 0, overwrite: bool = False, keep_tmp:bool=True):
    sample_name = fq_unit.sample
    mixed_fastq1, mixed_fastq2 = (fq_unit.paths.values())  # only works if the FASTQ has paired-end reads in two separate files
    mixed_fastq1=str(mixed_fastq1)
    mixed_fastq2=str(mixed_fastq2)
    print(f"csv: {refs_file_csv}")
    report_folder+="/"
    """
    makes dictionary of sequence objects
    """
    temp_ws = report_folder + "/" + sample_name + "_demultiplex_folders_and_files/"
    # final_sample_folder=temp_ws+"sample_fqs/"

    # print(temp_ws)

    print(repr(temp_ws))
    os.makedirs(temp_ws, exist_ok=True)
    # print("tempworkspace: ",temp_ws)
    """
    all the little stuff gets stored per sequence here, at least temporarily 
    """
    seq_data_folder = temp_ws + "sequence_data/"

    os.makedirs(seq_data_folder, exist_ok=True)
    sequence_objects = make_sequence_objects_from_csv(
        input_csv=refs_file_csv,
        barcode_start=barcode_start,
        barcode_end=barcode_end,
        fasta=fasta,
        fastq1_path=mixed_fastq1,
        fastq2_path=mixed_fastq2,
        paired=True,
        workspace=seq_data_folder)

    # print("\nhere:\n"+(sequence_objects))
    # print("workspace: ",vars(sequence_objects["3042-O-flank_1=hp1-DB"]))
    # demultiplex_workspace=demulti_workspace#"demultiplexed_sequences/"
    """
    makes a super fastq for memory efficent access to fastq reads

    checks values of ids to verify that the two fastq provided contain matching paired end reads
    """
    super_fq1 = super_fastq(mixed_fastq1, split, super_dir=temp_ws)
    super_set1 = super_fq1.split_fastq(True, False)

    super_fq2 = super_fastq(mixed_fastq2, split, super_dir=temp_ws)
    super_set2 = super_fq2.split_fastq(True, False)
    # print(super_set2)
    c = len(super_set1.difference(super_set2))

    super_set1.clear()
    super_set2.clear()

    if c != 0:
        raise Exception("Fastq ids do not match, please verify that the file")

    """
    runs grep in parallel 
    """
    if (parallel):
        parallel_grepping(sequence_objects=sequence_objects, fwd_clips=clipped, rev_clips=rev_clipped,
                          index_tolerence=index_tolerance, delete_fastq=False, mismatches=mismatch_tolerence,
                          overwrite=overwrite)
    else:
        regular_grepping(sequence_objects=sequence_objects, fwd_clips=clipped, rev_clips=rev_clipped,
                         index_tolerence=index_tolerance, delete_fastq=False, mismatches=mismatch_tolerence,
                         overwrite=overwrite)
    # print(cum+plus)

    """
    filters reads based on weather or not they map to multiple constructs 
    default is to delete multigrepped reads

    """
    unioned_sets_dictionary = finds_multigrepped_reads(sequence_objects=sequence_objects,
                                                       demultiplex_workspace=seq_data_folder)
    # debug#print(unioned_sets_dictionary.keys())

    # multigrepped_read_analysis()

    """
    organizes reads for writing
    uses super_fastq object to write fastqs in efficent way
    """

    fq1_paths = super_fq1.super_write_fastqs(unioned_sets_dictionary, report_folder, 1, sequence_objects,
                                             sample_name=sample_name)
    # super_fq1.destroy_temp_data()

    fq2_paths = super_fq2.super_write_fastqs(unioned_sets_dictionary, report_folder, 2, sequence_objects,
                                             sample_name=sample_name)
    # super_fq2.destroy_temp_data()

    """
    makes report on what amount of reads were found in which stage
    """
    # sequence_objects:dict,fq1:str,fq2:str,working_directory:str,unioned_sets:dict)
    print("creating report!!!")
    create_report(sequence_objects, mixed_fastq1, mixed_fastq2, report_folder, unioned_sets_dictionary,sample_name)

    return (), (), (report_folder + sample_name + "/",)

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
