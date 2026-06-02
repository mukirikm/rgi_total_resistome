import re
import warnings

from Bio import BiopythonWarning
from Bio.Seq import Seq
import traceback

from app.Base import BaseModel
from app.settings import *

warnings.simplefilter("ignore", BiopythonWarning)

class MutationsModule(BaseModel):
    """
    Class for mutation searches (SNVs, frameshifts, coSNPs, etc.).
    """

    def __init__ (self):
       pass
			
    def __repr__(self):
        """
        Returns Mutation class full object.
        """         
        return "Mutation({}".format(self.__dict__)

    def single_resistance_variant(self, detection_mode,  snp_dict_list, hsp_query, 
                                  hsp_sbjct_start, hsp_sbjct, orf_info, query_def,
                                  pred_genes_dict_prot=None, sub_prot_dict=None, hsp_query_start=None, 
                                  hsp_query_end=None, hsp_sbjct_end=None, real_qry_length=None, real_sbjct_length=None, strand=None) : 
        """
        Searches for SNVs in sequences.
        """

        for eachs in snp_dict_list:
            srv_output = {}
            srv_output["query_def"] = query_def

            pos = eachs["position"]
            ori = eachs["original"]
            chan = eachs["change"]
            
            if detection_mode == "PVM":
                if hsp_sbjct_start < pos and (hsp_sbjct_start + real_sbjct_length) > pos:
                    orf_protein_sequence = ""

                    if pred_genes_dict_prot:
                        if orf_info.strip() in pred_genes_dict_prot.keys():
                            orf_protein_sequence = pred_genes_dict_prot[orf_info.decode()].strip("*")
                            srv_output["eachs"] = eachs
                            srv_output["orf_protein_sequence"] = orf_protein_sequence
                            srv_output["chan"] = chan
                        else:
                            orf_protein_sequence = pred_genes_dict_prot[orf_info.decode()[:orf_info.decode().index(' # ')]].strip("*")
                            srv_output["eachs"] = eachs
                            srv_output["orf_protein_sequence"] = orf_protein_sequence
                            srv_output["chan"] = chan

                    if sub_prot_dict:
                        orf_protein_sequence = str(sub_prot_dict[orf_info.decode().split(" ")[0]])
                        srv_output["eachs"] = eachs
                        srv_output["orf_protein_sequence"] = orf_protein_sequence
                        srv_output["chan"] = chan 

                    # wildtype
                    wildtype = str(
                        hsp_sbjct[pos - hsp_sbjct_start + self.find_num_dash(hsp_sbjct, (pos-hsp_sbjct_start))])
                    srv_output["wildtype"] = wildtype

                    # Report ONLY if the SNPs are present
                    qry = int(
                        pos) - hsp_sbjct_start + self.find_num_dash(hsp_sbjct, (int(pos) - hsp_sbjct_start))
                    srv_output["qry"] = qry

                    # check for Var2
                    if str(chan) == "Var":
                        # update to the change
                        chan = str(
                            hsp_query[pos - hsp_sbjct_start + self.find_num_dash(hsp_sbjct, (pos-hsp_sbjct_start))])
                        
                        if hsp_query[qry] == chan and chan != wildtype:
                            # update eachs@change
                            eachs["change"] = chan
                        else:
                            # change same as wildtype, don't report
                            chan = ""
                            # eachs["change"] = chan

                    if hsp_query[qry] == chan: # if the amino acid at our specific position in the query sequence is the same as the NEW aa (same SNP change has occured)
                        query_snps = {}

                        # get position of mutation in the query sequence
                        d = int(
                            pos) - hsp_sbjct_start - self.find_num_dash(hsp_query, (int(pos) - hsp_sbjct_start))
                        query_snps = {
                            "original": ori, "change": chan ,"position": d+1}
                        # logger.debug("query_snp on frame {} {}".format(hsp.frame, json.dumps(query_snps, indent=2)))

                        srv_output["query_snps"] = query_snps   
                        srv_output["has_snp"] = True             
                        yield srv_output
                    else:
                        srv_output = None
                        srv_output = {"query_def": query_def,
                                      "has_snp": False}
                        yield srv_output

            if detection_mode == "POM":
                if hsp_sbjct_start < int(pos) and (hsp_sbjct_start + real_qry_length) > int(pos):
                    """Checks if there is a mutation."""
                    # logger.debug("Mutation check")
                    qry = int(
                        pos) - hsp_sbjct_start + self.find_num_dash(hsp_sbjct, (int(pos) - hsp_sbjct_start))
                    sbj = int(
                        pos) - hsp_sbjct_start + self.find_num_dash(hsp_sbjct, (int(pos) - hsp_sbjct_start))
                    
                    # bounds check before query indexing (making sure that for hsp_quer[qry], qry is neither negative nor greater than the actual length of the sequence)
                    if not (0 <= qry < len(hsp_query)):
                        srv_output["has_snp"] = False
                        yield srv_output

                    elif hsp_query[qry] == chan:
                        # logger.debug("Mutation detected")
                        srv_output["eachs"] = eachs
                        srv_output["has_snp"] = True
                        yield srv_output
                    else:
                        srv_output["has_snp"] = False
                        yield srv_output

            if detection_mode == "RGV":
                srv_output["eachs"] = eachs
                if hsp_query_start < pos and (hsp_query_start + real_qry_length) > pos:
                    # Report ONLY if the SNPs are present
                    qry = int(
                        pos) - hsp_query_start + self.find_num_dash(hsp_query, (int(pos) - hsp_query_start))
                    sbj = int(
                        pos) - hsp_query_start + self.find_num_dash(hsp_query, (int(pos) - hsp_query_start))

                    if hsp_sbjct[sbj].lower() == chan.lower():
                        query_snps = {}
                        logger.info(
                            "hsp.query_start: {}".format(hsp_query_start))
                        logger.info(
                            "hsp.query_end: {}".format(hsp_query_end))
                        logger.info(
                            "hsp.sbjct_start: {}".format(hsp_sbjct_start))
                        logger.info(
                            "hsp.sbjct_end: {}".format(hsp_sbjct_end))
                        d = 0
                        if strand == "+":
                            d = int(
                                pos) - hsp_query_start - self.find_num_dash(hsp_sbjct, (int(pos) - hsp_query_start))
                            # d = -1*(hsp.query_start - pos)
                            query_snps = {
                                "original": hsp_query[d], "change": hsp_sbjct[d], "position": (d + 1)}
                            srv_output["query_snps"] = query_snps
                            yield srv_output
                        else:
                            d = int(
                                pos) - hsp_query_start - self.find_num_dash(hsp_sbjct, (int(pos) - hsp_query_start))
                            # d = -1*(hsp.query_start - pos)
                            query_snps = {
                                "original": hsp_query[d], "change": hsp_sbjct[d], "position": (d + 1)}
                            srv_output["query_snps"] = query_snps
                            yield srv_output

                        logger.info(
                            "position in the query (3'->5') : {}".format(d))

                        # logger.info("query_snp on frame {} {}".format(hsp.frame, json.dumps(query_snps, indent=2)))

    def frameshift(self, hsp_query, hsp_sbjct, card_dna_ref, query_def, param_type, fs_dict_list=[]): 
        """
        Searches for frameshifts in sequences.
        """
        
        fs_result_prelim = {}
    
        fs_curated_list_reg = []
        fs_denovo_list_reg = []

        fs_curated_list_validation = []
        fs_denovo_list_validation = []

        fs_curated_result_HGVS = []
        fs_denovo_result_HGVS = []

        # for deletions
        qry_codon_pos = 0

        # for insertions
        sbjct_codon_pos = 0
                        
        split_ref = re.findall('.'*3, card_dna_ref)
        
        """for nucleotide deletions """            
        ### isolate the position of the gap(s) and the affected codon(s)
        if "-" in hsp_query:
            ## split the query sequence into a list of codons
            split_qry = re.findall('.'*3, hsp_query)
            stripped_qry = hsp_query.replace("-", "")

            ## translate the query sequence into a protein (seq stripped of gaps because Seq hates them)
            translated_stripped_qry = str(Seq(stripped_qry).translate(table=11))

            ## iterate through query codon list, find gaps, note position, and grab all relevant information
            for qry_codons in split_qry:
                if "-" in qry_codons and len(''.join(split_qry)) % 3 == 0:    
                    qry_codon_pos += 1  # in a biological context, codons do not start "indexing" at 0; they start at 1

                    if qry_codon_pos <= len(translated_stripped_qry):
                        aa_pos, affected_codon, corr_aa, translated_codon = self.single_fs(qry_codon_pos, translated_stripped_qry, split_ref)
                        fs_ter = self.termination(translated_stripped_qry, aa_pos)

                        if fs_dict_list:  # if there are curated frameshifts in the alignment title (PVM, POM)
                            for eachfs in fs_dict_list:
                                if eachfs["original_aa"] == translated_codon and eachfs["aa_position"] == aa_pos:
                                    # frameshift found is added to 3 lists in 3 different ways
                                    fs_curated_list_reg.append(f"{translated_codon}{aa_pos}{corr_aa}") ## e.g., A15A
                                    fs_curated_list_validation.append(f"{translated_codon}{aa_pos}fs") ## e.g., A15fs
                                    fs_curated_result_HGVS.append(f"{translated_codon}{aa_pos}{corr_aa}fsTer{fs_ter}") ## e.g., A15AfsTer9
                            
                            if fs_curated_list_reg:
                                for _ in fs_curated_list_reg:
                                    if f"{translated_codon}{aa_pos}{corr_aa}" not in fs_curated_list_reg and f"{translated_codon}{aa_pos}{corr_aa}" not in fs_denovo_list_reg:
                                        fs_denovo_list_reg.append(f"{translated_codon}{aa_pos}{corr_aa}") ## e.g., A15A
                                        fs_denovo_list_validation.append(f"{translated_codon}{aa_pos}fs") ## e.g., A15fs
                                        fs_denovo_result_HGVS.append(f"{translated_codon}{aa_pos}{corr_aa}fsTer{fs_ter}") ## e.g., A15AfsTer9
                            else:
                                if f"{translated_codon}{aa_pos}{corr_aa}" not in fs_curated_list_reg:
                                    fs_denovo_list_reg.append(f"{translated_codon}{aa_pos}{corr_aa}") ## e.g., A15A
                                    fs_denovo_list_validation.append(f"{translated_codon}{aa_pos}fs") ## e.g., A15fs
                                    fs_denovo_result_HGVS.append(f"{translated_codon}{aa_pos}{corr_aa}fsTer{fs_ter}") ## e.g., A15AfsTer9
                        else:  # homologs or models that don't have curated frameshifts
                            if f"{translated_codon}{aa_pos}{corr_aa}" not in fs_denovo_list_reg:
                                fs_denovo_list_reg.append(f"{translated_codon}{aa_pos}{corr_aa}") ## e.g., A15A
                                fs_denovo_list_validation.append(f"{translated_codon}{aa_pos}fs") ## e.g., A15fs
                                fs_denovo_result_HGVS.append(f"{translated_codon}{aa_pos}{corr_aa}fsTer{fs_ter}") ## e.g., A15AfsTer9
                    else:
                        pass

                ## for any other nucleotide in the sequence DO NOT COMMENT OUT
                else:
                    qry_codon_pos += 1
            
        """for nucleotide insertions"""
        ### isolate the position of the gap(s) and the affected codon(s)
        if "-" in hsp_sbjct:
            ## split the subject sequence into a list of codons
            split_sbjct = re.findall('.'*3, hsp_sbjct)
            stripped_sbjct = hsp_sbjct.replace("-", "")
            
            ## translate the subject sequence into a protein (seq stripped of gaps because Seq hates them)
            translated_stripped_sbjct = str(Seq(stripped_sbjct).translate(table=11, gap="-"))

            ## iterate through subject codon list, find gaps, note position, and grab all relevant information
            for sbjct_codons in split_sbjct:
                if "-" in sbjct_codons and len(''.join(split_sbjct)) % 3 == 0:
                    sbjct_codon_pos += 1
                    
                    if sbjct_codon_pos <= len(translated_stripped_sbjct):
                        aa_pos, affected_codon, corr_aa, translated_codon = self.single_fs(sbjct_codon_pos, translated_stripped_sbjct, split_ref)
                        fs_ter = self.termination(translated_stripped_sbjct, aa_pos)

                        if fs_dict_list:
                            for eachfs in fs_dict_list:
                                if eachfs["original_aa"] == translated_codon and eachfs["aa_position"] == aa_pos:
                                    fs_curated_list_reg.append(f"{translated_codon}{aa_pos}{corr_aa}") ## e.g., A15A
                                    fs_curated_list_validation.append(f"{translated_codon}{aa_pos}fs") ## e.g., A15fs
                                    fs_curated_result_HGVS.append(f"{translated_codon}{aa_pos}{corr_aa}fsTer{fs_ter}") ## e.g., A15AfsTer9
                                
                            if fs_curated_list_reg:
                                for _ in fs_curated_list_reg:
                                    if f"{translated_codon}{aa_pos}{corr_aa}" not in fs_curated_list_reg:
                                        fs_denovo_list_reg.append(f"{translated_codon}{aa_pos}{corr_aa}") ## e.g., A15A
                                        fs_denovo_list_validation.append(f"{translated_codon}{aa_pos}fs") ## e.g., A15fs
                                        fs_denovo_result_HGVS.append(f"{translated_codon}{aa_pos}{corr_aa}fsTer{fs_ter}") ## e.g., A15AfsTer9
                            else:
                                if f"{translated_codon}{aa_pos}{corr_aa}" not in fs_curated_list_reg:
                                    fs_denovo_list_reg.append(f"{translated_codon}{aa_pos}{corr_aa}") ## e.g., A15A
                                    fs_denovo_list_validation.append(f"{translated_codon}{aa_pos}fs") ## e.g., A15fs
                                    fs_denovo_result_HGVS.append(f"{translated_codon}{aa_pos}{corr_aa}fsTer{fs_ter}") ## e.g., A15AfsTer9
                        else:
                            if f"{translated_codon}{aa_pos}{corr_aa}" not in fs_denovo_list_reg:
                                fs_denovo_list_reg.append(f"{translated_codon}{aa_pos}{corr_aa}") ## e.g., A15A
                                fs_denovo_list_validation.append(f"{translated_codon}{aa_pos}fs") ## e.g., A15fs
                                fs_denovo_result_HGVS.append(f"{translated_codon}{aa_pos}{corr_aa}fsTer{fs_ter}") ## e.g., A15AfsTer9
                    else:
                        pass
                else:
                    sbjct_codon_pos += 1   

        """
        frameshift output (PVM, POM, PHM)
        """
        if fs_curated_result_HGVS or fs_denovo_result_HGVS:
            fs_result_prelim["query_def"] = str(query_def)
            fs_result_prelim["mutations"] = {"type": param_type}

            # you can change the output syntax here
            if fs_curated_result_HGVS:
                fs_result_prelim["mutations"]["curated"] = fs_curated_result_HGVS
            if fs_denovo_result_HGVS:
                fs_result_prelim["mutations"]["de_novo"] = fs_denovo_result_HGVS
        elif not fs_curated_result_HGVS and not fs_denovo_result_HGVS:
            return None
        
        return fs_result_prelim            

    def indel(self, hsp_query, hsp_sbjct, card_dna_ref, query_def, insert_type="", del_type="", curated_in_list=[], curated_del_list=[]):
        """
        Searches for insertions and deletions in sequences.
        WIP: separate indels by param_type? for now, indels are indels in the RGI output
        WIP: indels that cancel each other out (e.g., 1 ins/1 del)
        """
        
        # for deletions
        deletion = {}
        qry_codon_pos = 0

        # for insertions
        sbjct_codon_pos = 0

        # indel_curated_list_reg = []
        # indel_denovo_list_reg = []

        # indel_curated_list_validation = []
        # indel_denovo_list_validation = []

        indel_curated_result_HGVS = []
        indel_denovo_result_HGVS = []

        indel_result_prelim = {}
                        
        split_ref = re.findall('.'*3, card_dna_ref)
        split_sbjct = re.findall('.'*3, hsp_sbjct)
        split_qry = re.findall('.'*3, hsp_query)  # to grab our actual inserted stretch of sequence

        """ deletions """            
        ### isolate the position of the deletion and the affected codons
        if "-" in hsp_query:
            ## split the query sequence into a list of codons
            stripped_qry = hsp_query.replace("-", "")

            ## translate the query sequence into a protein (seq stripped of gaps because Seq hates them)
            translated_stripped_qry = str(Seq(stripped_qry).translate(table=11))

            ## iterate through query codon list, find gaps + flanks, and note position
            while qry_codon_pos < len(split_qry):  # cannot be <= here because # of items and # of indeces differ, 
                                                         # so split_qry[3] when there are 3 items (0,1,2) will fail
                deletion = {}
                qry_codons = split_qry[qry_codon_pos]

                if "-" in qry_codons:
                    del_gap_count = 0
                    qry_beginning_flank = split_qry[qry_codon_pos-1]
                    deletion[qry_codon_pos-1] = qry_beginning_flank

                    while qry_codon_pos < len(split_qry) and "-" in split_qry[qry_codon_pos]:
                        qry_current_codon = split_qry[qry_codon_pos] # snapshot of the current codon
                        deletion[qry_codon_pos] = qry_current_codon
                        del_gap_count += qry_current_codon.count("-")                        

                        qry_codon_pos += 1 # updates our index to the NEXT codon after successfully identifying a gap

                    if qry_codon_pos < len(split_qry):
                        qry_ending_flank = split_qry[qry_codon_pos]
                        deletion[qry_codon_pos] = qry_ending_flank
                        
                        if del_gap_count % 3 == 0:  # checking that our deletion is clean codons and doesn't shift the frame
                            del_result = self.indel_translator(deletion, split_ref, translated_stripped_qry, indel_type = "deletion")
                        else:
                            return None

                        if del_result is not None and curated_del_list:
                            result_range = range(del_result["first_pos"], del_result["last_pos"] + 1)

                            for curated_del in curated_del_list:
                                print("====================================================")
                                print("curated deletion:\n",curated_del,"\n*********************************\n")
                                if curated_del["pos2"] == "n/a":  # format 1 
                                                                  # if curated_del["deleted"] != "n/a"? do we need that?
                                    if del_result["first_aa"] == curated_del["aa1"] and del_result["first_pos"] == curated_del["pos1"]:
                                        indel_curated_result_HGVS.append(curated_del["full_indel"])
                                        print("positions match! curated deletion here!:", curated_del, "and also", del_result,"\n")
                                        print("====================================================\n")
                                    else:
                                        if del_result["deletion"] not in indel_denovo_result_HGVS:
                                            indel_denovo_result_HGVS.append(del_result["deletion"]) ## e.g., A15A
                                            print("this doesn't match! check if it's already in the HGVS list in case it's denovo:",del_result,"\n")
                                            print("====================================================\n")
                                else:  # format 2
                                    del_range = range(curated_del["pos1"], curated_del["pos2"] + 1)
                                    if all(n in result_range for n in del_range):  # if the positions of the indel found are within the range of the 
                                                                                   # curated indel
                                        indel_curated_result_HGVS.append(curated_del["full_indel"])
                                        print("curated indel found within the range!:",del_result, "and also", curated_del)
                                        if del_result["first_pos"] != curated_del["pos1"]:
                                            indel_denovo_result_HGVS.append(del_result["deletion"]) ## e.g., A15A
                                            print("..... buuuuut positions don't match! de novo:",del_result,"\n")
                                    else:
                                        if del_result["deletion"] not in indel_denovo_result_HGVS:
                                            indel_denovo_result_HGVS.append(del_result["deletion"]) ## e.g., A15A
                                            print("this doesn't match! check if it's already in the HGVS list in case it's denovo:",del_result,"\n")
                                            print("====================================================\n")
                        elif del_result is not None and not curated_del_list:
                            if del_result["deletion"] not in indel_denovo_result_HGVS:
                                indel_denovo_result_HGVS.append(del_result["deletion"]) ## e.g., A15A
                                print("this doesn't match! check if it's already in the HGVS list in case it's denovo:",del_result,"\n")
                                print("====================================================\n")
                    else:
                        qry_ending_flank = None

                ## for any other nucleotide in the sequence DO NOT COMMENT OUT  
                else:
                    qry_codon_pos += 1
            
        """ insertions """
        ### isolate the position of the insertion and the affected codons (plus flanking codons)
        if "-" in hsp_sbjct:
            ## split the subject sequence into a list of codons
            stripped_sbjct = hsp_sbjct.replace("-", "")
            
            ## translate the subject sequence into a protein (seq stripped of gaps because Seq hates them)
            translated_stripped_sbjct = str(Seq(stripped_sbjct).translate(table=11, gap="-"))

            ## iterate through query codon list, find gaps + flanks, and note position
            while sbjct_codon_pos < len(split_sbjct):  # cannot be <= here because # of items and # of indeces differ, 
                                                         # so split_qry[3] when there are 3 items (0,1,2) will fail
                insertion = {}
                sbjct_codons = split_sbjct[sbjct_codon_pos]

                if "-" in sbjct_codons:
                    insert_gap_count = 0
                    sbjct_beginning_flank = split_sbjct[sbjct_codon_pos-1]
                    insertion[sbjct_codon_pos-1] = sbjct_beginning_flank

                    while sbjct_codon_pos < len(split_sbjct) and "-" in split_sbjct[sbjct_codon_pos]:
                        sbjct_current_codon = split_sbjct[sbjct_codon_pos] # snapshot of the current codon
                        insertion[sbjct_codon_pos] = sbjct_current_codon
                        insert_gap_count += sbjct_current_codon.count("-")                        

                        sbjct_codon_pos += 1 # updates our index to the NEXT codon after successfully identifying a gap

                    if sbjct_codon_pos < len(split_sbjct):
                        sbjct_ending_flank = split_sbjct[sbjct_codon_pos]
                        insertion[sbjct_codon_pos] = sbjct_ending_flank

                        insertion_slice = split_qry[next(iter(insertion)) + 1:next(reversed(insertion))]

                        if insert_gap_count % 3 == 0:  # checking that our insertion is clean codons and doesn't shift the frame
                            in_result = self.indel_translator(insertion, split_ref, translated_stripped_sbjct, insertion_slice=insertion_slice, indel_type = "insertion")
                        else:
                            return None

                        if in_result is not None and curated_in_list:
                            result_range = range(in_result["first_pos"], in_result["last_pos"] + 1)

                            for curated_in in curated_in_list:
                                print("====================================================")
                                print("curated insertion:\n",curated_in,"\n*********************************\n")
                                if curated_in["pos2"] == "n/a":  # format 1 
                                                                  # if curated_in["inserted"] != "n/a"? do we need that?
                                    if in_result["first_aa"] == curated_in["aa1"] and in_result["first_pos"] == curated_in["pos1"]:
                                        indel_curated_result_HGVS.append(curated_in["full_indel"])
                                        print("positions match! curated insertion here!:", curated_in, "and also", in_result,"\n")
                                        print("====================================================\n")
                                    else:
                                        if in_result["insertion"] not in indel_denovo_result_HGVS:
                                            indel_denovo_result_HGVS.append(in_result["insertion"]) ## e.g., A15A
                                            print("this doesn't match! check if it's already in the HGVS list in case it's denovo:",in_result,"\n")
                                            print("====================================================\n")
                                else:  # format 2
                                    del_range = range(curated_in["pos1"], curated_in["pos2"] + 1)
                                    if all(n in result_range for n in del_range):  # if the positions of the indel found are within the range of the 
                                                                                   # curated indel
                                        indel_curated_result_HGVS.append(curated_in["full_indel"])
                                        print("curated indel found within the range!:",in_result, "and also", curated_in)
                                        if in_result["first_pos"] != curated_in["pos1"]:
                                            indel_denovo_result_HGVS.append(in_result["insertion"]) ## e.g., A15A
                                            print("..... buuuuut positions don't match! de novo:",in_result,"\n")
                                    else:
                                        if in_result["insertion"] not in indel_denovo_result_HGVS:
                                            indel_denovo_result_HGVS.append(in_result["insertion"]) ## e.g., A15A
                                            print("this doesn't match! check if it's already in the HGVS list in case it's denovo:",in_result,"\n")
                                            print("====================================================\n")
                        elif in_result is not None and not curated_in_list:
                            if in_result["insertion"] not in indel_denovo_result_HGVS:
                                indel_denovo_result_HGVS.append(in_result["insertion"]) ## e.g., A15A
                                print("this doesn't match! check if it's already in the HGVS list in case it's denovo:",in_result,"\n")
                                print("====================================================\n")
                    else:
                        sbjct_ending_flank = None

                ## for any other nucleotide in the sequence DO NOT COMMENT OUT  
                else:
                    sbjct_codon_pos += 1
        
        if indel_curated_result_HGVS or indel_denovo_result_HGVS:
            indel_result_prelim["query_def"] = str(query_def)
            indel_result_prelim["mutations"] = {"type": "indel mutation from peptide sequence"}

            # you can change the output syntax here
            if indel_curated_result_HGVS:
                indel_result_prelim["mutations"]["curated"] = indel_curated_result_HGVS
            if indel_denovo_result_HGVS:
                indel_result_prelim["mutations"]["de_novo"] = indel_denovo_result_HGVS
        elif not indel_curated_result_HGVS and not indel_denovo_result_HGVS:
            return None

        return indel_result_prelim            

    def single_fs(self, codon_count, translated_stripped_seq, split_ref):
        aa_pos = codon_count
        affected_codon = split_ref[aa_pos - 1]
        corr_aa = translated_stripped_seq[aa_pos - 1]  # index starts at 0
        translated_codon = str(Seq(affected_codon).translate(table=11))

        return aa_pos, affected_codon, corr_aa, translated_codon  # aa_pos is just codon_count... fix that
    
    def termination(self, translated_stripped_seq, aa_pos):
        aa_count = 0

        ## locating the frameshift in the translated protein (entire seq. chunk + position of termination)
        for aa in translated_stripped_seq[aa_pos - 1:]: # index starts at 0
            if aa == "*":
                break
            else:
                aa_count += 1
        
        return aa_count + 1
    
    def indel_translator(self, indel, split_ref, translated_stripped_seq, insertion_slice=[], indel_type=None):
        unpacked_indel = list(indel.items())

        # validating the positions in our indel to make sure nothing is out of bounds (i've learned my lesson)
        max_pos = max([pos for pos, codon in unpacked_indel])  # finding the max position (our upper bound)

        if max_pos - 1 >= len(translated_stripped_seq):
            return None

        del_codons = ""

        if indel_type == "insertion":
            in_dict = {}
            in_slice = str(Seq(insertion_slice[0]).translate(table=11))

            for i, (position, codon) in enumerate(unpacked_indel):  # we don't actually access codon... but you never know when you'll need it? :^)
                if i == 0:
                    affected_codon = split_ref[position]
                    original_aa = str(Seq(affected_codon).translate(table=11))

                    beginning_flank = original_aa
                    beginning_pos = position + 1  # adjusts the position so it isn't just the index of the string

                    in_dict["first_aa"] = beginning_flank
                    in_dict["first_pos"] = beginning_pos 
                elif i == 1:
                    middle_pos = position + 1

                elif i == len(unpacked_indel) - 1:
                    affected_codon = split_ref[position - 1]
                    original_aa = str(Seq(affected_codon).translate(table=11))

                    end_flank = original_aa
                    end_pos = position + 1

                    in_dict["last_aa"] = end_flank
                    in_dict["last_pos"] = end_pos

                else:  # everything in between the sandwich
                    new_aa = translated_stripped_seq[position - 1] # index starts at 0
            
            in_dict["insertion"] = f"{beginning_flank}{beginning_pos}_{end_flank}{middle_pos}ins{in_slice}"

            return in_dict
        else:
            del_dict = {}
            for i, (position, codon) in enumerate(unpacked_indel):
                if i == 1:
                    affected_codon = split_ref[position]  # looking at a list index position here
                    original_aa = str(Seq(affected_codon).translate(table=11))

                    first_aa = original_aa  # deletions don't have flanks--the first aa reported is where the deletion starts (index 1, not 0)
                    first_pos = position + 1  # looking at biological position here

                    del_dict["first_aa"] = first_aa
                    del_dict["first_pos"] = first_pos 
                elif i == len(unpacked_indel) - 2:
                    affected_codon = split_ref[position]
                    original_aa = str(Seq(affected_codon).translate(table=11))

                    last_aa = original_aa
                    last_pos = position + 1
                    del_dict["last_aa"] = last_aa
                    del_dict["last_pos"] = last_pos
                else: 
                    new_aa = translated_stripped_seq[position - 1]

                    del_codons += new_aa
            
            del_dict["deletion"] = f"{first_aa}{first_pos}_{last_aa}{last_pos}del{del_codons}"

            return del_dict

    def consolidate_mutations(self, input_type, hit_id, model_type, srv=None, other_mutations=[], phm=None, hsp_bitscore=None, pass_val=None):
        """
        Consolidates the results from all mutation functions and passes the output back into RGI's detection modules (PHM, PVM, POM, RGV).
        """

        has_snp = srv.get("has_snp", False) if srv is not None else False  # the PHM does not generate srv, so this avoids an AttributeError
        passes_eval = float(hsp_bitscore) >= float(pass_val)
        merged_mutations = {
            "query_def": "",
            "curated_mutations": None,
            "de_novo_mutations": None
        }

        # protein input
        if input_type == "protein":
            if srv is None and phm:
                return []
            elif srv is not None and has_snp:  # CASE 1: if the input is a protein there won't be a BLASTN xml generated
                return [srv]

        # nucleotide input
        elif input_type == "contig":
            if not other_mutations:
                if srv is not None and has_snp:
                    return [srv] # CASE 1: only SNP
                else:
                    return []  
            else:
                for mutations in other_mutations:
                    merged_mutations["query_def"] = mutations["query_def"]
                    mutation_id = mutations["query_def"].split()[0]

                    # building out our nested dictionary of mutations lists (in subdictionaries...)
                    # we want rich merged output
                    if mutations["mutations"].get("curated"):
                        curated_entry = {
                            mutations["mutations"]["type"]: mutations["mutations"]["curated"],
                        }
                        merged_mutations["curated_mutations"] = curated_entry

                    if mutations["mutations"].get("de_novo"):
                        de_novo_entry = {
                            mutations["mutations"]["type"]: mutations["mutations"]["de_novo"],
                        }
                        merged_mutations["de_novo_mutations"] = de_novo_entry

                has_curated_mutation = merged_mutations.get("curated_mutations") is not None
                has_denovo_mutation = merged_mutations.get("de_novo_mutations") is not None
                has_other_mutations = has_curated_mutation or has_denovo_mutation

                if srv is not None:
                    srv_id = srv["query_def"].split()[0]
                    
                    if mutation_id in srv_id:
                        # protein variant frameshift/SNP search
                        if model_type == "PVM":
                            if has_snp and has_other_mutations:  # CASE 2: SNP + other mutations
                                return [srv | merged_mutations]
                                    
                            elif not has_snp and has_other_mutations:  # CASE 3: other mutations found; SNPs were not found in any HSPs with SNP in the alignment title
                                if passes_eval and has_curated_mutation:  # Strict alignments
                                    merged_mutations["query_def"] += hit_id
                                    return [merged_mutations]
                                elif not passes_eval and has_curated_mutation:  # Loose alignments
                                    merged_mutations["query_def"] += hit_id
                                    return [merged_mutations]

                        # protein overexpression frameshift/SNP search
                        if model_type == "POM":
                            srv_id = srv["query_def"].split()[0]
                        
                            if has_snp and has_other_mutations:
                                return [srv | merged_mutations]
                                    
                            elif not has_snp and has_other_mutations:  # CASE 3: other mutations found; SNPs were not found in any HSPs with SNP in the alignment title
                                if passes_eval and has_curated_mutation:  # Strict alignments
                                    merged_mutations["query_def"] += hit_id
                                    return [merged_mutations]
                                elif not passes_eval and has_curated_mutation:  # Loose alignments
                                    merged_mutations["query_def"] += hit_id
                                    return [merged_mutations]
                                    
                # protein homolog frameshift search               
                if model_type == "PHM":
                    phm_id = phm["query_def"].split()[0]
                
                    if mutation_id in phm_id:
                        if passes_eval and has_denovo_mutation :
                            # Perfect PHMs will not have frameshifts (etc.) in them, so there's no need to add unique support for them here
                            merged_mutations["query_def"] += hit_id
                            return [merged_mutations]
                        elif not passes_eval and has_denovo_mutation:
                            merged_mutations["query_def"] += hit_id
                            return [merged_mutations]