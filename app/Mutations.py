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

    def frameshift(self, hsp_query, hsp_sbjct, card_dna_ref, query_def, hsp_sbjct_start=1, param_type=None, fs_dict_list=None): 
        """
        Searches for frameshifts in nucleotide sequences.

        ** Code optimized with Codex **
        """
        
        fs_result_prelim = {}
    
        fs_curated_list_reg = []
        fs_denovo_list_reg = []

        fs_curated_list_validation = []
        fs_denovo_list_validation = []

        fs_curated_result_HGVS = []
        fs_denovo_result_HGVS = []

        if fs_dict_list is None:
            fs_dict_list = []

        translated_qry = str(Seq(hsp_query.replace("-", "")).translate(table=11))  # Seq() hates dashed gaps, so we stripped them
        hsp_start_codon = (hsp_sbjct_start - 1) // 3

        # print(f"query_def: {query_def}")

        for raw_event in self.gap_events(hsp_query, hsp_sbjct, hsp_sbjct_start):

            if not raw_event["is_frameshift"] or raw_event["reference_nucl_index"] is None:
                continue

            # print(f"gap event\n{raw_event}\n")
            event, is_curated = self.gap_resolver(raw_event, card_dna_ref, fs_dict_list)

            aa_pos = (event["reference_nucl_index"] // 3) + 1  # converting our zero-based nucleotide start position to one-based amino acid space
            codon_start = (aa_pos - 1) * 3  # zero-based
            qry_aa_index = (aa_pos - 1) - hsp_start_codon  # zero-based

            if qry_aa_index < 0 or qry_aa_index >= len(translated_qry):
                continue

            translated_ref_codon = str(Seq(card_dna_ref[codon_start:codon_start + 3]).translate(table=11))
            new_aa = translated_qry[qry_aa_index]

            fs_termination = self.termination(translated_qry, qry_aa_index + 1)
            fs_reg = f"{translated_ref_codon}{aa_pos}{new_aa}"
            fs_validation = f"{translated_ref_codon}{aa_pos}fs"
            fs_hgvs = f"{fs_reg}fsTer{fs_termination}"

            # print(f"regular degular: {fs_reg}\nvalidation style: {fs_validation}\nHGVS syntax: {fs_hgvs}\n\n")

            if is_curated:
                fs_curated_list_reg.append(fs_reg)  # e.g., A151A
                fs_curated_list_validation.append(fs_validation)  # e.g., A15fs
                fs_curated_result_HGVS.append(fs_hgvs)  # e.g., A15AfsTer9
            else:
                fs_denovo_list_reg.append(fs_reg)
                fs_denovo_list_validation.append(fs_validation)
                fs_denovo_result_HGVS.append(fs_hgvs)

        fs_curated_list_validation = list(dict.fromkeys(fs_curated_list_validation))  
        fs_denovo_list_validation = list(dict.fromkeys(fs_denovo_list_validation))
        fs_curated_result_HGVS = list(dict.fromkeys(fs_curated_result_HGVS))
        fs_denovo_result_HGVS = list(dict.fromkeys(fs_denovo_result_HGVS))

        """
        frameshift output (PVM, POM, PHM)
        """
        if fs_curated_result_HGVS or fs_denovo_result_HGVS:
            fs_result_prelim["query_def"] = str(query_def)

            if param_type:
                fs_result_prelim["mutations"] = {"type": param_type}
            else:
                fs_result_prelim["mutations"] = {"type": "frameshift mutation"}

            # you can change the output syntax here
            if fs_curated_list_validation:
                fs_result_prelim["mutations"]["curated"] = fs_curated_list_validation
            if fs_denovo_result_HGVS:
                fs_result_prelim["mutations"]["de_novo"] = fs_denovo_list_validation

        elif not fs_curated_result_HGVS and not fs_denovo_result_HGVS:
            return None

        # print(f"{fs_result_prelim}\n")
        return fs_result_prelim    

    def indel(self, hsp_query, hsp_sbjct, card_dna_ref, query_def, insert_type="", del_type="", curated_in_list=None, curated_del_list=None):
        """
        Searches for insertions and deletions in sequences.
        WIP: indels that cancel each other out (e.g., 1 ins/1 del)
        """
        
        # for deletions
        deletion = {}
        qry_codon_pos = 0

        # for insertions
        sbjct_codon_pos = 0

        if curated_in_list is None:
            curated_in_list = []
        if curated_del_list is None:
            curated_del_list = []

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
                                if curated_del["pos2"] == "n/a":  # format 1 
                                                                  # if curated_del["deleted"] != "n/a"? do we need that?
                                    if del_result["first_aa"] == curated_del["aa1"] and del_result["first_pos"] == curated_del["pos1"]:
                                        indel_curated_result_HGVS.append(curated_del["full_indel"])
                                    else:
                                        if del_result["deletion"] not in indel_denovo_result_HGVS:
                                            indel_denovo_result_HGVS.append(del_result["deletion"]) ## e.g., A15A
                                else:  # format 2
                                    del_range = range(curated_del["pos1"], curated_del["pos2"] + 1)
                                    if all(n in result_range for n in del_range):  # if the positions of the indel found are within the range of the 
                                                                                   # curated indel
                                        indel_curated_result_HGVS.append(curated_del["full_indel"])
                                        if del_result["first_pos"] != curated_del["pos1"]:
                                            indel_denovo_result_HGVS.append(del_result["deletion"]) ## e.g., A15A
                                    else:
                                        if del_result["deletion"] not in indel_denovo_result_HGVS:
                                            indel_denovo_result_HGVS.append(del_result["deletion"]) ## e.g., A15A
                        elif del_result is not None and not curated_del_list:
                            if del_result["deletion"] not in indel_denovo_result_HGVS:
                                indel_denovo_result_HGVS.append(del_result["deletion"]) ## e.g., A15A
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
                                if curated_in["pos2"] == "n/a":  # format 1 
                                                                  # if curated_in["inserted"] != "n/a"? do we need that?
                                    if in_result["first_aa"] == curated_in["aa1"] and in_result["first_pos"] == curated_in["pos1"]:
                                        indel_curated_result_HGVS.append(curated_in["full_indel"])
                                    else:
                                        if in_result["insertion"] not in indel_denovo_result_HGVS:
                                            indel_denovo_result_HGVS.append(in_result["insertion"]) ## e.g., A15A
                                else:  # format 2
                                    del_range = range(curated_in["pos1"], curated_in["pos2"] + 1)
                                    if all(n in result_range for n in del_range):  # if the positions of the indel found are within the range of the 
                                                                                   # curated indel
                                        indel_curated_result_HGVS.append(curated_in["full_indel"])
                                        if in_result["first_pos"] != curated_in["pos1"]:
                                            indel_denovo_result_HGVS.append(in_result["insertion"]) ## e.g., A15A
                                    else:
                                        if in_result["insertion"] not in indel_denovo_result_HGVS:
                                            indel_denovo_result_HGVS.append(in_result["insertion"]) ## e.g., A15A
                        elif in_result is not None and not curated_in_list:
                            if in_result["insertion"] not in indel_denovo_result_HGVS:
                                indel_denovo_result_HGVS.append(in_result["insertion"]) ## e.g., A15A
                    else:
                        sbjct_ending_flank = None

                ## for any other nucleotide in the sequence DO NOT COMMENT OUT  
                else:
                    sbjct_codon_pos += 1
        
        indel_curated_result_HGVS = list(dict.fromkeys(indel_curated_result_HGVS)) # bandaid solution to dedupe mutations...
        indel_denovo_result_HGVS = list(dict.fromkeys(indel_denovo_result_HGVS))

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

    def nonsense(self, hsp_query, hsp_sbjct, card_dna_ref, query_def, param_type=None, ns_dict_list=None):
        """
        Searches for nonsense mutations in sequences.
        """
        
        ns_result_prelim = {}
    
        ns_curated_result_HGVS = []
        ns_denovo_result_HGVS = []

        if ns_dict_list is None:
            ns_dict_list = []

        split_ref = re.findall('.'*3, card_dna_ref)
        split_qry = re.findall('.'*3, hsp_query)  # to grab our actual inserted stretch of sequence

        stop_codons = ["UAA", "UAG", "UGA", "TAA", "TAG", "TGA"]

        ## test print statements
        # print(f"curated nonsense mutations:{ns_dict_list}\nsplit ref: {split_ref}\nsplit query: {split_qry}\n")

        for stop_pos, qry_codon in enumerate(split_qry, start=1):
            if qry_codon not in stop_codons:
                continue

            ref_index = stop_pos - 1

            # the aligned query and CARD reference can have different lengths
            # skip stops that cannot be mapped instead of indexing past split_ref
            if not 0 <= ref_index < len(split_ref):
                continue

            affected_aa = str(Seq(split_ref[ref_index]).translate(table=11))

            # a stop already present in CARD is the expected termination codon, not a newly introduced nonsense mutation
            if affected_aa == "*":
                continue

            nonsense_mutation = f"{affected_aa}{stop_pos}Ter"

            if ns_dict_list:
                for eachns in ns_dict_list:
                    if eachns["original_aa"] == affected_aa and eachns["aa_position"] == stop_pos:
                        ns_curated_result_HGVS.append(nonsense_mutation)

                if nonsense_mutation not in ns_curated_result_HGVS:
                    ns_denovo_result_HGVS.append(nonsense_mutation)
            elif nonsense_mutation not in ns_denovo_result_HGVS:
                ns_denovo_result_HGVS.append(nonsense_mutation)

        ns_curated_result_HGVS = list(dict.fromkeys(ns_curated_result_HGVS)) # bandaid solution to dedupe mutations...
        ns_denovo_result_HGVS = list(dict.fromkeys(ns_denovo_result_HGVS))

        """
        nonsense output (PVM, POM, PHM)
        """
        if ns_curated_result_HGVS or ns_denovo_result_HGVS:
            ns_result_prelim["query_def"] = str(query_def)
            if param_type:
                ns_result_prelim["mutations"] = {"type": param_type}
            else:
                ns_result_prelim["mutations"] = {"type": "nonsense mutation"}

            if ns_curated_result_HGVS:
                ns_result_prelim["mutations"]["curated"] = ns_curated_result_HGVS
            if ns_denovo_result_HGVS:
                ns_result_prelim["mutations"]["de_novo"] = ns_denovo_result_HGVS
        elif not ns_curated_result_HGVS and not ns_denovo_result_HGVS:
            return None
        
        return ns_result_prelim

    def gap_events(self, hsp_query, hsp_sbjct, hsp_sbjct_start):
        """
        Tracks non-triplet gap events (not full codon indels) relative to CARD reference nucleotide coordinates

        ** Code optimized with Codex **
        """
        if len(hsp_query) != len(hsp_sbjct):
            raise ValueError("Aligned query and subject must have equal lengths.")

        reference_index = hsp_sbjct_start - 1 # our subject: the CARD reference (zero-based)
        alignment_index = 0  # our pointer as we go character by character through a sequence string (zero-based)

        # print("\n====================\n"
        #     f"hsp_sbjct:\n{hsp_sbjct}\n"
        #     f"hsp_query:\n{hsp_query}\n"
        #     f"hsp_sbjct_start: {hsp_sbjct_start}\n"
        #     f"reference_index: {reference_index}\n"
        # )

        while alignment_index < len(hsp_query):  # while we aren't over the total length of our (potentially) gappy query
            query_base = hsp_query[alignment_index]
            subject_base = hsp_sbjct[alignment_index]

            """nucleotide deletions"""
            if query_base == "-" and subject_base != "-":  # if a gap (deletion) is found, the game is on
                gap_event_start = reference_index
                gap_event_bases = []

                while (  # once we've tracked a gap, keep going until it ends
                    alignment_index < len(hsp_query)
                    and hsp_query[alignment_index] == "-"
                    and hsp_sbjct[alignment_index] != "-"
                       ):
                    gap_event_bases.append(hsp_sbjct[alignment_index])
                    reference_index += 1
                    alignment_index += 1

                is_frameshift = len(gap_event_bases) % 3 != 0

                yield {
                    "type": "deletion",
                    "is_frameshift": is_frameshift,
                    "reference_nucl_index": gap_event_start,  # deletions are deletions (see insertion note for explanation)
                    "length": len(gap_event_bases),
                    "sequence": "".join(gap_event_bases),
                    }
                continue

            """nucleotide insertions"""
            if subject_base == "-" and query_base != "-":
                gap_event_bases = []

                while (  # game starts
                    alignment_index < len(hsp_query)
                    and hsp_sbjct[alignment_index] == "-"
                    and hsp_query[alignment_index] != "-"
                       ):
                    gap_event_bases.append(hsp_query[alignment_index])
                    alignment_index += 1

                is_frameshift = len(gap_event_bases) % 3 != 0  # we want to start counting insertions AFTER a reference base has been "consumed"
                                                               # note: the code doesn't yet handle insertions before the start of the reference sequence
                has_anchor = reference_index > 0  # lets us know if the insertion happens before or after the start of the reference sequence
                reference_nucl_index = reference_index - 1 if has_anchor else None  # insertions are reported relative to the preceding 
                                                                                    # reference base/codon because they happen in between bases

                yield {
                    "type": "insertion",
                    "is_frameshift": is_frameshift,
                    "reference_nucl_index": reference_nucl_index,
                    "length": len(gap_event_bases),
                    "sequence": "".join(gap_event_bases),
                    }
                continue

            if subject_base != "-":
                reference_index += 1

            alignment_index += 1

    def equivalent_gap_events(self, event, card_dna_ref):
        """
        Searches from the start of a gap event outwards and returns every equivalent alignment position for a gap event found in stretches of repeated sequence.

        ** Code optimized with Codex **
        """
        # print(f"card DNA reference\n{card_dna_ref}\n")
        candidates = [event]  # the intial gap event is the only candidate to start -- we intialized the list with our baseline gap event
        # print(f"initial gap event: {event}")

        """deletions"""
        if event["type"] == "deletion":
            event_start = event["reference_nucl_index"]
            event_length = event["length"]

            ## left outward search
            left_start = event_start

            # print(f"event sequence: {event["sequence"]}\nleft boundary: {card_dna_ref[left_start - 1]} and right boundary {card_dna_ref[left_start + event_length - 1]}")

            while(left_start > 0 
                  and card_dna_ref[left_start - 1] == card_dna_ref[left_start + event_length - 1]
                  ):
                left_start -= 1  # shift the index left
                candidate = event.copy()  # making a copy of the viable candidate we just found
                candidate["reference_nucl_index"] = left_start  # replacing its index with the shifted index
                candidates.append(candidate)

            ## right outward search
            right_start = event_start

            while(right_start + event_length < len(card_dna_ref)
                  and card_dna_ref[right_start] == card_dna_ref[right_start + event_length]
                  ):
                right_start += 1  # shift the index right
                candidate = event.copy()
                candidate["reference_nucl_index"] = right_start
                # print(f"\nright search candidate\n{candidate}\n")
                candidates.append(candidate)

            """insertions"""
        elif event["type"] == "insertion":
            boundary = event["reference_nucl_index"] + 1
            inserted_sequence = event["sequence"]
            # print(f"\ninsertion boundary: {boundary}\ninserted sequence: {inserted_sequence}\n")

            ## left outward search
            left_boundary = boundary
            left_sequence = inserted_sequence

            while(left_boundary > 0
                  and left_sequence
                  and card_dna_ref[left_boundary - 1] == left_sequence[-1]
                  ):
                left_sequence = (card_dna_ref[left_boundary - 1] + left_sequence[:-1])  # rotating the sequence 
                left_boundary -= 1

                candidate = event.copy()
                candidate["reference_nucl_index"] = left_boundary - 1
                candidate["sequence"] = left_sequence
                candidates.append(candidate)

            ##right outward search
            right_boundary = boundary
            right_sequence = inserted_sequence      
            while (right_boundary < len(card_dna_ref)
                   and right_sequence
                   and card_dna_ref[right_boundary] == right_sequence[0]
                   ):
                right_sequence = (right_sequence[1:] + card_dna_ref[right_boundary])
                right_boundary += 1

                candidate = event.copy()
                candidate["reference_nucl_index"] = right_boundary - 1
                candidate["sequence"] = right_sequence
                candidates.append(candidate)

        return candidates      

    def gap_resolver(self, event, card_dna_ref, fs_dict_list):
        """
        Looks at all gap event candidates, sees if there's a match to a CARD curated frameshift/indel, and chooses the curated option if possible.
        
        ** Code optimized with Codex **
        """
        curated_candidates = []

        for candidate in self.equivalent_gap_events(event, card_dna_ref):
            aa_pos = (candidate["reference_nucl_index"] // 3 + 1)
            codon_start_index = (aa_pos - 1) * 3
            reference_codon = card_dna_ref[codon_start_index:codon_start_index + 3]
            original_aa = str(Seq(reference_codon).translate(table=11))

            if any(
                eachfs["original_aa"] == original_aa
                and eachfs["aa_position"] == aa_pos
                for eachfs in fs_dict_list
                ):
                curated_candidates.append(candidate)

        if curated_candidates:
            return min(
                curated_candidates,
                key=lambda candidate: abs(candidate["reference_nucl_index"] - event["reference_nucl_index"])
                # super funky, but here we're trying to see which curated candidate is closest to the initial gap event
                ), True

        return event, False

    def single_fs(self, aa_pos, translated_stripped_seq, split_ref):
        affected_codon = split_ref[aa_pos - 1]
        corr_aa = translated_stripped_seq[aa_pos - 1]  # index starts at 0
        translated_codon = str(Seq(affected_codon).translate(table=11))

        return aa_pos, affected_codon, corr_aa, translated_codon
    
    def termination(self, translated_seq, aa_pos):
        aa_count = 0

        ## locating the frameshift in the translated protein (entire seq. chunk + position of termination)
        for aa in translated_seq[aa_pos - 1:]: # index starts at 0
            if aa == "*":
                break
            else:
                aa_count += 1
        
        return aa_count + 1
    
    def indel_translator(self, indel, split_ref, translated_stripped_seq, insertion_slice=None, indel_type=None):            
        unpacked_indel = list(indel.items())

        if insertion_slice is None:
            insertion_slice = []

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
                if i == len(unpacked_indel) - 2:
                    affected_codon = split_ref[position]
                    original_aa = str(Seq(affected_codon).translate(table=11))

                    last_aa = original_aa
                    last_pos = position + 1
                    del_dict["last_aa"] = last_aa
                    del_dict["last_pos"] = last_pos
                if i != 1 and i != len(unpacked_indel) - 2: 
                    del_codons += translated_stripped_seq[position - 1]
            
            del_dict["deletion"] = f"{first_aa}{first_pos}_{last_aa}{last_pos}del{del_codons}"

            return del_dict

    def consolidate_mutations(self, input_type, hit_id, model_type, srv=None, other_mutations=None, phm=None, hsp_bitscore=None, pass_val=None):
        """
        Consolidates the results from all mutation functions and passes the output back into RGI's detection modules (PHM, PVM, POM, RGV).
        """

        has_snp = srv.get("has_snp", False) if srv is not None else False  # the PHM does not generate srv, so this avoids an AttributeError
        passes_eval = float(hsp_bitscore) >= float(pass_val)
        merged_mutations = {
            "query_def": "",
            "curated_mutations": {},
            "de_novo_mutations": {}
        }

        if other_mutations is None:
            other_mutations = []

        # protein input
        # CASE 1: if the input is a protein there won't be a BLASTN xml generated
        if input_type == "protein":
            if srv is None:
                return []
            
            if model_type == "POM":
                return [srv]
            
            if has_snp:
                return [srv]
            
            return []

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
                    mutation_type = mutations["mutations"]["type"]

                    # building out our nested dictionary of mutations lists (in subdictionaries...)
                    # we want rich merged output
                    if mutations["mutations"].get("curated"):
                        merged_mutations["curated_mutations"].setdefault(mutation_type, [])
                        merged_mutations["curated_mutations"][mutation_type].extend(
                            mutations["mutations"]["curated"]
                        )

                    if mutations["mutations"].get("de_novo"):
                        merged_mutations["de_novo_mutations"].setdefault(mutation_type, [])
                        merged_mutations["de_novo_mutations"][mutation_type].extend(
                            mutations["mutations"]["de_novo"]
                        )

                for bucket in ["curated_mutations", "de_novo_mutations"]:
                    if merged_mutations[bucket]:
                        for mutation_type, values in merged_mutations[bucket].items():
                            merged_mutations[bucket][mutation_type] = list(dict.fromkeys(values))
                    else:
                        merged_mutations[bucket] = None

                if not merged_mutations["curated_mutations"]:
                    merged_mutations["curated_mutations"] = None

                if not merged_mutations["de_novo_mutations"]:
                    merged_mutations["de_novo_mutations"] = None

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
                                elif passes_eval and has_denovo_mutation:  # Strict alignments w de novo mutations
                                    merged_mutations["query_def"] += hit_id
                                    return [merged_mutations]
                                elif not passes_eval and has_denovo_mutation:  # Loose alignments w de novo mutations
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
                                elif passes_eval and has_denovo_mutation:  # Strict alignments w de novo mutations
                                    merged_mutations["query_def"] += hit_id
                                    return [merged_mutations]
                                elif not passes_eval and has_denovo_mutation:  # Loose alignments w de novo mutations
                                    merged_mutations["query_def"] += hit_id
                                    return [merged_mutations]
                                elif not passes_eval and has_curated_mutation:  # Loose alignments
                                    merged_mutations["query_def"] += hit_id
                                    return [merged_mutations]
                                    
                # protein homolog frameshift search               
                if model_type == "PHM":
                    phm_id = phm["query_def"].split()[0]
                
                    if mutation_id in phm_id:
                        if passes_eval and has_denovo_mutation:  # Strict alignments
                            # Perfect PHMs will not have frameshifts (etc.) in them, so there's no need to add unique support for them here
                            merged_mutations["query_def"] += hit_id
                            return [merged_mutations]
                        elif not passes_eval and has_denovo_mutation:  # Loose alignments (always de novo with PHM)
                            merged_mutations["query_def"] += hit_id
                            return [merged_mutations]
