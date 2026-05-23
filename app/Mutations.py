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

                    if hsp_query[qry] == chan:
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

    def frameshift(self, hsp_query, hsp_sbjct, card_dna_ref, query_def, fs_dict_list=[]): 
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
        
        if len(fs_dict_list) != 0:
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
                    if "-" in qry_codons:
                        qry_codon_pos += 1  # in a biological context, codons do not start "indexing" at 0; they start at 1

                        if qry_codon_pos <= len(translated_stripped_qry):
                            aa_pos, affected_codon, corr_aa, translated_codon = self.single_fs(qry_codon_pos, translated_stripped_qry, split_ref)
                            fs_ter = self.termination(translated_stripped_qry, aa_pos)

                            for eachfs in fs_dict_list:
                                if eachfs["original_aa"] == translated_codon and eachfs["aa_position"] == aa_pos:
                                    # frameshift found is added to 3 lists in 3 different ways
                                    fs_curated_list_reg.append(f"{translated_codon}{aa_pos}{corr_aa}") ## e.g., A15A
                                    fs_curated_list_validation.append(f"{translated_codon}{aa_pos}fs") ## e.g., A15fs
                                    fs_curated_result_HGVS.append(f"{translated_codon}{aa_pos}{corr_aa}fsTer{fs_ter}") ## e.g., A15AfsTer9
                            
                            if len(fs_curated_list_reg) != 0:
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
                    if "-" in sbjct_codons:
                        sbjct_codon_pos += 1
                        
                        if sbjct_codon_pos <= len(translated_stripped_sbjct):
                            aa_pos, affected_codon, corr_aa, translated_codon = self.single_fs(sbjct_codon_pos, translated_stripped_sbjct, split_ref)
                            fs_ter = self.termination(translated_stripped_sbjct, aa_pos)

                            for eachfs in fs_dict_list:
                                if eachfs["original_aa"] == translated_codon and eachfs["aa_position"] == aa_pos:
                                    # frameshift found is added to 3 lists in 3 different ways
                                    fs_curated_list_reg.append(f"{translated_codon}{aa_pos}{corr_aa}") ## e.g., A15A
                                    fs_curated_list_validation.append(f"{translated_codon}{aa_pos}fs") ## e.g., A15fs
                                    fs_curated_result_HGVS.append(f"{translated_codon}{aa_pos}{corr_aa}fsTer{fs_ter}") ## e.g., A15AfsTer9
                                
                            if len(fs_curated_list_reg) != 0:
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
                        sbjct_codon_pos += 1   

            """
            frameshift output (PVM & POM)
            """
            if len(fs_curated_result_HGVS) > 0 or len(fs_denovo_result_HGVS) > 0:
                fs_result_prelim["query_def"] = str(query_def)
                fs_result_prelim["has_fs"] = True

                # you can change the output syntax here
                if len(fs_curated_result_HGVS) > 0:
                    fs_result_prelim["curated_fs"] = fs_curated_result_HGVS
                if len(fs_denovo_result_HGVS) > 0:
                    fs_result_prelim["denovo_fs"] = fs_denovo_result_HGVS

            return fs_result_prelim
        
        ## homologs
        else: 
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
                    if "-" in qry_codons:
                        qry_codon_pos += 1

                        if qry_codon_pos <= len(translated_stripped_qry):
                            aa_pos, affected_codon, corr_aa, translated_codon = self.single_fs(qry_codon_pos, translated_stripped_qry, split_ref)
                            fs_ter = self.termination(translated_stripped_qry, aa_pos)

                            if f"{translated_codon}{aa_pos}{corr_aa}" not in fs_denovo_list_reg:
                                # frameshift found is added to 3 lists in 3 different ways
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
                    if "-" in sbjct_codons:
                        sbjct_codon_pos += 1
                        
                        if sbjct_codon_pos <= len(translated_stripped_sbjct):
                            aa_pos, affected_codon, corr_aa, translated_codon = self.single_fs(sbjct_codon_pos, translated_stripped_sbjct, split_ref)
                            fs_ter = self.termination(translated_stripped_sbjct, aa_pos)

                            if f"{translated_codon}{aa_pos}{corr_aa}" not in fs_denovo_list_reg:
                                # frameshift found is added to 3 lists in 3 different ways
                                fs_denovo_list_reg.append(f"{translated_codon}{aa_pos}{corr_aa}") ## e.g., A15A
                                fs_denovo_list_validation.append(f"{translated_codon}{aa_pos}fs") ## e.g., A15fs
                                fs_denovo_result_HGVS.append(f"{translated_codon}{aa_pos}{corr_aa}fsTer{fs_ter}") ## e.g., A15AfsTer9
                        else:
                            pass
                    else:
                        sbjct_codon_pos += 1   

            """
            frameshift output (PHM)
            """
            if len(fs_denovo_result_HGVS) > 0:
                # you can change the output syntax here
                fs_result_prelim["query_def"] = str(query_def)
                fs_result_prelim["denovo_fs"] = fs_denovo_result_HGVS
                fs_result_prelim["has_fs"] = True
            else:
                fs_result_prelim["query_def"] = str(query_def)
                fs_result_prelim["has_fs"] = False

            return fs_result_prelim

    def indel(self, hsp_query, hsp_sbjct, card_dna_ref, query_def):
        """
        Searches for insertions and deletions in sequences.
        Indels must be 3 nucleotides or greater to be detected. Anything else is detected by def_frameshift(). 
        
        Note: indels < 3 nucleotides, and handeled by def_frameshift(), do not necessarily result in a frameshift.

        https://hgvs-nomenclature.org/stable/recommendations/protein/insertion/
        Insertion: a sequence change between the translation initiation (start) and termination (stop) codon where, compared to the reference sequence, one or more amino acids are inserted, which is not a frameshift and where the insertion is not a copy of a sequence immediately N-terminal (5').
        """
        
        # for deletions
        qry_codon_count = 0

        # for insertions
        sbjct_codon_count = 0

        # indel_curated_list_reg = []
        indel_denovo_list_reg = []

        # indel_curated_list_validation = []
        indel_denovo_list_validation = []

        # indel_curated_result_HGVS = []
        indel_denovo_result_HGVS = []

        indel_result_prelim = {}
                        
        split_ref = re.findall('.'*3, card_dna_ref)

        """ deletions """            
        ### isolate the position of the indel and the affected codons (plus flanking codons)
        if "-" in hsp_query:
            ## split the query sequence into a list of codons
            split_qry = re.findall('.'*3, hsp_query)
            stripped_qry = hsp_query.replace("-", "")

            ## translate the query sequence into a protein (seq stripped of gaps because Seq hates them)
            translated_stripped_qry = str(Seq(stripped_qry).translate(table=11))

            ## iterate through query codon list, find gaps + flanks, and note position
            while qry_codon_count < len(split_qry):  # cannot be <= here because # of items and # of indeces differ, 
                                                     # so split_qry[3] when there are 3 items (0,1,2) will fail
                deletion = {}
                qry_codons = split_qry[qry_codon_count]

                if "-" in qry_codons:
                    beginning_flank = split_qry[qry_codon_count-1]
                    deletion[qry_codon_count-1] = beginning_flank

                    while qry_codon_count <= len(split_qry) and "-" in split_qry[qry_codon_count]:  # that <= may have to be a <, but i can't think straight right now
                        current_codon = split_qry[qry_codon_count]  # snapshot of the current codon
                        deletion[qry_codon_count] = current_codon
                        qry_codon_count += 1  # updates our index to the NEXT codon after successfully identifying a gap
                        
                    if qry_codon_count < len(split_qry):
                        ending_flank = split_qry[qry_codon_count]
                        deletion[qry_codon_count] = ending_flank
                        indel_denovo_result_HGVS.append(self.indel_translator(deletion, split_ref, translated_stripped_qry, indel_type = "deletion"))
                    else:
                        ending_flank = None

                ## for any other nucleotide in the sequence DO NOT COMMENT OUT  
                else:
                    qry_codon_count += 1

        """ insertions """
        ### isolate the position of the indel and the affected codons (plus flanking codons)
        if "-" in hsp_sbjct:
            ## split the subject sequence into a list of codons
            split_sbjct = re.findall('.'*3, hsp_sbjct)
            stripped_sbjct = hsp_sbjct.replace("-", "")
            
            ## translate the subject sequence into a protein (seq stripped of gaps because Seq hates them)
            translated_stripped_sbjct = str(Seq(stripped_sbjct).translate(table=11, gap="-"))

            ## iterate through query codon list, find gaps + flanks, and note position
            while sbjct_codon_count < len(split_sbjct):
                insertion = {}
                sbjct_codons = split_sbjct[sbjct_codon_count]

                if "-" in sbjct_codons:
                    sbjct_beginning_flank = split_sbjct[sbjct_codon_count-1]
                    insertion[sbjct_codon_count-1] = sbjct_beginning_flank

                    while sbjct_codon_count <= len(split_sbjct) and "-" in split_sbjct[sbjct_codon_count]:
                        sbjct_current_codon = split_sbjct[sbjct_codon_count]
                        insertion[sbjct_codon_count] = sbjct_current_codon
                        sbjct_codon_count += 1
                    if sbjct_codon_count < len(split_sbjct):
                        sbjct_ending_flank = split_sbjct[sbjct_codon_count]
                        insertion[sbjct_codon_count] = sbjct_ending_flank
                        indel_denovo_result_HGVS.append(self.indel_translator(insertion, split_ref, translated_stripped_sbjct, indel_type = "insertion"))
                    else:
                        sbjct_ending_flank = None

                ## for any other nucleotide in the sequence DO NOT COMMENT OUT  
                else:
                    sbjct_codon_count += 1

        if len(indel_denovo_result_HGVS) > 0:
            # you can change the output syntax here
            indel_result_prelim["query_def"] = str(query_def)
            indel_result_prelim["denovo_fs"] = indel_denovo_result_HGVS
            indel_result_prelim["has_indel"] = True
        else:
            indel_result_prelim["query_def"] = str(query_def)
            indel_result_prelim["has_indel"] = False

        return indel_result_prelim

    def single_fs(self, codon_count, translated_stripped_seq, split_ref):
        aa_pos = codon_count
        affected_codon = split_ref[aa_pos - 1]
        corr_aa = translated_stripped_seq[aa_pos - 1] # index starts at 0
        translated_codon = str(Seq(affected_codon).translate(table=11))

        return aa_pos, affected_codon, corr_aa, translated_codon ## aa_pos is just codon_count... fix that
    
    def termination(self, translated_stripped_seq, aa_pos):
        aa_count = 0

        ## locating the frameshift in the translated protein (entire seq. chunk + position of termination)
        for aa in translated_stripped_seq[aa_pos - 1:]: # index starts at 0
            if aa == "*":
                break
            else:
                aa_count += 1
        
        return aa_count + 1
    
    def indel_translator(self, indel, split_ref, translated_stripped_seq, indel_type=None):
        unpacked_indel = list(indel.items())
        inordel_codons = ""

        for i, (position, codon) in enumerate(unpacked_indel):  # we don't actually access codon... but you never know when you'll need it? :^)
            if i == 0:
                affected_codon = split_ref[position - 1]
                original_aa = str(Seq(affected_codon).translate(table=11))

                beginning_pos = position
                beginning_flank = original_aa
            elif i == len(unpacked_indel) - 1:
                affected_codon = split_ref[position - 1]
                original_aa = str(Seq(affected_codon).translate(table=11))

                end_pos = position
                end_flank = original_aa
            else:  # everything in between the sandwich
                new_aa = translated_stripped_seq[position - 1] # index starts at 0

                inordel_codons += new_aa

        if indel_type == "deletion":
            indel = f"{beginning_flank}{beginning_pos}_{end_flank}{end_pos}del{inordel_codons}"
            return indel
        elif indel_type == "insertion":
            indel = f"{beginning_flank}{beginning_pos}_{end_flank}{end_pos}ins{inordel_codons}"
            return indel

    def consolidate_mutations(self, input_type, hit_id, model_type, srv=None, fs=None, phm=None, hsp_bitscore=None, pass_val=None):
        """
        Consolidates the results from all mutation functions and passes the output back into RGI's detection modules (PHM, PVM, POM,    RGV).
        """

        has_snp = srv.get("has_snp", False) if srv is not None else False  # the PHM does not generate srv, so this avoids an AttributeError

        # protein input
        if input_type == "protein" and srv is not None:
            if has_snp:  # CASE 1: if the input is a protein there won't be a BLASTN xml generated
                return [srv]

            return []

        # nucleotide input
        else:
            if not fs:
                return []
            
            for fs_hit in fs:
                fs_id = fs_hit["query_def"].split()[0]
                has_curated_fs = "curated_fs" in fs_hit
                has_denovo_fs = "denovo_fs" in fs_hit  ## PHMs will only have de novo frameshifts
                                                       ## (nothing is curated for them right now)
                                                                                
                passes_eval = float(hsp_bitscore) >= float(pass_val)

                # protein variant frameshift/SNP search
                if model_type == "pvm":
                    srv_id = srv["query_def"].split()[0]

                    if has_snp:  # If SNPs were found
                        if fs_id in srv_id:
                            if not (has_curated_fs or has_denovo_fs):  # CASE 1 (only SNP)
                                return [srv]
                            else:  # CASE 2 (SNP and frameshift)
                                return [srv | fs_hit]
                            
                    else:  # CASE 3: frameshift found; SNPs were not found in any HSPs with SNP in the alignment title
                        if not has_snp:
                            if fs_id in srv_id:
                                if passes_eval and has_curated_fs:  # Strict alignments
                                    fs_hit["query_def"] += hit_id
                                    return [fs_hit]
                                elif not passes_eval and has_curated_fs:  # Loose alignments
                                    fs_hit["query_def"] += hit_id
                                    return [fs_hit]
                            return []
                        
                # protein overexpression frameshift/SNP search
                if model_type == "pom":
                    srv_id = srv["query_def"].split()[0]
                
                    if has_snp:  # If SNPs were found
                        if fs_id in srv_id:
                            if not (has_curated_fs or has_denovo_fs):  # CASE 1 (only SNP)
                                return [srv]
                            else:  # CASE 2 (SNP and frameshift)
                                return [srv | fs_hit]
                            
                    else:  # CASE 3: frameshift found; SNPs were not found in any HSPs with SNP in the alignment title
                        if not has_snp:
                            if fs_id in srv_id:
                                if passes_eval and has_curated_fs:  # Strict alignments
                                    fs_hit["query_def"] += hit_id
                                    return [fs_hit]
                                elif not passes_eval and has_curated_fs:  # Loose alignments
                                    fs_hit["query_def"] += hit_id
                                    return [fs_hit]
                            return []
                                
                # protein homolog frameshift search               
                if model_type == "phm":
                    phm_id = phm["query_def"].split()[0]
                
                    if fs_id in phm_id and has_denovo_fs and passes_eval:
                        # Perfect PHMs will not have frameshifts in them, so there's no need to add unique support for them here
                        fs_hit["query_def"] += hit_id
                        return [fs_hit]
                    return []