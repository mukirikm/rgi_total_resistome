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
    
        # fs_dict_list = []

        fs_curated_list_reg = []
        fs_denovo_list_reg = []

        fs_curated_result_HGVS = []
        fs_denovo_result_HGVS = []

        # for deletions
        qry_codon_count = 0

        # for insertions
        sbjct_codon_count = 0
                        
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
                        qry_codon_count += 1 # index starts at 1 not 0

                        if qry_codon_count <= len(translated_stripped_qry):
                            aa_pos, affected_codon, corr_aa, translated_codon = self.single_fs(qry_codon_count, translated_stripped_qry, split_ref)
                            fs_ter = self.termination(translated_stripped_qry, aa_pos)

                            for eachfs in fs_dict_list:
                                # print(eachfs)
                                if eachfs["original_aa"] == translated_codon and eachfs["aa_position"] == aa_pos:
                                    # logger.info("curated del fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                    # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": qry_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                    # fs_curated_list_reg.append("%s%s%s" % (translated_codon, aa_pos, corr_aa)) ## e.g., A15A
                                    fs_curated_list_reg.append("%s%sfs" % (translated_codon, aa_pos)) ## e.g., A15fs
                                    fs_curated_result_HGVS.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter)) ## e.g., A15AfsTer9
                            
                            if len(fs_curated_list_reg) != 0:
                                for _ in fs_curated_list_reg:
                                    if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list_reg and ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_denovo_list_reg:
                                        # logger.info("novel del fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                        # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": qry_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                        # fs_denovo_list_reg.append("%s%s%s" % (translated_codon, aa_pos, corr_aa)) ## e.g., A15A
                                        fs_denovo_list_reg.append("%s%sfs" % (translated_codon, aa_pos)) ## e.g., A15fs
                                        fs_denovo_result_HGVS.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter)) ## e.g., A15AfsTer9
                            else:
                                if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list_reg:
                                    # logger.info("novel del fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                    # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": qry_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                    # fs_denovo_list_reg.append("%s%s%s" % (translated_codon, aa_pos, corr_aa)) ## e.g., A15A
                                    fs_denovo_list_reg.append("%s%sfs" % (translated_codon, aa_pos)) ## e.g., A15fs
                                    fs_denovo_result_HGVS.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter)) ## e.g., A15AfsTer9

                    ## for any other nucleotide in the sequence DO NOT COMMENT OUT
                    else:
                        qry_codon_count += 1
                
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
                        sbjct_codon_count += 1 # index starts at 1 not 0
                        
                        if sbjct_codon_count <= len(translated_stripped_sbjct):
                            aa_pos, affected_codon, corr_aa, translated_codon = self.single_fs(sbjct_codon_count, translated_stripped_sbjct, split_ref)
                            fs_ter = self.termination(translated_stripped_sbjct, aa_pos)

                            for eachfs in fs_dict_list:
                                # print(eachfs)
                                if eachfs["original_aa"] == translated_codon and eachfs["aa_position"] == aa_pos:
                                    # logger.info("curated ins fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                    # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": sbjct_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                    # fs_curated_list_reg.append("%s%s%s" % (translated_codon, aa_pos, corr_aa)) ## e.g., A15A
                                    fs_curated_list_reg.append("%s%sfs" % (translated_codon, aa_pos)) ## e.g., A15fs
                                    fs_curated_result_HGVS.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter)) ## e.g., A15AfsTer9
                                
                            if len(fs_curated_list_reg) != 0:
                                for _ in fs_curated_list_reg:
                                    if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list_reg:
                                        # logger.info("novel ins fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                        # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": sbjct_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                        # fs_denovo_list_reg.append("%s%s%s" % (translated_codon, aa_pos, corr_aa)) ## e.g., A15A
                                        fs_denovo_list_reg.append("%s%sfs" % (translated_codon, aa_pos)) ## e.g., A15fs
                                        fs_denovo_result_HGVS.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter)) ## e.g., A15AfsTer9
                            else:
                                if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list_reg:
                                    # logger.info("novel ins fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                    # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": sbjct_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                    # fs_denovo_list_reg.append("%s%s%s" % (translated_codon, aa_pos, corr_aa)) ## e.g., A15A
                                    fs_denovo_list_reg.append("%s%sfs" % (translated_codon, aa_pos)) ## e.g., A15fs
                                    fs_denovo_result_HGVS.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter)) ## e.g., A15AfsTer9
                    else:
                        sbjct_codon_count += 1   

            fs_result_prelim["query_def"] = str(query_def)

            if len(fs_curated_result_HGVS) > 0 or len(fs_denovo_result_HGVS) > 0:
                fs_result_prelim["query_def"] = str(query_def)

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
                        qry_codon_count += 1 # index starts at 1 not 0

                        if qry_codon_count <= len(translated_stripped_qry):
                            aa_pos, affected_codon, corr_aa, translated_codon = self.single_fs(qry_codon_count, translated_stripped_qry, split_ref)
                            fs_ter = self.termination(translated_stripped_qry, aa_pos)

                            for _ in fs_curated_list_reg:
                                if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list_reg and ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_denovo_list_reg:
                                    # logger.info("novel del fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                    # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": qry_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                    # fs_denovo_list_reg.append("%s%s%s" % (translated_codon, aa_pos, corr_aa)) ## e.g., A15A
                                    fs_denovo_list_reg.append("%s%sfs" % (translated_codon, aa_pos)) ## e.g., A15fs
                                    fs_denovo_result_HGVS.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter)) ## e.g., A15AfsTer9
                        else:
                            if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list_reg:
                                # logger.info("novel del fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": qry_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                # fs_denovo_list_reg.append("%s%s%s" % (translated_codon, aa_pos, corr_aa)) ## e.g., A15A
                                fs_denovo_list_reg.append("%s%sfs" % (translated_codon, aa_pos)) ## e.g., A15fs
                                fs_denovo_result_HGVS.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter)) ## e.g., A15AfsTer9

                    ## for any other nucleotide in the sequence DO NOT COMMENT OUT
                    else:
                        qry_codon_count += 1

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
                        sbjct_codon_count += 1 # index starts at 1 not 0
                        
                        if sbjct_codon_count <= len(translated_stripped_sbjct):
                            aa_pos, affected_codon, corr_aa, translated_codon = self.single_fs(sbjct_codon_count, translated_stripped_sbjct, split_ref)
                            fs_ter = self.termination(translated_stripped_sbjct, aa_pos)

                        if len(fs_curated_list_reg) != 0:
                            for _ in fs_curated_list_reg:
                                if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list_reg:
                                    # logger.info("novel ins fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                    # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": sbjct_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                    # fs_denovo_list_reg.append("%s%s%s" % (translated_codon, aa_pos, corr_aa)) ## e.g., A15A
                                    fs_denovo_list_reg.append("%s%sfs" % (translated_codon, aa_pos)) ## e.g., A15fs
                                    fs_denovo_result_HGVS.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter)) ## e.g., A15AfsTer9
                        else:
                            if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list_reg:
                                # logger.info("novel ins fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": sbjct_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                # fs_denovo_list_reg.append("%s%s%s" % (translated_codon, aa_pos, corr_aa)) ## e.g., A15A
                                fs_denovo_list_reg.append("%s%sfs" % (translated_codon, aa_pos)) ## e.g., A15fs
                                fs_denovo_result_HGVS.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter)) ## e.g., A15AfsTer9
                    else:
                        sbjct_codon_count += 1   

            if len(fs_denovo_result_HGVS) > 0:
                fs_result_prelim["query_def"] = str(query_def)
                fs_result_prelim["denovo_fs"] = fs_denovo_result_HGVS
                fs_result_prelim["has_fs"] = True
            else:
                fs_result_prelim["query_def"] = str(query_def)
                fs_result_prelim["has_fs"] = False

            return fs_result_prelim

    def single_fs(self, codon_count, translated_stripped_seq, split_ref):
        aa_pos = codon_count
        affected_codon = split_ref[aa_pos - 1]
        # print("affected codon:", affected_codon)
        # print("trans stripped seq @ aapos:", translated_stripped_seq[aa_pos])
        # print("aa pos:", aa_pos)
        corr_aa = translated_stripped_seq[aa_pos - 1] # index starts at 0
        translated_codon = str(Seq(affected_codon).translate(table=11))

        return aa_pos, affected_codon, corr_aa, translated_codon
    
    def termination(self, translated_stripped_seq, aa_pos):
        aa_count = 0

        ## locating the frameshift in the translated protein (entire seq. chunk + position of termination)
        for aa in translated_stripped_seq[aa_pos - 1:]: # index starts at 0
            if aa == "*":
                break
            else:
                aa_count += 1
        
        return aa_count + 1

    def consolidate_mutations(self, input_type, hit_id, srv=None, fs=None, phm=None, hsp_bitscore=None, pass_val=None):
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
                if srv:
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

                            return[]
                                
                # protein homolog frameshift search               
                if phm:
                    phm_id = phm["query_def"].split()[0]
                
                    if fs_id in phm_id and has_denovo_fs and passes_eval:
                        # Perfect PHMs will not have frameshifts in them, so there's no need to add unique support for them here
                        fs_hit["query_def"] += hit_id
                        return [fs_hit]
                    
                    return []