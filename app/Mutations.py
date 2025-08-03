import re
import warnings

from Bio import BiopythonWarning
from Bio.Seq import Seq

from app.Base import BaseModel
from app.settings import *

warnings.simplefilter("ignore", BiopythonWarning)

class MutationsModule(BaseModel):
    """
    Class for mutation searches (SNVs, frameshifts, coSNPs, etc.).
    """

    def __init__ (self):
        self.srv_result = None
        self.fs_result = None
			
    def __repr__(self):
        """
        Returns Mutation class full object.
        """         
        return "Mutation({}".format(self.__dict__)

    def single_resistance_variant(self, predicted_genes_dict_protein, submitted_proteins_dict, snpl, real_sbjct_length, hsp_query, hsp_sbjct_start, hsp_sbjct, orf_info): 
        """
        Searches for SNVs in sequences.
        """
    
        snp_dict_list = []
        self.srv_result = []

        for each_snp in snpl:
            # print(each_snp)
            # snp_dict_list.append({"original": each_snp[0], "change": each_snp[-1], "position": int(each_snp[1:-1])})
            position = int(
                ''.join(filter(str.isdigit, each_snp)))

            original_change = (each_snp.split(
                ''.join(filter(str.isdigit, each_snp))))

            snp_dict_list.append(
                {"original": original_change[0], "change": original_change[-1], "position": position})
            
            # print(snp_dict_list)
            
        # print("-------------------------------- SNP DICTIONARY/LIST ORIGINAL ----------------------------------------------------")
        # for snp_cluster in snp_dict_list:
        #     print(snp_cluster)
        # print("-------------------------------- SNP DICTIONARY/LIST ORIGINAL ----------------------------------------------------")

        for eachs in snp_dict_list:
            # print(eachs)
            # print(eachs, "and", hit_id)
            srv_output = {}

            pos = eachs["position"]
            ori = eachs["original"]
            chan = eachs["change"]
            # print("normal:", chan)

            # wildtype
            # wildtype = str(
            #     hsp.sbjct[pos - hsp.sbjct_start + self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))])

            # check for Var
            if str(chan) == "Var":
                # update to the change
                chan = str(
                    hsp_query[pos - hsp_sbjct_start + self.find_num_dash(hsp_sbjct, (pos-hsp_sbjct_start))])
                # update eachs
                eachs["change"] = chan
                # print("var:", chan)

            if hsp_sbjct_start < pos and (hsp_sbjct_start + real_sbjct_length) > pos:
                orf_protein_sequence = ""

                # if predicted_genes_dict:
                #     if orf_info.strip() in predicted_genes_dict.keys():
                #         orf_protein_sequence = str(Seq(predicted_genes_dict[orf_info.decode()], generic_dna).translate(table=11)).strip("*")
                #     else:
                #         orf_protein_sequence = str(Seq(predicted_genes_dict[orf_info.decode()[:orf_info.decode().index(' # ')]], generic_dna).translate(table=11)).strip("*")

                if predicted_genes_dict_protein:
                    if orf_info.strip() in predicted_genes_dict_protein.keys():
                        orf_protein_sequence = predicted_genes_dict_protein[orf_info.decode()].strip("*")
                        srv_output["eachs"] = eachs
                        srv_output["orf_protein_sequence"] = orf_protein_sequence
                        srv_output["chan"] = chan
                        # print(orf_protein_sequence)
                    else:
                        orf_protein_sequence = predicted_genes_dict_protein[orf_info.decode()[:orf_info.decode().index(' # ')]].strip("*")
                        srv_output["eachs"] = eachs
                        srv_output["orf_protein_sequence"] = orf_protein_sequence
                        srv_output["chan"] = chan
                        # print(orf_protein_sequence)

                if submitted_proteins_dict:
                    # print("debug:", eachs, "and", hit_id)
                    # print(submitted_proteins_dict)
                    # print()
                    orf_protein_sequence = str(submitted_proteins_dict[orf_info.decode().split(" ")[0]])
                    srv_output["eachs"] = eachs
                    srv_output["orf_protein_sequence"] = orf_protein_sequence
                    srv_output["chan"] = chan 
                   # print(orf_protein_sequence)

                # logger.info("mutation | Model:"+str(model_id) + " | pos:" +str(pos) +" | change: "+str(hsp.query[pos - hsp.sbjct_start + \
                # 			self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(chan) + " AND wildtype: " + str(hsp.sbjct[pos - hsp.sbjct_start \
                # 			+self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(ori))

                # Report ONLY if the SNPs are present
                qry = int(pos) - hsp_sbjct_start + self.find_num_dash(hsp_sbjct, (int(pos) - hsp_sbjct_start))
                srv_output["qry"] = qry

                sbj = int(pos) - hsp_sbjct_start + self.find_num_dash(hsp_sbjct, (int(pos) - hsp_sbjct_start))
                # print(hsp_query[qry], "+", chan, ":", qry, ">", sbj)

                if hsp_query[qry] == chan: # if the amino acid at our specific position in the query sequence is the same as the NEW aa (same SNP change has occured)
                    # print(eachs)
                    query_snps = {}
                    # logger.debug("mutation | Model:"+str(model_id) + " | pos:" +str(pos) +" | change: "+str(hsp.query[pos - hsp.sbjct_start + \
                    # 		self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(chan) + " AND wildtype: " + str(hsp.sbjct[pos - hsp.sbjct_start \
                    # 		+self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(ori))

                    # get position of mutation in the query sequence
                    d = int(pos) - hsp_sbjct_start - self.find_num_dash(hsp_query, (int(pos) - hsp_sbjct_start))
                    # print(pos, hsp_sbjct_start, hsp_query)
                    # print(d, pos, hsp_query[qry], qry)
                    query_snps = {"original": ori, "change": chan ,"position": d+1}
                    # print(query_snps)
                    # logger.debug("query_snp on frame {} {}".format(hsp.frame, json.dumps(query_snps, indent=2)))

                    srv_output["query_snps"] = query_snps
     
                    self.srv_result.append(srv_output)

        # if len(self.srv_result) > 0: # this doesn't work...
        return self.srv_result

    def frameshift(self, fsl, hsp_query, hsp_sbjct, card_dna_ref): 
        """
        Searches for frameshifts in sequences.
        """

        self.fs_result = {}
    
        fs_dict_list = []

        fs_curated_list = []
        fs_denovo_list = []

        fs_curated_result = []
        fs_denovo_result = []

        # for deletions
        qry_codon_count = 0

        # for insertions
        sbjct_codon_count = 0

        ## grabbing curated frameshifts from blast XML (change to CARD JSON as input later?)
        for each_fs in fsl:
            if each_fs != None:
                # print(each_fs)
                position = int(
                    ''.join(filter(str.isdigit, each_fs)))

                original = (each_fs.split(
                    ''.join(filter(str.isdigit, each_fs))))
                
                fs_dict_list.append(
                    {"original_aa": original[0], "aa_position": position})
                
        split_ref = re.findall('.'*3, card_dna_ref)
        # print("split CARD reference below")
        # print(split_ref)
            
        """for nucleotide deletions """            
        ### isolate the position of the gap(s) and the affected codon(s)
        if "-" in hsp_query:
            ## split the query sequence into a list of codons
            split_qry = re.findall('.'*3, hsp_query)
            # print("\nsplit qry below")
            # print(split_qry)
            stripped_qry = hsp_query.replace("-", "")

            ## translate the query sequence into a protein (seq stripped of gaps because Seq hates them)
            translated_stripped_qry = str(Seq(stripped_qry).translate(table=11))
            # print("translated stripped qry below")
            # print(translated_stripped_qry)

            ## iterate through query codon list, find gaps, note position, and grab all relevant information
            for qry_codons in split_qry:
                if "-" in qry_codons:
                    qry_codon_count += 1 # index starts at 1 not 0

                    aa_pos, affected_codon, corr_aa, translated_codon = self.single_fs(qry_codon_count, translated_stripped_qry, split_ref)
                    fs_ter = self.termination(translated_stripped_qry, aa_pos)

                    for eachfs in fs_dict_list:
                        # print("current curated frameshift:", eachfs)
                        if eachfs["original_aa"] == translated_codon and eachfs["aa_position"] == aa_pos:
                            logger.info("curated del fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                            logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": qry_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                            fs_curated_list.append("%s%s%s" % (translated_codon, aa_pos, corr_aa))
                            fs_curated_result.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter))
                        else:
                            pass
                    
                    if len(fs_curated_list) != 0:
                        for curatedfs in fs_curated_list:
                            if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list and ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_denovo_list:
                                logger.info("novel del fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": qry_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                fs_denovo_list.append("%s%s%s" % (translated_codon, aa_pos, corr_aa))
                                fs_denovo_result.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter)) 
                    else:
                        if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list:
                            logger.info("novel del fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                            logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": qry_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                            fs_denovo_result.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter))

                ## for any other nucleotide in the sequence DO NOT COMMENT OUT
                else:
                    qry_codon_count += 1   
            
        """for nucleotide insertions"""
        ### isolate the position of the gap(s) and the affected codon(s)
        if "-" in hsp_sbjct:
            ## split the subject sequence into a list of codons
            split_sbjct = re.findall('.'*3, hsp_sbjct)
            # print("\nsplit sbjct below")
            # print(split_sbjct)
            stripped_sbjct = hsp_sbjct.replace("-", "")
            
            ## translate the subject sequence into a protein (seq stripped of gaps because Seq hates them)
            translated_stripped_sbjct = str(Seq(stripped_sbjct).translate(table=11, gap="-"))
            # print("translated stripped sb")
            # print(translated_stripped_sbjct)

            ## iterate through subject codon list, find gaps, note position, and grab all relevant information
            for sbjct_codons in split_sbjct:
                if "-" in sbjct_codons:
                    sbjct_codon_count += 1 # index starts at 1 not 0

                    aa_pos, affected_codon, corr_aa, translated_codon = self.single_fs(sbjct_codon_count, translated_stripped_sbjct, split_ref)
                    fs_ter = self.termination(translated_stripped_sbjct, aa_pos)

                    for eachfs in fs_dict_list:
                        # print("current curated frameshift:", eachfs)
                        if eachfs["original_aa"] == translated_codon and eachfs["aa_position"] == aa_pos:
                            logger.info("curated ins fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                            logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": sbjct_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                            fs_curated_list.append("%s%s%s" % (translated_codon, aa_pos, corr_aa))
                            fs_curated_result.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter))
                        else:
                            pass
                        
                    if len(fs_curated_list) != 0:
                        for curatedfs in fs_curated_list:
                            if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list:
                                logger.info("novel ins fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": sbjct_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                fs_denovo_list.append("%s%s%s" % (translated_codon, aa_pos, corr_aa))
                                fs_denovo_result.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter))
                    else:
                        if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list:
                            logger.info("novel ins fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                            logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": sbjct_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                            fs_denovo_result.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter))

                else:
                    sbjct_codon_count += 1   

        if len(fs_curated_result) > 0:
            self.fs_result["curated_fs"] = fs_curated_result
        if len(fs_denovo_result) > 0:
            self.fs_result["denovo_fs"] = fs_denovo_result
            
        return self.fs_result

    def single_fs(self, codon_count, translated_stripped_seq, split_ref):
        aa_pos = codon_count
        affected_codon = split_ref[aa_pos - 1]
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

    def consolidate_mutations(self):
        try:
            if self.srv_result:
                for srv_result_dict in self.srv_result:
                    return [srv_result_dict | self.fs_result]
            else:
                return self.fs_result
        except AttributeError: # if the input is a protein there won't be a BLASTN xml generated, thus, no frameshift output
             return self.srv_result