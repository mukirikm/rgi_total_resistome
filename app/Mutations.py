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
        self.srv_result = None
        self.fs_result = None

        self.hit_def_blastp = None
        self.hit_def_blastn = None
			
    def __repr__(self):
        """
        Returns Mutation class full object.
        """         
        return "Mutation({}".format(self.__dict__)

    def single_resistance_variant(self, predicted_genes_dict_protein, submitted_proteins_dict, snpl, 
                                  real_sbjct_length, hsp_query, hsp_sbjct_start, hsp_sbjct, orf_info, hit_def_blastp): 
        """
        Searches for SNVs in sequences.
        """

        snp_dict_list = []
        self.srv_result = None
        self.hit_def_blastp = hit_def_blastp
        # self.srv_result = [] # if we don't intialize this here, only the last SNP is considered
        # print("snp list\n", snpl)

        for each_snp in snpl:
            # print(each_snp)
            # snp_dict_list.append({"original": each_snp[0], "change": each_snp[-1], "position": int(each_snp[1:-1])})
            position = int(
                ''.join(filter(str.isdigit, each_snp)))

            original_change = (each_snp.split(
                ''.join(filter(str.isdigit, each_snp))))

            snp_dict_list.append(
                {"original": original_change[0], "change": original_change[-1], "position": position})
            
        # print("snp list=============")
        # print(snp_dict_list)
        # print()
            
        for eachs in snp_dict_list:
            # print(eachs)
            # print(eachs, "and", hit_id)
            srv_output = {}

            pos = eachs["position"]
            ori = eachs["original"]
            chan = eachs["change"]
            # print("normal:", chan)
            # print(pos, ori, chan)

            # print((hsp_sbjct_start), (hsp_sbjct_start + real_sbjct_length), pos)
            if hsp_sbjct_start < pos and (hsp_sbjct_start + real_sbjct_length) > pos:
                orf_protein_sequence = ""

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

                # wildtype
                wildtype = str(
                    hsp_sbjct[pos - hsp_sbjct_start + self.find_num_dash(hsp_sbjct, (pos-hsp_sbjct_start))])
                srv_output["wildtype"] = wildtype
                # print("wildtype:", wildtype)

                # Report ONLY if the SNPs are present
                qry = int(pos) - hsp_sbjct_start + self.find_num_dash(hsp_sbjct, (int(pos) - hsp_sbjct_start))
                srv_output["qry"] = qry
                # print("qry:", qry)
                # print(str(chan))
                # print(hsp_query[qry], chan, wildtype)

                # check for Var
                if str(chan) == "Var":
                    # print(eachs)
                    # update to the change
                    chan = str(
                        hsp_query[pos - hsp_sbjct_start + self.find_num_dash(hsp_sbjct, (pos-hsp_sbjct_start))])
                    # print(hsp_query[qry], chan, wildtype)
                    
                    if hsp_query[qry] == chan and chan != wildtype:
                        # update eachs@change
                        eachs["change"] = chan

                    else:
                        # change same as wildtype, don't report
                        chan = ""
                        eachs["change"] = chan

                if hsp_query[qry] == chan: # if the amino acid at our specific position in the query sequence is the same as the NEW aa (same SNP change has occured)
                    # print(eachs)
                    query_snps = {}

                    # get position of mutation in the query sequence
                    d = int(
                        pos) - hsp_sbjct_start - self.find_num_dash(hsp_query, (int(pos) - hsp_sbjct_start))
                    # print(pos, hsp_sbjct_start, hsp_query)
                    # print(d, pos, hsp_query[qry], qry)
                    query_snps = {
                        "original": ori, "change": chan ,"position": d+1}
                    # print(query_snps)
                    # logger.debug("query_snp on frame {} {}".format(hsp.frame, json.dumps(query_snps, indent=2)))

                    srv_output["query_snps"] = query_snps
                    # print(query_snps)

                    # print(srv_output)
                    # self.srv_result.append(srv_output)
                    self.srv_result = srv_output

                    # print("srv output")
                    # print(self.srv_result,"\n")
                    return self.srv_result, self.hit_def_blastp
        
    def frameshift(self, fsl, hsp_query, hsp_sbjct, card_dna_ref, hit_def_blastn): 
        """
        Searches for frameshifts in sequences.
        """

        fs_result_prelim = {}
        self.fs_result = None
        self.hit_def_blastn = hit_def_blastn
        # self.fs_result = []
    
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
            # print(type(each_fs), each_fs)
            position = int(
                ''.join(filter(str.isdigit, each_fs)))

            original = (each_fs.split(
                ''.join(filter(str.isdigit, each_fs))))
            
            fs_dict_list.append(
                {"original_aa": original[0], "aa_position": position})
                
        # print(fs_dict_list)
        # print(len(fs_dict_list))
                
        split_ref = re.findall('.'*3, card_dna_ref)
        # print("split CARD reference below")
        # print(split_ref)
         
        if len(fs_dict_list) != 0:
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
                # print("qry:", len(translated_stripped_qry))
                # print("translated stripped qry below")
                # print(translated_stripped_qry)

                ## iterate through query codon list, find gaps, note position, and grab all relevant information
                for qry_codons in split_qry:
                    if "-" in qry_codons:
                        qry_codon_count += 1 # index starts at 1 not 0

                        if qry_codon_count >= len(translated_stripped_qry):
                            aa_pos, affected_codon, corr_aa, translated_codon = self.single_fs(qry_codon_count, translated_stripped_qry, split_ref)
                            fs_ter = self.termination(translated_stripped_qry, aa_pos)

                            for eachfs in fs_dict_list:
                                # print("current curated frameshift:", eachfs)
                                if eachfs["original_aa"] == translated_codon and eachfs["aa_position"] == aa_pos:
                                    # logger.info("curated del fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                    # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": qry_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                    fs_curated_list.append("%s%s%s" % (translated_codon, aa_pos, corr_aa))
                                    fs_curated_result.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter))
                                else:
                                    pass
                            
                            if len(fs_curated_list) != 0:
                                for curatedfs in fs_curated_list:
                                    if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list and ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_denovo_list:
                                        # logger.info("novel del fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                        # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": qry_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                        fs_denovo_list.append("%s%s%s" % (translated_codon, aa_pos, corr_aa))
                                        fs_denovo_result.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter)) 
                            else:
                                if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list:
                                    # logger.info("novel del fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                    # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": qry_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                    fs_denovo_result.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter))

                        ## for any other nucleotide in the sequence DO NOT COMMENT OUT
                        else:
                            qry_codon_count += 1
                
            # try:
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
                # print("sbjct:", len(translated_stripped_sbjct))
                # print("translated stripped sb")
                # print(translated_stripped_sbjct)

                ## iterate through subject codon list, find gaps, note position, and grab all relevant information
                for sbjct_codons in split_sbjct:
                    if "-" in sbjct_codons:
                        sbjct_codon_count += 1 # index starts at 1 not 0
                        # print(sbjct_codon_count)
                        
                        if sbjct_codon_count <= len(translated_stripped_sbjct):
                            aa_pos, affected_codon, corr_aa, translated_codon = self.single_fs(sbjct_codon_count, translated_stripped_sbjct, split_ref)
                            fs_ter = self.termination(translated_stripped_sbjct, aa_pos)

                            for eachfs in fs_dict_list:
                                # print("current curated frameshift:", eachfs)
                                if eachfs["original_aa"] == translated_codon and eachfs["aa_position"] == aa_pos:
                                    # logger.info("curated ins fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                    # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": sbjct_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                    fs_curated_list.append("%s%s%s" % (translated_codon, aa_pos, corr_aa))
                                    fs_curated_result.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter))
                                else:
                                    pass
                                
                            if len(fs_curated_list) != 0:
                                for curatedfs in fs_curated_list:
                                    if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list:
                                        # logger.info("novel ins fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                        # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": sbjct_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                        fs_denovo_list.append("%s%s%s" % (translated_codon, aa_pos, corr_aa))
                                        fs_denovo_result.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter))
                            else:
                                if ("%s%s%s" % (translated_codon, aa_pos, corr_aa)) not in fs_curated_list:
                                    # logger.info("novel ins fs found: %s%s%s" % (translated_codon, aa_pos, corr_aa))
                                    # logger.info({"affected_codon": affected_codon, "translated_ref_aa": translated_codon, "ref_nucl_position": sbjct_codon_count, "aa_position": aa_pos, "new_aa": corr_aa, "stop_position": fs_ter})
                                    fs_denovo_result.append("%s%s%sfsTer%s" % (translated_codon, aa_pos, corr_aa, fs_ter))
                        else:
                            sbjct_codon_count += 1   
                            # print(sbjct_codon_count)
            # except Exception as e:
            #     traceback.print_exc()

            if len(fs_curated_result) > 0:
                fs_result_prelim["curated_fs"] = fs_curated_result
            if len(fs_denovo_result) > 0:
                fs_result_prelim["denovo_fs"] = fs_denovo_result

            # self.fs_result.append(fs_result_prelim)
            self.fs_result = fs_result_prelim
            # print("fs output")
            # print(self.fs_result,"\n")
            return self.fs_result, self.hit_def_blastn
            # print(fs_result_prelim)
                
            # print(self.fs_result)
        # else:
        #     # self.fs_result = []
        #     # print("fs output")
        #     # print(self.fs_result,"\n")
        #     self.fs_result = None
        #     return self.fs_result,  hit_def_blastn

    def single_fs(self, codon_count, translated_stripped_seq, split_ref):
        aa_pos = codon_count
        affected_codon = split_ref[aa_pos - 1]
        # print(affected_codon)
        # print("trans stripped seq @ aapos:", translated_stripped_seq[aa_pos])
        # print("aa pos:", aa_pos)
        corr_aa = translated_stripped_seq[aa_pos - 1] # index starts at 0
        # print(corr_aa)
        # print()
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

    def consolidate_mutations(self, input_type):
        print("blastn hit def:\n", self.hit_def_blastn)
        print("blastp hit def:\n", self.hit_def_blastp)
        if input_type == "protein": # if the input is a protein there won't be a BLASTN xml generated, thus, no frameshift output
            logger.info("no frameshift result generated (protein input); only SNV result exists")
            return [self.srv_result]
        else:
            # print("fs_resutlt:\n", self.fs_result)
            # print("srv_resutlt:\n", self.srv_result)
            if self.srv_result != None and self.fs_result != None:
                logger.info("SNV and frameshift result exists")
                # print([self.srv_result | self.fs_result])
                return [self.srv_result | self.fs_result]
                # for srv_result_unpacked in self.srv_result:
                #         print([srv_result_unpacked | self.fs_result])
                #     return [srv_result_unpacked | self.fs_result]
            
            elif self.fs_result != None and self.srv_result == None:
                logger.info("no frameshift result; only SNV result exists")
                return [self.srv_result]
            
            elif self.srv_result != None and self.fs_result == None:
                logger.info("no SNV result; only frameshift result exists")
                return [self.fs_result]