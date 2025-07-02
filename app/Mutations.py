from app.Base import BaseModel
from app.settings import *
import re
from Bio.Seq import Seq

class MutationsModule(BaseModel):
    """
    Class for mutation searches (SNVs, frameshifts, coSNPs, etc.).
    """
			
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
        srv_result =[]

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

                if hsp_query[qry] == chan:
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
            # print(srv_output)
            # print("============================= NEXT SNP =============================")
                    
        # print(srv_output, "\n", qry, chan)
            srv_result.append(srv_output)

        return srv_result

    def frameshift(self, fsl, real_sbjct_length, hsp_query, hsp_sbjct_start, hsp_sbjct, hit_id, midline, card_prot_ref): 
        """
        Searches for frameshifts in sequences.
        """
    
        fs_dict_list = []
        fs_result = []

        fs_discovery_list = []
        fs_curated_list = []

        # for deletions
        qry_nucl_count = 0
        qry_gap_count = 0
        aa_count = 0

        # for insertions
        sbjct_nucl_count = 0
        sbjct_gap_count = 0

        ## grabbing curated frameshifts from blast XML (change to CARD JSON as input later?)
        for each_fs in fsl:
            # print(each_fs)
            # print(fsl)
            if each_fs != None:
                # print(each_fs)
                position = int(
                    ''.join(filter(str.isdigit, each_fs)))

                original = (each_fs.split(
                    ''.join(filter(str.isdigit, each_fs))))
                
                fs_dict_list.append(
                    {"original_aa": original[0], "aa_position": position})
                
                # unused values
                """({"curated_nucl_position": round(position*3), "mutation_type": original[-1]})"""
                
        print("========= OG frameshift dictionary list ==========")              
        print(fs_dict_list)
            
        """ 
        for nucleotide deletions 
        """            
        ## "discovery" -- finding de novo frameshifts
        ### counting the number of gaps, isolate the position of the gap and the affected codon (works for a single nucleotide deletion)
        if "-" in hsp_query:
            ## translate the query sequence into a protein with gaps (aka as is)
            print("========= translated protein (with gaps) ==========")  
            stripped_qry = hsp_query.replace("-", "")

            # split_sbjct = re.findall('.'*3, hsp_sbjct)

            translated_stripped_qry = str(Seq(stripped_qry).translate(table=11, gap="-"))
            print(translated_stripped_qry)
            print("========= all frameshift information from query sequence ==========") 

            ## iterate through query sequence, find gaps, note position, and grab all relevant information
            for qry_nucl in hsp_query:
                if qry_nucl == "-":
                    qry_nucl_count += 1 # index starts at 1 not 0
                    print(qry_nucl, qry_nucl_count)
                    qry_gap_count += 1

                    affected_sbjct_codon = hsp_sbjct[qry_nucl_count-3:qry_nucl_count] # grabbing the codon (from frameshift - 3)
                    translated_sbjct_codon = str(Seq(affected_sbjct_codon).translate(table=11))
                    affected_sbjct_nucl = hsp_sbjct[qry_nucl_count-1] # accounting for index (subject indexing starts at 0 in this case)

                    ## locating the frameshift in the translated protein (entire seq. chunk + position of termination)
                    sbjct_aa_pos = qry_nucl_count//3
                    corr_sbjct_aa = translated_stripped_qry[(sbjct_aa_pos)-1] # index starts at 0 here

                    for aa in translated_stripped_qry[sbjct_aa_pos-1:]: # index starts at 0 in this case
                        if aa == "*":
                            break
                        else:
                            aa_count += 1
                            # print(aa)
                    fs_ter = aa_count + 1

                    ## final desired frameshift output format
                    print("%s%s%sfsTer%s" % (translated_sbjct_codon, sbjct_aa_pos, corr_sbjct_aa, fs_ter))

                    for eachfs in fs_dict_list:
                        print("current curated frameshift:", eachfs)
                        if eachfs["original_aa"] == translated_sbjct_codon and eachfs["aa_position"] == sbjct_aa_pos:
                            print("curated fs found: %s%s%s" % (translated_sbjct_codon, sbjct_aa_pos, corr_sbjct_aa))
                            fs_curated_list.append("%s%s%s" % (translated_sbjct_codon, sbjct_aa_pos, corr_sbjct_aa))
                        else:
                            pass
                        
                    for curatedfs in fs_curated_list:
                        if ("%s%s%s" % (translated_sbjct_codon, sbjct_aa_pos, corr_sbjct_aa)) not in fs_curated_list:
                            print("novel fs found: %s%s%s" % (translated_sbjct_codon, sbjct_aa_pos, corr_sbjct_aa))
                            fs_discovery_list.append("%s%s%s" % (translated_sbjct_codon, sbjct_aa_pos, corr_sbjct_aa))
                                  
                    """({"affected_sbjct_nucleotide": affected_sbjct_nucl, "affected_sbjct_codon": affected_sbjct_codon,
                                          "translated_sbjct_aa": translated_sbjct_codon, "sbjct_nucl_position": qry_nucl_count, "sbjct_aa_position": sbjct_aa_pos, "new_aa": corr_sbjct_aa, "stop_position": fs_ter})"""
                else:
                    qry_nucl_count += 1   
                    # print(qry_nucl, qry_nucl_count)
            
        """
        ========================================================================================================================
                                                                INSERTIONS -- TBA
        ========================================================================================================================
        """

        """
        ========================================================================================================================
                                                                INSERTIONS -- TBA
        ========================================================================================================================
        """       
        print("========= curated and discovery framshift lists ==========") 
        print(fs_curated_list)
        print(fs_discovery_list)
            
        # return fs_result
