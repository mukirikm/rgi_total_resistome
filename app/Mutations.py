from app.Base import BaseModel
from app.settings import *
from Bio.Blast import NCBIXML

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
        srv_output = {}

        for each_snp in snpl:
            snp_dict_list.append({"original": each_snp[0], "change": each_snp[-1], "position": int(each_snp[1:-1])})
        
        for eachs in snp_dict_list:
            pos = eachs["position"]
            ori = eachs["original"]
            chan = eachs["change"]

            srv_output["eachs"] = eachs

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
                        srv_output["orf_protein_sequence"] = orf_protein_sequence
                    else:
                        orf_protein_sequence = predicted_genes_dict_protein[orf_info.decode()[:orf_info.decode().index(' # ')]].strip("*")
                        srv_output["orf_protein_sequence"] = orf_protein_sequence

                if submitted_proteins_dict:
                    orf_protein_sequence = str(submitted_proteins_dict[orf_info.decode().split(" ")[0]])
                    srv_output["orf_protein_sequence"] = orf_protein_sequence

                # logger.info("mutation | Model:"+str(model_id) + " | pos:" +str(pos) +" | change: "+str(hsp.query[pos - hsp.sbjct_start + \
                # 			self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(chan) + " AND wildtype: " + str(hsp.sbjct[pos - hsp.sbjct_start \
                # 			+self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(ori))

                # Report ONLY if the SNPs are present
                qry = int(pos) - hsp_sbjct_start + self.find_num_dash(hsp_sbjct, (int(pos) - hsp_sbjct_start))
                sbj = int(pos) - hsp_sbjct_start + self.find_num_dash(hsp_sbjct, (int(pos) - hsp_sbjct_start))

                if hsp_query[qry] == chan:
                    query_snps = {}
                    # logger.debug("mutation | Model:"+str(model_id) + " | pos:" +str(pos) +" | change: "+str(hsp.query[pos - hsp.sbjct_start + \
                    # 		self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(chan) + " AND wildtype: " + str(hsp.sbjct[pos - hsp.sbjct_start \
                    # 		+self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(ori))

                    # get position of mutation in the query sequence
                    d = int(pos) - hsp_sbjct_start - self.find_num_dash(hsp_query, (int(pos) - hsp_sbjct_start))
                    query_snps = {"original": ori, "change": chan ,"position": d+1}
                    # logger.debug("query_snp on frame {} {}".format(hsp.frame, json.dumps(query_snps, indent=2)))

                    srv_output["query_snps"] = query_snps
            return srv_output

# def frameshift(self):
#     pass