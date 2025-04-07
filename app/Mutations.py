from app.Base import BaseModel
from app.settings import *

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
                    # print(d)
                    query_snps = {"original": ori, "change": chan ,"position": d+1}
                    # print(query_snps)
                    # logger.debug("query_snp on frame {} {}".format(hsp.frame, json.dumps(query_snps, indent=2)))

                    srv_output["query_snps"] = query_snps
            # print(srv_output)
            # print("============================= NEXT SNP =============================")
                    
        # print(srv_output, "\n", qry, chan)
            srv_result.append(srv_output)

        return srv_result

    def frameshift(self, predicted_genes_dict_protein, submitted_proteins_dict, fsl, real_sbjct_length, hsp_query, hsp_sbjct_start, hsp_sbjct, orf_info, hit_id): 
        """
        Searches for frameshifts in sequences.
        """
    
        fs_dict_list = []
        fs_result = []

        # print(predicted_genes_dict_protein)
        # print(submitted_proteins_dict)

        for each_fs in fsl:
            if each_fs != None:
                # print(each_fs)
                # snp_dict_list.append({"original": each_snp[0], "change": each_snp[-1], "position": int(each_snp[1:-1])})
                position = int(
                    ''.join(filter(str.isdigit, each_fs)))

                original = (each_fs.split(
                    ''.join(filter(str.isdigit, each_fs))))
                
                fs_dict_list.append(
                    {"original": original[0], "position": position})
                
                # print(fs_dict_list)

        for eachfs in fs_dict_list:
            # print(each_fs)
            # print(eachfs)
            # print(eachfs, "and", hit_id)
            fs_output = {}

            pos = eachfs["position"]
            ori = eachfs["original"]
            # chan = eachfs["change"]
            # print("normal:", chan)

            # print(hsp_sbjct_start, pos, (hsp_sbjct_start + real_sbjct_length))

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
                        fs_output["eachs"] = eachfs
                        fs_output["orf_protein_sequence"] = orf_protein_sequence
                        # fs_output["chan"] = chan
                        # print(orf_protein_sequence)
                    else:
                        orf_protein_sequence = predicted_genes_dict_protein[orf_info.decode()[:orf_info.decode().index(' # ')]].strip("*")
                        fs_output["eachs"] = eachfs
                        fs_output["orf_protein_sequence"] = orf_protein_sequence
                        # fs_output["chan"] = chan
                        # print(orf_protein_sequence)
                    
                if submitted_proteins_dict:
                    orf_protein_sequence = str(submitted_proteins_dict[orf_info.decode().split(" ")[0]])
                    fs_output["eachs"] = eachfs
                    fs_output["orf_protein_sequence"] = orf_protein_sequence
                    # fs_output["chan"] = chan 
                    print(orf_protein_sequence)

        # return fs_result
