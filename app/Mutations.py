from app.Base import BaseModel
from app.settings import *
# from app.VM import Variant
from Bio.Blast import NCBIXML

class MutationsModule(BaseModel):
    """
    Class for mutation searches (SNVs, frameshifts, coSNPs, etc.).
    """
	

    def __init__(self, xml_file, input_type, input_sequence, working_directory, align_title=""):
        # self.xml_file = xml_file
        self.input_type = input_type
        self.input_sequence = input_sequence
        self.working_directory = working_directory
        self.xml_file = xml_file

        # self.align_title = MutationsModule.baseline()
        self.align_title = align_title
        # self.orf_info = MutationsModule.baseline(self)
        # self.orf_from = MutationsModule.baseline(self)

        # self.model_type_id = MutationsModule.baseline(self)
        # print("__innit__")
			
    def __repr__(self):
        """
        Returns Mutation class full object.
        """         
        return "Mutation({}".format(self.__dict__)

    # @staticmethod
    # def baseline(self):
    #     with open(self.xml_file, 'r') as result_handle:
    #         blast_records = NCBIXML.parse(result_handle)
    #         for blast_record in blast_records:
    #             perfect = {}
    #             strict = {}
    #             loose = {}
    #             for self.alignment in blast_record.alignments:
    #                 self.align_title = self.alignment.title
    #                 self.orf_info = blast_record.query.encode('ascii','replace')

    #                 # return self.alignment
                    
    #                 # c = 0
    #                 # barc = 0
    #                 # for eachc in self.orf_info:
    #                 #     if barc >= 6:
    #                 #         break
    #                 #     elif eachc == '|':
    #                 #         barc += 1
    #                 #         c += 1
    #                 #     else:
    #                 #         c += 1
    #                 # self.orf_from = self.orf_info[c:]

    #                 model_type_id = self.extract_nth_bar(self.align_title, 0)
    #                 # # logger.info("model_type_id: {} ".format(model_type_id))
    #                 # space_pos = self.align_title.index(' ')
    #                 # hit_id = self.align_title[0:space_pos]
    #                 # hit_id = hit_id.encode('ascii','replace')
    #                 # model_descrpt = self.align_title[self.align_title.index(' ')+1:]
    #                 # underscore_in_MD = model_descrpt.index('_')
    #                 # self.model_id = model_descrpt[0:underscore_in_MD]
    #                 # seq_in_model = model_descrpt[underscore_in_MD+1: model_descrpt.index(' ')]
    #                 # pass_value = self.extract_nth_bar(alignment.title, 1)
    #                 # # logger.info("pass_value: {}".format(pass_value))

    #                 if model_type_id == 40293:
    #                     test_modelid = self.single_resistance_variant()
    #                     return test_modelid, self.align_title
    
    def indiv_snps(self, snp_dict_list):
        """
        Spits out individual SNPs in input data.
        """
        for eachs in snp_dict_list:
            return eachs


    def single_resistance_variant(self, predicted_genes_dict_protein, submitted_proteins_dict, snp_dict_list, real_query_length, real_sbjct_length, hsp_query, hsp_sbjct_start, hsp_sbjct, hsp_frame, orf_info): 
        """
        Searches for SNVs in sequences.
        """
        # predicted_genes_dict = {}
        # print("snp dict list\n", snp_dict_list,
        #       "\nsubmitted protein dict\n", submitted_proteins_dict,
        #       "\nreal query length\n", real_query_length, 
        #       "\nreal subject length\n", real_sbjct_length,  
        #       "\nhsp query\n", hsp_query,
        #       "\nhsp subject start\n", hsp_sbjct_start,  
        #       "\nhsp subject\n", hsp_sbjct,  
        #       "\nhsp frame\n", hsp_frame, 
        #       "\norf info\n", orf_info)
        
        for eachs in snp_dict_list:
            pos = eachs["position"]
            ori = eachs["original"]
            chan = eachs["change"]

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
                        yield orf_protein_sequence
                        # print(orf_protein_sequence)
                    else:
                        orf_protein_sequence = predicted_genes_dict_protein[orf_info.decode()[:orf_info.decode().index(' # ')]].strip("*")
                        yield orf_protein_sequence
                        # print(orf_protein_sequence)

                if submitted_proteins_dict:
                    orf_protein_sequence = str(submitted_proteins_dict[orf_info.decode().split(" ")[0]])
                    yield orf_protein_sequence
                    # print(orf_protein_sequence)

                # logger.info("mutation | Model:"+str(model_id) + " | pos:" +str(pos) +" | change: "+str(hsp.query[pos - hsp.sbjct_start + \
                # 			self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(chan) + " AND wildtype: " + str(hsp.sbjct[pos - hsp.sbjct_start \
                # 			+self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(ori))

                # Report ONLY if the SNPs are present
                qry = int(pos) - hsp_sbjct_start + self.find_num_dash(hsp_sbjct, (int(pos) - hsp_sbjct_start))
                sbj = int(pos) - hsp_sbjct_start + self.find_num_dash(hsp_sbjct, (int(pos) - hsp_sbjct_start))

                # print(qry, sbj)

                # print(chan)
                # print(hsp_query[qry])

                if hsp_query[qry] == chan:
                    # print("hsp query:", hsp_query[qry], "AND", "chan:", chan)
                    query_snps = {}
                    # logger.debug("mutation | Model:"+str(model_id) + " | pos:" +str(pos) +" | change: "+str(hsp.query[pos - hsp.sbjct_start + \
                    # 		self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(chan) + " AND wildtype: " + str(hsp.sbjct[pos - hsp.sbjct_start \
                    # 		+self.find_num_dash(hsp.sbjct, (pos-hsp.sbjct_start))]) + "=" + str(ori))

                    # get position of mutation in the query sequence
                    d = int(pos) - hsp_sbjct_start - self.find_num_dash(hsp_query, (int(pos) - hsp_sbjct_start))
                    # print("hsp.sbjct_start: ", hsp.sbjct_start)
                    query_snps = {"original": ori, "change": chan ,"position": d+1}
                    # logger.debug("query_snp on frame {} {}".format(hsp.frame, json.dumps(query_snps, indent=2)))

                    yield query_snps
            #         # print(query_snps)
            # return eachs

# # def frameshift(self):
# #     pass